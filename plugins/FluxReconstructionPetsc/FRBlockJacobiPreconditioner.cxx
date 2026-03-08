// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Petsc/PetscHeaders.hh" // must come before any header

#include "FluxReconstructionPetsc/FRBlockJacobiPreconditioner.hh"
#include "Petsc/PetscLSSData.hh"
#include "FluxReconstructionPetsc/FluxReconstructionPetsc.hh"

#include "Framework/CFL.hh"
#include "Framework/MeshData.hh"
#include "Framework/PhysicalModel.hh"
#include "Framework/SpaceMethod.hh"
#include "Framework/SpaceMethodData.hh"
#include "Framework/SubSystemStatus.hh"
#include "Framework/MethodStrategyProvider.hh"
#include "MathTools/MatrixInverter.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Petsc {

    extern PetscErrorCode FRBlockJacobiPcApply(PC pc, Vec X, Vec Y);

//////////////////////////////////////////////////////////////////////////////

using COOLFluiD::FluxReconstructionMethod::FluxReconstructionPetscModule;

MethodStrategyProvider<FRBlockJacobiPreconditioner,
                       PetscLSSData,
                       ShellPreconditioner,
                       FluxReconstructionPetscModule>
FRBlockJacobiPreconditionerProvider("FRBlockJacobi");

//////////////////////////////////////////////////////////////////////////////

FRBlockJacobiPreconditioner::FRBlockJacobiPreconditioner(const std::string& name) :
  ShellPreconditioner(name),
  socket_updateCoeff("updateCoeff"),
  socket_volumes("volumes", false),
  _pcc(),
  _inverter()
{
}

//////////////////////////////////////////////////////////////////////////////

FRBlockJacobiPreconditioner::~FRBlockJacobiPreconditioner()
{
}

//////////////////////////////////////////////////////////////////////////////

void FRBlockJacobiPreconditioner::setPreconditioner()
{
  CFLog(VERBOSE, "FRBlockJacobiPreconditioner::setPreconditioner()\n");

  _pcc.pJFC = getMethodData().getJFContext();
  _pcc.nEqs = getMethodData().getNbSysEquations();

  DataHandle<State*, GLOBAL> states = _pcc.pJFC->states->getDataHandle();
  const CFuint nbStates = states.size();

  // build state local ID -> updatable state index mapping
  std::vector<CFint> stateToUpIdx(nbStates, -1);
  CFuint nbUpdatable = 0;
  for (CFuint i = 0; i < nbStates; ++i)
  {
    if (states[i]->isParUpdatable())
    {
      stateToUpIdx[i] = nbUpdatable++;
    }
  }

  // build cell-to-updatable-state-index mapping
  SafePtr<TopologicalRegionSet> cells =
    MeshDataStack::getActive()->getTrs("InnerCells");
  const CFuint nbCells = cells->getLocalNbGeoEnts();

  _pcc.cellToUpStateIdx.clear();
  _pcc.cellToUpStateIdx.reserve(nbCells);
  _pcc.nUpdatableCells = 0;

  CFuint maxSolPtsPerCell = 0;

  for (CFuint iCell = 0; iCell < nbCells; ++iCell)
  {
    const CFuint nbStatesInCell = cells->getNbStatesInGeo(iCell);

    // check if this cell has at least one updatable state
    bool isUpdatable = false;
    for (CFuint iSt = 0; iSt < nbStatesInCell; ++iSt)
    {
      const CFuint stateID = cells->getStateID(iCell, iSt);
      if (stateToUpIdx[stateID] >= 0)
      {
        isUpdatable = true;
        break;
      }
    }

    if (isUpdatable)
    {
      std::vector<CFuint> upIndices;
      upIndices.reserve(nbStatesInCell);
      for (CFuint iSt = 0; iSt < nbStatesInCell; ++iSt)
      {
        const CFuint stateID = cells->getStateID(iCell, iSt);
        if (stateToUpIdx[stateID] >= 0)
        {
          upIndices.push_back(static_cast<CFuint>(stateToUpIdx[stateID]));
        }
      }
      if (upIndices.size() > maxSolPtsPerCell)
      {
        maxSolPtsPerCell = upIndices.size();
      }
      _pcc.cellToUpStateIdx.push_back(upIndices);
      _pcc.nUpdatableCells++;
    }
  }

  // allocate inverse element blocks
  _pcc.invElemBlocks.resize(_pcc.nUpdatableCells);
  const CFuint nEqs = _pcc.nEqs;
  CFuint totalBlockMem = 0;
  for (CFuint i = 0; i < _pcc.nUpdatableCells; ++i)
  {
    const CFuint nSolPts = _pcc.cellToUpStateIdx[i].size();
    const CFuint blockSize = nSolPts * nEqs;
    _pcc.invElemBlocks[i].resize(blockSize, blockSize);
    _pcc.invElemBlocks[i] = 0.0;
    totalBlockMem += blockSize * blockSize * sizeof(CFreal);
  }

  // create matrix inverter for the largest element block
  const CFuint maxBlockSize = maxSolPtsPerCell * nEqs;
  _inverter.reset(MatrixInverter::create(maxBlockSize, false));

  CFLog(INFO, "FRBlockJacobiPreconditioner: " << _pcc.nUpdatableCells
    << " updatable cells, " << nEqs << " equations, max "
    << maxSolPtsPerCell << " sol pts/cell, block size "
    << maxBlockSize << "x" << maxBlockSize
    << ", total block memory: " << totalBlockMem / (1024.0*1024.0) << " MB\n");

  // register the shell preconditioner with PETSc
  CF_CHKERRCONTINUE(PCShellSetContext(
    _pcc.pJFC->petscData->getPreconditioner(), &_pcc));
  CF_CHKERRCONTINUE(PCShellSetApply(
    _pcc.pJFC->petscData->getPreconditioner(), FRBlockJacobiPcApply));
}

//////////////////////////////////////////////////////////////////////////////

void FRBlockJacobiPreconditioner::computeBeforeSolving()
{
  CFLog(VERBOSE, "FRBlockJacobiPreconditioner::computeBeforeSolving() => START\n");

  SafePtr<SpaceMethod> spaceMtd = _pcc.pJFC->spaceMethod;
  SafePtr<SpaceMethodData> spaceData = spaceMtd->getSpaceMethodData();

  // ---- Step 1: Save flags ----
  const bool savedDoComputeJacob = spaceData->doComputeJacobian();
  const bool savedFillPrecondMat = spaceData->fillPreconditionerMatrix();

  // ---- Step 2: Backup updateCoeff and rhs ----
  DataHandle<CFreal> updateCoeff = _pcc.pJFC->updateCoeff->getDataHandle();
  DataHandle<CFreal> rhs = _pcc.pJFC->rhs->getDataHandle();
  DataHandle<State*, GLOBAL> states = _pcc.pJFC->states->getDataHandle();
  const CFuint nbStates = states.size();
  const CFuint nEqs = _pcc.nEqs;
  const CFuint rhsSize = nbStates * nEqs;

  RealVector bkpUpdateCoeff(nbStates);
  for (CFuint i = 0; i < nbStates; ++i)
  {
    bkpUpdateCoeff[i] = updateCoeff[i];
  }

  // Backup socket_rhs: computeSpaceResidual/computeTimeResidual will overwrite it.
  // Although computeJFMat recomputes rhs before reading it, the backup ensures
  // no stale data persists if the code flow changes in the future.
  RealVector bkpRhs(rhsSize);
  for (CFuint i = 0; i < rhsSize; ++i)
  {
    bkpRhs[i] = rhs[i];
  }

  // ---- Step 3: Enable Jacobian assembly into preconditioner matrix ----
  spaceData->setComputeJacobianFlag(true);
  spaceData->setFillPreconditionerMatrix(true);

  // ---- Step 4: Zero the preconditioner matrix ----
  PetscMatrix& precondMat = _pcc.pJFC->petscData->getPreconditionerMatrix();
  precondMat.resetToZeroEntries();

  // ---- Step 5: Reset updateCoeff and compute residuals ----
  // This fills the preconditioner matrix with the full Jacobian
  updateCoeff = 0.0;
  spaceMtd->computeSpaceResidual(1.0);
  spaceMtd->computeTimeResidual(1.0);

  // ---- Step 6: Finalize preconditioner matrix assembly ----
  precondMat.finalAssembly();

  // ---- Step 7: Restore flags ----
  spaceData->setComputeJacobianFlag(savedDoComputeJacob);
  spaceData->setFillPreconditionerMatrix(savedFillPrecondMat);

  // ---- Step 8: Restore updateCoeff and rhs ----
  for (CFuint i = 0; i < nbStates; ++i)
  {
    updateCoeff[i] = bkpUpdateCoeff[i];
  }
  for (CFuint i = 0; i < rhsSize; ++i)
  {
    rhs[i] = bkpRhs[i];
  }

  // ---- Step 9: Extract and invert per-cell diagonal blocks ----

  // FRBlockJacobi operates on a SeqBAIJ matrix (useBlockPreconditionerMatrix=true),
  // which uses 0-based local indexing. The cellToUpStateIdx already stores
  // locally-updatable indices (0-based), so we use upIdx * nEqs + iEq directly.
  // Note: the old code used upStatesGlobalIDs (MPI-global indices), which works
  // in serial (global==local) but fails in parallel where global IDs have a
  // rank-specific offset that doesn't match SeqBAIJ's local indexing.

  // Track current inverter size to handle mixed meshes (different nSolPts per cell)
  CFuint currentInverterSize = 0;

  for (CFuint iCell = 0; iCell < _pcc.nUpdatableCells; ++iCell)
  {
    const std::vector<CFuint>& upStateIndices = _pcc.cellToUpStateIdx[iCell];
    const CFuint nSolPts = upStateIndices.size();
    const CFuint blockSize = nSolPts * nEqs;

    if (nSolPts == 0) continue;

    // Build row/col index arrays for MatGetValues using local 0-based indices
    std::vector<CFint> matIndices(blockSize);
    for (CFuint iSol = 0; iSol < nSolPts; ++iSol)
    {
      const CFuint upIdx = upStateIndices[iSol];
      for (CFuint iEq = 0; iEq < nEqs; ++iEq)
      {
        matIndices[iSol * nEqs + iEq] = static_cast<CFint>(upIdx * nEqs + iEq);
      }
    }

    // Extract the diagonal block from the preconditioner matrix
    // Result is in row-major order: block[i*blockSize + j]
    RealMatrix elemBlock(blockSize, blockSize);
    precondMat.getValues(blockSize, &matIndices[0],
                         blockSize, &matIndices[0],
                         &elemBlock[0]);

    // Recreate inverter if block size differs (mixed meshes with different nSolPts)
    if (blockSize != currentInverterSize)
    {
      _inverter.reset(MatrixInverter::create(blockSize, false));
      currentInverterSize = blockSize;
    }

    // Invert the element block
    _inverter->invert(elemBlock, _pcc.invElemBlocks[iCell]);
  }

  CFLog(VERBOSE, "FRBlockJacobiPreconditioner::computeBeforeSolving() => DONE\n");
}

//////////////////////////////////////////////////////////////////////////////

void FRBlockJacobiPreconditioner::computeAfterSolving()
{
  // reset inverse element blocks
  for (CFuint i = 0; i < _pcc.nUpdatableCells; ++i)
  {
    _pcc.invElemBlocks[i] = 0.0;
  }
}

//////////////////////////////////////////////////////////////////////////////

PetscErrorCode FRBlockJacobiPcApply(PC pc, Vec X, Vec Y)
{
  // X is the input vector (Krylov vector to be preconditioned)
  // Y is the output vector (preconditioned result)
  //
  // Element block-Jacobi: for each cell, gather the cell's DOFs from X,
  // multiply by the inverted element diagonal block, scatter to Y.

  const CFreal* x;
  CFreal* y;

  PetscFunctionBegin;

  void* ctx;
  CF_CHKERRCONTINUE(PCShellGetContext(pc, &ctx));

  FRBlockJacobiPcJFContext* pcContext = (FRBlockJacobiPcJFContext*)(ctx);

  const CFuint nEqs = pcContext->nEqs;
  const CFuint nCells = pcContext->nUpdatableCells;

  CF_CHKERRCONTINUE(VecGetArrayRead(X, &x));
  CF_CHKERRCONTINUE(VecGetArray(Y, &y));

  // initialize Y to zero
  PetscInt ySize;
  CF_CHKERRCONTINUE(VecGetLocalSize(Y, &ySize));
  for (PetscInt i = 0; i < ySize; ++i)
  {
    y[i] = 0.0;
  }

  for (CFuint iCell = 0; iCell < nCells; ++iCell)
  {
    const std::vector<CFuint>& upStateIndices = pcContext->cellToUpStateIdx[iCell];
    const CFuint nSolPts = upStateIndices.size();

    if (nSolPts == 0) continue;

    const CFuint blockSize = nSolPts * nEqs;
    const RealMatrix& invBlock = pcContext->invElemBlocks[iCell];

    // Gather: build local x_block from the PETSc vector
    // Apply: y_block = invBlock * x_block
    // Scatter: write y_block back to PETSc vector
    for (CFuint iRow = 0; iRow < blockSize; ++iRow)
    {
      const CFuint iSol = iRow / nEqs;
      const CFuint iEq  = iRow % nEqs;
      const CFuint yIdx = upStateIndices[iSol] * nEqs + iEq;
      cf_assert(yIdx < (CFuint)ySize);

      CFreal val = 0.0;
      for (CFuint jCol = 0; jCol < blockSize; ++jCol)
      {
        const CFuint jSol = jCol / nEqs;
        const CFuint jEq  = jCol % nEqs;
        const CFuint xIdx = upStateIndices[jSol] * nEqs + jEq;

        val += invBlock(iRow, jCol) * x[xIdx];
      }
      y[yIdx] = val;
    }
  }

  CF_CHKERRCONTINUE(VecRestoreArrayRead(X, &x));
  CF_CHKERRCONTINUE(VecRestoreArray(Y, &y));

  PetscFunctionReturn(0);
}

//////////////////////////////////////////////////////////////////////////////

vector<SafePtr<BaseDataSocketSink> > FRBlockJacobiPreconditioner::needsSockets()
{
  vector<SafePtr<BaseDataSocketSink> > result = ShellPreconditioner::needsSockets();
  result.push_back(&socket_updateCoeff);
  result.push_back(&socket_volumes);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace Petsc

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
