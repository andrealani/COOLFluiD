// Copyright (C) 2016 KU Leuven, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Petsc/PetscHeaders.hh" // must come before any header

#include "FluxReconstructionPetsc/FRP0Preconditioner.hh"
#include "Petsc/PetscLSSData.hh"
#include "FluxReconstructionPetsc/FluxReconstructionPetsc.hh"
#include "FluxReconstructionMethod/FluxReconstructionSolverData.hh"

#include "Framework/CFL.hh"
#include "Framework/MeshData.hh"
#include "Framework/PhysicalModel.hh"
#include "Framework/SpaceMethod.hh"
#include "Framework/SpaceMethodData.hh"
#include "Framework/SubSystemStatus.hh"
#include "Framework/MethodStrategyProvider.hh"
#include "Common/PE.hh"
#include "MathTools/MatrixInverter.hh"

#include <map>
#include <cmath>

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Petsc {

    extern PetscErrorCode FRP0PcApply(PC pc, Vec X, Vec Y);

//////////////////////////////////////////////////////////////////////////////

using COOLFluiD::FluxReconstructionMethod::FluxReconstructionPetscModule;

MethodStrategyProvider<FRP0Preconditioner,
                       PetscLSSData,
                       ShellPreconditioner,
                       FluxReconstructionPetscModule>
FRP0PreconditionerProvider("FRP0Precond");

//////////////////////////////////////////////////////////////////////////////

void FRP0Preconditioner::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption<CFreal>("SmootherOmega",
    "Relaxation parameter for the smoother (default 1.0)");
  options.addConfigOption<std::string>("BlockMode",
    "Smoother block mode: ElementBlock (default) or PointBlock (per-DOF)");
  options.addConfigOption<bool>("DirectBlocks",
    "Compute element blocks directly without PETSc matrix (default true)");
  options.addConfigOption<std::string>("CoarseSolveType",
    "P0 coarse solve type: FaceCoupled (Galerkin-projected sparse matrix + direct LU) "
    "or BlockDiag (element-local block inversion, default)");
  options.addConfigOption<bool>("Multiplicative",
    "Use multiplicative (defect-based) coarse correction instead of additive (default false)");
}

//////////////////////////////////////////////////////////////////////////////

FRP0Preconditioner::FRP0Preconditioner(const std::string& name) :
  ShellPreconditioner(name),
  socket_updateCoeff("updateCoeff"),
  socket_volumes("volumes", false),
  _pcc(),
  _inverter(),
  _p0Inverter(),
  _smootherOmega(1.0),
  _blockMode("ElementBlock"),
  _directBlocks(true),
  _coarseSolveType("BlockDiag"),
  _multiplicative(false)
{
  addConfigOptionsTo(this);
  setParameter("SmootherOmega", &_smootherOmega);
  setParameter("BlockMode", &_blockMode);
  setParameter("DirectBlocks", &_directBlocks);
  setParameter("CoarseSolveType", &_coarseSolveType);
  setParameter("Multiplicative", &_multiplicative);
}

//////////////////////////////////////////////////////////////////////////////

FRP0Preconditioner::~FRP0Preconditioner()
{
  if (_pcc.useP0FaceCoupled)
  {
    if (_pcc.p0Ksp) { KSPDestroy(&_pcc.p0Ksp); _pcc.p0Ksp = CFNULL; }
    if (_pcc.p0Mat) { MatDestroy(&_pcc.p0Mat); _pcc.p0Mat = CFNULL; }
    if (_pcc.p0RhsVec) { VecDestroy(&_pcc.p0RhsVec); _pcc.p0RhsVec = CFNULL; }
    if (_pcc.p0SolVec) { VecDestroy(&_pcc.p0SolVec); _pcc.p0SolVec = CFNULL; }
  }
  if (_pcc.defectVec) { VecDestroy(&_pcc.defectVec); _pcc.defectVec = CFNULL; }
}

//////////////////////////////////////////////////////////////////////////////

void FRP0Preconditioner::setPreconditioner()
{
  CFLog(VERBOSE, "FRP0Preconditioner::setPreconditioner()\n");

  _pcc.pJFC = getMethodData().getJFContext();
  _pcc.nEqs = getMethodData().getNbSysEquations();
  _pcc.nDim = PhysicalModelStack::getActive()->getDim();

  // Parse block mode
  _pcc.useElementBlocks = (_blockMode == "ElementBlock");
  CFLog(INFO, "FRP0Preconditioner: BlockMode = " << _blockMode
    << " (useElementBlocks = " << _pcc.useElementBlocks << ")"
    << ", DirectBlocks = " << _directBlocks << "\n");

  DataHandle<State*, GLOBAL> states = _pcc.pJFC->states->getDataHandle();
  const CFuint nbStates = states.size();

  // Build state local ID -> updatable state index mapping
  std::vector<CFint> stateToUpIdx(nbStates, -1);
  CFuint nbUpdatable = 0;
  for (CFuint i = 0; i < nbStates; ++i)
  {
    if (states[i]->isParUpdatable())
    {
      stateToUpIdx[i] = nbUpdatable++;
    }
  }

  // Build cell-to-updatable-state-index mapping
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
    const CFuint firstStateID = cells->getStateID(iCell, 0);

    if (!states[firstStateID]->isParUpdatable()) continue;

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
      maxSolPtsPerCell = upIndices.size();

    _pcc.cellToUpStateIdx.push_back(upIndices);
    _pcc.nUpdatableCells++;
  }

  _pcc.nSolPtsPerCell = maxSolPtsPerCell;
  _pcc.smootherOmega = _smootherOmega;

  const CFuint nEqs = _pcc.nEqs;
  const CFuint nCells = _pcc.nUpdatableCells;

  // Allocate smoother blocks
  CFuint totalBlockMem = 0;
  if (_pcc.useElementBlocks)
  {
    // ElementBlock: full (nSolPts*nEqs)^2 per cell
    _pcc.invElemBlocks.resize(nCells);
    for (CFuint i = 0; i < nCells; ++i)
    {
      const CFuint nSolPts = _pcc.cellToUpStateIdx[i].size();
      const CFuint blockSize = nSolPts * nEqs;
      _pcc.invElemBlocks[i].resize(blockSize, blockSize);
      _pcc.invElemBlocks[i] = 0.0;
      totalBlockMem += blockSize * blockSize * sizeof(CFreal);
    }

    // Create inverter for largest element block
    const CFuint maxBlockSize = maxSolPtsPerCell * nEqs;
    _inverter.reset(MatrixInverter::create(maxBlockSize, false));
  }
  else
  {
    // PointBlock: per-DOF (nEqs x nEqs) blocks
    _pcc.smootherInvBlocks.resize(nbUpdatable);
    for (CFuint i = 0; i < nbUpdatable; ++i)
    {
      _pcc.smootherInvBlocks[i].resize(nEqs, nEqs, 0.0);
    }
    totalBlockMem = nbUpdatable * nEqs * nEqs * sizeof(CFreal);

    // Create inverter for nEqs x nEqs blocks
    _inverter.reset(MatrixInverter::create(nEqs, false));
  }

  // ---- P0 coarse level setup ----
  _pcc.useP0FaceCoupled = (_coarseSolveType == "FaceCoupled");
  _pcc.p0Residual.resize(nCells * nEqs, 0.0);
  _pcc.p0Correction.resize(nCells * nEqs, 0.0);

  if (_pcc.useP0FaceCoupled)
  {
    // ---- Build cell-to-updatable-cell-index mapping ----
    _pcc.cellToUpdIdx.assign(nbCells, -1);
    CFuint updIdx = 0;
    for (CFuint iCell = 0; iCell < nbCells; ++iCell)
    {
      const CFuint firstStateID = cells->getStateID(iCell, 0);
      if (states[firstStateID]->isParUpdatable())
      {
        _pcc.cellToUpdIdx[iCell] = updIdx++;
      }
    }
    cf_assert(updIdx == nCells);

    // ---- Build cell adjacency from node-sharing ----
    // Two cells are face-neighbors if they share >= nDim nodes
    const CFuint nDim = _pcc.nDim;

    // Step 1: find max node ID and build node-to-cell mapping
    CFuint maxNodeID = 0;
    for (CFuint iCell = 0; iCell < nbCells; ++iCell)
    {
      for (CFuint iNode = 0; iNode < cells->getNbNodesInGeo(iCell); ++iNode)
      {
        const CFuint nodeID = cells->getNodeID(iCell, iNode);
        if (nodeID > maxNodeID) maxNodeID = nodeID;
      }
    }

    std::vector<std::vector<CFuint> > nodeToCells(maxNodeID + 1);
    for (CFuint iCell = 0; iCell < nbCells; ++iCell)
    {
      for (CFuint iNode = 0; iNode < cells->getNbNodesInGeo(iCell); ++iNode)
      {
        nodeToCells[cells->getNodeID(iCell, iNode)].push_back(iCell);
      }
    }

    // Step 2: for each updatable cell, find face neighbors
    _pcc.cellNeighbors.resize(nCells);
    CFuint updCell = 0;
    for (CFuint iCell = 0; iCell < nbCells; ++iCell)
    {
      if (_pcc.cellToUpdIdx[iCell] < 0) continue;

      // Count shared nodes with every other cell
      std::map<CFuint, CFuint> sharedCount;
      for (CFuint iNode = 0; iNode < cells->getNbNodesInGeo(iCell); ++iNode)
      {
        const CFuint nodeID = cells->getNodeID(iCell, iNode);
        for (CFuint j = 0; j < nodeToCells[nodeID].size(); ++j)
        {
          const CFuint otherCell = nodeToCells[nodeID][j];
          if (otherCell != iCell)
          {
            sharedCount[otherCell]++;
          }
        }
      }

      // Face neighbor = shares >= nDim nodes AND is updatable
      for (std::map<CFuint, CFuint>::iterator it = sharedCount.begin();
           it != sharedCount.end(); ++it)
      {
        if (it->second >= nDim && _pcc.cellToUpdIdx[it->first] >= 0)
        {
          _pcc.cellNeighbors[updCell].push_back(
            static_cast<CFuint>(_pcc.cellToUpdIdx[it->first]));
        }
      }
      updCell++;
    }

    // ---- Create P0 sparse matrix (SeqBAIJ with nEqs block size) ----
    // Compute per-row nnz (number of block entries per row)
    std::vector<PetscInt> p0nnz(nCells);
    CFuint totalP0Nnz = 0;
    for (CFuint i = 0; i < nCells; ++i)
    {
      p0nnz[i] = 1 + static_cast<PetscInt>(_pcc.cellNeighbors[i].size()); // diag + neighbors
      totalP0Nnz += p0nnz[i];
    }

    const PetscInt p0Size = static_cast<PetscInt>(nCells * nEqs);
    CF_CHKERRCONTINUE(MatCreateSeqBAIJ(PETSC_COMM_SELF,
      static_cast<PetscInt>(nEqs), p0Size, p0Size,
      0, &p0nnz[0], &_pcc.p0Mat));
    CF_CHKERRCONTINUE(MatSetOption(_pcc.p0Mat, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE));

    // ---- Create P0 vectors ----
    CF_CHKERRCONTINUE(VecCreateSeq(PETSC_COMM_SELF, p0Size, &_pcc.p0RhsVec));
    CF_CHKERRCONTINUE(VecCreateSeq(PETSC_COMM_SELF, p0Size, &_pcc.p0SolVec));

    // ---- Create inner KSP (direct LU solve) ----
    CF_CHKERRCONTINUE(KSPCreate(PETSC_COMM_SELF, &_pcc.p0Ksp));
    CF_CHKERRCONTINUE(KSPSetType(_pcc.p0Ksp, KSPPREONLY));

    PC p0pc;
    CF_CHKERRCONTINUE(KSPGetPC(_pcc.p0Ksp, &p0pc));
    CF_CHKERRCONTINUE(PCSetType(p0pc, PCLU));

    const CFuint p0MatMem = totalP0Nnz * nEqs * nEqs * sizeof(CFreal);
    totalBlockMem += p0MatMem;

    CFLog(INFO, "FRP0Preconditioner: P0 face-coupled coarse solve (direct LU)"
      << ", P0 matrix: " << nCells << "x" << nCells << " blocks of "
      << nEqs << "x" << nEqs
      << ", avg neighbors/cell: "
      << (nCells > 0 ? (CFreal)(totalP0Nnz - nCells) / nCells : 0)
      << ", P0 matrix memory: " << p0MatMem / (1024.0*1024.0) << " MB\n");
  }
  else
  {
    // ---- BlockDiag mode: allocate P0 inverse blocks ----
    _pcc.p0InvBlocks.resize(nCells);
    for (CFuint i = 0; i < nCells; ++i)
    {
      _pcc.p0InvBlocks[i].resize(nEqs, nEqs, 0.0);
    }
    totalBlockMem += nCells * nEqs * nEqs * sizeof(CFreal);

    // Create P0 inverter (nEqs x nEqs)
    _p0Inverter.reset(MatrixInverter::create(nEqs, false));
  }

  CFLog(INFO, "FRP0Preconditioner: " << nCells
    << " updatable cells, " << nEqs << " equations, "
    << maxSolPtsPerCell << " sol pts/cell, omega = " << _smootherOmega
    << ", CoarseSolveType = " << _coarseSolveType
    << ", total block memory: " << totalBlockMem / (1024.0*1024.0) << " MB\n");

  // Multiplicative mode: defect Vec allocated lazily in FRP0PcApply (needs matching parallel layout)
  _pcc.useMultiplicative = _multiplicative;
  if (_pcc.useMultiplicative)
  {
    CFLog(INFO, "FRP0Preconditioner: multiplicative mode enabled\n");
  }

  // Register the shell preconditioner with PETSc
  CF_CHKERRCONTINUE(PCShellSetContext(
    _pcc.pJFC->petscData->getPreconditioner(), &_pcc));
  CF_CHKERRCONTINUE(PCShellSetApply(
    _pcc.pJFC->petscData->getPreconditioner(), FRP0PcApply));
}

//////////////////////////////////////////////////////////////////////////////

void FRP0Preconditioner::computeBeforeSolving()
{
  CFLog(VERBOSE, "FRP0Preconditioner::computeBeforeSolving() => assembling Jacobian + extracting blocks"
    << " (DirectBlocks = " << _directBlocks << ")\n");

  const CFuint nEqs = _pcc.nEqs;
  const CFuint nCells = _pcc.nUpdatableCells;

  // Pass omega to context
  _pcc.smootherOmega = _smootherOmega;

  SafePtr<SpaceMethod> spaceMtd = _pcc.pJFC->spaceMethod;
  SafePtr<SpaceMethodData> spaceData = spaceMtd->getSpaceMethodData();

  // Downcast to FR solver data for block assembly methods
  using FluxReconstructionMethod::FluxReconstructionSolverData;
  SafePtr<FluxReconstructionSolverData> frData =
    spaceData.d_castTo<FluxReconstructionSolverData>();

  // Save flags
  const bool savedDoComputeJacob = spaceData->doComputeJacobian();
  const bool savedFillPrecondMat = spaceData->fillPreconditionerMatrix();

  // Backup updateCoeff and rhs
  DataHandle<CFreal> updateCoeff = _pcc.pJFC->updateCoeff->getDataHandle();
  DataHandle<CFreal> rhs = _pcc.pJFC->rhs->getDataHandle();
  DataHandle<State*, GLOBAL> states = _pcc.pJFC->states->getDataHandle();
  const CFuint nbStates = states.size();
  const CFuint rhsSize = nbStates * nEqs;

  RealVector bkpUpdateCoeff(nbStates);
  for (CFuint i = 0; i < nbStates; ++i)
  {
    bkpUpdateCoeff[i] = updateCoeff[i];
  }

  RealVector bkpRhs(rhsSize);
  for (CFuint i = 0; i < rhsSize; ++i)
  {
    bkpRhs[i] = rhs[i];
  }

  if (_directBlocks)
  {
    // ======== DIRECT BLOCKS MODE: bypass PETSc matrix entirely ========

    SafePtr<TopologicalRegionSet> cells =
      MeshDataStack::getActive()->getTrs("InnerCells");
    const CFuint nbAllCells = cells->getLocalNbGeoEnts();

    // P0 off-diagonal map (if FaceCoupled mode)
    std::map<std::pair<CFuint,CFuint>, RealMatrix> p0OffDiagMap;
    if (_pcc.useP0FaceCoupled)
    {
      frData->setP0OffDiagBlocks(&p0OffDiagMap);
    }

    // Mode-specific storage (only one is populated)
    std::vector<RealMatrix> cellBlocks;                    // ElementBlock mode
    std::vector<std::vector<RealMatrix>> pointBlocks;      // PointBlock mode

    if (_pcc.useElementBlocks)
    {
      // Allocate full element-diagonal blocks (temporary, freed after extraction)
      cellBlocks.resize(nbAllCells);
      for (CFuint i = 0; i < nbAllCells; ++i)
      {
        const CFuint nSolPts = cells->getNbStatesInGeo(i);
        const CFuint blockSize = nSolPts * nEqs;
        cellBlocks[i].resize(blockSize, blockSize, 0.0);
      }

      updateCoeff = 0.0;
      frData->computeDiagBlocks(cellBlocks, &(*spaceMtd));
    }
    else
    {
      // Allocate compact per-DOF blocks (temporary, freed after extraction)
      pointBlocks.resize(nbAllCells);
      for (CFuint i = 0; i < nbAllCells; ++i)
      {
        const CFuint nSolPts = cells->getNbStatesInGeo(i);
        pointBlocks[i].resize(nSolPts);
        for (CFuint j = 0; j < nSolPts; ++j)
          pointBlocks[i][j].resize(nEqs, nEqs, 0.0);
      }

      updateCoeff = 0.0;
      frData->computePointDiagBlocks(pointBlocks, &(*spaceMtd));
    }

    // Common: clear P0 off-diag pointer, zero P0 matrix, restore backup
    if (_pcc.useP0FaceCoupled) frData->setP0OffDiagBlocks(CFNULL);

    if (_pcc.useP0FaceCoupled)
    {
      CF_CHKERRCONTINUE(MatZeroEntries(_pcc.p0Mat));
    }

    for (CFuint i = 0; i < nbStates; ++i)
      updateCoeff[i] = bkpUpdateCoeff[i];
    for (CFuint i = 0; i < rhsSize; ++i)
      rhs[i] = bkpRhs[i];

    // ---- Extract, project, and invert blocks ----
    CFuint updCellIdx = 0;
    RealMatrix p0Block(nEqs, nEqs);
    CFuint currentInverterSize = 0;

    for (CFuint iCell = 0; iCell < nbAllCells; ++iCell)
    {
      const CFuint firstStateID = cells->getStateID(iCell, 0);
      if (!states[firstStateID]->isParUpdatable()) continue;

      const CFuint nSolPts = cells->getNbStatesInGeo(iCell);
      const CFuint blockSize = nSolPts * nEqs;

      if (nSolPts == 0) { updCellIdx++; continue; }

      if (_pcc.useElementBlocks)
      {
        // ElementBlock: extract from full element-diagonal blocks
        RealMatrix& elemBlock = cellBlocks[iCell];

        // Compute P0 block via Galerkin projection
        p0Block = 0.0;
        for (CFuint i = 0; i < nSolPts; ++i)
          for (CFuint j = 0; j < nSolPts; ++j)
            for (CFuint e = 0; e < nEqs; ++e)
              for (CFuint f = 0; f < nEqs; ++f)
                p0Block(e, f) += elemBlock(i * nEqs + e, j * nEqs + f);

        const CFreal invNSol = 1.0 / (CFreal)nSolPts;
        for (CFuint e = 0; e < nEqs; ++e)
          for (CFuint f = 0; f < nEqs; ++f)
            p0Block(e, f) *= invNSol;

        // Invert element block for smoother
        if (blockSize != currentInverterSize)
        {
          _inverter.reset(MatrixInverter::create(blockSize, false));
          currentInverterSize = blockSize;
        }
        {
          CFreal blockNorm = 0.0;
          const CFuint bsz = elemBlock.size();
          for (CFuint _k = 0; _k < bsz; ++_k) {
            const CFreal av = std::abs(elemBlock.ptr()[_k]);
            if (av > blockNorm) blockNorm = av;
          }
          RealMatrix& invBlock = _pcc.invElemBlocks[updCellIdx];
          if (blockNorm < 1e-14) {
            invBlock = 0.0;
            for (CFuint k = 0; k < blockSize; ++k) invBlock(k, k) = 1.0;
          } else {
            _inverter->invert(elemBlock, invBlock);
          }
        }
      }
      else
      {
        // PointBlock: extract from compact per-DOF blocks (no full element block needed)
        const std::vector<CFuint>& upStateIndices = _pcc.cellToUpStateIdx[updCellIdx];
        p0Block = 0.0;

        for (CFuint iSol = 0; iSol < nSolPts; ++iSol)
        {
          RealMatrix& dofBlock = pointBlocks[iCell][iSol];

          // Accumulate for P0
          for (CFuint e = 0; e < nEqs; ++e)
            for (CFuint f = 0; f < nEqs; ++f)
              p0Block(e, f) += dofBlock(e, f);

          // Invert DOF block for smoother
          const CFuint upIdx = upStateIndices[iSol];
          {
            CFreal blockNorm = 0.0;
            for (CFuint _k = 0; _k < dofBlock.size(); ++_k) {
              const CFreal av = std::abs(dofBlock.ptr()[_k]);
              if (av > blockNorm) blockNorm = av;
            }
            RealMatrix& invBlock = _pcc.smootherInvBlocks[upIdx];
            if (blockNorm < 1e-14) {
              invBlock = 0.0;
              for (CFuint k = 0; k < nEqs; ++k) invBlock(k, k) = 1.0;
            } else {
              _inverter->invert(dofBlock, invBlock);
            }
          }
        }

        const CFreal invNSol = 1.0 / (CFreal)nSolPts;
        for (CFuint e = 0; e < nEqs; ++e)
          for (CFuint f = 0; f < nEqs; ++f)
            p0Block(e, f) *= invNSol;
      }

      // ---- P0 coarse level: store or invert the projected block ----
      if (_pcc.useP0FaceCoupled)
      {
        PetscInt blockRow = static_cast<PetscInt>(updCellIdx);
        CF_CHKERRCONTINUE(MatSetValuesBlocked(_pcc.p0Mat,
          1, &blockRow, 1, &blockRow, &p0Block[0], INSERT_VALUES));
      }
      else
      {
        {
          CFreal blockNorm = 0.0;
          for (CFuint _k = 0; _k < p0Block.size(); ++_k) {
            const CFreal av = std::abs(p0Block.ptr()[_k]);
            if (av > blockNorm) blockNorm = av;
          }
          RealMatrix& invBlock = _pcc.p0InvBlocks[updCellIdx];
          if (blockNorm < 1e-14) {
            invBlock = 0.0;
            for (CFuint k = 0; k < nEqs; ++k) invBlock(k, k) = 1.0;
          } else {
            _p0Inverter->invert(p0Block, invBlock);
          }
        }
      }

      updCellIdx++;
    }

    // Insert actual off-diagonal blocks from Galerkin-projected face cross-blocks
    if (_pcc.useP0FaceCoupled)
    {
      for (auto it = p0OffDiagMap.begin(); it != p0OffDiagMap.end(); ++it)
      {
        const CFuint cellI = it->first.first;   // TRS-local cell index
        const CFuint cellJ = it->first.second;  // TRS-local cell index

        // Map TRS-local cell indices to updatable cell indices
        if (_pcc.cellToUpdIdx[cellI] < 0 || _pcc.cellToUpdIdx[cellJ] < 0) continue;

        PetscInt rowIdx = static_cast<PetscInt>(_pcc.cellToUpdIdx[cellI]);
        PetscInt colIdx = static_cast<PetscInt>(_pcc.cellToUpdIdx[cellJ]);
        CF_CHKERRCONTINUE(MatSetValuesBlocked(_pcc.p0Mat,
          1, &rowIdx, 1, &colIdx, &(it->second[0]), INSERT_VALUES));
      }
    }

    // Assemble P0 sparse matrix for FaceCoupled mode
    if (_pcc.useP0FaceCoupled)
    {
      CF_CHKERRCONTINUE(MatAssemblyBegin(_pcc.p0Mat, MAT_FINAL_ASSEMBLY));
      CF_CHKERRCONTINUE(MatAssemblyEnd(_pcc.p0Mat, MAT_FINAL_ASSEMBLY));
      CF_CHKERRCONTINUE(KSPSetOperators(_pcc.p0Ksp, _pcc.p0Mat, _pcc.p0Mat));
    }
  }
  else
  {
    // ======== PETSc MATRIX MODE: assemble into preconditioner matrix ========

    // Allocate P0 off-diagonal map for capturing face cross-blocks
    std::map<std::pair<CFuint,CFuint>, RealMatrix> p0OffDiagMap;
    if (_pcc.useP0FaceCoupled)
    {
      frData->setP0OffDiagBlocks(&p0OffDiagMap);
    }

    // Enable Jacobian assembly into preconditioner matrix
    spaceData->setComputeJacobianFlag(true);
    spaceData->setFillPreconditionerMatrix(true);
    spaceData->setLinearResidualMode(true);

    // Zero the preconditioner matrix
    PetscMatrix& precondMat = _pcc.pJFC->petscData->getPreconditionerMatrix();
    precondMat.resetToZeroEntries();

    // Reset updateCoeff and compute residuals
    updateCoeff = 0.0;
    spaceMtd->computeSpaceResidual(1.0);
    spaceMtd->computeTimeResidual(1.0);

    // Finalize preconditioner matrix assembly
    precondMat.finalAssembly();

    // Restore flags
    spaceData->setComputeJacobianFlag(savedDoComputeJacob);
    spaceData->setFillPreconditionerMatrix(savedFillPrecondMat);
    spaceData->setLinearResidualMode(false);

    // Unregister off-diagonal map
    frData->setP0OffDiagBlocks(CFNULL);

    // Zero the P0 sparse matrix for FaceCoupled mode
    if (_pcc.useP0FaceCoupled)
    {
      CF_CHKERRCONTINUE(MatZeroEntries(_pcc.p0Mat));
    }

    // Restore updateCoeff and rhs
    for (CFuint i = 0; i < nbStates; ++i)
      updateCoeff[i] = bkpUpdateCoeff[i];
    for (CFuint i = 0; i < rhsSize; ++i)
      rhs[i] = bkpRhs[i];

    // Extract and invert blocks from PETSc matrix
    const std::vector<CFint>& globalIDs = _pcc.pJFC->upStatesGlobalIDs;
    RealMatrix p0Block(nEqs, nEqs);
    CFuint currentInverterSize = 0;

    if (_pcc.useElementBlocks)
    {
      for (CFuint iCell = 0; iCell < nCells; ++iCell)
      {
        const std::vector<CFuint>& upStateIndices = _pcc.cellToUpStateIdx[iCell];
        const CFuint nSolPts = upStateIndices.size();
        const CFuint blockSize = nSolPts * nEqs;

        if (nSolPts == 0) continue;

        // Build index array for MatGetValues
        std::vector<CFint> matIndices(blockSize);
        for (CFuint iSol = 0; iSol < nSolPts; ++iSol)
        {
          const CFuint upIdx = upStateIndices[iSol];
          for (CFuint iEq = 0; iEq < nEqs; ++iEq)
            matIndices[iSol * nEqs + iEq] = globalIDs[upIdx * nEqs + iEq];
        }

        // Extract the full element diagonal block
        RealMatrix elemBlock(blockSize, blockSize);
        precondMat.getValues(blockSize, &matIndices[0],
                             blockSize, &matIndices[0],
                             &elemBlock[0]);

        // Compute P0 block via Galerkin projection
        p0Block = 0.0;
        for (CFuint i = 0; i < nSolPts; ++i)
          for (CFuint j = 0; j < nSolPts; ++j)
            for (CFuint e = 0; e < nEqs; ++e)
              for (CFuint f = 0; f < nEqs; ++f)
                p0Block(e, f) += elemBlock(i * nEqs + e, j * nEqs + f);

        const CFreal invNSol = 1.0 / (CFreal)nSolPts;
        for (CFuint e = 0; e < nEqs; ++e)
          for (CFuint f = 0; f < nEqs; ++f)
            p0Block(e, f) *= invNSol;

        if (blockSize != currentInverterSize)
        {
          _inverter.reset(MatrixInverter::create(blockSize, false));
          currentInverterSize = blockSize;
        }
        _inverter->invert(elemBlock, _pcc.invElemBlocks[iCell]);

        // P0 coarse level
        if (_pcc.useP0FaceCoupled)
        {
          PetscInt blockRow = static_cast<PetscInt>(iCell);
          CF_CHKERRCONTINUE(MatSetValuesBlocked(_pcc.p0Mat,
            1, &blockRow, 1, &blockRow, &p0Block[0], INSERT_VALUES));

          // Off-diagonal blocks are filled below from actual face cross-blocks
        }
        else
        {
          _p0Inverter->invert(p0Block, _pcc.p0InvBlocks[iCell]);
        }
      }
    }
    else
    {
      RealMatrix dofBlock(nEqs, nEqs);

      for (CFuint iCell = 0; iCell < nCells; ++iCell)
      {
        const std::vector<CFuint>& upStateIndices = _pcc.cellToUpStateIdx[iCell];
        const CFuint nSolPts = upStateIndices.size();

        if (nSolPts == 0) continue;

        p0Block = 0.0;

        std::vector<CFint> dofIndices(nEqs);
        for (CFuint iSol = 0; iSol < nSolPts; ++iSol)
        {
          const CFuint upIdx = upStateIndices[iSol];
          for (CFuint iEq = 0; iEq < nEqs; ++iEq)
            dofIndices[iEq] = globalIDs[upIdx * nEqs + iEq];

          precondMat.getValues(nEqs, &dofIndices[0],
                               nEqs, &dofIndices[0],
                               &dofBlock[0]);

          for (CFuint e = 0; e < nEqs; ++e)
            for (CFuint f = 0; f < nEqs; ++f)
              p0Block(e, f) += dofBlock(e, f);

          _inverter->invert(dofBlock, _pcc.smootherInvBlocks[upIdx]);
        }

        const CFreal invNSol = 1.0 / (CFreal)nSolPts;
        for (CFuint e = 0; e < nEqs; ++e)
          for (CFuint f = 0; f < nEqs; ++f)
            p0Block(e, f) *= invNSol;

        // P0 coarse level
        if (_pcc.useP0FaceCoupled)
        {
          PetscInt blockRow = static_cast<PetscInt>(iCell);
          CF_CHKERRCONTINUE(MatSetValuesBlocked(_pcc.p0Mat,
            1, &blockRow, 1, &blockRow, &p0Block[0], INSERT_VALUES));

          // Off-diagonal blocks are filled below from actual face cross-blocks
        }
        else
        {
          _p0Inverter->invert(p0Block, _pcc.p0InvBlocks[iCell]);
        }
      }
    }

    // Insert actual off-diagonal blocks from Galerkin-projected face cross-blocks
    if (_pcc.useP0FaceCoupled)
    {
      for (auto it = p0OffDiagMap.begin(); it != p0OffDiagMap.end(); ++it)
      {
        const CFuint cellI = it->first.first;
        const CFuint cellJ = it->first.second;

        if (_pcc.cellToUpdIdx[cellI] < 0 || _pcc.cellToUpdIdx[cellJ] < 0) continue;

        PetscInt rowIdx = static_cast<PetscInt>(_pcc.cellToUpdIdx[cellI]);
        PetscInt colIdx = static_cast<PetscInt>(_pcc.cellToUpdIdx[cellJ]);
        CF_CHKERRCONTINUE(MatSetValuesBlocked(_pcc.p0Mat,
          1, &rowIdx, 1, &colIdx, &(it->second[0]), INSERT_VALUES));
      }
    }

    // Assemble P0 sparse matrix for FaceCoupled mode
    if (_pcc.useP0FaceCoupled)
    {
      CF_CHKERRCONTINUE(MatAssemblyBegin(_pcc.p0Mat, MAT_FINAL_ASSEMBLY));
      CF_CHKERRCONTINUE(MatAssemblyEnd(_pcc.p0Mat, MAT_FINAL_ASSEMBLY));
      CF_CHKERRCONTINUE(KSPSetOperators(_pcc.p0Ksp, _pcc.p0Mat, _pcc.p0Mat));
    }
  }

  CFLog(VERBOSE, "FRP0Preconditioner::computeBeforeSolving() => DONE\n");
}

//////////////////////////////////////////////////////////////////////////////

void FRP0Preconditioner::computeAfterSolving()
{
  // Blocks are fully overwritten on the next computeBeforeSolving() call —
  // no need to zero them here.
}

//////////////////////////////////////////////////////////////////////////////

PetscErrorCode FRP0PcApply(PC pc, Vec X, Vec Y)
{
  // Two-level p-multigrid preconditioner (multiplicative or additive).
  //
  // Multiplicative (default):
  //   y = omega * B^{-1} * X           (smoother)
  //   d = X - A * y                     (defect via MatMFFD)
  //   Y = y + P * D0^{-1} * R * d      (coarse correction on defect)
  //
  // Additive:
  //   Y = omega * B^{-1} * X + P * D0^{-1} * R * X
  //
  // B^{-1}: inverted element-diagonal blocks (ElementBlock or PointBlock)
  // D0^{-1}: Galerkin-projected P0 operator (BlockDiag or FaceCoupled)
  // R: restriction (uniform averaging to cell mean)
  // P: prolongation (constant injection to all DOFs)

  PetscFunctionBegin;

  void* ctx;
  CF_CHKERRCONTINUE(PCShellGetContext(pc, &ctx));
  FRP0PcJFContext* pcContext = (FRP0PcJFContext*)(ctx);

  const CFuint nEqs = pcContext->nEqs;
  const CFuint nCells = pcContext->nUpdatableCells;
  const CFreal omega = pcContext->smootherOmega;
  const bool useElemBlocks = pcContext->useElementBlocks;

  // Get raw arrays
  const CFreal* x;
  CFreal* y;
  CF_CHKERRCONTINUE(VecGetArrayRead(X, &x));
  CF_CHKERRCONTINUE(VecGetArray(Y, &y));

  PetscInt ySize;
  CF_CHKERRCONTINUE(VecGetLocalSize(Y, &ySize));

  // Initialize Y to zero
  for (PetscInt i = 0; i < ySize; ++i)
    y[i] = 0.0;

  // Zero P0 buffers
  const CFuint p0Size = nCells * nEqs;
  for (CFuint i = 0; i < p0Size; ++i)
  {
    pcContext->p0Residual[i] = 0.0;
    pcContext->p0Correction[i] = 0.0;
  }

  // ---- SMOOTHER + RESTRICTION (combined cell loop) ----
  for (CFuint iCell = 0; iCell < nCells; ++iCell)
  {
    const std::vector<CFuint>& upIdx = pcContext->cellToUpStateIdx[iCell];
    const CFuint cellNSol = upIdx.size();
    if (cellNSol == 0) continue;

    if (useElemBlocks)
    {
      // ElementBlock smoother: y_block = omega * invElemBlock * x_block
      const CFuint blockSize = cellNSol * nEqs;
      const RealMatrix& invBlock = pcContext->invElemBlocks[iCell];

      for (CFuint iRow = 0; iRow < blockSize; ++iRow)
      {
        const CFuint iSol = iRow / nEqs;
        const CFuint iEq  = iRow % nEqs;
        const CFuint yIdx = upIdx[iSol] * nEqs + iEq;

        CFreal val = 0.0;
        for (CFuint jCol = 0; jCol < blockSize; ++jCol)
        {
          const CFuint jSol = jCol / nEqs;
          const CFuint jEq  = jCol % nEqs;
          val += invBlock(iRow, jCol) * x[upIdx[jSol] * nEqs + jEq];
        }
        y[yIdx] = omega * val;
      }
    }
    else
    {
      // PointBlock smoother: per-DOF nEqs×nEqs mat-vec
      for (CFuint k = 0; k < cellNSol; ++k)
      {
        const RealMatrix& invBlk = pcContext->smootherInvBlocks[upIdx[k]];
        const CFuint baseIdx = upIdx[k] * nEqs;
        for (CFuint e = 0; e < nEqs; ++e)
        {
          CFreal val = 0.0;
          for (CFuint f = 0; f < nEqs; ++f)
            val += invBlk(e, f) * x[baseIdx + f];
          y[baseIdx + e] = omega * val;
        }
      }
    }

    // Restriction: p0rhs[iCell] = (1/nSolPts) * sum_k src[upIdx[k]]
    // In additive mode, src = x (raw input).
    // In multiplicative mode, src = defect (computed below after this loop).
    if (!pcContext->useMultiplicative)
    {
      const CFreal weight = 1.0 / (CFreal)cellNSol;
      for (CFuint e = 0; e < nEqs; ++e)
      {
        CFreal sum = 0.0;
        for (CFuint k = 0; k < cellNSol; ++k)
          sum += x[upIdx[k] * nEqs + e];
        pcContext->p0Residual[iCell * nEqs + e] = weight * sum;
      }
    }
  }

  // ---- MULTIPLICATIVE DEFECT + RESTRICTION ----
  if (pcContext->useMultiplicative)
  {
    // Release arrays before PETSc Vec operations
    CF_CHKERRCONTINUE(VecRestoreArrayRead(X, &x));
    CF_CHKERRCONTINUE(VecRestoreArray(Y, &y));

    // Lazy-allocate defect Vec with same parallel layout as X
    if (!pcContext->defectVec)
    {
      CF_CHKERRCONTINUE(VecDuplicate(X, &pcContext->defectVec));
    }

    // defect = X - A*Y (what the smoother missed)
    Mat Amat;
    CF_CHKERRCONTINUE(PCGetOperators(pc, &Amat, CFNULL));
    CF_CHKERRCONTINUE(MatMult(Amat, Y, pcContext->defectVec));
    // defectVec = X - defectVec
    CF_CHKERRCONTINUE(VecAYPX(pcContext->defectVec, -1.0, X));

    // Re-acquire arrays for restriction
    const CFreal* defect;
    CF_CHKERRCONTINUE(VecGetArrayRead(pcContext->defectVec, &defect));
    CF_CHKERRCONTINUE(VecGetArrayRead(X, &x));
    CF_CHKERRCONTINUE(VecGetArray(Y, &y));

    // Restrict the defect to P0
    for (CFuint iCell = 0; iCell < nCells; ++iCell)
    {
      const std::vector<CFuint>& upIdx = pcContext->cellToUpStateIdx[iCell];
      const CFuint cellNSol = upIdx.size();
      if (cellNSol == 0) continue;

      const CFreal weight = 1.0 / (CFreal)cellNSol;
      for (CFuint e = 0; e < nEqs; ++e)
      {
        CFreal sum = 0.0;
        for (CFuint k = 0; k < cellNSol; ++k)
          sum += defect[upIdx[k] * nEqs + e];
        pcContext->p0Residual[iCell * nEqs + e] = weight * sum;
      }
    }

    CF_CHKERRCONTINUE(VecRestoreArrayRead(pcContext->defectVec, &defect));
  }

  // ---- P0 COARSE SOLVE ----
  if (pcContext->useP0FaceCoupled)
  {
    // Direct LU solve: p0Mat * p0Sol = p0Rhs

    // Copy restricted residual into PETSc Vec
    CFreal* p0rhs;
    CF_CHKERRCONTINUE(VecGetArray(pcContext->p0RhsVec, &p0rhs));
    for (CFuint i = 0; i < p0Size; ++i)
      p0rhs[i] = pcContext->p0Residual[i];
    CF_CHKERRCONTINUE(VecRestoreArray(pcContext->p0RhsVec, &p0rhs));

    // Solve
    CF_CHKERRCONTINUE(KSPSolve(pcContext->p0Ksp, pcContext->p0RhsVec,
                                pcContext->p0SolVec));

    // Copy solution back to p0Correction buffer
    const CFreal* p0sol;
    CF_CHKERRCONTINUE(VecGetArrayRead(pcContext->p0SolVec, &p0sol));
    for (CFuint i = 0; i < p0Size; ++i)
      pcContext->p0Correction[i] = p0sol[i];
    CF_CHKERRCONTINUE(VecRestoreArrayRead(pcContext->p0SolVec, &p0sol));
  }
  else
  {
    // Block-diagonal P0 solve: p0corr = p0InvBlock * p0rhs (per-cell nEqs×nEqs)
    for (CFuint iCell = 0; iCell < nCells; ++iCell)
    {
      const RealMatrix& p0Inv = pcContext->p0InvBlocks[iCell];
      for (CFuint e = 0; e < nEqs; ++e)
      {
        CFreal val = 0.0;
        for (CFuint f = 0; f < nEqs; ++f)
          val += p0Inv(e, f) * pcContext->p0Residual[iCell * nEqs + f];
        pcContext->p0Correction[iCell * nEqs + e] = val;
      }
    }
  }

  // ---- PROLONGATION: Y[upIdx[k]] += p0corr (constant injection to all DOFs) ----
  for (CFuint iCell = 0; iCell < nCells; ++iCell)
  {
    const std::vector<CFuint>& upIdx = pcContext->cellToUpStateIdx[iCell];
    const CFuint cellNSol = upIdx.size();
    if (cellNSol == 0) continue;

    for (CFuint k = 0; k < cellNSol; ++k)
    {
      for (CFuint e = 0; e < nEqs; ++e)
      {
        y[upIdx[k] * nEqs + e] += pcContext->p0Correction[iCell * nEqs + e];
      }
    }
  }

  CF_CHKERRCONTINUE(VecRestoreArrayRead(X, &x));
  CF_CHKERRCONTINUE(VecRestoreArray(Y, &y));

  PetscFunctionReturn(0);
}

//////////////////////////////////////////////////////////////////////////////

vector<SafePtr<BaseDataSocketSink> > FRP0Preconditioner::needsSockets()
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
