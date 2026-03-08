// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Petsc/PetscHeaders.hh" // must come before any header

#include "FluxReconstructionPetsc/FRJacobAssembler.hh"
#include "Petsc/PetscLSSData.hh"
#include "FluxReconstructionPetsc/FluxReconstructionPetsc.hh"

#include "Framework/MeshData.hh"
#include "Framework/PhysicalModel.hh"
#include "Framework/SpaceMethod.hh"
#include "Framework/SpaceMethodData.hh"
#include "Framework/MethodStrategyProvider.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Petsc {

//////////////////////////////////////////////////////////////////////////////

using COOLFluiD::FluxReconstructionMethod::FluxReconstructionPetscModule;

MethodStrategyProvider<FRJacobAssembler,
                       PetscLSSData,
                       ShellPreconditioner,
                       FluxReconstructionPetscModule>
FRJacobAssemblerProvider("FRJacobAssembler");

//////////////////////////////////////////////////////////////////////////////

FRJacobAssembler::FRJacobAssembler(const std::string& name) :
  ShellPreconditioner(name),
  socket_updateCoeff("updateCoeff")
{
}

//////////////////////////////////////////////////////////////////////////////

FRJacobAssembler::~FRJacobAssembler()
{
}

//////////////////////////////////////////////////////////////////////////////

void FRJacobAssembler::setPreconditioner()
{
  // No-op: we do NOT call PCShellSetApply() here.
  // This preserves PETSc's built-in PC type (PCILU, PCASM, etc.)
  // set in BaseSetup::setKSP() from the user's PCType config option.
  CFLog(INFO, "FRJacobAssembler::setPreconditioner() => no-op "
    "(PETSc built-in PC will be used)\n");
}

//////////////////////////////////////////////////////////////////////////////

void FRJacobAssembler::computeBeforeSolving()
{
  CFLog(INFO, "FRJacobAssembler::computeBeforeSolving() => assembling Jacobian\n");

  // FRJacobAssembler relies on PETSc's native PCILU operating on precondMat.
  // Bug #4: Requires DifferentPreconditionerMatrix=true to create a real matrix.
  cf_assert(getMethodData().getDifferentPreconditionerMatrix() &&
    "FRJacobAssembler requires DifferentPreconditionerMatrix = true in CFcase");
  // Bug #5: Requires the ParBAIJ path (!useBlockPreconditionerMatrix) so that
  // KSPSetOperators passes precondMat (not shellMat) as the PC matrix.
  // With useBlockPreconditionerMatrix=true, ParJFSolveSys routes
  // KSPSetOperators(ksp, shellMat, shellMat), which makes PETSc try
  // PCILU on a MATSHELL => crash.
  cf_assert(!getMethodData().useBlockPreconditionerMatrix() &&
    "FRJacobAssembler requires UseBlockPreconditioner = false (needs ParBAIJ path for PETSc PCILU)");

  JFContext* jfc = getMethodData().getJFContext();
  SafePtr<SpaceMethod> spaceMtd = jfc->spaceMethod;
  SafePtr<SpaceMethodData> spaceData = spaceMtd->getSpaceMethodData();

  // ---- Step 1: Save flags ----
  const bool savedDoComputeJacob = spaceData->doComputeJacobian();
  const bool savedFillPrecondMat = spaceData->fillPreconditionerMatrix();

  // ---- Step 2: Backup updateCoeff and rhs ----
  DataHandle<CFreal> updateCoeff = jfc->updateCoeff->getDataHandle();
  DataHandle<CFreal> rhs = jfc->rhs->getDataHandle();
  DataHandle<State*, GLOBAL> states = jfc->states->getDataHandle();
  const CFuint nbStates = states.size();
  const CFuint nEqs = getMethodData().getNbSysEquations();
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

  // ---- Step 3: Enable Jacobian assembly into preconditioner matrix ----
  spaceData->setComputeJacobianFlag(true);
  spaceData->setFillPreconditionerMatrix(true);

  // ---- Step 4: Zero the preconditioner matrix ----
  PetscMatrix& precondMat = jfc->petscData->getPreconditionerMatrix();
  precondMat.resetToZeroEntries();

  // ---- Step 5: Reset updateCoeff and compute residuals ----
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

  CFLog(INFO, "FRJacobAssembler::computeBeforeSolving() => DONE\n");
}

//////////////////////////////////////////////////////////////////////////////

void FRJacobAssembler::computeAfterSolving()
{
  // No-op: the preconditioner matrix is managed by PETSc
}

//////////////////////////////////////////////////////////////////////////////

vector<SafePtr<BaseDataSocketSink> > FRJacobAssembler::needsSockets()
{
  vector<SafePtr<BaseDataSocketSink> > result = ShellPreconditioner::needsSockets();
  result.push_back(&socket_updateCoeff);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace Petsc

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
