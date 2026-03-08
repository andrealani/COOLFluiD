// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Petsc/PetscHeaders.hh" // must come before any header

#include "Framework/MeshData.hh"
#include "Common/CFLog.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/SubSystemStatus.hh"
#include "Framework/State.hh"
#include "Framework/SpaceMethod.hh"

#include "Petsc/Petsc.hh"
#include "Petsc/ParJFSolveSys.hh"
#include "Petsc/LUSGSPreconditioner.hh"
#include "Petsc/DPLURPreconditioner.hh"
#include "Petsc/TridiagPreconditioner.hh"
#include "Common/PE.hh"
#include "Environment/CFEnvVars.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Environment;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace Petsc {

// Forward declaration of MFFD residual callback
PetscErrorCode computeResidualForMFFD(void* ctx, Vec U, Vec F);

extern PetscErrorCode LUSGSPcApply(void *ctx, Vec X, Vec Y);
extern PetscErrorCode DPLURPcApply(void *ctx, Vec X, Vec Y);
extern PetscErrorCode TridiagPcApply(void *ctx, Vec X, Vec Y);

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<ParJFSolveSys, PetscLSSData, PetscModule>
parJFSolveSysProvider("ParJFSolveSys");

MethodCommandProvider<ParJFSolveSys, PetscLSSData, PetscModule>
seqJFSolveSysProvider("SeqJFSolveSys");

//////////////////////////////////////////////////////////////////////////////

ParJFSolveSys::ParJFSolveSys(const string& name) :
  StdParSolveSys(name),
  socket_updateCoeff("updateCoeff")
{
}

//////////////////////////////////////////////////////////////////////////////

ParJFSolveSys::~ParJFSolveSys()
{
}

//////////////////////////////////////////////////////////////////////////////

void ParJFSolveSys::execute()
{
  CFAUTOTRACE;

  JFContext* jfc = getMethodData().getJFContext();

  DataHandle<State*, GLOBAL> states = jfc->states->getDataHandle();
  const CFuint nbEqs = states[0]->size();
  const CFuint nbStates = states.size();

  DataHandle<CFreal> rhs = jfc->rhs->getDataHandle();

  // loop over states - doing backup of state and rhs vectors
  for(CFuint i = 0; i < nbStates; ++i)
  {
    const CFuint iTimesEq = i*nbEqs;
    for(CFuint j = 0; j < nbEqs; ++j)
    {
      const CFuint iTimesj = iTimesEq + j;
      jfc->bkpStates[iTimesj] = (*states[i])[j];
    }
  }

  cf_assert(_upLocalIDs.size() == _upStatesGlobalIDs.size());
  const CFuint vecSize = _upLocalIDs.size();

  PetscMatrix& mat = getMethodData().getMatrix();
  PetscMatrix& precondMat = getMethodData().getPreconditionerMatrix();

  PetscVector& rhsVec = getMethodData().getRhsVector();
  PetscVector& solVec = getMethodData().getSolVector();

  KSP& ksp = getMethodData().getKSP();

  // the rhs is copied into the PetscVector for the rhs
  // _upStatesGlobalIDs[i] is different from _upLocalIDs[i]
  // in case multiple LSS are used
  for (CFuint i = 0; i < vecSize; ++i) {
    // cout << "nbSysEq = " << nbEqs << ", upLocalIDs = " << _upLocalIDs[i]
    // 	 << ", upStatesGlobalIDs = " << _upStatesGlobalIDs[i] << endl;
    rhsVec.setValue(_upStatesGlobalIDs[i], rhs[_upLocalIDs[i]]);
  }

  // assemble the rhs vector
  rhsVec.assembly();

  if(jfc->differentPreconditionerMatrix) {
    getMethodData().getShellPreconditioner()->computeBeforeSolving();

    if(!getMethodData().useBlockPreconditionerMatrix()) {
      //cout << "\n\n\n Setting up J-F different preconditioner matrix ParBAIJ with Petsc preconditioner \n\n\n";
      precondMat.finalAssembly();

#if PETSC_VERSION_MINOR >= 6
      CF_CHKERRCONTINUE(KSPSetOperators(ksp, mat.getMat(), precondMat.getMat()));
#else
      CF_CHKERRCONTINUE(KSPSetOperators(ksp, mat.getMat(), precondMat.getMat(), DIFFERENT_NONZERO_PATTERN));
#endif
    }
    else {
      //cout << "\n\n\n Setting up J-F different preconditioner matrix with Shell preconditioner \n\n\n";
      precondMat.finalAssembly();
#if PETSC_VERSION_MINOR >= 6
      CF_CHKERRCONTINUE(KSPSetOperators(ksp, mat.getMat(), mat.getMat()));
#else
      CF_CHKERRCONTINUE(KSPSetOperators(ksp, mat.getMat(), mat.getMat(), DIFFERENT_NONZERO_PATTERN));
#endif
    }
  }
  else {
    //cout << "\n\n\n Setting up J-F with Shell preconditioner or without preconditioner\n\n\n" ;
#if PETSC_VERSION_MINOR >= 6
    CF_CHKERRCONTINUE(KSPSetOperators(ksp, mat.getMat(), mat.getMat()));
#else
    CF_CHKERRCONTINUE(KSPSetOperators(ksp, mat.getMat(), mat.getMat(), DIFFERENT_NONZERO_PATTERN));
#endif
    getMethodData().getShellPreconditioner()->computeBeforeSolving();
  }

  CF_CHKERRCONTINUE(KSPSetUp(ksp));

  // --- MatMFFDSetBase: set the linearization point (U, F(U)) ---
  // Pack current states into stateVec for PETSc's MatMFFD
  // Option 1: pass F=NULL, PETSc evaluates F(U) via callback on first MatMult
  // (1 extra residual eval per Newton step — negligible vs 60+ KSP evals)
  for (CFuint i = 0; i < vecSize; ++i) {
    const CFint localID = _upLocalIDs[i];
    const CFuint stateIdx = localID / nbEqs;
    const CFuint eqIdx = localID % nbEqs;
    CF_CHKERRCONTINUE(VecSetValue(jfc->stateVec, _upStatesGlobalIDs[i],
      (*states[stateIdx])[eqIdx], INSERT_VALUES));
  }
  CF_CHKERRCONTINUE(VecAssemblyBegin(jfc->stateVec));
  CF_CHKERRCONTINUE(VecAssemblyEnd(jfc->stateVec));
  CF_CHKERRCONTINUE(MatMFFDSetBase(mat.getMat(), jfc->stateVec, PETSC_NULL));

  // --- Eisenstat-Walker adaptive KSP tolerance ---
  if (jfc->useEisenstatWalker) {
    // Reset EW state at the start of each new time step / outer iteration
    const CFuint currTimeStep = SubSystemStatusStack::getActive()->getNbIter();
    if (currTimeStep != jfc->ewLastTimeStep) {
      jfc->prevNonlinResNorm = -1.0;
      jfc->prevEta = 0.5;
      jfc->ewLastTimeStep = currTimeStep;
    }

    // Get current nonlinear residual norm from the rhs vector
    PetscReal currResNorm;
    CF_CHKERRCONTINUE(VecNorm(rhsVec.getVec(), NORM_2, &currResNorm));

    CFreal eta;
    if (jfc->prevNonlinResNorm > 0.0) {
      // Eisenstat-Walker Choice 2: eta = gamma * (||F_k|| / ||F_{k-1}||)^alpha
      eta = jfc->ewGamma * std::pow(currResNorm / jfc->prevNonlinResNorm, jfc->ewAlpha);

      // Safeguard: prevent eta from decreasing too rapidly
      // eta_k >= gamma * eta_{k-1}^alpha  (Eisenstat & Walker 1996, Section 2)
      const CFreal etaSafeguard = jfc->ewGamma * std::pow(jfc->prevEta, jfc->ewAlpha);
      eta = std::max(eta, etaSafeguard);

      // Clamp to [ewMinEta, ewMaxEta]
      eta = std::max(jfc->ewMinEta, std::min(eta, jfc->ewMaxEta));
    }
    else {
      eta = 0.5;  // first Newton iteration: moderate tolerance
    }

    CF_CHKERRCONTINUE(KSPSetTolerances(ksp, eta, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT));
    CFLog(INFO, "JFNK Eisenstat-Walker: eta = " << eta
      << ", ||F|| = " << currResNorm << "\n");

    jfc->prevNonlinResNorm = currResNorm;
    jfc->prevEta = eta;
  }

  CF_CHKERRCONTINUE(KSPSolve(ksp, rhsVec.getVec(), solVec.getVec()));
  CFint iter = 0;
  CF_CHKERRCONTINUE(KSPGetIterationNumber(ksp, &iter));

  PetscReal kspResNorm;
  CF_CHKERRCONTINUE(KSPGetResidualNorm(ksp, &kspResNorm));
  KSPConvergedReason reason;
  CF_CHKERRCONTINUE(KSPGetConvergedReason(ksp, &reason));
  CFLog(INFO, "KSP convergence reached at iteration: " << iter
    << ", residual norm: " << kspResNorm
    << ", reason: " << reason << "\n");

  getMethodData().getShellPreconditioner()->computeAfterSolving();

  solVec.copy(&rhs[0], &_upLocalIDs[0], vecSize);
}

//////////////////////////////////////////////////////////////////////////////

void ParJFSolveSys::setup()
{
  CFAUTOTRACE;

  StdParSolveSys::setup();

  JFContext* jfc = getMethodData().getJFContext();
  jfc->upLocalIDs = _upLocalIDs;
  jfc->upStatesGlobalIDs = _upStatesGlobalIDs;
  jfc->updateCoeff = &socket_updateCoeff;

  const CFuint nbStates = socket_states.getDataHandle().size();
  jfc->bkpUpdateCoeff.resize(nbStates);

  // In JF mode, Jacobian assembly must be disabled: the main matrix is a
  // MatMFFD that does not support MatSetValues.  The Jacob commands
  // (ConvRHSJacob, StdTimeRHSJacob, etc.) check doComputeJacobian() and
  // skip all matrix operations when the flag is false, while still computing
  // the RHS residual.  NewtonIterator resets the flag each sub-iteration
  // via setComputeJacobianFlag(getDoComputeJacobFlag()), so users of
  // NewtonIterator should set DoComputeJacobian = false in their CFcase.
  jfc->spaceMethod->setComputeJacobianFlag(false);

  // Register the MFFD residual evaluation callback.
  // PETSc's MatMFFD calls this function to evaluate F(U) at perturbed states;
  // it then computes the finite difference (F(U+hv) - F(U))/h internally.
  PetscMatrix& mat = getMethodData().getMatrix();
  CF_CHKERRCONTINUE(MatMFFDSetFunction(mat.getMat(),
    (PetscErrorCode (*)(void*, Vec, Vec))computeResidualForMFFD, (void*)(jfc)));

  CFLog(INFO, "ParJFSolveSys::setup() => MatMFFD callback registered\n");

  // FGMRES recommendation when using variable preconditioner
  if (jfc->differentPreconditionerMatrix) {
    KSPType currentType;
    CF_CHKERRCONTINUE(KSPGetType(getMethodData().getKSP(), &currentType));
    if (std::string(currentType) != std::string(KSPFGMRES)) {
      CFLog(WARN, "JFNK with DifferentPreconditionerMatrix: FGMRES is recommended "
        "for variable preconditioners (inner KSP). Current KSP type: " << currentType
        << ". Set KSPType = KSPFGMRES in CFcase if using ILU or nested KSP.\n");
    } else {
      // FGMRES is inherently right-preconditioned. Set explicitly for clarity.
      CF_CHKERRCONTINUE(KSPSetPCSide(getMethodData().getKSP(), PC_RIGHT));
      CFLog(INFO, "JFNK: FGMRES with right preconditioning (required for variable PC)\n");
    }
  }

  // set the preconditioner for later use
  //if(!jfc->differentPreconditionerMatrix) {
    getMethodData().getShellPreconditioner()->setPreconditioner();
  //}
}

//////////////////////////////////////////////////////////////////////////////

/// MatMFFD residual callback: evaluates F(U) = -R(U) at the given state U.
/// PETSc calls this function during MatMult to evaluate the operator at
/// perturbed states. PETSc handles the finite differencing internally:
///   A·v ≈ (F(U+hv) - F(U)) / h
/// where A = -dR/dU is the Jacobian operator.
///
/// Sign convention: COOLFluiD's rhs = R(U). We need dF/dU = A = -dR/dU,
/// so F(U) = -R(U) = -rhs.
PetscErrorCode computeResidualForMFFD(void* ctx, Vec U, Vec F)
{
  PetscFunctionBeginUser;

  JFContext* jfc = (JFContext*)(ctx);
  DataHandle<State*, GLOBAL> states = jfc->states->getDataHandle();
  DataHandle<CFreal> rhs = jfc->rhs->getDataHandle();
  DataHandle<CFreal> updateCoeff = jfc->updateCoeff->getDataHandle();

  const CFuint nbEqs = states[0]->size();
  const CFuint nbStates = states.size();

  // --- Step 1: Copy Vec U (perturbed state from PETSc) into COOLFluiD states ---
  const CFreal* uArray;
  CF_CHKERRCONTINUE(VecGetArrayRead(U, &uArray));

  CFuint idx = 0;
  for (CFuint i = 0; i < nbStates; ++i) {
    if (states[i]->isParUpdatable()) {
      const CFuint idxTimesEq = idx * nbEqs;
      State& currState = *states[i];
      for (CFuint j = 0; j < nbEqs; ++j) {
        currState[j] = uArray[idxTimesEq + j];
      }
      idx++;
    }
    // Backup updateCoeff for all states (including ghosts)
    jfc->bkpUpdateCoeff[i] = updateCoeff[i];
  }
  CF_CHKERRCONTINUE(VecRestoreArrayRead(U, &uArray));

  // --- Step 2: Synchronize ghost states ---
  if (CFEnv::getInstance().getVars()->SyncAlgo != "Old") {
    states.synchronize();
  }
  else {
    states.beginSync();
    states.endSync();
  }

  // --- Step 3: Evaluate residual R(U) ---
  jfc->spaceMethod->setComputeJacobianFlag(false);
  updateCoeff = 0.0;
  jfc->spaceMethod->computeSpaceResidual(1.0);
  jfc->spaceMethod->computeTimeResidual(1.0);

  // --- Step 4: Copy F = -R(U) = -rhs into output Vec F ---
  // Sign: F = -rhs so that dF/dU = -dR/dU = A (the Jacobian operator)
  const CFuint vecSize = jfc->upLocalIDs.size();
  for (CFuint i = 0; i < vecSize; ++i) {
    CF_CHKERRCONTINUE(VecSetValue(F, jfc->upStatesGlobalIDs[i],
      -rhs[jfc->upLocalIDs[i]], INSERT_VALUES));
  }
  CF_CHKERRCONTINUE(VecAssemblyBegin(F));
  CF_CHKERRCONTINUE(VecAssemblyEnd(F));

  // --- Step 5: Restore states and updateCoeff from backup ---
  // This is critical: the preconditioner may read COOLFluiD states between
  // MatMult calls, so states must be restored to the base point after each
  // MFFD function evaluation.
  const RealVector& bkpStates = jfc->bkpStates;
  for (CFuint i = 0; i < nbStates; ++i) {
    const CFuint iTimesEq = i * nbEqs;
    State& currState = *states[i];
    for (CFuint j = 0; j < nbEqs; ++j) {
      currState[j] = bkpStates[iTimesEq + j];
    }
    updateCoeff[i] = jfc->bkpUpdateCoeff[i];
  }

  PetscFunctionReturn(0);
}

//////////////////////////////////////////////////////////////////////////////

vector<SafePtr<BaseDataSocketSink> > ParJFSolveSys::needsSockets()
{
  vector<SafePtr<BaseDataSocketSink> > result = StdParSolveSys::needsSockets();
  result.push_back(&socket_updateCoeff);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace Petsc

} // namespace COOLFluiD
