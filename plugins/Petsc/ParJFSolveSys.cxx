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

extern PetscErrorCode computeJFMat(Mat petscMat, Vec x, Vec y);
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

  mat.setJFFunction((void (*)(void))computeJFMat);

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
    if(!getMethodData().useBlockPreconditionerMatrix()) {
      //cout << "\n\n\n Setting up J-F different preconditioner matrix ParBAIJ with Petsc preconditioner \n\n\n";
      precondMat.finalAssembly();
      
#if PETSC_VERSION_MINOR==6 || PETSC_VERSION_MINOR==7 || PETSC_VERSION_MINOR==9 || PETSC_VERSION_MINOR==11 || PETSC_VERSION_MINOR==12 || PETSC_VERSION_MINOR==15 || PETSC_VERSION_MINOR==18
      CF_CHKERRCONTINUE(KSPSetOperators(ksp, mat.getMat(), precondMat.getMat()));
#else
      CF_CHKERRCONTINUE(KSPSetOperators(ksp, mat.getMat(), precondMat.getMat(), DIFFERENT_NONZERO_PATTERN));
#endif	   
    }
    else {
      //cout << "\n\n\n Setting up J-F different preconditioner matrix with Shell preconditioner \n\n\n";
      precondMat.finalAssembly();
#if PETSC_VERSION_MINOR==6 || PETSC_VERSION_MINOR==7 || PETSC_VERSION_MINOR==9 || PETSC_VERSION_MINOR==11 || PETSC_VERSION_MINOR==12 || PETSC_VERSION_MINOR==15 || PETSC_VERSION_MINOR==18
      CF_CHKERRCONTINUE(KSPSetOperators(ksp, mat.getMat(), mat.getMat()));
#else
      CF_CHKERRCONTINUE(KSPSetOperators(ksp, mat.getMat(), mat.getMat(), DIFFERENT_NONZERO_PATTERN));
#endif
      getMethodData().getShellPreconditioner()->computeBeforeSolving();
    }
  }
  else {
    //cout << "\n\n\n Setting up J-F with Shell preconditioner or without preconditioner\n\n\n" ;
#if PETSC_VERSION_MINOR==6 || PETSC_VERSION_MINOR==7 || PETSC_VERSION_MINOR==9 || PETSC_VERSION_MINOR==11 || PETSC_VERSION_MINOR==12 || PETSC_VERSION_MINOR==15 || PETSC_VERSION_MINOR==18
    CF_CHKERRCONTINUE(KSPSetOperators(ksp, mat.getMat(), mat.getMat()));
#else
    CF_CHKERRCONTINUE(KSPSetOperators(ksp, mat.getMat(), mat.getMat(), DIFFERENT_NONZERO_PATTERN));
#endif
    getMethodData().getShellPreconditioner()->computeBeforeSolving();
  }
  
  CF_CHKERRCONTINUE(KSPSetUp(ksp));

  void* ctx;
  mat.getJFContext(&ctx);
	
	// Testing stuff
	/*
	/// begin
	//ierr = KSPSetInitialGuessNonzero(ksp, PETSC_TRUE); CHKERRCONTINUE(ierr);
	// getting the LU-SGS preconditioner context
	void* pcContext;
	PC _pc;
	KSPGetPC(ksp, &_pc);
	PCShellGetContext(_pc, &pcContext);
	//cout << "after get LUSGSPc-context" << endl;
	/// DP-LUR
	//DPLURPcApply((void*)pcContext, rhsVec.getVec(), solVec.getVec());
	/// LUSGS
	//LUSGSPcApply((void*)pcContext, rhsVec.getVec(), solVec.getVec());
	//cout << "after LUSGSPc-apply" << endl;
	TridiagPcApply((void*)pcContext, rhsVec.getVec(), solVec.getVec());
	/// end
	*/
	
  CF_CHKERRCONTINUE(KSPSolve(ksp, rhsVec.getVec(), solVec.getVec()));
  CFint iter = 0;
  CF_CHKERRCONTINUE(KSPGetIterationNumber(ksp, &iter));

  CFLog(INFO, "KSP convergence reached at iteration: " << iter << "\n");

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

  // set the preconditioner for later use
  //if(!jfc->differentPreconditionerMatrix) {
    getMethodData().getShellPreconditioner()->setPreconditioner();
  //}
}

//////////////////////////////////////////////////////////////////////////////

PetscErrorCode computeJFMat(Mat petscMat, Vec x, Vec y)
{
  void* ctx;

  CF_CHKERRCONTINUE(MatShellGetContext(petscMat, &ctx));

  JFContext* jfc = (JFContext*)(ctx);
  DataHandle<State*, GLOBAL> states = jfc->states->getDataHandle(); 

  const CFuint nbEqs = states[0]->size();
  const CFuint nbStates = states.size();

  DataHandle<CFreal> rhs = jfc->rhs->getDataHandle(); 
  DataHandle<CFreal> updateCoeff = jfc->updateCoeff->getDataHandle();
  SafePtr<PetscVector> rhsVec = jfc->rhsVec;

  CFreal* statesArray;
  CF_CHKERRCONTINUE(VecGetArray(x, &statesArray));

  // loop over states - doing backup of states vector
  RealVector& bkpStates = jfc->bkpStates; 
  RealVector& bkpUpdateCoeff = jfc->bkpUpdateCoeff;

  const CFreal eps = jfc->eps;

  CFuint idx = 0;

  CFreal* rhsArray;
  CF_CHKERRCONTINUE(VecGetArray(rhsVec->getVec(), &rhsArray));

  for(CFuint i = 0; i < nbStates; ++i) 
  {
    if (states[i]->isParUpdatable()) 
    {
      const CFuint idxTimesEq = idx*nbEqs;
      const CFuint iTimesEq = i*nbEqs;
      State& currState = *states[i];
      for(CFuint j = 0; j < nbEqs; ++j) {
        // states = U + eps*delta_U, U = bkpStates, delta_U = statesArray
        currState[j] = bkpStates[iTimesEq + j] + eps*statesArray[idxTimesEq + j];
      }
      idx++;
    }

    // back up the update coefficient
    bkpUpdateCoeff[i] = updateCoeff[i];
  }

  // syncronize the states after the modifications
  if (CFEnv::getInstance().getVars()->SyncAlgo != "Old") {
    states.synchronize();
  } 
  else {
    states.beginSync();
    states.endSync();
  }
  
  // computation of F(U + eps*delta_U)
  jfc->spaceMethod->setComputeJacobianFlag(false);

  // reset to 0. the updateCoeff storage
  updateCoeff = 0.0;
  jfc->spaceMethod->computeSpaceResidual(1.0);
  jfc->spaceMethod->computeTimeResidual(1.0);

  if (!jfc->jfApprox2ndOrder) { // if 1st order of J-F approximation was chosen
    // final assignment into PetscVec y
    const CFuint vecSize = jfc->upLocalIDs.size();
    const CFreal invEps = 1.0/eps;
    for(CFuint i = 0; i < vecSize; ++i) {
      const CFreal Fv = (rhsArray[i] - rhs[jfc->upLocalIDs[i]])*invEps;
      VecSetValue(y, jfc->upStatesGlobalIDs[i], Fv, INSERT_VALUES);
    }
  }
  else { // if 2nd order of J-F approximation was chosen
    // assignment of RHS(U + eps*deltaU) into PetscVec y
    const CFuint vecSize = jfc->upLocalIDs.size();
    for(CFuint i = 0; i < vecSize; ++i) {
      VecSetValue(y, jfc->upStatesGlobalIDs[i], -rhs[jfc->upLocalIDs[i]], INSERT_VALUES);
    }

    CF_CHKERRCONTINUE(VecAssemblyBegin(y));
    CF_CHKERRCONTINUE(VecAssemblyEnd(y));

    idx = 0;
    for(CFuint i = 0; i < nbStates; ++i) {
      if (states[i]->isParUpdatable()) {
        const CFuint idxTimesEq = idx*nbEqs;
        const CFuint iTimesEq = i*nbEqs;
        State& currState = *states[i];
        for(CFuint j = 0; j < nbEqs; ++j) {
          // states = U - eps*delta_U, U = bkpStates, delta_U = statesArray
          currState[j] = bkpStates[iTimesEq + j] - eps*statesArray[idxTimesEq + j];
        }
        idx++;
      }
    }

    // syncronize the states after the modifications
    if (CFEnv::getInstance().getVars()->SyncAlgo != "Old") {
      states.synchronize();
    } 
    else {
      states.beginSync();
      states.endSync();
    }

    // computation of F(U + eps*delta_U)
    jfc->spaceMethod->setComputeJacobianFlag(false);

    // reset to 0. the updateCoeff storage
    updateCoeff = 0.0;
    jfc->spaceMethod->computeSpaceResidual(1.0);
    jfc->spaceMethod->computeTimeResidual(1.0);

    // assignment of RHS(U - eps*deltaU) to y vector
    for(CFuint i = 0; i < vecSize; ++i) {
      VecSetValue(y, jfc->upStatesGlobalIDs[i], rhs[jfc->upLocalIDs[i]], ADD_VALUES);
    }

    CF_CHKERRCONTINUE(VecAssemblyBegin(y));
    CF_CHKERRCONTINUE(VecAssemblyEnd(y));

    // multiplying with 1.0/2eps
    const CFreal inv2Eps = 0.5/eps;
    CF_CHKERRCONTINUE(VecScale (y, inv2Eps));
  }

  CF_CHKERRCONTINUE(VecAssemblyBegin(y));
  CF_CHKERRCONTINUE(VecAssemblyEnd(y));

  CF_CHKERRCONTINUE(VecRestoreArray(x, &statesArray));
  CF_CHKERRCONTINUE(VecRestoreArray(rhsVec->getVec(), &rhsArray));

  /// AL: the restoring of the states should be done here otherwise the states
  ///     the preconditioning is performed starting from wrong states
  for(CFuint i = 0; i < nbStates; ++i)
  {
    const CFuint iTimesEq = i*nbEqs;
    State& currState = *states[i];

    for(CFuint j = 0; j < nbEqs; ++j) {
      currState[j] = bkpStates[iTimesEq + j];
    }
    // restore the update coefficient
    updateCoeff[i] = bkpUpdateCoeff[i];
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

//  LocalWords:  rhsArrayVecstatesArray
