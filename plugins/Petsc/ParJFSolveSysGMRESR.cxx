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
#include "Framework/LSSIdxMapping.hh"
#include "Framework/State.hh"
#include "Framework/SpaceMethod.hh"

#include "Petsc/Petsc.hh"
#include "Petsc/ParJFSolveSysGMRESR.hh"
#include "Common/PE.hh"
#include <fstream>
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

extern PetscErrorCode computeJFMatGMRESR(Mat petscMat, Vec x, Vec y);
extern PetscErrorCode GMRESR(KSP*, PetscMatrix*, Vec*, Vec*);

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<ParJFSolveSysGMRESR, PetscLSSData, PetscModule> 
	parJFSolveSysGMRESRProvider("ParJFSolveSysGMRESR");

MethodCommandProvider<ParJFSolveSysGMRESR, PetscLSSData, PetscModule> 
	seqJFSolveSysGMRESRProvider("SeqJFSolveSysGMRESR");

//////////////////////////////////////////////////////////////////////////////

ParJFSolveSysGMRESR::ParJFSolveSysGMRESR(const string& name) :
	StdParSolveSys(name)
{
}

//////////////////////////////////////////////////////////////////////////////

ParJFSolveSysGMRESR::~ParJFSolveSysGMRESR()
{
}

//////////////////////////////////////////////////////////////////////////////

void ParJFSolveSysGMRESR::execute()
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
	mat.setJFFunction((void (*)(void))computeJFMatGMRESR);
	
	
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

#if PETSC_VERSION_MINOR==6 || PETSC_VERSION_MINOR==7 || PETSC_VERSION_MINOR==9 || PETSC_VERSION_MINOR==11 || PETSC_VERSION_MINOR==12 || PETSC_VERSION_MINOR==15 || PETSC_VERSION_MINOR==18
	CFuint ierr = KSPSetOperators(ksp,mat.getMat(), mat.getMat());
#else
	CFuint ierr = KSPSetOperators(ksp,mat.getMat(), mat.getMat(),DIFFERENT_NONZERO_PATTERN);
#endif	
	ierr = KSPSetUp(ksp);
	CHKERRCONTINUE(ierr);
	
	void* ctx;// = &temp_pointer1;
	mat.getJFContext(&ctx);
	
	//ierr = KSPSolve(ksp, rhsVec.getVec(), solVec.getVec());
	//CHKERRCONTINUE(ierr);
	//cout << "\n\n\n\nSolve - GMRESR - start\n\n\n\n";
	//ierr = KSPSetInitialGuessNonzero(ksp, PETSC_FALSE);
	//CHKERRCONTINUE(ierr);
	
	Vec B;
	VecDuplicate(rhsVec.getVec(), &B);
	VecCopy(rhsVec.getVec(), B);
	VecAssemblyBegin(B);
	VecAssemblyEnd(B);
	Vec X;
	VecDuplicate(solVec.getVec(), &X);
	VecCopy(solVec.getVec(), X);
	VecAssemblyBegin(X);
	VecAssemblyEnd(X);
	
	ierr = GMRESR(&ksp, &mat, &B, &X);
	CHKERRCONTINUE(ierr);
	
	solVec.copyVec(X);

	// v3.0
	// VecDestroy(B);
	// VecDestroy(X);
	
	// v3.2
	VecDestroy(&B);
	VecDestroy(&X);
	
   // restore original states before updating
	for(CFuint i = 0; i < nbStates; ++i)
	{
		const CFuint iTimesEq = i*nbEqs;
		for(CFuint j = 0; j < nbEqs; ++j)
		{
	   // states = original unperturbed states 
			(*states[i])[j] = jfc->bkpStates[iTimesEq + j];
		}
	}
	
	CFint iter = 0;
	ierr = KSPGetIterationNumber(ksp, &iter);
	CHKERRCONTINUE(ierr);
	
	//CFLog(INFO, "KSP convergence reached at iteration: " << iter << "\n");
	solVec.copy(&rhs[0], &_upLocalIDs[0], vecSize);
}

//////////////////////////////////////////////////////////////////////////////

void ParJFSolveSysGMRESR::setup()
{
	CFAUTOTRACE;
	
	StdParSolveSys::setup();
	
	JFContext* jfc = getMethodData().getJFContext();
	jfc->upLocalIDs = _upLocalIDs;
	jfc->upStatesGlobalIDs = _upStatesGlobalIDs;
}

//////////////////////////////////////////////////////////////////////////////

PetscErrorCode computeJFMatGMRESR(Mat petscMat, Vec x, Vec y)
{
	void* ctx;
	
	CF_CHKERRCONTINUE(MatShellGetContext(petscMat, &ctx));
	
	JFContext* jfc = (JFContext*)(ctx);
	DataHandle<State*, GLOBAL> states = jfc->states->getDataHandle(); 
	
	const CFuint nbEqs = states[0]->size();
	const CFuint nbStates = states.size();
	
	DataHandle<CFreal> rhs = jfc->rhs->getDataHandle();
	SafePtr<PetscVector> rhsVec = jfc->rhsVec;
	
	CFreal* statesArray;
	CF_CHKERRCONTINUE(VecGetArray(x, &statesArray));
	
  // loop over states - doing backup of states vector
	RealVector& bkpStates = jfc->bkpStates;
	const CFreal eps = jfc->eps;
	CFuint idx = 0;
	
  // to remove -----------------------------------------
  //  string fileName = "fileIDs_" + to_str<CFuint>(PE::GetPE().GetRank());
  // ofstream fout(fileName.c_str());
  // ---------------------------------------------------
	
	CFreal* rhsArray;
	CF_CHKERRCONTINUE(VecGetArray(rhsVec->getVec(), &rhsArray));
	
	for(CFuint i = 0; i < nbStates; ++i) 
	{
		if (states[i]->isParUpdatable()) 
		{
			const CFuint idxTimesEq = idx*nbEqs;
			const CFuint iTimesEq = i*nbEqs;
			
			for(CFuint j = 0; j < nbEqs; ++j)
			{
// 	fout << "(i,j)= "<< i << " " << j << ", iTj= " 
// 	      << iTimesEq + j << ",idxTj = " << idxTimesEq + j 
// 	      <<  ", gID = " <<  jfc->upStatesGlobalIDs[idxTimesEq + j]<< "\n";
				
	// states = U + eps*delta_U, U = bkpStates, delta_U = statesArray
				(*states[i])[j] = bkpStates[iTimesEq + j] + eps*statesArray[idxTimesEq + j];
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
	jfc->spaceMethod->computeSpaceResidual(1.0);
	jfc->spaceMethod->computeTimeResidual(1.0);
	
  // final assignment into PetscVec y
	const CFreal invEps = 1.0/eps;
	const CFuint vecSize = jfc->upLocalIDs.size();
	for(CFuint i = 0; i < vecSize; ++i)
	{
		const CFreal Fv = (rhsArray[i] - rhs[jfc->upLocalIDs[i]])*invEps;
		VecSetValue(y, jfc->upStatesGlobalIDs[i], Fv, INSERT_VALUES);
	}
	
	CF_CHKERRCONTINUE(VecAssemblyBegin(y));
	CF_CHKERRCONTINUE(VecAssemblyEnd(y));
	
	CF_CHKERRCONTINUE(VecRestoreArray(x, &statesArray));
	CF_CHKERRCONTINUE(VecRestoreArray(rhsVec->getVec(), &rhsArray));
	
	PetscFunctionReturn(0);
}

//////////////////////////////////////////////////////////////////////////////

PetscErrorCode GMRESR(KSP* ksp, PetscMatrix* A, Vec* B, Vec* X)
{
	CFint ierr;
	//cout << "Vec B " << endl;
	//ierr = VecView(*B, PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
/// GMRESR - BEGIN
	Vec R;
	CFint gmresRMaxIter = 20;
	//Vec Z;
	Vec** C = new Vec*[gmresRMaxIter];
	Vec** U = new Vec*[gmresRMaxIter];
	
	PetscReal normR;
	PetscReal normC;
	PetscReal scalarProdCR;
	PetscReal alpha;
	
	PetscFunctionBegin;
	ierr = VecDuplicate(*X, &R); CHKERRQ(ierr);
	
	MatMult(A->getMat(), *X, R); // R = Ax
	VecAYPX(R, -1, *B); // R = B - Ax
	
	//cout << "GMRESR - begin" << endl;
/// GMRESR - algorithm - begin
	PetscInt k = -1;
	ierr = VecNorm(R, NORM_2, &normR); CHKERRQ(ierr);
	while((normR > 1.e-15 || k == -1) && k < (gmresRMaxIter-1))
	{
		k++;
	// creating new vector C[k], U[k]
		C[k] = new Vec;
		U[k] = new Vec;
		ierr = VecDuplicate(R, &C[k][0]); CHKERRQ(ierr);
		ierr = VecDuplicate(R, &U[k][0]); CHKERRQ(ierr);
		ierr = VecCopy(R, U[k][0]); CHKERRQ(ierr); // inital guess for U = R - temporary vector
		ierr = KSPSolve(*ksp, R, U[k][0]); CHKERRQ(ierr);
		
		MatMult(A->getMat(), U[k][0], C[k][0]); // C = AZ
		
		for(int i = 0; i < k; i++)
		{
			VecDot(C[i][0], C[k][0], &alpha); // alpha = C[i]*C[k] (scalar product)
			VecAXPY(C[k][0], -alpha, C[i][0]); // C[k] = C[k] - alpha*C[i]
			VecAXPY(U[k][0], -alpha, U[i][0]); // U[k] = U[k] - alpha*U[i]
		}
	// dividing by norm of C[k]
		ierr = VecNorm(C[k][0], NORM_2, &normC); CHKERRQ(ierr);
		VecScale(C[k][0], 1.0/normC);
		VecScale(U[k][0], 1.0/normC);
		
		VecDot(C[k][0], R, &scalarProdCR);
		VecAXPY(*X, scalarProdCR, U[k][0]); // X = X + scalarProdCR*U[k]
		VecAXPY(R, -scalarProdCR, C[k][0]); //  R = R - scalarProdCR*C[k]
		
		ierr = VecNorm(R, NORM_2, &normR); CHKERRQ(ierr);
	};
	cout << "GMRESR iter = " << k << " res = " << normR << endl;
/// GMRESR - algorithm - end
	// cleaning memory
	for(int i = 0; i <= k; i++)
	{
	  VecDestroy(&U[i][0]);
	  delete U[i];
	  VecDestroy(&C[i][0]);
	  delete C[i];
	}
	delete U;
	delete C;
	VecDestroy(&R);
	
	//KSPSolve(*ksp, *B, *X);
	//cout << "GMRESR iter = " << k << " residual = " << normR << endl;
	//ierr = VecView(*X, PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
	/// GMRESR - END
	
	//cout << "Vec X " << endl;
	//ierr = VecView(*X, PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
	//cout << "Vec R " << endl;
	//ierr = VecView(R, PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
	//cout << "Vec B " << endl;
	//ierr = VecView(*B, PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
	//ierr = KSPSolve(*ksp, *B, *X); CHKERRQ(ierr);
	PetscFunctionReturn(0);
}


} // namespace Petsc

} // namespace COOLFluiD

//  LocalWords:  rhsArrayVecstatesArray
