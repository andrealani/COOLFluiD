// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Petsc/DPLURPreconditioner.hh"
#include "Petsc/PetscLSSData.hh"
#include "Petsc/Petsc.hh"

#include "Framework/SpaceMethod.hh"
#include "Framework/SpaceMethodData.hh"
#include "Framework/GlobalJacobianSparsity.hh"
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

    extern PetscErrorCode DPLURPcApply(PC pc, Vec X, Vec Y);

//////////////////////////////////////////////////////////////////////////////

MethodStrategyProvider<DPLURPreconditioner,
                       PetscLSSData,
                       ShellPreconditioner,
                       PetscModule>
DPLURPreconditionerProvider("DPLUR");

//////////////////////////////////////////////////////////////////////////////

void DPLURPreconditioner::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< CFreal >("omega","Relaxation constant (0 < omega < 1 - underrelaxation; 1 < omega < 2 - overrelaxation)");
  options.addConfigOption< CFuint >("nbSweeps", "Number of sweeps in the DP-LUR method");
}

//////////////////////////////////////////////////////////////////////////////

DPLURPreconditioner::DPLURPreconditioner(const std::string& name) :
  ShellPreconditioner(name),
  socket_diagMatrices("diagMatrices"),
  socket_upLocalIDsAll("upLocalIDsAll"),
  _pcc(),
  _inverter(CFNULL),
  _omega(),
  _nbSweeps()
{
  addConfigOptionsTo(this);

  _omega = 1.0;
  setParameter("omega", &_omega);

  _nbSweeps = 6;
  setParameter("nbSweeps", &_nbSweeps);
}

//////////////////////////////////////////////////////////////////////////////

DPLURPreconditioner::~DPLURPreconditioner()
{
}

//////////////////////////////////////////////////////////////////////////////

void DPLURPreconditioner::setPreconditioner()
{
  // getting the JFContext pointer into DPLURPcJFContext pcc
	_pcc.pJFC = getMethodData().getJFContext();
	_pcc.diagMatrices = &socket_diagMatrices;
	_pcc.upLocalIDsAll = &socket_upLocalIDsAll;
	_pcc.omega = _omega;
	_pcc.nbSweeps = _nbSweeps;
	
	SelfRegistPtr<GlobalJacobianSparsity> sparsity =
	getMethodData().getCollaborator<SpaceMethod>()->createJacobianSparsity();
	sparsity->computeMatrixPattern(*_pcc.pJFC->states, _pcc.stateNeighbors);
	
	_inverter.reset(MatrixInverter::create(getMethodData().getNbSysEquations(), false));
	
	// pointer to the SpaceMethod
	SafePtr<SpaceMethod> spaceMethod = _pcc.pJFC->spaceMethod;
	SafePtr<SpaceMethodData::PreconditionerData> pData = spaceMethod->getSpaceMethodData()->getPreconditionerData();
	pData->result.resize(getMethodData().getNbSysEquations());
	DataHandle<State*, GLOBAL> states = _pcc.pJFC->states->getDataHandle();
	
  // here set the shell preconditioner
  CF_CHKERRCONTINUE(PCShellSetContext(_pcc.pJFC->petscData->getPreconditioner(), &_pcc));
  CF_CHKERRCONTINUE(PCShellSetApply(_pcc.pJFC->petscData->getPreconditioner(), DPLURPcApply));
}

//////////////////////////////////////////////////////////////////////////////

void DPLURPreconditioner::computeBeforeSolving()
{
  DataHandle<CFreal> diagMatrices = socket_diagMatrices.getDataHandle();
  DataHandle<State*, GLOBAL> states = _pcc.pJFC->states->getDataHandle();

  const CFuint nbEqs = getMethodData().getNbSysEquations();
  const CFuint nbEqs2 = nbEqs*nbEqs;
  const CFuint nbUpdatableStates = diagMatrices.size()/nbEqs2;
  
  RealMatrix invMat(nbEqs, nbEqs);
  RealMatrix matIter(nbEqs, nbEqs, &diagMatrices[0]);
  
  for (CFuint i = 0; i < nbUpdatableStates; ++i) {
    // the matrix is used as an iterator here
    matIter.wrap(nbEqs, nbEqs, &diagMatrices[i*nbEqs2]);
    _inverter->invert(matIter,invMat);
    
    for (CFuint m = 0; m < nbEqs2; ++m) {
      matIter[m] = invMat[m];
    }
  }
}
    
//////////////////////////////////////////////////////////////////////////////

void DPLURPreconditioner::computeAfterSolving() 
{
  // reset to 0 all the matrices
  DataHandle<CFreal> diagMatrices = socket_diagMatrices.getDataHandle(); 
  for (CFuint i = 0 ; i < diagMatrices.size(); ++i) {
    diagMatrices[i] = 0.0;
  }
}
    
//////////////////////////////////////////////////////////////////////////////

PetscErrorCode DPLURPcApply(PC pc, Vec X, Vec Y)
{
	// X is input vector - vector to be preconditioned
	// Y is output vector - preconditioned vector X
	// input vector X is the RHS in our formulation
	// ouput vector Y is the SOLUTION vector
	CFreal* x;
	CFreal* y;
	
	PetscFunctionBegin;
	
	// getting the DP-LUR preconditioner context
	 // v3.2
	void* ctx;  CF_CHKERRCONTINUE(PCShellGetContext(pc,&ctx));
	
	DPLURPcJFContext* pcContext = (DPLURPcJFContext*)(ctx);
	// getting data handle to states - states are stored in JFContext
	DataHandle<State*, GLOBAL> states = pcContext->pJFC->states->getDataHandle();
	// getting data handle to diagonal inverted matrices
	DataHandle<CFreal> diagMatInv = pcContext->diagMatrices->getDataHandle(); // inverted diagonal matrices
	DataHandle<CFint> upLocalIDsAll = pcContext->upLocalIDsAll->getDataHandle(); // inverted diagonal matrices
	
	// putting X vector (vector to be preconditioned) into array "rhs" and Y output preconditioned vector into array "y"
	CF_CHKERRCONTINUE(VecGetArray(X, &x));
	CF_CHKERRCONTINUE(VecGetArray(Y, &y));
	
	// getting number of equations
	// pointer to the SpaceMethod
	SafePtr<SpaceMethod> spaceMethod = pcContext->pJFC->spaceMethod;
	SafePtr<SpaceMethodData::PreconditionerData> pData = spaceMethod->getSpaceMethodData()->getPreconditionerData();
	
	const CFuint nbStates = states.size();
	const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
	const CFuint nbEqs2 = nbEqs*nbEqs;
	const CFreal eps = pcContext->pJFC->eps;
	const CFreal invEps = -1.0/eps;
	RealVector& bkpStates = pcContext->pJFC->bkpStates;
	
	ConnectivityTable<CFuint>& stateNeighbors = pcContext->stateNeighbors;
	RealVector& sumDeltaR = pData->result;
	
	RealVector tmpY(nbEqs, &y[0]);
	RealMatrix invMatIter(nbEqs, nbEqs, &diagMatInv[0]);
	
	pData->useAllStateIDs = true; // all neighbors are needed in DP-LUR
	pData->useBiggerStateIDs = false;
	
// some testing stuff
	//const CFreal maxRelax = 0.5;
	//const CFreal minRelax = 0.5;
	//const CFreal dRelax = (maxRelax - minRelax)/(nbSweeps - 1.0);
	
	//cout << "omega = " << pcContext->omega << endl;
	//cout << "nbSweeps = " << pcContext->nbSweeps << endl;
	
	const CFuint nbSweeps = pcContext->nbSweeps; // number of sweeps in DP-LUR
	const CFreal relax = pcContext->omega; // under/over relaxation factor
	
	// sweeping
	for(CFuint sweep = 0; sweep < nbSweeps; ++sweep) {
		for(CFuint i = 0; i < nbStates; ++i)
		{
			if (states[i]->isParUpdatable()) {
				const CFuint stateID = states[i]->getLocalID();
				pData->currentStateID = stateID; // state ID corresponding to the selected i state
				
				// Jacobian free approach delta_R[j] = (RHS(U[j] + eps*y[j]) - RHS(U[j]))/eps;  - delta_R[j] is an block vector
				
				// compute sum of sumDeltaR[i] = delta_R[j];
				sumDeltaR = 0.0; // putting zeroes into sumDeltaR array
				
				spaceMethod->computeSpaceRhsForStatesSet(1.0);
				sumDeltaR *= -1.0;
				// in sumDeltaR is now RHS(U)
				
				const CFuint nbNeighbors = stateNeighbors.nbCols(i);
				
				for(CFuint j = 0; j < nbNeighbors; ++j) {
					const CFuint stateNeighborID = stateNeighbors(i,j);
					
					if ((states[stateNeighborID]->isParUpdatable())) {
						const CFint updatableStateID = upLocalIDsAll[stateNeighborID];
						cf_assert(updatableStateID >= 0);
						
						const CFuint jTimesEqUp = updatableStateID*nbEqs;
						for(CFuint iEq = 0; iEq < nbEqs; ++iEq) {
							// states = U[j] + eps*y[j], U[j] = bkpStates
							(*states[stateNeighborID])[iEq] += eps*y[jTimesEqUp + iEq]; 
						}
					}
				}
				
				spaceMethod->computeSpaceRhsForStatesSet(1.0);
				// in sumDeltaR is now -RHS(U[j] + eps*y[j]) + RHS(U[j])
				
				// restoring original states
				for(CFuint j = 0; j < nbNeighbors; ++j) {
					const CFuint stateNeighborID = stateNeighbors(i,j);
					if ((states[stateNeighborID]->isParUpdatable())) {
					  //	const CFint updatableStateID = upLocalIDsAll[stateNeighborID];
					  //cf_assert(updatableStateID >= 0);
						
						const CFuint jTimesEq   = stateNeighborID*nbEqs;
						for(CFuint iEq = 0; iEq < nbEqs; ++iEq) {
							(*states[stateNeighborID])[iEq] = bkpStates[jTimesEq + iEq];
						}
					}
				}
				
				sumDeltaR *= invEps;
				// in sumDeltaR is now (RHS(U[j] - eps*y[j]) + RHS(U[j]))/eps
				
				const CFuint startX = upLocalIDsAll[stateID]*nbEqs;
				for(CFuint iEq = 0; iEq < nbEqs; ++iEq) {
					sumDeltaR[iEq] += x[startX + iEq];
				}
				
				invMatIter.wrap(nbEqs, nbEqs, &diagMatInv[i*nbEqs2]);
				tmpY.wrap(nbEqs, &y[startX]);
				tmpY = (1.0 - relax)*tmpY + relax*(invMatIter*sumDeltaR);
			}
		}
	// some testing stuff
		//relax += dRelax;
	}
	
	// restoring of arrays X - vector to be preconditioned and Y - preconditioned vector
	CF_CHKERRCONTINUE(VecRestoreArray(X, &x));
	CF_CHKERRCONTINUE(VecRestoreArray(Y, &y));
	
	PetscFunctionReturn(0);
}

//////////////////////////////////////////////////////////////////////////////

vector<SafePtr<BaseDataSocketSink> > DPLURPreconditioner::needsSockets()
{
  vector<SafePtr<BaseDataSocketSink> > result = ShellPreconditioner::needsSockets();
  result.push_back(&socket_diagMatrices); 
  result.push_back(&socket_upLocalIDsAll);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace Petsc

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
