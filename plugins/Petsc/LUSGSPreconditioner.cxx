// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Petsc/LUSGSPreconditioner.hh"
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
      
      // v3.0
      // extern PetscErrorCode LUSGSPcApply(void *ctx, Vec X, Vec Y);
      // v3.2
      extern PetscErrorCode LUSGSPcApply(PC pc, Vec X, Vec Y);
      
//////////////////////////////////////////////////////////////////////////////

MethodStrategyProvider<LUSGSPreconditioner,
                       PetscLSSData,
                       ShellPreconditioner,
                       PetscModule>
LUSGSPreconditionerProvider("LUSGS");

//////////////////////////////////////////////////////////////////////////////

void LUSGSPreconditioner::defineConfigOptions(Config::OptionList& options)
{
	options.addConfigOption< CFreal >("omega","Relaxation constant (0 < omega < 1 - underrelaxation; 1 < omega < 2 - overrelaxation)");
}

//////////////////////////////////////////////////////////////////////////////

LUSGSPreconditioner::LUSGSPreconditioner(const std::string& name) :
  ShellPreconditioner(name),
  socket_diagMatrices("diagMatrices"),
  socket_upLocalIDsAll("upLocalIDsAll"),
  _pcc(),
  _omega(),
  _inverter(CFNULL)
{
  addConfigOptionsTo(this);

  _omega = 1.0;
  setParameter("omega", &_omega);
}

//////////////////////////////////////////////////////////////////////////////

LUSGSPreconditioner::~LUSGSPreconditioner()
{
}

//////////////////////////////////////////////////////////////////////////////

void LUSGSPreconditioner::setPreconditioner()
{
  // getting the JFContext pointer into LUSGSPcJFContext pcc
  _pcc.pJFC = getMethodData().getJFContext();
  _pcc.diagMatrices = &socket_diagMatrices;
  _pcc.upLocalIDsAll = &socket_upLocalIDsAll;
  _pcc.omega = _omega; //cout << "\n\n\n\n\nOMEGA = " << _omega << "\n\n\n" << endl;

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
  CF_CHKERRCONTINUE(PCShellSetApply(_pcc.pJFC->petscData->getPreconditioner(), LUSGSPcApply));
}

//////////////////////////////////////////////////////////////////////////////

void LUSGSPreconditioner::computeBeforeSolving()
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

void LUSGSPreconditioner::computeAfterSolving() 
{
  // reset to 0 all the matrices
  DataHandle<CFreal> diagMatrices = socket_diagMatrices.getDataHandle(); 
  for (CFuint i = 0 ; i < diagMatrices.size(); ++i) {
    diagMatrices[i] = 0.0;
  }
}

//////////////////////////////////////////////////////////////////////////////

// v3.0
// PetscErrorCode LUSGSPcApply(void *ctx, Vec X, Vec Y)
// v3.2
PetscErrorCode LUSGSPcApply(PC pc, Vec X, Vec Y)
{
  // X is input vector - vector to be preconditioned
  // Y is output vector - preconditioned vector X
  // input vector X is the RHS in our formulation
  // ouput vector Y is the SOLUTION vector
  CFreal* x;
  CFreal* y;
  PetscFunctionBegin;

  // getting the LU-SGS preconditioner context
  
  // v3.2
  void* ctx;  CF_CHKERRCONTINUE(PCShellGetContext(pc,&ctx));
  
  LUSGSPcJFContext* pcContext = (LUSGSPcJFContext*)(ctx);
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
  // relaxation factor
  const CFreal relax = pcContext->omega;

  const CFreal invEps = -(relax/eps);
  RealVector& bkpStates = pcContext->pJFC->bkpStates;

  ConnectivityTable<CFuint>& stateNeighbors = pcContext->stateNeighbors;
  RealVector& sumDeltaR = pData->result;

  RealVector tmpY(nbEqs, &y[0]);
  RealMatrix invMatIter(nbEqs, nbEqs, &diagMatInv[0]);

  // NumericalJacobian& numJacob = spaceMethod->getSpaceMethodData()->getNumericalJacobian();
  // const RealVector& refValues = PhysicalModelStack::getActive()->getImplementor()->getRefStateValues();
  //  RealVector eps(nbEqs);
  //   RealVector invEps(nbEqs);
  //   for (CFuint i = 0; i < nbEqs; ++i) {
  //     eps[i] = numJacob.computeEps(i, refValues[i]); 
  //     invEps[i] = 1./eps[i];
  //   }

  pData->useAllStateIDs = false;
  // first step of LU-SGS preconditioning - loop over all states
  pData->useBiggerStateIDs = false;

  for(CFuint i = 0; i < nbStates; ++i)
  {
    if (states[i]->isParUpdatable()) {
      const CFuint stateID = states[i]->getLocalID();

      pData->currentStateID = stateID; // state ID corresponding to the selected i state

      // at first - compute Jacobians*vector product of cells situated on the left side form diag
      // i.e. all delta_R[j] j < i - but must be neighbours of cell[i]

      // Jacobian free approach delta_R[j] = (RHS(U[j] + eps*y[j]) - RHS(U[j]))/eps;  - delta_R[j] is an block vector
      // compute sum of sumDeltaR[i] = delta_R[j];
      sumDeltaR = 0.0; // putting zeroes into sumDeltaR array

      //cout <<endl << "FORWARD(j<i)  (U_" << i << ") => ";

      spaceMethod->computeSpaceRhsForStatesSet(1.0);

      sumDeltaR *= -1.0;

      // there may be a problem with RHS (delta_R[j]) sign!!! - like before in Jacobian free approach
      const CFuint nbNeighbors = stateNeighbors.nbCols(i);

      //cout << endl << "neighbors(" << nbNeighbors << ") => ";

      for(CFuint j = 0; j < nbNeighbors; ++j) {
        const CFuint stateNeighborID = stateNeighbors(i,j);

        //cout << stateNeighborID << " ";

        if ((states[stateNeighborID]->isParUpdatable())) {
          if (stateNeighborID < stateID) {
            const CFint updatableStateID = upLocalIDsAll[stateNeighborID];
            cf_assert(updatableStateID >= 0);

            //cout << "neighborID = " << stateNeighborID << ", upStateID = " << updatableStateID << endl;
            const CFuint jTimesEqUp = updatableStateID*nbEqs;
            for(CFuint iEq = 0; iEq < nbEqs; ++iEq) {
              // states = U[j] + eps*y[j], U[j] = bkpStates
              (*states[stateNeighborID])[iEq] += eps*y[jTimesEqUp + iEq]; 
            }
          }
        }
      }

      //cout <<endl << "FORWARD(j<i)  (U_" << i << "+ dU) => ";
      spaceMethod->computeSpaceRhsForStatesSet(1.0);
      //cout <<endl;

      // restore the original states
      // there may be a problem with RHS (delta_R[j]) sign!!! - like before in Jacobian free approach
      for(CFuint j = 0; j < nbNeighbors; ++j) {
        const CFuint stateNeighborID = stateNeighbors(i,j);
        if ((states[stateNeighborID]->isParUpdatable())) {
          if (stateNeighborID < stateID) {
            // const CFint updatableStateID = upLocalIDsAll[stateNeighborID];
            cf_assert(upLocalIDsAll[stateNeighborID] >= 0);

            const CFuint jTimesEq   = stateNeighborID*nbEqs;
            for(CFuint iEq = 0; iEq < nbEqs; ++iEq) {
              (*states[stateNeighborID])[iEq] = bkpStates[jTimesEq + iEq];
            }
            //cout << "states[" << stateNeighborID <<"] = " << (*states[stateNeighborID]) << endl;
          }
        }
      }

      sumDeltaR *= invEps;

      const CFuint startX = upLocalIDsAll[stateID]*nbEqs;
      for(CFuint iEq = 0; iEq < nbEqs; ++iEq) {
        sumDeltaR[iEq] += x[startX + iEq];
      }

      invMatIter.wrap(nbEqs, nbEqs, &diagMatInv[i*nbEqs2]);
      tmpY.wrap(nbEqs, &y[startX]);
      tmpY = invMatIter*sumDeltaR;
    }
  }
  // second step of LU-SGS preconditioning
  pData->useBiggerStateIDs = true;

  // here you have to use a CFint because otherwise the condition >=0 is ALWAYS satisfied
  for(CFint i = nbStates - 1; i >= 0; --i) {
    if (states[i]->isParUpdatable()) {
      const CFuint stateID = states[i]->getLocalID();
      pData->currentStateID = stateID; // state ID corresponding to the selected i state

      sumDeltaR = 0.0; // putting zeroes into sumDeltaR array

      //cout <<endl << "BACKWARD(j>i)  (U_" << i << ") => ";
      spaceMethod->computeSpaceRhsForStatesSet(1.0);
      sumDeltaR *= -1.0;

      // there may be a problem with RHS (delta_R[j]) sign!!! - like before in Jacobian free approach
      const CFuint nbNeighbors = stateNeighbors.nbCols(i);
      //cout << endl << "neighbors(" << nbNeighbors << ") => ";

      for(CFuint j = 0; j < nbNeighbors; ++j) {
        const CFuint stateNeighborID = stateNeighbors(i,j);

        //cout << stateNeighborID << " ";

        if ((states[stateNeighborID]->isParUpdatable())) {
          if (stateNeighborID > stateID) {
            const CFint updatableStateID = upLocalIDsAll[stateNeighborID];
            cf_assert(updatableStateID >= 0);

            // we have to say to space method that the vector of unknowns is y not deltaU
            // states = U + eps*y, U = bkpStates
            const CFuint jTimesEqUp = updatableStateID*nbEqs;

            for(CFuint iEq = 0; iEq < nbEqs; ++iEq) {
              // states = U[j] + eps*y[j], U[j] = bkpState   
              (*states[stateNeighborID])[iEq] += eps*y[jTimesEqUp + iEq];
            }
            //cout << "states[" << stateNeighborID <<"] = " << (*states[stateNeighborID]) << endl;
          }
        }
      }

      //cout <<endl << "BACKWARD(j>i)  (U_" << i << "+ dU) => ";
      spaceMethod->computeSpaceRhsForStatesSet(1.0);
      //cout <<endl;

      // restore the original states
      // there may be a problem with RHS (delta_R[j]) sign!!! - like before in Jacobian free approach
      for(CFuint j = 0; j < nbNeighbors; ++j) {
        const CFuint stateNeighborID = stateNeighbors(i,j);
        if ((states[stateNeighborID]->isParUpdatable())) {
          if (stateNeighborID > stateID) {
            //const CFint updatableStateID = upLocalIDsAll[stateNeighborID];
            cf_assert(upLocalIDsAll[stateNeighborID] >= 0);

            // we have to say to space method that the vector of unknowns is y not deltaU
            // states = U + eps*y, U = bkpStates
            const CFuint jTimesEq   = stateNeighborID*nbEqs;
            for(CFuint iEq = 0; iEq < nbEqs; ++iEq) {
              (*states[stateNeighborID])[iEq] = bkpStates[jTimesEq + iEq];
            }
            //cout << "states[" << stateNeighborID <<"] = " << (*states[stateNeighborID]) << endl;
          }
        }
      }

      sumDeltaR *= invEps;

      invMatIter.wrap(nbEqs, nbEqs, &diagMatInv[i*nbEqs2]);
      const CFuint startX = upLocalIDsAll[stateID]*nbEqs;
      tmpY.wrap(nbEqs, &y[startX]);
      tmpY += invMatIter*sumDeltaR;
    }
  }

  // restoring of arrays X - vector to be preconditioned and Y - preconditioned vector
  CF_CHKERRCONTINUE(VecRestoreArray(X, &x));
  CF_CHKERRCONTINUE(VecRestoreArray(Y, &y));

  PetscFunctionReturn(0);
}

//////////////////////////////////////////////////////////////////////////////

vector<SafePtr<BaseDataSocketSink> > LUSGSPreconditioner::needsSockets()
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

