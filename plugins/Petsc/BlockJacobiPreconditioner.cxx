// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Petsc/BlockJacobiPreconditioner.hh"
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

    extern PetscErrorCode BlockJacobiPcApply(PC pc, Vec X, Vec Y);
    
//////////////////////////////////////////////////////////////////////////////

MethodStrategyProvider<BlockJacobiPreconditioner,
                       PetscLSSData,
                       ShellPreconditioner,
                       PetscModule>
BlockJacobiPreconditionerProvider("BJacobi");

//////////////////////////////////////////////////////////////////////////////

BlockJacobiPreconditioner::BlockJacobiPreconditioner(const std::string& name) :
  ShellPreconditioner(name),
  socket_diagMatrices("diagMatrices"),
  socket_upLocalIDsAll("upLocalIDsAll"),
  _pcc(),
  _inverter(CFNULL)
{
}

//////////////////////////////////////////////////////////////////////////////

BlockJacobiPreconditioner::~BlockJacobiPreconditioner()
{
}

//////////////////////////////////////////////////////////////////////////////

void BlockJacobiPreconditioner::setPreconditioner()
{
  // getting the JFContext pointer into BlockJacobiPcJFContext pcc
  _pcc.pJFC = getMethodData().getJFContext();
  _pcc.diagMatrices = &socket_diagMatrices;
  _pcc.upLocalIDsAll = &socket_upLocalIDsAll;

  _inverter.reset(MatrixInverter::create(getMethodData().getNbSysEquations(), false));

  DataHandle<State*, GLOBAL> states = _pcc.pJFC->states->getDataHandle();

  CF_CHKERRCONTINUE(PCShellSetContext(_pcc.pJFC->petscData->getPreconditioner(), &_pcc));
  CF_CHKERRCONTINUE(PCShellSetApply(_pcc.pJFC->petscData->getPreconditioner(), BlockJacobiPcApply));
}

//////////////////////////////////////////////////////////////////////////////

void BlockJacobiPreconditioner::computeBeforeSolving() 
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

void BlockJacobiPreconditioner::computeAfterSolving() 
{
  // reset to 0 all the matrices
  DataHandle<CFreal> diagMatrices = socket_diagMatrices.getDataHandle(); 
  for (CFuint i =0 ; i < diagMatrices.size(); ++i) {
    diagMatrices[i] = 0.0;
  }
}

//////////////////////////////////////////////////////////////////////////////

PetscErrorCode BlockJacobiPcApply(PC pc, Vec X, Vec Y)
{
  // X is input vector - vector to be preconditioned
  // Y is output vector - preconditioned vector X
  // input vector X is the RHS in our formulation 
  // ouput vector Y is the SOLUTION vector
  CFreal* x;
  CFreal* y;

  PetscFunctionBegin;
  
  void* ctx;  CF_CHKERRCONTINUE(PCShellGetContext(pc,&ctx));
  
  BlockJacobiPcJFContext* pcContext = (BlockJacobiPcJFContext*)(ctx);
  
  // getting data handle to diagonal inverted matrices
  DataHandle<CFreal> diagMatInv   = pcContext->diagMatrices->getDataHandle(); // inverted diagonal matrices
  DataHandle<CFint> upLocalIDsAll = pcContext->upLocalIDsAll->getDataHandle(); // updatable state IDs (-1 i not updatable)

  // putting X vector (vector to be preconditioned) into array "rhs" and Y output preconditioned vector into array "y"
  CF_CHKERRCONTINUE(VecGetArray(X, &x));
  CF_CHKERRCONTINUE(VecGetArray(Y, &y));

  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  const CFuint nbEqs2 = nbEqs*nbEqs;

  // getting nuber of equations
  DataHandle<State*, GLOBAL> states = pcContext->pJFC->states->getDataHandle();
  const CFint nbUpdatableStates = diagMatInv.size()/nbEqs2;

  RealVector tmpX(nbEqs, &x[0]);
  RealVector tmpY(nbEqs, &y[0]);
  RealMatrix invMatIter(nbEqs, nbEqs, &diagMatInv[0]);

  for(CFint i = 0; i < nbUpdatableStates; ++i)
  {
    const CFuint startIdx = i*nbEqs;
    tmpX.wrap(nbEqs,&x[startIdx]);
    tmpY.wrap(nbEqs,&y[startIdx]);
    invMatIter.wrap(nbEqs, nbEqs, &diagMatInv[i*nbEqs2]);
    
    tmpY = invMatIter*tmpX;
  }

  // restoring of arrays X - vector to be preconditioned and Y - preconditioned vector
  CF_CHKERRCONTINUE(VecRestoreArray(X, &x));
  CF_CHKERRCONTINUE(VecRestoreArray(Y, &y));

  PetscFunctionReturn(0);
}

//////////////////////////////////////////////////////////////////////////////

vector<SafePtr<BaseDataSocketSink> > BlockJacobiPreconditioner::needsSockets()
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
