// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Petsc/TridiagPreconditioner.hh"
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

    extern PetscErrorCode TridiagPcApply(PC pc, Vec X, Vec Y);

//////////////////////////////////////////////////////////////////////////////

MethodStrategyProvider<TridiagPreconditioner,
                       PetscLSSData,
                       ShellPreconditioner,
                       PetscModule>
TridiagPreconditionerProvider("Tridiag");

//////////////////////////////////////////////////////////////////////////////

TridiagPreconditioner::TridiagPreconditioner(const std::string& name) :
  ShellPreconditioner(name),
  socket_diagMatrices("diagMatrices"),
  socket_underDiagMatrices("underDiagMatrices"),
  socket_aboveDiagMatrices("aboveDiagMatrices"),
  socket_upLocalIDsAll("upLocalIDsAll"),
  socket_dplrIsFirstInLine("dplrIsFirstInLine"),
  socket_dplrToLocalIDs("dplrToLocalIDs"),
  socket_localToDplrIDs("localToDplrIDs"),
  socket_dplrCellInLine("dplrCellInLine"),
  _pcc(),
  _inverter(CFNULL)
{
}

//////////////////////////////////////////////////////////////////////////////

TridiagPreconditioner::~TridiagPreconditioner()
{
}

//////////////////////////////////////////////////////////////////////////////

void TridiagPreconditioner::setPreconditioner()
{
// cout << "TridiagPreconditioner::setPreconditioner() - begin" << endl;

  // getting the JFContext pointer into TridiagPcJFContext pcc
  // cout << "getting the JFContext pointer into TridiagPcJFContext pcc ";
  //_pcc.pJFC = getMethodData().getJFContext();
  // cout << "... done" << endl;

  // cout << "getting the diagMatrices data socket ";
  _pcc.diagMatrices = &socket_diagMatrices;
  _pcc.underDiagMatrices = &socket_underDiagMatrices;
  _pcc.aboveDiagMatrices = &socket_aboveDiagMatrices;
  // cout << "... done" << endl;

  // cout << "getting the upLocalIDsAll data socket ";
  _pcc.upLocalIDsAll = &socket_upLocalIDsAll;
  // cout << "... done" << endl;

  _pcc.dplrIsFirstInLine = &socket_dplrIsFirstInLine;
  _pcc.dplrToLocalIDs = &socket_dplrToLocalIDs;
  _pcc.localToDplrIDs = &socket_localToDplrIDs;
  _pcc.dplrCellInLine = &socket_dplrCellInLine;

  // cout << "getting the inverter.reset ";
  _inverter.reset(MatrixInverter::create(getMethodData().getNbSysEquations(), false));
  // cout << "... done" << endl;

  // cout << "getting the States data handle ";
  // DataHandle<State*, GLOBAL> states = _pcc.pJFC->states->getDataHandle();
  // cout << "... done" << endl;

  //CF_CHKERRCONTINUE(PCShellSetContext(_pcc.pJFC->petscData->getPreconditioner(), &_pcc)); /// debug
  //CF_CHKERRCONTINUE(PCShellSetApply(_pcc.pJFC->petscData->getPreconditioner(), TridiagPcApply)); /// debug

  _pcc.petscData = &getMethodData(); /// debug

  CF_CHKERRCONTINUE(PCShellSetContext(_pcc.petscData->getPreconditioner(), &_pcc)); /// debug
  CF_CHKERRCONTINUE(PCShellSetApply(_pcc.petscData->getPreconditioner(), TridiagPcApply)); /// debug
}

//////////////////////////////////////////////////////////////////////////////

void TridiagPreconditioner::computeBeforeSolving()
{
  DataHandle<CFreal> diagMatrices = socket_diagMatrices.getDataHandle();
  DataHandle<CFreal> underDiagMatrices = socket_underDiagMatrices.getDataHandle();
  DataHandle<CFreal> aboveDiagMatrices = socket_aboveDiagMatrices.getDataHandle();

  //DataHandle<State*, GLOBAL> states = _pcc.pJFC->states->getDataHandle();

  const CFuint nbEqs = getMethodData().getNbSysEquations();
  const CFuint nbEqs2 = nbEqs*nbEqs;
  const CFuint nbUpdatableStates = diagMatrices.size()/nbEqs2;

  RealMatrix invMat(nbEqs, nbEqs);
  RealMatrix tmpMat(nbEqs, nbEqs);

  RealMatrix matIterDiag(nbEqs, nbEqs, &diagMatrices[0]);
  RealMatrix matIterUnderDiag(nbEqs, nbEqs, &underDiagMatrices[0]);
  RealMatrix matIterAboveDiag(nbEqs, nbEqs, &aboveDiagMatrices[0]);

  /*
  CFreal elem = 1.;
  matIterDiag = 0.;
  matIterUnderDiag = 0.;
  matIterAboveDiag = 0.;
  for(CFuint j = 0; j < nbEqs2; j += nbEqs+1) {
    matIterDiag[j] = elem; elem += 1.;
    matIterAboveDiag[j] = elem; elem += 1.;
  }
  for(CFint i = 1; i < nbUpdatableStates-1; i++) {
    matIterDiag.wrap(&diagMatrices[i*nbEqs2]);
    matIterUnderDiag.wrap(&underDiagMatrices[i*nbEqs2]);
    matIterAboveDiag.wrap(&aboveDiagMatrices[i*nbEqs2]);

    matIterDiag = 0.;
    matIterUnderDiag = 0.;
    matIterAboveDiag = 0.;

    for(CFuint j = 0; j < nbEqs2; j += nbEqs+1) {
      matIterUnderDiag[j] = elem; elem += 1.;
      matIterDiag[j] = elem; elem += 1.;
      matIterAboveDiag[j] = elem; elem += 1.;
    }
    // << "Diag = " << matIterDiag << endl;
    //cout << "Under-Diag = " << matIterUnderDiag << endl;
    //cout << "Above-Diag = " << matIterAboveDiag << endl;
  }
  matIterDiag.wrap(&diagMatrices[(nbUpdatableStates-1)*nbEqs2]);
  matIterUnderDiag.wrap(&underDiagMatrices[(nbUpdatableStates-1)*nbEqs2]);
  matIterAboveDiag.wrap(&aboveDiagMatrices[(nbUpdatableStates-1)*nbEqs2]);

  matIterDiag = 0.;
  matIterUnderDiag = 0.;
  matIterAboveDiag = 0.;

  for(CFuint j = 0; j < nbEqs2; j += nbEqs+1) {
    matIterUnderDiag[j] = elem; elem += 1.;
    matIterDiag[j] = elem; elem += 1.;
  }
  */

  for(CFuint i = 0; i < nbUpdatableStates; i++) {
    matIterDiag.wrap(nbEqs, nbEqs, &diagMatrices[i*nbEqs2]);
    matIterUnderDiag.wrap(nbEqs, nbEqs,&underDiagMatrices[i*nbEqs2]);
    matIterAboveDiag.wrap(nbEqs, nbEqs,&aboveDiagMatrices[i*nbEqs2]);

    // cout << "Diag_" << i << " = " << matIterDiag << endl;
    // cout << "Under-Diag_" << i << " = " << matIterUnderDiag << endl;
    // cout << "Above-Diag_" << i << " = " << matIterAboveDiag << endl;
  }


  // c - under-diag jacobians
  // a - diag jacobians
  // b - above diag jacobians
  matIterDiag.wrap(nbEqs, nbEqs,&diagMatrices[0]);
  matIterUnderDiag.wrap(nbEqs, nbEqs,&underDiagMatrices[0]);
  matIterAboveDiag.wrap(nbEqs, nbEqs,&aboveDiagMatrices[0]);
  // a_0 = a_0^{-1}
  _inverter->invert(matIterDiag, invMat);
  for (CFuint m = 0; m < nbEqs2; ++m) {
    matIterDiag[m] = invMat[m];
  }
  //cout << "nu0 = " << matIterDiag << endl;

  // b_0 = -a_0^{-1}*b_0 -> b_0 == mju_0
  tmpMat = matIterDiag*matIterAboveDiag;
  tmpMat *= -1.0;
  for (CFuint m = 0; m < nbEqs2; ++m) {
    matIterAboveDiag[m] = tmpMat[m];
  }
  //cout << "mju0 = " << matIterAboveDiag << endl;

  for (CFuint i = 1; i < nbUpdatableStates; ++i) {

    // the matrix is used as an iterator here
    // a_{i}
    matIterDiag.wrap(nbEqs, nbEqs,&diagMatrices[i*nbEqs2]);
    //cout << "diag" << i << " = " << matIterDiag << endl;

    // b_{i-1} == mju_{i-1}
    matIterAboveDiag.wrap(nbEqs, nbEqs,&aboveDiagMatrices[(i-1)*nbEqs2]);
    //cout << "mju" << i-1 << " = " << matIterAboveDiag << endl;

    // c_{i}
    matIterUnderDiag.wrap(nbEqs, nbEqs,&underDiagMatrices[i*nbEqs2]);
    //cout << "c" << i << " = " << matIterUnderDiag << endl;

    // tmpMat = c_{i}*b_{i-1} - because b_{i-1} = mju_{i-1}
    tmpMat = matIterUnderDiag*matIterAboveDiag;
    //cout << "c" << i << "*mju" << i-1 << " = " << tmpMat << endl;

    // tmpMat += a_{i}
    tmpMat += matIterDiag;
    //cout << "c" << i << "*mju" << i-1 << " + a" << i << " = " << tmpMat << endl;

    // inversion of tmpMat (tmpMat = c_{i}*mju_{i-1} + a_{i})
    _inverter->invert(tmpMat, invMat);
    // assigning of the new a_{i} jacobian-block a_{i} := (c_{i}*mju_{i-1} + a_{i})^{-1}
    for (CFuint m = 0; m < nbEqs2; ++m) {
      matIterDiag[m] = invMat[m];
    }
    //cout << "nu" << i << " = " << matIterDiag << endl;

    // b_{i}
    matIterAboveDiag.wrap(nbEqs, nbEqs,&aboveDiagMatrices[i*nbEqs2]);
    // b_{i} = -a_{i}*b_{i} -> b_{i} = mju_{i}
    tmpMat = matIterDiag*matIterAboveDiag;
    tmpMat *= -1.0;
    for (CFuint m = 0; m < nbEqs2; ++m) {
      matIterAboveDiag[m] = tmpMat[m];
    }
    //cout << "mju" << i << " = " << matIterAboveDiag << endl;
  }
  // now in diag = a_{i} := (c_{i}*mju_{i-1} + a_{i})^{-1}
  // now in above-diag = b_{i} := mju_{i}
  // now in under-diag = c_{i} := c_{i}
}

//////////////////////////////////////////////////////////////////////////////

void TridiagPreconditioner::computeAfterSolving()
{
  // reset to 0 all the matrices
  DataHandle<CFreal> diagMatrices = socket_diagMatrices.getDataHandle();
  DataHandle<CFreal> underDiagMatrices = socket_underDiagMatrices.getDataHandle();
  DataHandle<CFreal> aboveDiagMatrices = socket_aboveDiagMatrices.getDataHandle();

  for (CFuint i = 0; i < diagMatrices.size(); ++i) {
    diagMatrices[i] = 0.0;
    underDiagMatrices[i] = 0.0;
    aboveDiagMatrices[i] = 0.0;
  }
}

//////////////////////////////////////////////////////////////////////////////

PetscErrorCode TridiagPcApply(PC pc, Vec X, Vec Y)
{
  // cout << "TridiagPcApply - begin" << endl;
  // X is input vector - vector to be preconditioned
  // Y is output vector - preconditioned vector X
  // input vector X is the RHS in our formulation 
  // ouput vector Y is the SOLUTION vector
  CFreal* x;
  CFreal* y;

  /// I have to do something with this - this is not the best implementation
  CFreal* rho;

  PetscFunctionBegin;
  
  // v3.2
  void* ctx;  CF_CHKERRCONTINUE(PCShellGetContext(pc,&ctx));
   
  TridiagPcJFContext* pcContext = (TridiagPcJFContext*)(ctx);

  // getting data handle to diagonal inverted matrices
  DataHandle<CFreal> diagMatrices = pcContext->diagMatrices->getDataHandle(); // (c_{i}*mju_{i} + a_{i})^(-1)
  DataHandle<CFreal> aboveDiagMatrices = pcContext->aboveDiagMatrices->getDataHandle(); // mju_{i}
  DataHandle<CFreal> underDiagMatrices = pcContext->underDiagMatrices->getDataHandle(); // c_{i}

  DataHandle<CFint> upLocalIDsAll = pcContext->upLocalIDsAll->getDataHandle(); // updatable state IDs (-1 i not updatable)

  // putting X vector (vector to be preconditioned) into array "rhs" and Y output preconditioned vector into array "y"
  CF_CHKERRCONTINUE(VecGetArray(X, &x));
  CF_CHKERRCONTINUE(VecGetArray(Y, &y));

  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  const CFuint nbEqs2 = nbEqs*nbEqs;

  // getting nuber of equations
  //DataHandle<State*, GLOBAL> states = pcContext->pJFC->states->getDataHandle();

  const CFint nbUpdatableStates = diagMatrices.size()/nbEqs2;

  /// I have to do something with this - this is not the best implementation
  rho = new CFreal[nbUpdatableStates*nbEqs];

  RealVector tmpX(nbEqs, &x[0]);
  RealVector tmpY(nbEqs, &y[0]);
  RealVector tmpRho(nbEqs, &rho[0]);
  RealVector tmpRhoPrev(nbEqs, &rho[0]);

  /*
  CFreal elem = 0.;
  CFuint count = 0;
  for(CFuint j = 0; j < nbEqs; j++) {
    elem = count+1 + count+2; count += 2;
    tmpX[j] = elem;
  }
  for(CFint i = 1; i < nbUpdatableStates-1; i++) {
    tmpX.wrap(&x[i*nbEqs]);
    for(CFuint j = 0; j < nbEqs; j++) {
      elem = count+1 + count+2 + count+3; count += 3;
      tmpX[j] = elem;
    }
  }
  tmpX.wrap(&x[(nbUpdatableStates-1)*nbEqs]);
  for(CFuint j = 0; j < nbEqs; j++) {
    elem = count+1 + count+2; count += 2;
    tmpX[j] = elem;
  }

  for(int i = 0; i < nbUpdatableStates; i++) {
    tmpX.wrap(&x[i*nbEqs]);
    cout << "RHS" << i << " = " << tmpX << endl;
  }
  */

  RealMatrix nu(nbEqs, nbEqs, &diagMatrices[0]);
  RealMatrix mju(nbEqs, nbEqs, &aboveDiagMatrices[0]);
  RealMatrix c(nbEqs, nbEqs, &underDiagMatrices[0]);

  // filling of rho_{i} = nu_{i}*(RHS_{i} - c_{i}*rho_{i-1})
  tmpX.wrap(nbEqs,&x[0]);
  tmpRho.wrap(nbEqs,&rho[0]);
  tmpRho = nu*tmpX;
  //cout << "i = " << 0 << " tmpRho = " << tmpRho << endl;
  for(CFint i = 1; i < nbUpdatableStates; ++i) {
    const CFuint startIdxVec = i*nbEqs;
    const CFuint startIdxMat = i*nbEqs2;

    tmpX.wrap(nbEqs,&x[startIdxVec]);
    tmpRho.wrap(nbEqs,&rho[startIdxVec]);
    tmpRhoPrev.wrap(nbEqs,&rho[startIdxVec-nbEqs]);

    c.wrap(nbEqs, nbEqs, &underDiagMatrices[startIdxMat]);
    nu.wrap(nbEqs, nbEqs, &diagMatrices[startIdxMat]);

    tmpRho = nu*(tmpX - c*tmpRhoPrev);
    tmpRho = c*tmpRhoPrev;
    tmpX -= tmpRho;
    tmpRho = nu*tmpX;
    //cout << "i = " << i << " tmpRho = " << tmpRho << endl;
    //cout << "i = " << i << " c = " << c << endl;
    //cout << "i = " << i << " nu = " << nu << endl;
    //mju.wrap(&aboveDiagMatrices[i*nbEqs2]); cout << "i = " << i << " mju = " << mju << endl;
    //cout << "i = " << i << " tmpX = " << tmpX << endl;
  }

  // backward step
  tmpY.wrap(nbEqs, &y[(nbUpdatableStates-1)*nbEqs]);
  tmpRho.wrap(nbEqs, &rho[(nbUpdatableStates-1)*nbEqs]);
  tmpY = tmpRho;
  for(CFint i = nbUpdatableStates-2; i >= 0; i--) {
    const CFuint startIdxVec = i*nbEqs;

    tmpY.wrap(nbEqs, &y[startIdxVec]);
    tmpX.wrap(nbEqs, &y[startIdxVec+nbEqs]);
    tmpRho.wrap(nbEqs, &rho[startIdxVec]);

    mju.wrap(nbEqs, nbEqs, &aboveDiagMatrices[i*nbEqs2]);

    tmpY = mju*tmpX;
    //cout << "i = " << i << " tmpX = " << tmpX << endl;
    tmpY += tmpRho;
    //cout << "i = " << i << " RESULT = " << tmpY << endl;
  }

  // restoring of arrays X - vector to be preconditioned and Y - preconditioned vector
  CF_CHKERRCONTINUE(VecRestoreArray(X, &x));
  CF_CHKERRCONTINUE(VecRestoreArray(Y, &y));

  /// I have to do something with this - this is not the best implementation
  delete rho;

// cout << "TridiagPcApply - begin" << endl;

  PetscFunctionReturn(0);
}

//////////////////////////////////////////////////////////////////////////////

vector<SafePtr<BaseDataSocketSink> > TridiagPreconditioner::needsSockets()
{
  vector<SafePtr<BaseDataSocketSink> > result = ShellPreconditioner::needsSockets();

  result.push_back(&socket_diagMatrices);
  result.push_back(&socket_underDiagMatrices);
  result.push_back(&socket_aboveDiagMatrices);

  result.push_back(&socket_upLocalIDsAll);

  result.push_back(&socket_dplrIsFirstInLine);
  result.push_back(&socket_dplrToLocalIDs);
  result.push_back(&socket_localToDplrIDs);
  result.push_back(&socket_dplrCellInLine);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

} // namespace Petsc

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
