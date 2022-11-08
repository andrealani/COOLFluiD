// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Petsc/ILUPreconditioner.hh"
#include "Petsc/PetscLSSData.hh"
#include "Petsc/Petsc.hh"

#include "Framework/SpaceMethodData.hh"
#include "Framework/MethodStrategyProvider.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace Petsc {

    extern PetscErrorCode ILUPcApply(PC pc, Vec X, Vec Y);

//////////////////////////////////////////////////////////////////////////////

MethodStrategyProvider<ILUPreconditioner,
                       PetscLSSData,
                       ShellPreconditioner,
                       PetscModule>
ILUPreconditionerProvider("ILU");

//////////////////////////////////////////////////////////////////////////////

void ILUPreconditioner::defineConfigOptions(Config::OptionList& options) { }

//////////////////////////////////////////////////////////////////////////////

ILUPreconditioner::ILUPreconditioner(const std::string& name) :
  ShellPreconditioner(name),
  _pcc()
{
  addConfigOptionsTo(this);
}

//////////////////////////////////////////////////////////////////////////////

ILUPreconditioner::~ILUPreconditioner()
{
  // release of PETSc index set
  // v3.0
  // CF_CHKERRCONTINUE(ISDestroy(_indexSet));
  // v3.2 
  CF_CHKERRCONTINUE(ISDestroy(&_indexSet));
  delete _indexSetArray;
}

//////////////////////////////////////////////////////////////////////////////

void ILUPreconditioner::setPreconditioner()
{
  // getting the JFContext pointer into ILUPcContext pcc
  _pcc.pJFC = getMethodData().getJFContext();

  PetscMatrix& precondMat = getMethodData().getPreconditionerMatrix();
  _pcc.precondMat = &precondMat;

  // here set the shell preconditioner
  CF_CHKERRCONTINUE(PCShellSetContext(_pcc.pJFC->petscData->getPreconditioner(), &_pcc));
  CF_CHKERRCONTINUE(PCShellSetApply(_pcc.pJFC->petscData->getPreconditioner(), ILUPcApply));

  _info.levels = 0; // it has to be set to 0, otherwise MatILUFactor() is not done in-place

  CFint rangeM;
  CFint rangeN;
  CF_CHKERRCONTINUE(MatGetOwnershipRange(_pcc.precondMat->getMat(), &rangeM, &rangeN)); // getting the range over processors

  const CFuint matSize = rangeN - rangeM;

  _indexSetArray = new CFint[matSize];
  for(CFuint i = 0; i < matSize; i++) _indexSetArray[i] = i;
  
  // old unsupported
  // CF_CHKERRCONTINUE(ISCreateGeneralWithArray(PETSC_COMM_SELF, matSize, _indexSetArray, &_indexSet));
  
  
  CF_CHKERRCONTINUE(ISCreateGeneral(PETSC_COMM_SELF, matSize, _indexSetArray, PETSC_USE_POINTER, &_indexSet));
  CF_CHKERRCONTINUE(ISSetPermutation(_indexSet));
}

//////////////////////////////////////////////////////////////////////////////

void ILUPreconditioner::computeBeforeSolving()
{
  CF_CHKERRCONTINUE(MatILUFactor(_pcc.precondMat->getMat(), _indexSet, _indexSet, &_info));
}

//////////////////////////////////////////////////////////////////////////////

void ILUPreconditioner::computeAfterSolving()
{
  MatSetUnfactored(_pcc.precondMat->getMat());
}

//////////////////////////////////////////////////////////////////////////////

PetscErrorCode ILUPcApply(PC pc, Vec X, Vec Y)
{
  // X is input vector - vector to be preconditioned
  // Y is output vector - preconditioned vector X
  PetscFunctionBegin;

  // v3.2
  void* ctx;  CF_CHKERRCONTINUE(PCShellGetContext(pc,&ctx));
    
  // getting the LU-SGS preconditioner context
  ILUPcContext* pcContext = (ILUPcContext*)(ctx);

  CFint localRangeM;
  CFint localRangeN;
  CF_CHKERRCONTINUE(MatGetOwnershipRange(pcContext->precondMat->getMat(), &localRangeM, &localRangeN));
  const CFint localSize = localRangeN - localRangeM;

  CFreal* xArray;
  Vec xLocal;
  CF_CHKERRCONTINUE(VecGetArray(X, &xArray));
#if PETSC_VERSION_MINOR==4 || PETSC_VERSION_MINOR==6 || PETSC_VERSION_MINOR==7 || PETSC_VERSION_MINOR==9 || PETSC_VERSION_MINOR==11 || PETSC_VERSION_MINOR==12 || PETSC_VERSION_MINOR==15 || PETSC_VERSION_MINOR==18
  CFuint blockSize = 1;
  CF_CHKERRCONTINUE(VecCreateSeqWithArray(PETSC_COMM_SELF, blockSize, localSize, xArray, &xLocal));
#else
  CF_CHKERRCONTINUE(VecCreateSeqWithArray(PETSC_COMM_SELF, localSize, xArray, &xLocal));
#endif

  CFreal* yArray;
  Vec yLocal;
  CF_CHKERRCONTINUE(VecGetArray(Y, &yArray));
#if PETSC_VERSION_MINOR==4 || PETSC_VERSION_MINOR==6 || PETSC_VERSION_MINOR==7 || PETSC_VERSION_MINOR==9 || PETSC_VERSION_MINOR==11 || PETSC_VERSION_MINOR==12 || PETSC_VERSION_MINOR==15 || PETSC_VERSION_MINOR==18
  CF_CHKERRCONTINUE(VecCreateSeqWithArray(PETSC_COMM_SELF, blockSize, localSize, yArray, &yLocal));
#else
  CF_CHKERRCONTINUE(VecCreateSeqWithArray(PETSC_COMM_SELF, localSize, yArray, &yLocal));
#endif
  CF_CHKERRCONTINUE(MatSolve(pcContext->precondMat->getMat(), xLocal, yLocal));

  CF_CHKERRCONTINUE(VecRestoreArray(X, &xArray));
  CF_CHKERRCONTINUE(VecRestoreArray(Y, &yArray));

  CF_CHKERRCONTINUE(VecDestroy(&xLocal));
  CF_CHKERRCONTINUE(VecDestroy(&yLocal));

  PetscFunctionReturn(0);
}

//////////////////////////////////////////////////////////////////////////////

vector<SafePtr<BaseDataSocketSink> > ILUPreconditioner::needsSockets()
{
  vector<SafePtr<BaseDataSocketSink> > result = ShellPreconditioner::needsSockets();

  return result;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace Petsc

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
