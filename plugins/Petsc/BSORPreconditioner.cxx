// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Petsc/BSORPreconditioner.hh"
#include "Petsc/PetscLSSData.hh"
#include "Petsc/Petsc.hh"

#include "petscis.h"  

#include "Framework/MethodStrategyProvider.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Petsc {
    
    extern PetscErrorCode BSORPcApply(PC pc, Vec x, Vec y);
    
//////////////////////////////////////////////////////////////////////////////

MethodStrategyProvider<BSORPreconditioner,
                       PetscLSSData,
                       ShellPreconditioner,
                       PetscModule>
BSORPreconditionerProvider("BSOR");

//////////////////////////////////////////////////////////////////////////////

void BSORPreconditioner::defineConfigOptions(Config::OptionList& options)
{
  
}

//////////////////////////////////////////////////////////////////////////////

BSORPreconditioner::BSORPreconditioner(const std::string& name) :
                     ShellPreconditioner(name),
                     _pcc()
{
  addConfigOptionsTo(this);
}

//////////////////////////////////////////////////////////////////////////////

BSORPreconditioner::~BSORPreconditioner()
{
}

//////////////////////////////////////////////////////////////////////////////

void BSORPreconditioner::setPreconditioner()
{
  // getting the MFContext pointer into BSORPcContext pcc
  _pcc.pMFC = getMethodData().getMFContext();
  //_pcc.relaxType = _relaxType;

  // here set the shell preconditioner
  CF_CHKERRCONTINUE(PCShellSetContext(_pcc.pMFC->petscData->getPreconditioner(), &_pcc));
  CF_CHKERRCONTINUE(PCShellSetApply(_pcc.pMFC->petscData->getPreconditioner(), BSORPcApply));
}

//////////////////////////////////////////////////////////////////////////////

PetscErrorCode BSORPcApply(PC pc, Vec x, Vec y)
{ 
  // v3.2
  void* ctx;  CF_CHKERRCONTINUE(PCShellGetContext(pc,&ctx));
  
  BSORPcContext* pcContext = (BSORPcContext*)(ctx);
  Mat mat = pcContext->pMFC->mat->getMat();

  // under/over relaxation factor - but current version of PETSc (3.0.0-p2) does not have implemented this functionality 
  // for this type of matrix (BAIJ - petsc block matrix)
  const CFreal omega = 1.0;

  static CFint localRangeM; // M - the global index of the first local row
  static CFint localRangeN; // N - one more than the global index of the last local row
  static CFint localSize;

  static Mat* subMat;

  CFreal *xArray;
  CFreal *yArray;

  Vec xLocal;
  Vec yLocal;

  // this can be static because it doesn't change in the whole computation
  static bool idxIsSet = false;
  static IS indexSet;

  PetscFunctionBegin;

  // if indexes of submatrices haven't been set yet - they are set in first call of this function
  if (!idxIsSet) {
    idxIsSet = true;

    CF_CHKERRCONTINUE(MatGetOwnershipRange(mat, &localRangeM, &localRangeN)); // getting the range over processors
    localSize = localRangeN - localRangeM;

    CFint* idx = new CFint[localSize];
    for (CFint i = 0; i < localSize; i++) {
      idx[i] = i+localRangeM;
    }

    // v3.2: which one to use 
    // typedef enum { PETSC_COPY_VALUES, PETSC_OWN_POINTER, PETSC_USE_POINTER} PetscCopyMode;

    CF_CHKERRCONTINUE(ISCreateGeneral(MPI_COMM_SELF, localSize, idx, PETSC_COPY_VALUES, &indexSet));
    delete idx;

    CF_CHKERRCONTINUE(MatGetSubMatrices(mat, 1, &indexSet, &indexSet, MAT_INITIAL_MATRIX, &subMat));
  }
  else {
    CF_CHKERRCONTINUE(MatGetSubMatrices(mat, 1, &indexSet, &indexSet, MAT_REUSE_MATRIX, &subMat));
  }

  CF_CHKERRCONTINUE(VecGetArray(x, &xArray));
#if PETSC_VERSION_MINOR==4 || PETSC_VERSION_MINOR==6 || PETSC_VERSION_MINOR==7 || PETSC_VERSION_MINOR==9 || PETSC_VERSION_MINOR==11 || PETSC_VERSION_MINOR==12 || PETSC_VERSION_MINOR==15 || PETSC_VERSION_MINOR==18
  CFuint blockSize = 1; 
  CF_CHKERRCONTINUE(VecCreateSeqWithArray(PETSC_COMM_SELF, blockSize, localSize, xArray, &xLocal));
#else
  CF_CHKERRCONTINUE(VecCreateSeqWithArray(PETSC_COMM_SELF, localSize, xArray, &xLocal));
#endif
  CF_CHKERRCONTINUE(VecGetArray(y, &yArray));
#if PETSC_VERSION_MINOR==4 || PETSC_VERSION_MINOR==6 || PETSC_VERSION_MINOR==7 || PETSC_VERSION_MINOR==9 || PETSC_VERSION_MINOR==11 || PETSC_VERSION_MINOR==12 || PETSC_VERSION_MINOR==15 || PETSC_VERSION_MINOR==18 
  CF_CHKERRCONTINUE(VecCreateSeqWithArray(PETSC_COMM_SELF, blockSize, localSize, yArray, &yLocal));
#else
  CF_CHKERRCONTINUE(VecCreateSeqWithArray(PETSC_COMM_SELF, localSize, yArray, &yLocal));
#endif

  //CF_CHKERRCONTINUE(MatGetSubMatrices(mat, 1, &indexSet, &indexSet, MAT_INITIAL_MATRIX, &subMat));
  
  // MatPBRelax removed in v3.1
  // CF_CHKERRCONTINUE(MatPBRelax(subMat[0], xLocal, omega, MatSORType(SOR_LOCAL_SYMMETRIC_SWEEP+SOR_ZERO_INITIAL_GUESS), 0., 1, 1, yLocal)); 
  CF_CHKERRCONTINUE(MatSOR(subMat[0], xLocal, omega, MatSORType(SOR_LOCAL_SYMMETRIC_SWEEP+SOR_ZERO_INITIAL_GUESS), 0., 1, 1, yLocal));
  //CF_CHKERRCONTINUE(MatDestroyMatrices(1, &subMat));
  
  CF_CHKERRCONTINUE(VecRestoreArray(x, &xArray));
  CF_CHKERRCONTINUE(VecRestoreArray(y, &yArray));
  
  // v3.0
  // CF_CHKERRCONTINUE(VecDestroy(xLocal));
  // CF_CHKERRCONTINUE(VecDestroy(yLocal));
  // v3.2
  CF_CHKERRCONTINUE(VecDestroy(&xLocal));
  CF_CHKERRCONTINUE(VecDestroy(&yLocal));
    
  PetscFunctionReturn(0);
}

//////////////////////////////////////////////////////////////////////////////

vector<SafePtr<BaseDataSocketSink> > BSORPreconditioner::needsSockets()
{
  vector<SafePtr<BaseDataSocketSink> > result = ShellPreconditioner::needsSockets();

  return result;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace Petsc

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

