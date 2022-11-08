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
#include "Petsc/ParMFSolveSys.hh"
#include "Common/PE.hh"
#include <fstream>

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Petsc {

    extern PetscErrorCode computeMFMat(Mat petscMat, Vec x, Vec y);

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<ParMFSolveSys, PetscLSSData, PetscModule> 
parMFSolveSysProvider("ParMFSolveSys");

MethodCommandProvider<ParMFSolveSys, PetscLSSData, PetscModule> 
seqMFSolveSysProvider("SeqMFSolveSys");

//////////////////////////////////////////////////////////////////////////////

ParMFSolveSys::ParMFSolveSys(const string& name) :
  StdParSolveSys(name)
{
}

//////////////////////////////////////////////////////////////////////////////

ParMFSolveSys::~ParMFSolveSys()
{
}

//////////////////////////////////////////////////////////////////////////////

void ParMFSolveSys::execute()
{
  CFAUTOTRACE;

  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();
  DataHandle< CFreal> rhs = socket_rhs.getDataHandle();

  cf_assert(_upLocalIDs.size() == _upStatesGlobalIDs.size());
  const CFuint vecSize = _upLocalIDs.size();

  PetscMatrix& precMat = getMethodData().getMatrix(); // preconditioner matrix
  precMat.finalAssembly(); // assemble the matrix

  PetscMatrix& jfMat = getMethodData().getJFMatrix(); // Jacobian-free matrix
  jfMat.setJFFunction((void (*)(void))computeMFMat);

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
  CFuint ierr = KSPSetOperators(ksp, jfMat.getMat(), precMat.getMat());
#else
  // once you use a Petsc matrix as preconditioner matrix, you can choose standard preconditioners from Petsc 
  CFuint ierr = KSPSetOperators(ksp, jfMat.getMat(), precMat.getMat(), DIFFERENT_NONZERO_PATTERN);
#endif
  
  ierr = KSPSetUp(ksp);
  CHKERRCONTINUE(ierr);
  precMat.getMat();

  ierr = KSPSolve(ksp, rhsVec.getVec(), solVec.getVec());
  CHKERRCONTINUE(ierr);

  CFint iter = 0;
  ierr = KSPGetIterationNumber(ksp, &iter);
  CHKERRCONTINUE(ierr);

  CFLog(INFO, "KSP convergence reached at iteration: " << iter << "\n");

  solVec.copy(&rhs[0], &_upLocalIDs[0], vecSize);
}

//////////////////////////////////////////////////////////////////////////////

void ParMFSolveSys::setup()
{
  CFAUTOTRACE;

  StdParSolveSys::setup();

  MFContext* mfc = getMethodData().getMFContext();
  PetscMatrix& mat = getMethodData().getMatrix(); // linear system matrix
  mfc->mat = &mat;

  // set the preconditioner for later use
  getMethodData().getShellPreconditioner()->setPreconditioner();
}

//////////////////////////////////////////////////////////////////////////////

PetscErrorCode computeMFMat(Mat petscMat, Vec x, Vec y)
{
  void* ctx;

  PetscFunctionBegin;

  CF_CHKERRCONTINUE(MatShellGetContext(petscMat, &ctx));

  MFContext* mfc = (MFContext*)(ctx);

  CF_CHKERRCONTINUE(MatMult(mfc->mat->getMat(), x, y));

  PetscFunctionReturn(0);
}

//////////////////////////////////////////////////////////////////////////////

vector<SafePtr<BaseDataSocketSink> > ParMFSolveSys::needsSockets()
{
  vector<SafePtr<BaseDataSocketSink> > result = StdParSolveSys::needsSockets();

  return result;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace Petsc

} // namespace COOLFluiD

//  LocalWords:  rhsArrayVecstatesArray
