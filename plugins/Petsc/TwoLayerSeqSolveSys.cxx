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

#include "Petsc/Petsc.hh"
#include "Petsc/TwoLayerSeqSolveSys.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace Petsc {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<TwoLayerSeqSolveSys, PetscLSSData, PetscModule> twoLayerSeqSolveSysProvider("TwoLayerSeqSolveSys");

//////////////////////////////////////////////////////////////////////////////

void TwoLayerSeqSolveSys::execute()
{
  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();
  DataHandle<CFreal> rhs = socket_rhs.getDataHandle();
  DataHandle<CFreal> interRhs = socket_interRhs.getDataHandle();


  cf_assert(_idx.size() == (rhs.size() + interRhs.size()));

  const CFuint nbStates = states.size();
  const CFuint nbEqs = getMethodData().getNbSysEquations();

  PetscMatrix& mat = getMethodData().getMatrix();
  PetscVector& rhsVec = getMethodData().getRhsVector();
  PetscVector& solVec = getMethodData().getSolVector();
  KSP& ksp = getMethodData().getKSP();

  // assemble the matrix
  mat.finalAssembly();

  // after the first iteration freeze the non zero structure
  // of the matrix
  if (SubSystemStatusStack::getActive()->getNbIter() == 1) {
    mat.freezeNonZeroStructure();
  }

  // the rhs is copied into the PetscVector for the rhs
  /// for the two layer system, the interRhs and rhs
  /// are copied alteratively
  for (CFuint i=0;i< nbStates ;++i){
    rhsVec.setValues(nbEqs, &_idx[2*i*nbEqs], &interRhs[i*nbEqs]);
    rhsVec.setValues(nbEqs, &_idx[(2*i*nbEqs)+nbEqs], &rhs[i*nbEqs]);
  }

  // assemble the rhs vector
  rhsVec.assembly();
  
#if PETSC_VERSION_MINOR==6 || PETSC_VERSION_MINOR==7 || PETSC_VERSION_MINOR==9 || PETSC_VERSION_MINOR==11 || PETSC_VERSION_MINOR==12 || PETSC_VERSION_MINOR==15 || PETSC_VERSION_MINOR==18
  CFuint ierr = KSPSetOperators(ksp, mat.getMat(), mat.getMat());
#else
  CFuint ierr = KSPSetOperators(ksp,
				mat.getMat(),
				mat.getMat(),
				DIFFERENT_NONZERO_PATTERN);
#endif
  CHKERRCONTINUE(ierr);

  // PetscLogInfo(mat.getMat(), "MATRIX memory allocation");

  ierr = KSPSetUp(ksp);
  CHKERRCONTINUE(ierr);

  ierr = KSPSolve(ksp, rhsVec.getVec(), solVec.getVec());
  CHKERRCONTINUE(ierr);

  CFint iter = 0;
  ierr = KSPGetIterationNumber(ksp, &iter);
  CHKERRCONTINUE(ierr);

  CFLog(INFO, "KSP convergence reached at iteration " <<
        iter << "\n");

  // copy solVec into the rhs
  /// @todo change this for returning separately the interRhs and rhs
  for (CFuint i=0;i < nbStates;++i){
    solVec.copyValues(&interRhs[i*nbEqs], _idx[2*i*nbEqs], nbEqs);
    solVec.copyValues(&rhs[i*nbEqs], _idx[(2*i*nbEqs)+nbEqs], nbEqs);
  }

}

//////////////////////////////////////////////////////////////////////////////

void TwoLayerSeqSolveSys::setup()
{
 DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();

  const CFuint vecSize =  2*states.size()*getMethodData().getNbSysEquations();
  // indexes for the insertion of elements in a PetscVector
  _idx.resize(vecSize);
  for (CFuint i = 0; i < vecSize; ++i) {
    _idx[i] = static_cast<CFint>(i);
  }
}

//////////////////////////////////////////////////////////////////////////////

vector<SafePtr<BaseDataSocketSink> > TwoLayerSeqSolveSys::needsSockets()
{
  vector<SafePtr<BaseDataSocketSink> > result;

  result.push_back(&socket_states);
  result.push_back(&socket_rhs);
  result.push_back(&socket_interRhs);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace Petsc

} // namespace COOLFluiD
