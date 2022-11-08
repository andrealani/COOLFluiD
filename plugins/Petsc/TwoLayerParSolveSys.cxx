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

#include "Petsc/Petsc.hh"
#include "Petsc/TwoLayerParSolveSys.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace Petsc {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<TwoLayerParSolveSys, PetscLSSData, PetscModule> twoLayerParSolveSysProvider("TwoLayerParSolveSys");

//////////////////////////////////////////////////////////////////////////////

TwoLayerParSolveSys::TwoLayerParSolveSys(const std::string& name) :
  PetscLSSCom(name),
  socket_states("states"),
  socket_rhs("rhs"),
  socket_interRhs("interRhs"),
  _upLocalIDs(),
  _upStatesGlobalIDs()
{
}

//////////////////////////////////////////////////////////////////////////////

TwoLayerParSolveSys::~TwoLayerParSolveSys()
{
}

//////////////////////////////////////////////////////////////////////////////

vector<SafePtr<BaseDataSocketSink> >
TwoLayerParSolveSys::needsSockets()
{
  vector<SafePtr<BaseDataSocketSink> > result;

  result.push_back(&socket_states);
  result.push_back(&socket_rhs);
  result.push_back(&socket_interRhs);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

void TwoLayerParSolveSys::execute()
{

  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();
  DataHandle<CFreal> rhs = socket_rhs.getDataHandle();
  DataHandle<CFreal> interRhs = socket_interRhs.getDataHandle();

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
  // For the Two Layer system, the interRhs and rhs
  // are copied alteratively
  // Entries whose indexes are == -1 (in _upStatesGlobalIDs)
  // will be ignored
  for (CFuint i=0;i< nbStates ;++i){
/*CFout << "State i: " << i << "  - isParUpdatable: " << _states[i]->isParUpdatable() << "\n" << CFendl;
CFout << "State i: " << i << "  - upStatesGlobalIDs: " << _upStatesGlobalIDs[2*i*nbEqs] << "\n" << CFendl;*/
      rhsVec.setValues(nbEqs, &_upStatesGlobalIDs[2*i*nbEqs], &interRhs[i*nbEqs]);
      rhsVec.setValues(nbEqs, &_upStatesGlobalIDs[(2*i*nbEqs)+nbEqs], &rhs[i*nbEqs]);
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

  rhs = 0.0;
  interRhs = 0.0;

  // copy solVec into the rhs
  /// @todo change this for returning separately the interRhs and rhs
  CFuint countUp = 0;
  for (CFuint i=0;i < nbStates;++i){
    if(states[i]->isParUpdatable()){
      solVec.copyValues(&interRhs[i*nbEqs], _upLocalIDs[2*countUp*nbEqs], nbEqs);
      solVec.copyValues(&rhs[i*nbEqs], _upLocalIDs[(2*countUp*nbEqs)+nbEqs], nbEqs);
//       CFout << "Cop rhs[" <<i << "] to " << _upLocalIDs[(2*countUp*nbEqs)+nbEqs] << "\n" << CFendl;
      ++countUp;
    }
  }
// CFout << "OUTPUT solVec " <<"\n" << CFendl;
// solVec.printToScreen();
// CFout << "OUTPUT rhs " <<"\n" << CFendl;
// for(CFuint i=0;i<_rhs.size();++i)
// {
// CFout << _rhs[i] << "\n" << CFendl;
// }

}

//////////////////////////////////////////////////////////////////////////////

void TwoLayerParSolveSys::setup()
{

  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();

  const CFuint nbEqs = getMethodData().getNbSysEquations();
  const CFuint nbStates = states.size();
  const CFuint vecSize =  2 * nbStates * nbEqs;

  // indexes for the insertion of elements in a PetscVector
  _upStatesGlobalIDs.resize(vecSize);

  const LSSIdxMapping& idxMapping = getMethodData().getLocalToGlobalMapping();

  // indexes are set to global ID for updatable states or
  // -1 for not updatable states
  for (CFuint i = 0; i < nbStates; ++i) {
    const CFint interStartID = 2 * i * nbEqs;
    const CFint interGlobalID = 2 * static_cast<CFint>(idxMapping.getColID(states[i]->getLocalID())) * nbEqs;

    const CFint startID = 2 * i * nbEqs + nbEqs;
    const CFint globalID = (2 * static_cast<CFint>(idxMapping.getColID(states[i]->getLocalID())) * nbEqs) + nbEqs;

    const bool isUpdatableState = states[i]->isParUpdatable();

    for (CFuint iEq = 0; iEq < nbEqs; ++iEq) {
      _upStatesGlobalIDs[interStartID + iEq] = (isUpdatableState) ? static_cast<CFint>(interGlobalID + iEq) : -1;
    }
    for (CFuint iEq = 0; iEq < nbEqs; ++iEq) {
      _upStatesGlobalIDs[startID + iEq] = (isUpdatableState) ? static_cast<CFint>(globalID + iEq) : -1;
    }
  }

/*  for (CFuint i = 0; i < _upStatesGlobalIDs.size(); ++i) {
CFout << i << "   -   " << _upStatesGlobalIDs[i] << "\n" << CFendl;
CFout << i << "   -   " << _upStatesGlobalIDs[i+1] << "\n" << CFendl;
}*/
  // resize _upGlobalIDs and _upLocalIDs
  CFuint nbUpdatableStates = 0;
  for (CFuint i = 0; i < nbStates; ++i) {
    if (states[i]->isParUpdatable()) {
      ++nbUpdatableStates;
    }
  }
  _upLocalIDs.reserve(2*nbUpdatableStates*nbEqs);

  // set the global and local IDs of the updatable states
CFuint countUp = 0;
  for (CFuint i = 0; i < nbStates; ++i) {
    if (states[i]->isParUpdatable()) {
      CFuint interLocalID = 2*countUp*nbEqs;
      CFuint localID      = 2*countUp*nbEqs + nbEqs;
      for (CFuint iEq = 0; iEq < nbEqs; ++iEq) {
        _upLocalIDs.push_back(static_cast<CFint>(interLocalID + iEq));
      }
      for (CFuint iEq = 0; iEq < nbEqs; ++iEq) {
        _upLocalIDs.push_back(static_cast<CFint>(localID + iEq));
      }
countUp++;
    }
  }

  cf_assert(nbUpdatableStates == _upLocalIDs.size()/(2*nbEqs));
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace Petsc

} // namespace COOLFluiD
