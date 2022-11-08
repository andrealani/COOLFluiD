// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Petsc/PetscHeaders.hh" // must come before any header

#include "Common/CFLog.hh"

#include "Framework/MeshData.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/SubSystemStatus.hh"
#include "Framework/LSSIdxMapping.hh"
#include "Framework/State.hh"
#include "Framework/PhysicalModel.hh"

#include "Petsc/Petsc.hh"
#include "Petsc/StdParSolveSys.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace Petsc {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<StdParSolveSys, PetscLSSData, PetscModule> stdParSolveSysProvider("StdParSolveSys");

MethodCommandProvider<StdParSolveSys, PetscLSSData, PetscModule> stdSeqSolveSysProvider("StdSeqSolveSys");

//////////////////////////////////////////////////////////////////////////////

StdParSolveSys::StdParSolveSys(const std::string& name) :
  PetscLSSCom(name),
  socket_states("states"),
  socket_nodes("nodes"),
  socket_rhs("rhs"), 
  _upLocalIDs(),
  _upStatesGlobalIDs()
{
}

//////////////////////////////////////////////////////////////////////////////

StdParSolveSys::~StdParSolveSys()
{
}

//////////////////////////////////////////////////////////////////////////////

void StdParSolveSys::execute()
{
  CFAUTOTRACE;
  
Stopwatch<WallTime> stopTimer;
stopTimer.start();

  DataHandle< CFreal> rhs = socket_rhs.getDataHandle();
  cf_assert(_upLocalIDs.size() == _upStatesGlobalIDs.size());
  const CFuint vecSize = _upLocalIDs.size();
  
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
  // _upStatesGlobalIDs[i] is different from _upLocalIDs[i]
  // in case multiple LSS are used
  for (CFuint i = 0; i < vecSize; ++i) {
    // cout << "nbSysEq = " << nbEqs << ", upLocalIDs = " << _upLocalIDs[i]
    // 	 << ", upStatesGlobalIDs = " << _upStatesGlobalIDs[i] << endl;
    rhsVec.setValue(_upStatesGlobalIDs[i], rhs[_upLocalIDs[i]]);
  }

  // assemble the rhs vector
  rhsVec.assembly();
  
  /*if (getMethodData().getNbSysEquations() == 9) {
    mat.printToFile("matPETSC.txt");
    rhsVec.printToFile("rhsPETSC.txt");
    exit(1);
    }*/
  
  const CFuint nbIter = SubSystemStatusStack::getActive()->getNbIter();
  if (getMethodData().getSaveRate() > 0) {
    if (getMethodData().isSaveSystemToFile() || (nbIter%getMethodData().getSaveRate() == 0)) { 
      const string mFile = "mat-iter" + StringOps::to_str(nbIter) + ".dat";
      mat.printToFile(mFile.c_str());
      // mat.printToScreen();
      
      const string vFile = "rhs-iter" + StringOps::to_str(nbIter) + ".dat";
      rhsVec.printToFile(vFile.c_str());
      // rhsVec.printToScreen();
    }
  }
 
  // reuse te preconditioner
#if PETSC_VERSION_MINOR==7 || PETSC_VERSION_MINOR==9 || PETSC_VERSION_MINOR==11 || PETSC_VERSION_MINOR==12 || PETSC_VERSION_MINOR==15 || PETSC_VERSION_MINOR==18
  PetscBool reusePC = ((nbIter-1)%getMethodData().getPreconditionerRate() == 0) ?
     PETSC_FALSE : PETSC_TRUE;
  CFLog(VERBOSE, "StdParSolveSys::execute() => reusePC [" << reusePC <<"]\n");
  CHKERRCONTINUE(KSPSetReusePreconditioner(ksp,reusePC));
  PC& pc = getMethodData().getPreconditioner();
  CHKERRCONTINUE(PCSetReusePreconditioner(pc,reusePC));
#endif
 
#if PETSC_VERSION_MINOR==6 || PETSC_VERSION_MINOR==7 || PETSC_VERSION_MINOR==9 || PETSC_VERSION_MINOR==11 || PETSC_VERSION_MINOR==12 || PETSC_VERSION_MINOR==15 || PETSC_VERSION_MINOR==18
  CFuint ierr = KSPSetOperators(ksp, mat.getMat(), mat.getMat());
#else
  CFuint ierr = KSPSetOperators
    (ksp, mat.getMat(), mat.getMat(),DIFFERENT_NONZERO_PATTERN);
#endif
  
  //This is to allow viewing the matrix structure in X windows
  //Works only if Petsc is compiled with X
  if(getMethodData().isShowMatrixStructure())
  {
    ierr = MatView(mat.getMat(), PETSC_VIEWER_DRAW_WORLD);
    CHKERRCONTINUE(ierr);
  }

  ierr = KSPSetUp(ksp);
  CHKERRCONTINUE(ierr);

#if PETSC_VERSION_MINOR==7 || PETSC_VERSION_MINOR==9 || PETSC_VERSION_MINOR==11 || PETSC_VERSION_MINOR==12 || PETSC_VERSION_MINOR==15 || PETSC_VERSION_MINOR==18
  if (getMethodData().getPCType()==PCBJACOBI && getMethodData().useGPU()) {
    PetscInt its,nlocal,first;
    KSP* subksp;
    PC   subpc;  
    PC& pc = getMethodData().getPreconditioner();
    PCBJacobiGetSubKSP(pc,&nlocal,&first,&subksp);
    for (CFuint i=0; i<nlocal; i++) {
      KSPGetPC(subksp[i],&subpc);
      PCSetType(subpc,PCILU);
    }
  }
#endif

  ierr = KSPSolve(ksp, rhsVec.getVec(), solVec.getVec());
  CHKERRCONTINUE(ierr);
  
  //solVec.printToFile("solPETSC.txt");
  
  CFint iter = 0;
  ierr = KSPGetIterationNumber(ksp, &iter);
  CHKERRCONTINUE(ierr);
  
  // Ask to stop the simulation if convergence is achieved at iteration 0 (i.e. LSS was not solved)
  if (iter == 0) {
    SubSystemStatusStack::getActive()->setStopSimulation(true);
  }
  
  if (nbIter%getMethodData().getKSPConvergenceShowRate() == 0) {
    CFLog(INFO, "KSP convergence reached at iteration: " << iter << "\n");
  }
  
  solVec.copy(&rhs[0], &_upLocalIDs[0], vecSize);
  CFLog(VERBOSE, "StdParSolveSys::execute() took " << stopTimer << "s\n");
}

//////////////////////////////////////////////////////////////////////////////

void StdParSolveSys::setup()
{
  CFAUTOTRACE;

  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();
  DataHandle < Framework::Node*, Framework::GLOBAL > nodes = socket_nodes.getDataHandle();
  const bool useNodeBased = getMethodData().useNodeBased();
  const CFuint nbStates = (!useNodeBased) ? states.size() : nodes.size();
  const CFuint nbEqs = getMethodData().getNbSysEquations();
  
  // count the number of updatable states
  CFuint nbUpdatableStates = 0;
  for (CFuint i = 0; i < nbStates; ++i) {
    if ((!useNodeBased && states[i]->isParUpdatable()) ||
	(useNodeBased && nodes[i]->isParUpdatable())) {
      ++nbUpdatableStates;
    }
  }

  const CFuint vecSize = nbUpdatableStates * nbEqs;
  // indexes for the insertion of elements in a PetscVector
  _upStatesGlobalIDs.reserve(vecSize);
  _upLocalIDs.reserve(vecSize);

  const LSSIdxMapping& idxMapping = getMethodData().getLocalToGlobalMapping();
  const CFuint totalNbEqs = PhysicalModelStack::getActive()->getNbEq();
  const std::valarray<bool>& maskArray = *getMethodData().getMaskArray();

  // indexes are set to global ID for updatable states
  for (CFuint i = 0; i < nbStates; ++i) {
    const CFuint localID = i*totalNbEqs;
    if ((!useNodeBased && states[i]->isParUpdatable()) ||
	(useNodeBased && nodes[i]->isParUpdatable())) {
      const CFuint sID = (!useNodeBased) ? states[i]->getLocalID() : nodes[i]->getLocalID();
      CFint globalID = static_cast<CFint>(idxMapping.getColID(sID))*nbEqs;
      
      for (CFuint iEq = 0; iEq < totalNbEqs; ++iEq) {
	if (maskArray[iEq]) {
	  _upStatesGlobalIDs.push_back(globalID++);
	  _upLocalIDs.push_back(static_cast<CFint>(localID + iEq));
	}
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

vector<SafePtr<BaseDataSocketSink> > StdParSolveSys::needsSockets()
{
  vector<SafePtr<BaseDataSocketSink> > result;
  
  result.push_back(&socket_states);
  result.push_back(&socket_nodes);
  result.push_back(&socket_rhs);
  
  return result;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace Petsc

} // namespace COOLFluiD
