// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Petsc/BaseSetup.hh" // must come first because includes PetscHeaders

#include "Framework/MeshData.hh"
#include "Framework/State.hh"
#include "Common/PE.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace Petsc {

//////////////////////////////////////////////////////////////////////////////

/**
 * Monitor PETSc solver convergence
 */
int monitorPetscConvergence(KSP ksp, CFint it, CFdouble rnorm, void *ctx)
{
  CFout << "PETSc Iter: " << it << " Norm: " << rnorm << "\n";
  return 0;
}
      
//////////////////////////////////////////////////////////////////////////////

BaseSetup::BaseSetup(const std::string& name) :
  PetscLSSCom(name),
  socket_bStatesNeighbors("bStatesNeighbors"),
  socket_states("states"),
  socket_nodes("nodes"),
  socket_rhs("rhs")
{
}

//////////////////////////////////////////////////////////////////////////////

BaseSetup::~BaseSetup()
{
}

//////////////////////////////////////////////////////////////////////////////

void BaseSetup::execute()
{
  CFAUTOTRACE;
  
  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();
  DataHandle < Framework::Node*, Framework::GLOBAL > nodes = socket_nodes.getDataHandle();
  const bool useNodeBased = getMethodData().useNodeBased();
  const CFuint nbStates = (!useNodeBased) ? states.size() : nodes.size();
  
  // set the index mapping (global IDs to global Petsc IDs)
  setIdxMapping();

  CFuint nbUpdatableStates = 0;
  for (CFuint i = 0; i < nbStates; ++i) {
    if ((!useNodeBased && states[i]->isParUpdatable()) ||
	(useNodeBased && nodes[i]->isParUpdatable())) {
      ++nbUpdatableStates;
    }
  }
  
  CFuint totalNbStates = 0;
  if (PE::GetPE().IsParallel()) {
    totalNbStates = (!useNodeBased) ? states.getGlobalSize() : nodes.getGlobalSize();
  }
  else {
    totalNbStates = nbStates;
  }

  // set the vectors
  setVectors(nbUpdatableStates, totalNbStates);

  // set the matrix
  setMatrix(nbUpdatableStates, totalNbStates);

  // set the Krylov Sub Space solver
  setKSP();
}

//////////////////////////////////////////////////////////////////////////////

void BaseSetup::setKSP()
{
  CFAUTOTRACE;

  // Creation of the linear solver context
  PC& pc = getMethodData().getPreconditioner();
  KSP& ksp = getMethodData().getKSP();
  const string nsp = getMethodData().getNamespace();
  
  int ierr = 0;
  ierr = KSPCreate(PE::GetPE().GetCommunicator(nsp),&ksp);
  CHKERRCONTINUE(ierr);
  
  // set the pointer to the preconditioner
  ierr = KSPGetPC(ksp, &pc);
  CHKERRCONTINUE(ierr);
  
  // set the type of preconditioner
  ierr = PCSetType(pc, getMethodData().getPCType());
  CHKERRCONTINUE(ierr);
  
  ierr = PCFactorSetLevels(pc, getMethodData().getILULevels());
  CHKERRCONTINUE(ierr);

  ierr = PCFactorSetMatOrderingType(pc, getMethodData().getMatOrderType());
  CHKERRCONTINUE(ierr);

  // set the type of Krylov subspace method
  ierr = KSPSetType(ksp, getMethodData().getKSPType());
  CHKERRCONTINUE(ierr);

  // set the number of krylov spaces
  ierr = KSPGMRESSetRestart(ksp, static_cast<int>(getMethodData().getNbKSP()));
  CHKERRCONTINUE(ierr);

  // set convergence tolerances
  ierr = KSPSetTolerances(ksp,
                          getMethodData().getRelativeTol(),
                          getMethodData().getAbsoluteTol(),
                          getMethodData().getDivergenceTol(),
                          getMethodData().getMaxIterations());
  CHKERRCONTINUE(ierr);

  // sets the monitor for the convergence
  if (getMethodData().isOutput()) {
    ierr = KSPMonitorSet(ksp,monitorPetscConvergence,CFNULL,CFNULL);
    CHKERRCONTINUE(ierr);
  }
}

//////////////////////////////////////////////////////////////////////////////

void BaseSetup::setVectors(const CFuint localSize,
			   const CFuint globalSize)
{
  CFAUTOTRACE;

  PetscVector& sol = getMethodData().getSolVector();
  PetscVector& rhs = getMethodData().getRhsVector();

  const CFuint nbEqs = getMethodData().getNbSysEquations();
  const string nsp = getMethodData().getNamespace();
  
  // create the solution and the rhs vectors
  rhs.setGPU(getMethodData().useGPU());
  rhs.create(PE::GetPE().GetCommunicator(nsp), localSize*nbEqs,
             globalSize*nbEqs,"rhs");
  
  sol.setGPU(getMethodData().useGPU());
  sol.create(PE::GetPE().GetCommunicator(nsp), localSize*nbEqs,
             globalSize*nbEqs,"Solution");
  
  // initialize the two vectors
  sol.initialize(PE::GetPE().GetCommunicator(nsp),1.0);
  rhs.initialize(PE::GetPE().GetCommunicator(nsp),0.0);
}
      
//////////////////////////////////////////////////////////////////////////////

vector<SafePtr<BaseDataSocketSink> > BaseSetup::needsSockets()
{
  vector<SafePtr<BaseDataSocketSink> > result;

  result.push_back(&socket_bStatesNeighbors);
  result.push_back(&socket_states);
  result.push_back(&socket_nodes);
  result.push_back(&socket_rhs);
  
  return result;
}

//////////////////////////////////////////////////////////////////////////////

} // namespace Petsc

} // namespace COOLFluiD
