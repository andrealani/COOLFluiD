// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Petsc/PetscHeaders.hh" // must come before any header

#include "Framework/MeshData.hh"
#include "Framework/GlobalJacobianSparsity.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/LSSIdxMapping.hh"
#include "Framework/State.hh"
#include "Framework/SpaceMethod.hh"
#include "Framework/MethodData.hh"

#include "Petsc/PetscIdxMapping.hh"
#include "Petsc/Petsc.hh"
#include "Petsc/NewParSetup.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace Petsc {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<NewParSetup, PetscLSSData, PetscModule> newParSetupProvider("NewParSetup");

//////////////////////////////////////////////////////////////////////////////

NewParSetup::NewParSetup(const std::string& name) :
  BaseSetup(name)
{
}

//////////////////////////////////////////////////////////////////////////////

NewParSetup::~NewParSetup()
{
}

//////////////////////////////////////////////////////////////////////////////

void NewParSetup::setIdxMapping()
{
  CFAUTOTRACE;

  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();
  DataHandle < Framework::Node*, Framework::GLOBAL > nodes = socket_nodes.getDataHandle();
  const bool useNodeBased = getMethodData().useNodeBased();
  const CFuint nbStates = (!useNodeBased) ? states.size() : nodes.size();
  
  std::valarray<CFuint> stateGlobalIDs(nbStates);
  std::valarray<bool> isGhost(nbStates);
  
  if (!useNodeBased) {
    // build a contiguos mapping
    states.buildContiguosGlobal();
    for (CFuint i = 0; i< nbStates; ++i) {
      stateGlobalIDs[i] = states.getContiguosID(states[i]->getLocalID());
      isGhost[i] = !states[i]->isParUpdatable();
    }
    // free contiguos mapping
    states.freeContiguosGlobal();
  }
  else {
    // build a contiguos mapping
    nodes.buildContiguosGlobal();
    for (CFuint i = 0; i< nbStates; ++i) {
      stateGlobalIDs[i] = nodes.getContiguosID(nodes[i]->getLocalID());
      isGhost[i] = !nodes[i]->isParUpdatable();
    }
    // free contiguos mapping
    nodes.freeContiguosGlobal();
  }
  
  getMethodData().getLocalToGlobalMapping().createMapping(stateGlobalIDs, isGhost);
}
      
//////////////////////////////////////////////////////////////////////////////

void NewParSetup::setMatrix(const CFuint localSize,
                            const CFuint globalSize)
{
  CFAUTOTRACE;
  
  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();
  DataHandle < Framework::Node*, Framework::GLOBAL > nodes = socket_nodes.getDataHandle();
  const bool useNodeBased = getMethodData().useNodeBased();
  const CFuint nbStates = (!useNodeBased) ? states.size() : nodes.size();
  const CFuint nbEqs = getMethodData().getNbSysEquations();
  
  CFLog(VERBOSE, "NewParSetup::setMatrix() => useNodeBased[" << useNodeBased << "], nbEqs[" << nbEqs << "]\n");
  
  // all non zero entries
  std::valarray<CFint> allNonZero(nbStates);
  allNonZero = 0;

  // ghost neighbors (out of diagonal submatrix non zero entries)
  std::valarray<CFint> outDiagNonZero(nbStates);
  outDiagNonZero = 0;

  SelfRegistPtr<GlobalJacobianSparsity> sparsity =
    getMethodData().getCollaborator<SpaceMethod>()->createJacobianSparsity();
  
  sparsity->setDataSockets(socket_states, socket_nodes, socket_bStatesNeighbors);
  (!useNodeBased) ? sparsity->computeNNz(allNonZero, outDiagNonZero) :
    sparsity->computeNNzNodeBased(allNonZero, outDiagNonZero);
  
  const CFuint localSizeUp = (!getMethodData().useGPU() && !getMethodData().useAIJ()) ? localSize : localSize*nbEqs; 
  std::valarray<CFint> allNonZeroUp((CFint)0,localSizeUp);
  std::valarray<CFint> outDiagNonZeroUp((CFint)0,localSizeUp);
  
  CFuint countUp = 0;
  if (!getMethodData().useGPU() && !getMethodData().useAIJ()) {
    for (CFuint i = 0; i < nbStates; ++i) {
      if ((!useNodeBased && states[i]->isParUpdatable()) ||
	  (useNodeBased && nodes[i]->isParUpdatable())) {
	allNonZeroUp[countUp]     = allNonZero[i];
	outDiagNonZeroUp[countUp] = outDiagNonZero[i];
	++countUp;
      }
    }
  }
  else {
    for (CFuint i = 0; i < nbStates; ++i) {
      if ((!useNodeBased && states[i]->isParUpdatable()) ||
	  (useNodeBased && nodes[i]->isParUpdatable())) {
	const CFuint start = i*nbEqs;
	for (CFuint iEq = 0; iEq < nbEqs; ++iEq) {
	  const CFuint idx = start+iEq;  
	  allNonZeroUp[countUp]     = allNonZero[idx]; 
	  outDiagNonZeroUp[countUp] = outDiagNonZero[idx];
	  ++countUp;
	}
      }
    }
    cf_assert(countUp == localSize*nbEqs);
  }
  
  PetscMatrix& mat = getMethodData().getMatrix();
  const string nsp = getMethodData().getNamespace();
  
  // create a parallel sparse matrix in block compressed row format
  mat.setGPU(getMethodData().useGPU());
  mat.setAIJ(getMethodData().useAIJ());
  mat.createParBAIJ(PE::GetPE().GetCommunicator(nsp),
                    nbEqs,
                    localSize*nbEqs,
                    localSize*nbEqs,
                    globalSize*nbEqs,
                    globalSize*nbEqs,
                    0, &allNonZeroUp[0],
                    0, &outDiagNonZeroUp[0],
                    "Jacobian");
}
      
//////////////////////////////////////////////////////////////////////////////

    } // namespace Petsc

} // namespace COOLFluiD
