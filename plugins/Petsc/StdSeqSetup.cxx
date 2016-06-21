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

#include "Petsc/Petsc.hh"
#include "Petsc/StdSeqSetup.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace Petsc {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<StdSeqSetup, PetscLSSData, PetscModule> stdSeqSetupProvider("StdSeqSetup");

MethodCommandProvider<StdSeqSetup, PetscLSSData, PetscModule> newSeqSetupProvider("NewSeqSetup");

//////////////////////////////////////////////////////////////////////////////

StdSeqSetup::StdSeqSetup(const std::string& name) :
  BaseSetup(name)
{
}

//////////////////////////////////////////////////////////////////////////////

StdSeqSetup::~StdSeqSetup()
{
}

//////////////////////////////////////////////////////////////////////////////

void StdSeqSetup::setIdxMapping()
{
  CFAUTOTRACE;

  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();
  DataHandle < Framework::Node*, Framework::GLOBAL > nodes = socket_nodes.getDataHandle();
  const bool useNodeBased = getMethodData().useNodeBased();
  const CFuint nbStates = (!useNodeBased) ? states.size() : nodes.size();
  
  vector<CFuint> stateGlobalIDs(nbStates);
  for (CFuint i = 0; i < nbStates; ++i) {
    stateGlobalIDs[i] = (!useNodeBased) ? states[i]->getLocalID() : nodes[i]->getLocalID();
  }
  
  // create the mapping global IDs to global Petsc IDs
  getMethodData().getLocalToGlobalMapping().createStdSequential(stateGlobalIDs);
}

//////////////////////////////////////////////////////////////////////////////

void StdSeqSetup::setMatrix(const CFuint localSize,
                            const CFuint globalSize)
{
  CFAUTOTRACE;

  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();
  DataHandle < Framework::Node*, Framework::GLOBAL > nodes = socket_nodes.getDataHandle();
  const bool useNodeBased = getMethodData().useNodeBased();
  const CFuint nbStates = (!useNodeBased) ? states.size() : nodes.size();
  const CFuint nbEqs = getMethodData().getNbSysEquations();
  
  CFLog(VERBOSE, "StdSeqSetup::setMatrix() => useNodeBased[" << useNodeBased << "], nbEqs[" << nbEqs << "]\n");
  
  // all non zero entries (block matrix format is not supported on GPU)
  std::valarray<CFint> allNonZero(nbStates);
  allNonZero = 0;
  
  // ghost neighbors (out of diagonal submatrix non zero entries
  // dummy, since this is sequential)
  std::valarray<CFint> outDiagNonZero(nbStates);
  outDiagNonZero = 0;

  SelfRegistPtr<GlobalJacobianSparsity> sparsity =
    getMethodData().getCollaborator<SpaceMethod>()->createJacobianSparsity();
  
  sparsity->setDataSockets(socket_states, socket_nodes, socket_bStatesNeighbors);
  (!useNodeBased) ? sparsity->computeNNz(allNonZero, outDiagNonZero) :
    sparsity->computeNNzNodeBased(allNonZero, outDiagNonZero);    
  
  // allocating the matrix
  PetscMatrix& mat = getMethodData().getMatrix();
  
  const CFuint blockSize   = nbEqs;
  const CFuint nbRows      = nbStates*nbEqs;
  const CFuint nbCols      = nbRows;
  const CFuint nbNonZeroBlocks = 0;
  
  // create a sequential sparse matrix in block compressed row
  // format
  mat.setGPU(getMethodData().useGPU());
  mat.setAIJ(getMethodData().useAIJ());
  mat.createSeqBAIJ(blockSize,
		    nbRows,
		    nbCols,
		    nbNonZeroBlocks,
		    &allNonZero[0],
		    "Jacobian");
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace Petsc

} // namespace COOLFluiD
