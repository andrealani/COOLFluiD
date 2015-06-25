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
#include "Petsc/TwoLayerParSetup.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace Petsc {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<TwoLayerParSetup, PetscLSSData, PetscModule> twoLayerParSetupProvider("TwoLayerParSetup");

//////////////////////////////////////////////////////////////////////////////

TwoLayerParSetup::TwoLayerParSetup(const std::string& name) :
  BaseSetup(name)
{
}

//////////////////////////////////////////////////////////////////////////////

TwoLayerParSetup::~TwoLayerParSetup()
{
}

//////////////////////////////////////////////////////////////////////////////

void TwoLayerParSetup::setIdxMapping()
{
  CFAUTOTRACE;

  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();
  const CFuint nbStates = states.size();

  // build a contiguos mapping
  states.buildContiguosGlobal();

  std::valarray<CFuint> stateGlobalIDs(nbStates);
  std::valarray<bool> isGhost(nbStates);
  for (CFuint i = 0; i< nbStates; ++i) {
    stateGlobalIDs[i] = states.getContiguosID(states[i]->getLocalID());
    isGhost[i] = !states[i]->isParUpdatable();
  }

  // free contiguos mapping
  states.freeContiguosGlobal();

  getMethodData().getLocalToGlobalMapping().createMapping(stateGlobalIDs, isGhost);
}

//////////////////////////////////////////////////////////////////////////////

void TwoLayerParSetup::setMatrix(const CFuint localSize,
                            const CFuint globalSize)
{
  CFAUTOTRACE;

  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();

  const CFuint nbStates = states.size();
  const CFuint nbEqs = getMethodData().getNbSysEquations();

  // all non zero entries
  std::valarray<CFint> allNonZero(nbStates);
  allNonZero = 0;

  // ghost neighbors (out of diagonal submatrix non zero entries)
  std::valarray<CFint> outDiagNonZero(nbStates);
  outDiagNonZero = 0;

  SelfRegistPtr<GlobalJacobianSparsity> sparsity =
    getMethodData().getCollaborator<SpaceMethod>()->createJacobianSparsity();

  sparsity->setDataSockets(socket_states, socket_nodes, socket_bStatesNeighbors);
  sparsity->computeNNz(allNonZero, outDiagNonZero);

  // the allNonZeroUp and outDiagNonZeroUp have to be double sized,
  // with twice the same info...in the correct order!!!!

//   for (CFuint i=0; i < nbStates;++i){
//     totalAllNonZero[2*i] = allNonZero[i];
//     totalAllNonZero[(2*i)+1] = allNonZero[i];
//   }

  std::valarray<CFint> allNonZeroUp(localSize*2);
  std::valarray<CFint> outDiagNonZeroUp(localSize*2);
  CFuint countUp = 0;
  for (CFuint i = 0; i < nbStates; ++i) {
    if (states[i]->isParUpdatable()) {
      allNonZeroUp[2*countUp]         = allNonZero[i];
      allNonZeroUp[(2*countUp) + 1]   = allNonZero[i];
      outDiagNonZeroUp[2*countUp]     = outDiagNonZero[i];
      outDiagNonZeroUp[(2*countUp)+1] = outDiagNonZero[i];
      ++countUp;
    }
  }

  PetscMatrix& mat = getMethodData().getMatrix();

  const CFuint blockSize     = 2 * nbEqs;
  const CFuint localNbRows   = 2 * localSize * nbEqs;
  const CFuint localNbCols   = localNbRows;
  const CFuint globalNbRows  = 2 * globalSize * nbEqs;
  const CFuint globalNbCols  = globalNbRows;

  const string nsp = getMethodData().getNamespace();
  
  // create a sequential sparse matrix in block compressed row
  // format
  mat.createParBAIJ(PE::GetPE().GetCommunicator(nsp),
                    blockSize,
                    localNbRows,
                    localNbCols,
                    globalNbRows,
                    globalNbCols,
                    0, &allNonZeroUp[0],
                    0, &outDiagNonZeroUp[0],
                    "Jacobian");
}

//////////////////////////////////////////////////////////////////////////////

void TwoLayerParSetup::setVectors(const CFuint localSize,
				  const CFuint globalSize)
{
  const CFuint nbEqs = getMethodData().getNbSysEquations();
  const CFuint localVecSize  = 2 * localSize  * nbEqs;
  const CFuint globalVecSize = 2 * globalSize * nbEqs;
  const string nsp = getMethodData().getNamespace();
 
  PetscVector& sol = getMethodData().getSolVector();
  PetscVector& rhs = getMethodData().getRhsVector();
  
  // create the solution and the rhs vectors
  rhs.create(PE::GetPE().GetCommunicator(nsp), localVecSize, globalVecSize,"Rhs");
  sol.create(PE::GetPE().GetCommunicator(nsp), localVecSize, globalVecSize,"Solution");
  
  // initialize the two vectors
  sol.initialize(PE::GetPE().GetCommunicator(nsp), 1.0);
  rhs.initialize(PE::GetPE().GetCommunicator(nsp), 0.0);
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace Petsc

} // namespace COOLFluiD
