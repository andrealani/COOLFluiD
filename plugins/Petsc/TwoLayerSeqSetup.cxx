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
#include "Petsc/TwoLayerSeqSetup.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace Petsc {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<TwoLayerSeqSetup, PetscLSSData, PetscModule> twoLayerSeqSetupProvider("TwoLayerSeqSetup");

//////////////////////////////////////////////////////////////////////////////

TwoLayerSeqSetup::TwoLayerSeqSetup(const std::string& name) :
  BaseSetup(name)
{
}

//////////////////////////////////////////////////////////////////////////////

TwoLayerSeqSetup::~TwoLayerSeqSetup()
{
}

//////////////////////////////////////////////////////////////////////////////

void TwoLayerSeqSetup::setIdxMapping()
{
  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();

  const CFuint nbStates = states.size();
  vector<CFuint> stateGlobalIDs(nbStates);
  for (CFuint i = 0; i < nbStates; ++i) {
    stateGlobalIDs[i] = states[i]->getLocalID();
  }

  // create the mapping global IDs to global Petsc IDs
  getMethodData().getLocalToGlobalMapping().createStdSequential(stateGlobalIDs);
}

//////////////////////////////////////////////////////////////////////////////

void TwoLayerSeqSetup::setMatrix(const CFuint localSize,
                            const CFuint globalSize)
{
  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();

  const CFuint nbStates = states.size();
  const CFuint nbEqs = getMethodData().getNbSysEquations();

  // all non zero entries
  std::valarray<CFint> allNonZero(nbStates);
  std::valarray<CFint> totalAllNonZero(2*nbStates);
  allNonZero = 0;
  totalAllNonZero = 0;

  // ghost neighbors (out of diagonal submatrix non zero entries)
  std::valarray<CFint> outDiagNonZero(nbStates);
  outDiagNonZero = 0;

  SelfRegistPtr<GlobalJacobianSparsity> sparsity =
    getMethodData().getCollaborator<SpaceMethod>()->createJacobianSparsity();

  sparsity->setDataSockets(socket_states, socket_nodes, socket_bStatesNeighbors);
  sparsity->computeNNz(allNonZero, outDiagNonZero);

  // the totalAllNonZero has to be double sized
  // with twice the same info ... in the correct order!
  for (CFuint i=0; i < nbStates;++i){
    totalAllNonZero[2*i] = allNonZero[i];
    totalAllNonZero[(2*i)+1] = allNonZero[i];
  }

  PetscMatrix& mat = getMethodData().getMatrix();

  const CFuint blockSize   = 2*nbEqs;
  const CFuint nbRows      = 2*nbStates*nbEqs;
  const CFuint nbCols      = nbRows;
  const CFuint nbNonZeroBlocks = 0;

  // create a sequential sparse matrix in block compressed row format
  mat.createSeqBAIJ(blockSize,
                    nbRows,
                    nbCols,
                    nbNonZeroBlocks,
                    &totalAllNonZero[0],
                    "Jacobian");
}

//////////////////////////////////////////////////////////////////////////////

void TwoLayerSeqSetup::setVectors(const CFuint localSize,
                             const CFuint globalSize)
{
  cf_assert(localSize == globalSize);
  const CFuint vecSize  = 2*localSize*getMethodData().getNbSysEquations();

  PetscVector& sol = getMethodData().getSolVector();
  PetscVector& rhs = getMethodData().getRhsVector();

  // create the solution and the rhs vectors
  sol.create(PETSC_COMM_WORLD, vecSize, vecSize, "Solution");
  rhs.create(PETSC_COMM_WORLD, vecSize, vecSize, "rhs");

  // initialize the two vectors
  sol.initialize(PETSC_COMM_WORLD, 1.0);
  rhs.initialize(PETSC_COMM_WORLD, 0.0);
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace Petsc

} // namespace COOLFluiD
