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
#include "Framework/SpaceMethodData.hh"
#include "Framework/MethodData.hh"

#include "Petsc/PetscIdxMapping.hh"
#include "Petsc/Petsc.hh"
#include "Petsc/ParJFSetup.hh"

#include "Common/PE.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace Petsc {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<ParJFSetup, PetscLSSData, PetscModule> 
parJFSetupProvider("ParJFSetup");

//////////////////////////////////////////////////////////////////////////////

ParJFSetup::ParJFSetup(const string& name) :
BaseSetup(name)
{
  addConfigOptionsTo(this);

  _epsilon = 1e-6;
  setParameter("Epsilon",&_epsilon);

  _jfApprox2ndOrder = false;
  setParameter("JFApprox2ndOrder", &_jfApprox2ndOrder);
}

//////////////////////////////////////////////////////////////////////////////

ParJFSetup::~ParJFSetup()
{
}

//////////////////////////////////////////////////////////////////////////////

void ParJFSetup::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< CFreal >("Epsilon","Epsilon for computing numerical derivative");

  options.addConfigOption< bool >("JFApprox2ndOrder", "2nd order of the Jacobian-free matrix vector product approximation (options: true/false)");
}

//////////////////////////////////////////////////////////////////////////////

void ParJFSetup::setIdxMapping()
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

  if (getMethodData().useBlockPreconditionerMatrix()) {
    std::valarray<CFuint> upLocalIDs(nbStates);
    CFuint nbUpdatables = 0;
    for (CFuint i = 0; i< nbStates; ++i) {
      if (states[i]->isParUpdatable()) {
        upLocalIDs[i] = nbUpdatables++;
      }
    }
    getMethodData().getLocalToLocallyUpdatableMapping().createMapping(upLocalIDs, isGhost);
  }
}

//////////////////////////////////////////////////////////////////////////////

void ParJFSetup::setMatrix(const CFuint localSize,
                           const CFuint globalSize)
{
  CFAUTOTRACE;

  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();
  DataHandle < CFreal > rhs  = socket_rhs.getDataHandle();

  const CFuint nbStates = states.size();
  const CFuint nbEqs = getMethodData().getNbSysEquations();
  PetscMatrix& mat = getMethodData().getMatrix();

  // set the data in the context 
  JFContext* ctx = getMethodData().getJFContext();

  ctx->petscData = &getMethodData();
  ctx->states = &socket_states;
  ctx->rhs = &socket_rhs;
  ctx->spaceMethod = getMethodData().getCollaborator<SpaceMethod>();
  ctx->eps = _epsilon;
  ctx->jfApprox2ndOrder = _jfApprox2ndOrder;
  ctx->differentPreconditionerMatrix = getMethodData().getDifferentPreconditionerMatrix();

  ctx->bkpStates.resize(nbStates*nbEqs);
  ctx->rhsVec = &getMethodData().getRhsVector();

  cout << Common::PE::GetPE().GetRank() << " matrix localSize = " << localSize << ", globalSize = " << globalSize << endl;

  // create a parallel Jacobian-Free matrix
  mat.createParJFMat(PETSC_COMM_WORLD,
                     localSize *nbEqs,
                     localSize *nbEqs,
                     globalSize *nbEqs,
                     globalSize *nbEqs,
                     (void*)(ctx),
                     "JFMatrix");

  // if one would like to use matrix preconditioner
  if(ctx->differentPreconditionerMatrix) {
    getMethodData().getCollaborator<SpaceMethod>()->getSpaceMethodData()->setFillPreconditionerMatrix(true);
    // all non zero entries
    std::valarray<CFint> allNonZero(nbStates);
    allNonZero = 0;

    // ghost neighbors (out of diagonal submatrix non zero entries)
    std::valarray<CFint> outDiagNonZero(nbStates);
    outDiagNonZero = 0;

    SelfRegistPtr<GlobalJacobianSparsity> sparsity = getMethodData().getCollaborator<SpaceMethod>()->createJacobianSparsity();
    
    sparsity->setDataSockets(socket_states, socket_nodes, socket_bStatesNeighbors);
    sparsity->computeNNz(allNonZero, outDiagNonZero);

    std::valarray<CFint> allNonZeroUp(localSize);
    CFuint countUp = 0;
    for (CFuint i = 0; i < nbStates; ++i) {
      if (states[i]->isParUpdatable()) {
        allNonZeroUp[countUp]     = allNonZero[i];
        ++countUp;
      }
    }
    //cout << "Number of states = " << nbStates << endl;
    PetscMatrix& precondMat = getMethodData().getPreconditionerMatrix();
    // create a parallel sparse matrix in block compressed row format
    if(!getMethodData().useBlockPreconditionerMatrix()) {
      cout << "Using different preconditioner with ParBAIJ matrix" << endl;
      precondMat.createParBAIJ(PETSC_COMM_WORLD,
                               nbEqs,
                               localSize*nbEqs,
                               localSize*nbEqs,
                               globalSize*nbEqs,
                               globalSize*nbEqs,
                               0, &allNonZeroUp[0],
                               0, CFNULL,
                               "PreconditionerMatrix");
    }
    else {
      cout << "Using Block Preconditioner with SeqBAIJ matrix" << endl;
      precondMat.createSeqBAIJ(nbEqs,
                               localSize*nbEqs,
                               localSize*nbEqs,
                               0, &allNonZeroUp[0],
                               "PreconditionerMatrix");
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void ParJFSetup::setVectors(const CFuint localSize, 
                            const CFuint globalSize)
{
  CFAUTOTRACE;

  PetscVector& sol = getMethodData().getSolVector();
  PetscVector& rhs = getMethodData().getRhsVector();

  const CFuint nbEqs = getMethodData().getNbSysEquations();

  // create the solution and the rhs vectors
  rhs.create(PETSC_COMM_WORLD, localSize*nbEqs,
             globalSize*nbEqs,"rhs");
  sol.create(PETSC_COMM_WORLD, localSize*nbEqs,
             globalSize*nbEqs,"Solution");

  cout << Common::PE::GetPE().GetRank() << " vector localSize = " << localSize << ", globalSize = " << globalSize << endl;

  // initialize the two vectors
  sol.initialize(PETSC_COMM_WORLD, 1.0);
  rhs.initialize(PETSC_COMM_WORLD, 0.0);
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace Petsc

} // namespace COOLFluiD
