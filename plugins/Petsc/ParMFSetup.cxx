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
#include "Petsc/ParMFSetup.hh"

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

MethodCommandProvider<ParMFSetup, PetscLSSData, PetscModule> 
parMFSetupProvider("ParMFSetup");

//////////////////////////////////////////////////////////////////////////////

ParMFSetup::ParMFSetup(const string& name) :
BaseSetup(name)
{
  addConfigOptionsTo(this);
}

//////////////////////////////////////////////////////////////////////////////

ParMFSetup::~ParMFSetup()
{
}

//////////////////////////////////////////////////////////////////////////////

void ParMFSetup::defineConfigOptions(Config::OptionList& options)
{
}

//////////////////////////////////////////////////////////////////////////////

void ParMFSetup::setIdxMapping()
{
  CFAUTOTRACE;

  DataHandle < Framework::State*, Framework::GLOBAL > states = 
    socket_states.getDataHandle();
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

void ParMFSetup::setMatrix(const CFuint localSize,
                           const CFuint globalSize)
{
  CFAUTOTRACE;

  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();
  DataHandle < CFreal > rhs  = socket_rhs.getDataHandle();

  const CFuint nbStates = states.size();
  const CFuint nbEqs = getMethodData().getNbSysEquations();
  PetscMatrix& jfMat = getMethodData().getJFMatrix(); // Jacobian-free matrix
  PetscMatrix& precMat = getMethodData().getMatrix(); // preconditioner matrix

  // set the data in the context 
  MFContext* ctx = getMethodData().getMFContext();
  
  ctx->petscData = &getMethodData();
  
  const std::string nsp = getMethodData().getNamespace();
  CFLog(VERBOSE, "ParMFSetup::setMatrix() => P" << Common::PE::GetPE().GetRank(nsp) << 
	" matrix localSize = " << localSize << ", globalSize = " << globalSize << "\n");
  
  // creates a parallel Jacobian-Free matrix
  jfMat.createParJFMat(PE::GetPE().GetCommunicator(nsp),
		       localSize *nbEqs,
		       localSize *nbEqs,
		       globalSize *nbEqs,
		       globalSize *nbEqs,
		       (void*)(ctx),
		       "JFMatrix");
  
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

  std::valarray<CFint> allNonZeroUp(localSize);
  std::valarray<CFint> outDiagNonZeroUp(localSize);
  CFuint countUp = 0;
  for (CFuint i = 0; i < nbStates; ++i) {
    if (states[i]->isParUpdatable()) {
      allNonZeroUp[countUp]     = allNonZero[i];
      outDiagNonZeroUp[countUp] = outDiagNonZero[i];
      ++countUp;
    }
  }

  // creates a parallel preconditioner sparse matrix in block compressed row format
  precMat.createParBAIJ(PE::GetPE().GetCommunicator(nsp),
			nbEqs,
			localSize*nbEqs,
			localSize*nbEqs,
			globalSize*nbEqs,
			globalSize*nbEqs,
			0, &allNonZeroUp[0],
			0, &outDiagNonZeroUp[0],
			"PreconditionerMatrix");
}

//////////////////////////////////////////////////////////////////////////////

void ParMFSetup::setVectors(const CFuint localSize,
                            const CFuint globalSize)
{
  CFAUTOTRACE;

  PetscVector& sol = getMethodData().getSolVector();
  PetscVector& rhs = getMethodData().getRhsVector();

  const CFuint nbEqs = getMethodData().getNbSysEquations();
  const std::string nsp = getMethodData().getNamespace();
  
  // create the solution and the rhs vectors
  rhs.create(PE::GetPE().GetCommunicator(nsp), localSize*nbEqs, globalSize*nbEqs,"rhs");
  sol.create(PE::GetPE().GetCommunicator(nsp), localSize*nbEqs, globalSize*nbEqs,"Solution");
  
  // initialize the two vectors
  sol.initialize(PE::GetPE().GetCommunicator(nsp),1.0);
  rhs.initialize(PE::GetPE().GetCommunicator(nsp),0.0);
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace Petsc

} // namespace COOLFluiD
