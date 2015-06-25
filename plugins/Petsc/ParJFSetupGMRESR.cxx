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
#include "Petsc/ParJFSetupGMRESR.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

namespace Petsc {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<ParJFSetupGMRESR, PetscLSSData, PetscModule> 
	parJFSetupGMRESRProvider("ParJFSetupGMRESR");

//////////////////////////////////////////////////////////////////////////////

ParJFSetupGMRESR::ParJFSetupGMRESR(const string& name) :
	BaseSetup(name)
{
	addConfigOptionsTo(this);
	
	_epsilon = 1e-6;
	setParameter("Epsilon",&_epsilon);
	
	//_gmresRestart = 10;
	//setParameter("GMRESRestart", &_gmresRestart);
	
	//_gmresRelTol = PETSC_DEFAULT;
	//setParameter("GMRESRelTol", &_gmresRelTol);
	
	//_gmresAbsTol = PETSC_DEFAULT;
	//setParameter("GMRESAbsTol", &_gmresAbsTol);
	
	//_gmresDivTol = PETSC_DEFAULT;
	//setParameter("GMRESDivTol", &_gmresDivTol);
	
	//_gmresPreconditioner = PCNONE;
	//setParameter("GMRESPreconditioner", &_gmresPreconditioner);
	
	_gmresrMaxIter = 20;
	setParameter("GMRESRMaxIter", &_gmresrMaxIter);
	
	_gmresrTol = 1e-5;
	setParameter("GMRESRTol", &_gmresrTol);
	
	_gmresrSubspaceDim = 0;
	setParameter("GMRESRSubspaceDim", &_gmresrSubspaceDim);
}

//////////////////////////////////////////////////////////////////////////////

ParJFSetupGMRESR::~ParJFSetupGMRESR()
{
}

//////////////////////////////////////////////////////////////////////////////

void ParJFSetupGMRESR::defineConfigOptions(Config::OptionList& options)
{
	options.addConfigOption< CFreal >("Epsilon","Epsilon for computing numerical derivative");
	//options.addConfigOption< CFint >("GMRESrestart","restart of GMRES method inside GMRESR solver");
	//options.addConfigOption< CFint >("GMRESRelTol","Relative tolerance of GMRES method inside GMRESR solver");
	//options.addConfigOption< CFint >("GMRESAbsTol","Absolute tolerance of GMRES method inside GMRESR solver");
	//options.addConfigOption< CFint >("GMRESDivTol","Divergence tolerance of GMRES method inside GMRESR solver");
	//options.addConfigOption< PCType >("GMRESPreconditioner","Type of preconditioner for GMRES method inside GMRESR solver");
	options.addConfigOption< CFint >("GMRESRMaxIter","Max number of GMRESR iterations (if it does not converge)");
	options.addConfigOption< CFreal >("GMRESRTol","Tolerance for GMRESR");
	options.addConfigOption< CFint >("GMRESRSubspaceDim","Maximum dimension of GMRESR solver subspace (0 = dimension of system)");
}

//////////////////////////////////////////////////////////////////////////////

void ParJFSetupGMRESR::setIdxMapping()
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

void ParJFSetupGMRESR::setMatrix(const CFuint localSize,
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
  
  ctx->states = &socket_states;
  ctx->rhs = &socket_rhs;
  ctx->spaceMethod = getMethodData().getCollaborator<SpaceMethod>();
  ctx->eps = _epsilon;
  ctx->bkpStates.resize(nbStates*nbEqs);
  ctx->rhsVec = &getMethodData().getRhsVector();
  
  const string nsp = getMethodData().getNamespace();
  // create a parallel Jacobian-Free matrix
  mat.createParJFMat(PE::GetPE().GetCommunicator(nsp),
		     localSize *nbEqs,
		     localSize *nbEqs,
		     globalSize *nbEqs,
		     globalSize *nbEqs,
		     (void*)(ctx),
		     "JFMatrix");
}

//////////////////////////////////////////////////////////////////////////////

void ParJFSetupGMRESR::setVectors(const CFuint localSize,
				  const CFuint globalSize)
{
  CFAUTOTRACE;
  
  PetscVector& sol = getMethodData().getSolVector();
  PetscVector& rhs = getMethodData().getRhsVector();
  
  const CFuint nbEqs = getMethodData().getNbSysEquations();
  const string nsp = getMethodData().getNamespace();
  
  // create the solution and the rhs vectors
  rhs.create(PE::GetPE().GetCommunicator(nsp), localSize*nbEqs,
	     globalSize*nbEqs,"rhs");
  sol.create(PE::GetPE().GetCommunicator(nsp), localSize*nbEqs,
	     globalSize*nbEqs,"Solution");
  
  // initialize the two vectors
  sol.initialize(PE::GetPE().GetCommunicator(nsp),1.0);
  rhs.initialize(PE::GetPE().GetCommunicator(nsp),0.0);
}

//////////////////////////////////////////////////////////////////////////////

} // namespace Petsc

} // namespace COOLFluiD

