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

  _useEisenstatWalker = true;
  setParameter("UseEisenstatWalker", &_useEisenstatWalker);
}

//////////////////////////////////////////////////////////////////////////////

ParJFSetup::~ParJFSetup()
{
}

//////////////////////////////////////////////////////////////////////////////

void ParJFSetup::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< CFreal >("Epsilon","DEPRECATED: Epsilon for numerical derivative (ignored with MatMFFD, which uses adaptive epsilon)");

  options.addConfigOption< bool >("JFApprox2ndOrder", "DEPRECATED: 2nd-order FD not supported with MatMFFD (ignored, logged as warning if true)");

  options.addConfigOption< bool >("UseEisenstatWalker", "Enable Eisenstat-Walker adaptive KSP tolerance (default true)");
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
  ctx->eps = _epsilon;  // kept for GMRESR variant and LUSGS/DPLUR preconditioners
  ctx->differentPreconditionerMatrix = getMethodData().getDifferentPreconditionerMatrix();

  // Eisenstat-Walker config
  ctx->useEisenstatWalker = _useEisenstatWalker;
  ctx->prevNonlinResNorm = -1.0;
  ctx->prevEta = 0.5;

  // DEPRECATED options — log warnings
  ctx->jfApprox2ndOrder = _jfApprox2ndOrder;
  if (_jfApprox2ndOrder) {
    CFLog(WARN, "ParJFSetup: JFApprox2ndOrder=true is DEPRECATED and ignored. "
      "MatMFFD uses adaptive 1st-order FD (DS formula) which is more robust.\n");
  }
  if (std::abs(_epsilon - 1e-6) > 1e-15) {
    CFLog(WARN, "ParJFSetup: User-specified Epsilon=" << _epsilon << " is DEPRECATED and ignored. "
      "MatMFFD computes epsilon adaptively via the DS formula.\n");
  }

  ctx->bkpStates.resize(nbStates*nbEqs);
  ctx->rhsVec = &getMethodData().getRhsVector();

  const std::string nsp = getMethodData().getNamespace();

  CFLog(VERBOSE, "ParJFSetup::setMatrix() => P" <<  Common::PE::GetPE().GetRank(nsp) << " matrix localSize = "
	<< localSize << ", globalSize = " << globalSize << "\n");

  // Create a parallel Matrix-Free Finite Difference (MFFD) matrix.
  // PETSc handles epsilon computation adaptively (DS formula by default).
  mat.createParMFFDMat(PE::GetPE().GetCommunicator(nsp),
                       localSize * nbEqs,
                       localSize * nbEqs,
                       globalSize * nbEqs,
                       globalSize * nbEqs,
                       "JFMatrix");

  // Allocate a PETSc Vec for packing current states (used by MatMFFDSetBase in execute())
  CF_CHKERRCONTINUE(VecDuplicate(getMethodData().getRhsVector().getVec(), &ctx->stateVec));

  CFLog(INFO, "ParJFSetup: Using PETSc MatMFFD with adaptive epsilon (DS formula)\n");
  if (_useEisenstatWalker) {
    CFLog(INFO, "ParJFSetup: Eisenstat-Walker adaptive KSP tolerance enabled "
      "(gamma=" << ctx->ewGamma << ", alpha=" << ctx->ewAlpha << ")\n");
  }

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
    std::valarray<CFint> outDiagNonZeroUp(localSize);
    CFuint countUp = 0;
    for (CFuint i = 0; i < nbStates; ++i) {
      if (states[i]->isParUpdatable()) {
        allNonZeroUp[countUp]     = allNonZero[i];
        outDiagNonZeroUp[countUp] = outDiagNonZero[i];
        ++countUp;
      }
    }

    PetscMatrix& precondMat = getMethodData().getPreconditionerMatrix();
    // create a parallel sparse matrix in block compressed row format
    if(!getMethodData().useBlockPreconditionerMatrix()) {
      CFLog(INFO, "ParJFSetup: Using different preconditioner with ParBAIJ matrix\n");
      precondMat.createParBAIJ(PE::GetPE().GetCommunicator(nsp),
                               nbEqs,
                               localSize*nbEqs,
                               localSize*nbEqs,
                               globalSize*nbEqs,
                               globalSize*nbEqs,
                               0, &allNonZeroUp[0],
                               0, &outDiagNonZeroUp[0],
                               "PreconditionerMatrix");
    }
    else {
      CFLog(INFO, "ParJFSetup: Using Block Preconditioner with SeqBAIJ matrix\n");
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
  const string nsp = getMethodData().getNamespace();
  
  // create the solution and the rhs vectors
  rhs.create(PE::GetPE().GetCommunicator(nsp), localSize*nbEqs,
             globalSize*nbEqs,"rhs");
  sol.create(PE::GetPE().GetCommunicator(nsp), localSize*nbEqs,
             globalSize*nbEqs,"Solution");
  
  CFLog(VERBOSE, "ParJFSetup::setVectors() => P" << Common::PE::GetPE().GetRank(nsp) << 
	" vector localSize = " << localSize << ", globalSize = " << globalSize << "\n");
  
  // initialize the two vectors
  sol.initialize(PE::GetPE().GetCommunicator(nsp), 1.0);
  rhs.initialize(PE::GetPE().GetCommunicator(nsp), 0.0);
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace Petsc

} // namespace COOLFluiD
