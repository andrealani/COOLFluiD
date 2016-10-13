// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Paralution/StdSetup.hh"
#include "Paralution/Paralution.hh"

#include "Framework/SpaceMethod.hh"
#include "Framework/GlobalJacobianSparsity.hh"
#include "Framework/MeshData.hh"
#include "Framework/State.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Common/PE.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace Paralution {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<StdSetup, ParalutionLSSData, ParalutionModule> 
stdSetupParalutionProvider("StdSetup");

//////////////////////////////////////////////////////////////////////////////

StdSetup::StdSetup(const std::string& name) :
  ParalutionLSSCom(name),
  socket_bStatesNeighbors("bStatesNeighbors"),
  socket_states("states"),
  socket_nodes("nodes"),
  socket_rhs("rhs")
{
}

//////////////////////////////////////////////////////////////////////////////

StdSetup::~StdSetup()
{
}

//////////////////////////////////////////////////////////////////////////////

void StdSetup::execute()
{
  CFAUTOTRACE;
  CFLog(NOTICE, "StdSetup::execute() \n");

 // paralution::init_paralution();
  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();
  DataHandle < Framework::Node*, Framework::GLOBAL > nodes = socket_nodes.getDataHandle();
  
  const bool useNodeBased = getMethodData().useNodeBased();
  const CFuint nbStates = (!useNodeBased) ? states.size() : nodes.size();
  
  // set the index mapping (global IDs to global Paralution IDs)
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
  // setVectors(nbUpdatableStates, totalNbStates);
  
  // set the matrix
   setMatrix(nbUpdatableStates, totalNbStates);

  // set the Krylov Sub Space solver
  // setKSP();
}

//////////////////////////////////////////////////////////////////////////////

void StdSetup::setKSP()
{
  CFAUTOTRACE;

  // Allocate and configure the LSS solver
}

//////////////////////////////////////////////////////////////////////////////

void StdSetup::setVectors(const CFuint localSize,
			   const CFuint globalSize)
{
  CFAUTOTRACE;
  CFLog(NOTICE, "StdSetup::setVectors() \n");
  // ParalutionVector& sol = getMethodData().getSolVector();
  // ParalutionVector& rhs = getMethodData().getRhsVector();

   const CFuint nbEqs = getMethodData().getNbSysEquations();
   const string nsp = getMethodData().getNamespace();
  
  // // create the solution and the rhs vectors
  // rhs.setGPU(getMethodData().useGPU());
  // rhs.create(PE::GetPE().GetCommunicator(nsp), localSize*nbEqs,
  //            globalSize*nbEqs,"rhs");
  
  // sol.setGPU(getMethodData().useGPU());
  // sol.create(PE::GetPE().GetCommunicator(nsp), localSize*nbEqs,
  //            globalSize*nbEqs,"Solution");
  
  // // initialize the two vectors
  // sol.initialize(PE::GetPE().GetCommunicator(nsp),1.0);
  // rhs.initialize(PE::GetPE().GetCommunicator(nsp),0.0);
}
      
//////////////////////////////////////////////////////////////////////////////

void StdSetup::setIdxMapping()
{
  CFAUTOTRACE;
  CFLog(NOTICE, "StdSetup::setIdxMapping() \n");
  //TODO ESTO ESTABA COMENTADO
   DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();
   DataHandle < Framework::Node*, Framework::GLOBAL > nodes = socket_nodes.getDataHandle();
   bool useNodeBased = getMethodData().useNodeBased();
   CFuint nbStates = (!useNodeBased) ? states.size() : nodes.size();
  
  vector<CFuint> stateGlobalIDs(nbStates);
  for (CFuint i = 0; i < nbStates; ++i) {
    stateGlobalIDs[i] = (!useNodeBased) ? states[i]->getLocalID() : nodes[i]->getLocalID();
  }
  
  // create the mapping global IDs to global Petsc IDs
  getMethodData().getLocalToGlobalMapping().createStdSequential(stateGlobalIDs);


//From ParJFSetup.cxx

//  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();
  nbStates = states.size();

  // build a contiguos mapping
  states.buildContiguosGlobal();

  std::valarray<CFuint> stateGlobalIDs2(nbStates);
  std::valarray<bool> isGhost(nbStates);
  for (CFuint i = 0; i< nbStates; ++i) {
    stateGlobalIDs2[i] = states.getContiguosID(states[i]->getLocalID());
    isGhost[i] = !states[i]->isParUpdatable();
  }

  // free contiguos mapping
  states.freeContiguosGlobal();

  getMethodData().getLocalToGlobalMapping().createMapping(stateGlobalIDs2, isGhost);

//  if (getMethodData().useBlockPreconditionerMatrix()) {
    std::valarray<CFuint> upLocalIDs(nbStates);
    CFuint nbUpdatables = 0;
    for (CFuint i = 0; i< nbStates; ++i) {
      if (states[i]->isParUpdatable()) {
        upLocalIDs[i] = nbUpdatables++;
      }
    }
    getMethodData().getLocalToLocallyUpdatableMapping().createMapping(upLocalIDs, isGhost);
//  }
// End


}

//////////////////////////////////////////////////////////////////////////////

void StdSetup::setMatrix(const CFuint localSize,
			 const CFuint globalSize)
{
  CFAUTOTRACE;

  CFAUTOTRACE;
  std::cout << "localsize: " << localSize << " globalsize: " << globalSize << "\n";
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
  ParalutionMatrix& mat = getMethodData().getMatrix();

for (CFint i=0; i<allNonZero.size(); i++){ if (allNonZero[i] != 0){ std::cout << i << " " << allNonZero[i] << "\n";}}
  std::cout << allNonZero.size() <<"\n"; 
 // const CFuint blockSize   = nbEqs;
 // const CFuint nbRows      = nbStates*nbEqs;
 // const CFuint nbCols      = nbRows;
 // const CFuint nbNonZeroBlocks = 0;
 // std::cout << blockSize << " " << nbRows << " " << nbCols << " " << nbNonZeroBlocks << "\n";
 
   mat.createCSR(allNonZero, nbEqs);

  // create a sequential sparse matrix in block compressed row
  // format
//  mat.setGPU(getMethodData().useGPU());
//  mat.setAIJ(getMethodData().useAIJ());
//  mat.createSeqBAIJ(blockSize,
//		    nbRows,
//		    nbCols,
//		    nbNonZeroBlocks,
//		    &allNonZero[0],
//		    "Jacobian");
}

//////////////////////////////////////////////////////////////////////////////

vector<SafePtr<BaseDataSocketSink> > StdSetup::needsSockets()
{
  vector<SafePtr<BaseDataSocketSink> > result;

  result.push_back(&socket_bStatesNeighbors);
  result.push_back(&socket_states);
  result.push_back(&socket_nodes);
  result.push_back(&socket_rhs);
  
  return result;
}

//////////////////////////////////////////////////////////////////////////////

} // namespace Paralution

} // namespace COOLFluiD
