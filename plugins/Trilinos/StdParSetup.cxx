// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Common/PE.hh"
#include "Framework/MeshData.hh"
#include "Framework/GlobalJacobianSparsity.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/LSSIdxMapping.hh"
#include "Framework/State.hh"
#include "Framework/SpaceMethod.hh"
#include "Framework/MethodData.hh"

#include "Trilinos/Trilinos.hh"
#include "Trilinos/TrilinosIdxMapping.hh"
#include "Trilinos/StdParSetup.hh"

#include "Epetra_MpiComm.h"
#include <fstream>

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace Trilinos {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<StdParSetup, TrilinosLSSData, TrilinosModule> stdParSetupProvider("StdParSetup");

//////////////////////////////////////////////////////////////////////////////

StdParSetup::StdParSetup(const std::string& name) :
  TrilinosLSSCom(name),
  socket_bStatesNeighbors("bStatesNeighbors"),
  socket_states("states"),
  socket_nodes("nodes")
{
}

//////////////////////////////////////////////////////////////////////////////

StdParSetup::~StdParSetup()
{
}

//////////////////////////////////////////////////////////////////////////////

void StdParSetup::execute()
{
  CFTRACEBEGIN;

  CFAUTOTRACE;

  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();

// unused //  const CFuint nbStates = states.size();

  // set the Epetra communicator
  setEpetraComm();

  // set the index mapping (global IDs to global Trilinos IDs)
  setMapping();

  // set the vectors
  setVectors();

  // set the matrix
  setMatrix();

  CFTRACEEND;
}

//////////////////////////////////////////////////////////////////////////////

void StdParSetup::setEpetraComm()
{
#ifndef CF_HAVE_MPI
  CFerr << "ABORT StdParSetup: recompile coolfluid and/or trilinos with mpi-compilers" << "\n";
  exit(1):
#endif // CF_HAVE_MPI

  Epetra_MpiComm *comm = new Epetra_MpiComm(
    Common::PE::GetPE().GetCommunicator() );    // deletion --> in UnSetup
  getMethodData().setEpetraComm(comm);
}

//////////////////////////////////////////////////////////////////////////////

void StdParSetup::setMapping()
{

  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();

  const CFuint nbStates = states.size(); // states 'seen' by this processor (owned + ghost)
  const CFuint nbEqs    = getMethodData().getNbSysEquations();

  vector<int> stateGlobalIDs(nbStates);
  vector<bool> isNonLocalRow(nbStates);

  for (CFuint localID = 0; localID < nbStates; ++localID) {
    stateGlobalIDs[localID] = states[localID]->getGlobalID();
    isNonLocalRow[localID] = !states[localID]->isParUpdatable();
  }

  // TrilinosIdxMapping
  TrilinosIdxMapping idxMap(stateGlobalIDs, isNonLocalRow);

  // configuration of the mapping:
  //   coolfluid local state ID --> trilinos global state ID
  getMethodData().getLocalToGlobalMapping().createMapping(idxMap);

  // determining the global IDs of the locally OWNED unknowns
  vector<int> myGlobalUnknowns(nbStates * nbEqs); // locally owned unknowns (NOT the states, and EXcluding the ghost unknowns), worst case allocation
  int nbMyUnknowns = 0;
  for (CFuint s=0; s<nbStates; s++) {
    if (states[s]->isParUpdatable()) {
      for (CFuint e=0; e<nbEqs; e++) {
        myGlobalUnknowns[nbMyUnknowns] = nbEqs*stateGlobalIDs[s] + e;
        nbMyUnknowns++;
      }
    }
  }

  const std::string nsp = getMethodData().getNamespace();
  MPI_Comm comm = PE::GetPE().GetCommunicator(nsp);
  
// FOR TESTING PURPOSES
//  char filename[10];
//  int rank;
//  MPI_Comm_rank(comm, &rank);
//  sprintf(&filename[0], "stateID%d", rank);
//  ofstream file(filename);
//  file << "nbEqs = " << nbEqs << "\n";
//  for (CFuint s=0; s<nbStates; s++) {
//    file << states[s]->getGlobalID() << " : ";
//    if (states[s]->isParUpdatable()) file << "local" << "\n";
//    else file << "ghost" << "\n";
//  }
//  file.close();

// determining the total number of unknowns --> not needed: calculated by Epetra_Map constructor (nbGlobal = -1)
//  const CFuint totalNbStates = MeshDataStack::getActive()->getDataStorage()->
//    getData<CFreal>("statedata").getGlobalSize();
//  int nbGlobalUnknowns = totalNbStates * nbEqs;

  // construction of the Epetra map --> deletion in StdParUnsetup
  Epetra_Comm *comm = getMethodData().getEpetraComm();
  Epetra_Map *map = new Epetra_Map(-1, nbMyUnknowns, &myGlobalUnknowns[0], 0, *comm);
  getMethodData().setEpetraMap(map);

/*  CFLog(VERBOSE, "***********StdParSetup::setMapping() begin : nbProcessors = " << Common::PE::GetPE().GetProcessorCount() << "\n" );
  //CFLog(VERBOSE, "***********StdParSetup::setMapping()\n");

  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();

  const CFuint nbStates = states.size();
  const CFuint nbEqs    = getMethodData().getNbSysEquations();

  DataHandle<CFreal, GLOBAL> statedata = MeshDataStack::getActive()->getDataStorage()->
    getData<CFreal>("statedata");

  /// @todo remove this
MPI_Barrier(comm);
  CFLog(VERBOSE, "***********StdParSetup::setMapping() 1 : nbProcessors = " << Common::PE::GetPE().GetProcessorCount() << "\n" );
  vector<CFuint> nbUpdatablesPerRank = statedata.getNbUpdatablesPerRank();

MPI_Barrier(comm);
  CFLog(VERBOSE, "***********StdParSetup::setMapping() 2 : nbProcessors = " << Common::PE::GetPE().GetProcessorCount() << "\n" );
  const vector< vector<CFuint> >& updatablesPerRank =
      statedata.setAndGetUpdatablesListPerRank();

  // store all the state global IDs
MPI_Barrier(comm);
  CFLog(VERBOSE, "***********StdParSetup::setMapping() 3 : nbProcessors = " << Common::PE::GetPE().GetProcessorCount() << "\n" );
  vector<CFuint> stateGlobalIDs(nbStates);
  vector<bool> isGhostVec(nbStates);
  for (CFuint i = 0; i< nbStates; ++i) {
    stateGlobalIDs[i] = states[i]->getGlobalID();
    isGhostVec[i] = !states[i]->isParUpdatable();
  }

MPI_Barrier(comm);
  CFLog(VERBOSE, "***********StdParSetup::setMapping() 4 : nbProcessors = " << Common::PE::GetPE().GetProcessorCount() << "\n" );
  TrilinosIdxMapping trilinosIdxMapping(updatablesPerRank,
                                  stateGlobalIDs,
                                  isGhostVec);

  // create the mapping global IDs to global Trilinos IDs
MPI_Barrier(comm);
  CFLog(VERBOSE, "***********StdParSetup::setMapping() 5 : nbProcessors = " << Common::PE::GetPE().GetProcessorCount() << "\n" );
  LSSIdxMapping::getInstance().createParallel(trilinosIdxMapping);

MPI_Barrier(comm);
  CFLog(VERBOSE, "***********StdParSetup::setMapping() 6 : nbProcessors = " << Common::PE::GetPE().GetProcessorCount() << "\n" );
  int nbMyUnknowns = 0;
  int *myGlobalUnknowns = new int[nbEqs * nbStates];
  for (CFuint i = 0; i< nbStates; ++i) {
    if (isGhostVec[i]) {
      int trilinosID = LSSIdxMapping::getInstance().getID(i);
      for (CFuint e=0; e<nbEqs; e++) {
        myGlobalUnknowns[nbMyUnknowns] = trilinosID*nbEqs + e;
        nbMyUnknowns++;
      }
    }
  }

  // construction of the Epetra map --> deletion in StdParUnsetup
  CFLog(VERBOSE, "***********StdParSetup::setMapping() 7 : nbProcessors = " << Common::PE::GetPE().GetProcessorCount() << "\n" );
  Epetra_Comm *comm = getMethodData().getEpetraComm();
  Epetra_Map *map = new Epetra_Map(-1, nbMyUnknowns, myGlobalUnknowns, 0, *comm);
  getMethodData().setEpetraMap(map);

  delete[] myGlobalUnknowns;

  CFLog(VERBOSE, "***********StdParSetup::setMapping() end : nbProcessors = " << Common::PE::GetPE().GetProcessorCount() << "\n" );
*/
}

//////////////////////////////////////////////////////////////////////////////

void StdParSetup::setMatrix()
{
  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();

  const CFuint nbStates = states.size(); // locally owned + ghost states
  const CFuint nbEqs = getMethodData().getNbSysEquations();

  // all non zero entries
  std::valarray<CFint> allNonZero(nbStates);
  allNonZero = 0;

  // ghost neighbors (out of diagonal submatrix non zero entries)
  std::valarray<CFint> outDiagNonZero(nbStates);
  outDiagNonZero = 0;

  // dimensions
  Epetra_Map *map = getMethodData().getEpetraMap();
  int nbMyUnknowns = map->NumMyElements();
  int nbGlobalUnknowns = map->NumGlobalElements();

  SafePtr<TopologicalRegionSet> elements =
    MeshDataStack::getActive()->getTrs("InnerCells");

  SelfRegistPtr<GlobalJacobianSparsity> sparsity =
    getMethodData().getCollaborator<SpaceMethod>()->createJacobianSparsity();

  sparsity->setDataSockets(socket_states,socket_nodes,  socket_bStatesNeighbors);
  sparsity->computeNNz(allNonZero, outDiagNonZero);

  int *dnnz = new int[nbStates]; // BLOCK diagonal entries !! worst case size
  int *onnz = new int[nbStates]; // BLOCK diagonal entries !! worst case size
  int nbMyStates = 0;
  for (CFuint i = 0; i < nbStates; ++i) {
    if (states[i]->isParUpdatable()) {
      dnnz[nbMyStates] = allNonZero[i];
      onnz[nbMyStates] = outDiagNonZero[i];
      ++nbMyStates;
    }
  }
  //CFout << "nbMyStates = " << nbMyStates << "\n";
  //CFout << "nbEqs = " << nbEqs << "\n";
  //CFout << "nbMyUnknowns = " << nbMyUnknowns << "\n";
  cf_assert((int)(nbMyStates*nbEqs) == (int)(nbMyUnknowns));
  
  const std::string nsp = getMethodData().getNamespace();
  MPI_Comm comm = PE::GetPE().GetCommunicator(nsp);
  
  // allocation of the matrix
  TrilinosMatrix* mat = getMethodData().getMatrix();
  mat->setEpetraMap(map);
  mat->createParBAIJ(comm,
		     nbEqs,
		     nbMyUnknowns,
		     nbMyUnknowns,
		     nbGlobalUnknowns,
		     nbGlobalUnknowns,
		     0, dnnz,
		     0, onnz,
		     "Jacobian");
  
  // deallocation
  delete[] dnnz;
  delete[] onnz;
}

//////////////////////////////////////////////////////////////////////////////

void StdParSetup::setVectors()
{
  // dimensions
  Epetra_Map *map = getMethodData().getEpetraMap();
  int nbMyUnknowns = map->NumMyElements();
  int nbGlobalUnknowns = map->NumGlobalElements();
  
  const std::string nsp = getMethodData().getNamespace();
  MPI_Comm comm = PE::GetPE().GetCommunicator(nsp);
  
  // solution vector
  TrilinosVector* sol = getMethodData().getSolVector();
  sol->setEpetraMap(map);
  sol->create(comm, nbMyUnknowns, nbGlobalUnknowns, "solution");
  sol->initialize(comm, 0.0);

  // rhs vector
  TrilinosVector* rhs = getMethodData().getRhsVector();
  rhs->setEpetraMap(map);
  rhs->create(comm, nbMyUnknowns, nbGlobalUnknowns, "rhs");
  rhs->initialize(comm, 0.0);
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
StdParSetup::needsSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSink> > result;

  result.push_back(&socket_bStatesNeighbors);
  result.push_back(&socket_states);
  result.push_back(&socket_nodes);
  
  return result;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace Trilinos

  } // namespace Numerics

} // namespace COOLFluiD
