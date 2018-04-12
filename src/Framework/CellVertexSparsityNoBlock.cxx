// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include <iostream>

#include "Environment/ObjectProvider.hh"

#include "Common/BadValueException.hh"
#include "Common/NotImplementedException.hh"

#include "Framework/Framework.hh"
#include "Framework/MeshData.hh"
#include "Framework/CellVertexSparsityNoBlock.hh"
#include "Framework/Storage.hh"
#include "Framework/PartitionerPeriodicTools.hh"
#include "Framework/PhysicalModel.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<CellVertexSparsityNoBlock,
                            GlobalJacobianSparsity,
                            FrameworkLib>
aCellVertexSparsityNoBlockProvider("CellVertexNoBlock");

//////////////////////////////////////////////////////////////////////////////

CellVertexSparsityNoBlock::CellVertexSparsityNoBlock()
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

CellVertexSparsityNoBlock::~CellVertexSparsityNoBlock()
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

void CellVertexSparsityNoBlock::initPeriodicData(std::vector<int>& crossNodes)
{
  CFAUTOTRACE;

  // AL: this implementation may be inconsistent with AIJ format
  
  /* int ndim= (int)PhysicalModelStack::getActive()->getDim();

  // read periodic file and bail out silent if not present
  std::string name0("FILE_NOT_EXISTS"),name1("FILE_NOT_EXISTS");
  std::vector<int> idx0(0),idx1(0);
  std::vector<CFreal> coord0(0),coord1(0);
  PartitionerPeriodicTools::readPeriodicInfo(ndim,name0,idx0,coord0,name1,idx1,coord1);
  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();
  CFuint nstates=states.size();
  crossNodes.resize(nstates);
  crossNodes.assign(nstates,-1);
  if (name0=="FILE_NOT_EXISTS") return;

  // get the trs pair
  SafePtr<TopologicalRegionSet> trs0=MeshDataStack::getActive()->getTrs(name0);
  if (trs0.isNull()) throw BadValueException(FromHere(),"Trs0 with name '" + name0 + "' could not be found.");
  SafePtr<TopologicalRegionSet> trs1=MeshDataStack::getActive()->getTrs(name1);
  if (trs1.isNull()) throw BadValueException(FromHere(),"Trs1 with name '" + name1 + "' could not be found.");

  // replacing global indices with local, based on a sorted node lookup
  std::vector<PartitionerPeriodicTools::PeriodicInfoItem> p0 =
    PartitionerPeriodicTools::fillPeriodicInfoItemVector(ndim, idx0, coord0);
  std::vector<PartitionerPeriodicTools::PeriodicInfoItem> p1 = 
    PartitionerPeriodicTools::fillPeriodicInfoItemVector(ndim, idx1, coord1);
  std::vector<PartitionerPeriodicTools::PeriodicInfoItem> l0 =
    PartitionerPeriodicTools::fillPeriodicInfoItemVector(trs0,&states[0],false);
  std::vector<PartitionerPeriodicTools::PeriodicInfoItem> l1 =
    PartitionerPeriodicTools::fillPeriodicInfoItemVector(trs1,&states[0],false);
  std::vector<PartitionerPeriodicTools::PeriodicInfoItem> common0 = 
    PartitionerPeriodicTools::findCommonNodes(l0,p0);
  std::vector<PartitionerPeriodicTools::PeriodicInfoItem> common1 = 
    PartitionerPeriodicTools::findCommonNodes(l1,p1);

//  std::cout << "PERPERPER: common0  " << common0.size() << "\n" << std::flush;
//  std::cout << "PERPERPER: common1  " << common1.size() << "\n" << std::flush;

  // match the pairs vice versa
  std::vector<CFreal> delta(ndim);
  if (p1.size()!=0) // if p1 is empty, p0 is also empty
  {
    for (int i=0; i<ndim; ++i) delta[i]=p0[0].crd[i]-p1[0].crd[i];
    for (int i=0; i<(const int)common1.size(); ++i)
      for (int j=0; j<(const int)ndim; ++j)
        common1[i].crd[j]+=delta[j];
  }
  std::vector<PartitionerPeriodicTools::PeriodicInfoItem> cross0 = 
    PartitionerPeriodicTools::findCommonNodes(common1,common0);
  std::vector<PartitionerPeriodicTools::PeriodicInfoItem> cross1 = 
    PartitionerPeriodicTools::findCommonNodes(common0,common1);
  if (cross0.size()!=cross1.size()) throw BadValueException(FromHere(),"Not the same number of nodes.");
  
  //  std::cout << "PERPERPER: cross0  " << cross0.size() << "\n" << std::flush;
  //  std::cout << "PERPERPER: cross1  " << cross1.size() << "\n" << std::flush;
  
  // fill crossnodes
  // crossNodes: -1 if not periodic node, othervise its the local idx of its periodic pair
  for (int i=0; i<(const int)cross0.size(); ++i) {
    crossNodes[cross0[i].idx]=cross1[i].idx;
    crossNodes[cross1[i].idx]=cross0[i].idx;
    }*/
  
}

//////////////////////////////////////////////////////////////////////////////

void CellVertexSparsityNoBlock::computeNNz
(std::valarray<CFint>& nnzOut, std::valarray<CFint>& ghostNnzOut)
{
  CFAUTOTRACE;

  cf_assert(socket_states.isConnected());
  cf_assert(socket_bStatesNeighbors.isConnected());

  DataHandle < Framework::State*, Framework::GLOBAL > states =
    socket_states.getDataHandle();
  DataHandle<std::valarray< State* > > bStatesNeighbors =
    socket_bStatesNeighbors.getDataHandle();

  cf_assert(nnzOut.size() == ghostNnzOut.size());

  // count cells on the diagonal
  valarray<CFint> nnz(nnzOut.size());
  valarray<CFint> ghostNnz(nnzOut.size());
  
  std::vector<std::string> tags;
  tags.push_back ( "inner" );
  tags.push_back ( "cell" );

  std::vector< Common::SafePtr<TopologicalRegionSet> > innerCellsList = MeshDataStack::getActive()->getFilteredTrsList(tags);
  const CFuint nbStates = states.size();
  std::valarray< bool > isBoundaryState(false, nbStates);
  computeBoundaryStatesFlag(isBoundaryState);
  
  vector< vector< CFuint > > neighborCells(nbStates);
  vector< vector< Common::SafePtr<TopologicalRegionSet> > > neighborCellsTrs(nbStates);

  // detect all the neighbor Cells for each State
  const CFuint nbGroups = innerCellsList.size();
  for(CFuint iGroup=0; iGroup < nbGroups; ++iGroup){
    Common::SafePtr<TopologicalRegionSet> currTrs = innerCellsList[iGroup];
    SafePtr< ConnectivityTable< CFuint > > cellStates =
      MeshDataStack::getActive()->getConnectivity("cellStates_" + currTrs->getName());

    const CFuint nbCells = cellStates->nbRows();
    for (CFuint iCell = 0; iCell < nbCells; ++iCell) {
      const CFuint nbStatesInCell = cellStates->nbCols(iCell);
      for (CFuint iState = 0; iState < nbStatesInCell; ++iState) {
        const CFuint stateID = (*cellStates)(iCell,iState);
        cf_assert(stateID < nbStates);
        neighborCells[stateID].push_back(iCell);
        neighborCellsTrs[stateID].push_back(currTrs);
      }
    }
  }

  /*std::vector<int> crossNodes;
    initPeriodicData(crossNodes);*/
  
  typedef set< State*, less< State* > > SetOfNeighbors;

  for (CFuint iState = 0; iState < nbStates; ++iState)
  {
    const CFuint nbNeighborCells = neighborCells[iState].size();
    SetOfNeighbors* neighborStates = new SetOfNeighbors();
    SetOfNeighbors* ghostNeighStates = new SetOfNeighbors();

    for (CFuint iCell = 0; iCell < nbNeighborCells; ++iCell) {
      const CFuint cellTrsID = neighborCells[iState][iCell];
      Common::SafePtr<TopologicalRegionSet> cellTrs = neighborCellsTrs[iState][iCell];
      SafePtr< ConnectivityTable< CFuint > > cellStates =
        MeshDataStack::getActive()->getConnectivity("cellStates_" + cellTrs->getName());

      const CFuint nbStatesInCell = cellStates->nbCols(cellTrsID);
      for (CFuint jState = 0; jState < nbStatesInCell; ++jState) {
        // store all the neighbor states + the considered state itself
        State *const currState = states[(*cellStates)(cellTrsID,jState)];
        neighborStates->insert(SetOfNeighbors::value_type(currState));
        if (!currState->isParUpdatable()) {
          // if the state is not updatable store it in ghostNeighStates
          ghostNeighStates->insert(SetOfNeighbors::value_type(currState));
        }
      }

      // if there is periodicity, add the counter-pairs connectivity
      /*if (crossNodes[iState]!=-1) {
        const CFuint nbNeighborCells = neighborCells[crossNodes[iState]].size();
	
        for (CFuint iCell = 0; iCell < nbNeighborCells; ++iCell) {
          const CFuint cellTrsID = neighborCells[crossNodes[iState]][iCell];
          Common::SafePtr<TopologicalRegionSet> cellTrs = neighborCellsTrs[crossNodes[iState]][iCell];
          SafePtr< ConnectivityTable< CFuint > > cellStates =
            MeshDataStack::getActive()->getConnectivity("cellStates_" + cellTrs->getName());
	  
          const CFuint nbStatesInCell = cellStates->nbCols(cellTrsID);
          for (CFuint jState = 0; jState < nbStatesInCell; ++jState) {
            // store all the neighbor states + the considered state itself
            State *const currState = states[(*cellStates)(cellTrsID,jState)];
            neighborStates->insert(SetOfNeighbors::value_type(currState));
            if (!currState->isParUpdatable()) {
              // if the state is not updatable store it in ghostNeighStates
              ghostNeighStates->insert(SetOfNeighbors::value_type(currState));
            }
          }
	}
	}*/
    }
    
    // number of neighbor states == number of non zero entries in the matrix
    const CFuint nbNeighStates = neighborStates->size();
    nnz[iState] = nbNeighStates;
    ghostNnz[iState] = ghostNeighStates->size();

    // if the State is on the boundary store all the neighbor States
    if (isBoundaryState[iState]) {
      bStatesNeighbors[iState].resize(nbNeighStates);
      SetOfNeighbors::const_iterator itn;
      CFuint in = 0;
      for (itn = neighborStates->begin();
           itn != neighborStates->end();
           ++itn, ++in) {
        bStatesNeighbors[iState][in] = *itn;
      }
    }

    deletePtr( neighborStates );
    deletePtr( ghostNeighStates );
  }

  // check periodics
/*  for(CFuint i=0; i<nbStates; ++i) if (crossNodes[i]!=-1)
  {
//    std::cout << "nods: "<< i << " " << crossNodes[i] << " nnz: " << nnz[i] << " " << nnz[crossNodes[i]] << " ghostnnz: " << ghostNnz[i] << " " << ghostNnz[crossNodes[i]] << " bsn: " << bStatesNeighbors[i].size() << " " << bStatesNeighbors[crossNodes[i]].size() << " : ";
    if (bStatesNeighbors[i].size()!=bStatesNeighbors[crossNodes[i]].size()) throw BadValueException(FromHere(),"bStateNeighbors size does not match for periodic pairs.");
    for (CFuint j=0; j<bStatesNeighbors[i].size(); ++j)
      if (bStatesNeighbors[i][j]->getGlobalID()!=bStatesNeighbors[crossNodes[i]][j]->getGlobalID())
        throw BadValueException(FromHere(),"Mismatching node IDs in nieghbors of periodic pairs");
//    for (CFuint j=0; j<bStatesNeighbors[i].size(); ++j)
//      std::cout << " (" << bStatesNeighbors[i][j]->getLocalID() << " " << bStatesNeighbors[crossNodes[i]][j]->getLocalID() << ")";
//    std::cout << "\n" << std::flush;
}*/
//  std::cout << "NNZOUT SUM:      " << nnz.sum() << "\n" << std::flush;
//  std::cout << "NNZOUTGHOST SUM: " << ghostNnz.sum() << "\n" << std::flush;
//  sleep(2);

  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  const CFuint outSize = nnz.size()*nbEqs;
  nnzOut.resize(outSize);
  ghostNnzOut.resize(outSize);
  
  for (CFuint iState = 0; iState < nbStates; ++iState) {
    const CFuint startRow = iState*nbEqs;
    // CFLog(INFO, "CellVertexSparsityNoBlock::computeNNz() => (nnz, ghostNnz)["
    //<< iState << "] = (" << nnz[iState] << ", " <<ghostNnz[iState] << ")\n");
    
    for (CFuint iEq = 0; iEq < nbEqs; ++iEq) {
      nnzOut[startRow+iEq] = nnz[iState]*nbEqs;
      ghostNnzOut[startRow+iEq] = ghostNnz[iState]*nbEqs;
    }
  }  
}

//////////////////////////////////////////////////////////////////////////////

void CellVertexSparsityNoBlock::computeMatrixPattern( std::valarray<CFint>& nnz,
						      std::valarray<CFint>& ghostNnz,
						      vector< vector<CFuint> >& matrixPattern)
{
  CFAUTOTRACE;

  throw NotImplementedException
    (FromHere(), "CellVertexSparsityNoBlock::computeMatrixPattern()");
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
