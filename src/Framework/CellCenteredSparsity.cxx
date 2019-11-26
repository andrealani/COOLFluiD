// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include <fstream>

#include "Common/NotImplementedException.hh"

#include "Environment/ObjectProvider.hh"

#include "Framework/CellCenteredSparsity.hh"
#include "Framework/MeshData.hh"
#include "Framework/MapGeoToTrsAndIdx.hh"
#include "Framework/State.hh"
#include "Framework/CFSide.hh"
#include "Framework/Framework.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

 namespace Framework {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<CellCenteredSparsity,
               GlobalJacobianSparsity,
               FrameworkLib>
aCellCenteredSparsityProvider("CellCentered");

//////////////////////////////////////////////////////////////////////////////

CellCenteredSparsity::CellCenteredSparsity() : GlobalJacobianSparsity()
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

CellCenteredSparsity::~CellCenteredSparsity()
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

void CellCenteredSparsity::computeNNz(std::valarray<CFint>& nnz,
                                      std::valarray<CFint>& ghostNnz)
{
  CFAUTOTRACE;

  cf_assert(nnz.size() == ghostNnz.size());

  cf_assert(socket_states.isConnected());

  // get states datahandle
  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();

  // get InnerCells and InnerFaces TRSs
  SafePtr<TopologicalRegionSet> innerCells =
    MeshDataStack::getActive()->getTrs("InnerCells");

  SafePtr<TopologicalRegionSet> innerFaces =
    MeshDataStack::getActive()->getTrs("InnerFaces");

  // get connectivities
  Common::SafePtr< Common::ConnectivityTable<CFuint> > faceToCells =
    MeshDataStack::getActive()->getConnectivity("InnerFaces-Faces2Cells");

  SafePtr< ConnectivityTable< CFuint > > cellStates =
    MeshDataStack::getActive()->getConnectivity("cellStates_InnerCells");

  // number of cells
  const CFuint nbCells = cellStates->nbRows();
  cf_assert(nbCells == innerCells->getLocalNbGeoEnts());

  // loop on cells and add contributions
  for (CFuint iCell = 0; iCell < nbCells; ++iCell)
  {
    // number of states in this cell
    const CFuint nbStates = cellStates->nbCols(iCell);
    cf_assert(nbStates > 0);

    // add number of non-zero contributions to each state in cell
    const CFuint firstStateID = (*cellStates)(iCell,0);
    if (states[firstStateID]->isParUpdatable()) // if the first cell-state is parUpdatable, they all should be
    {
      for (CFuint iState = 0; iState < nbStates; ++iState)
      {
        const CFuint stateID = (*cellStates)(iCell,iState);
        cf_assert(stateID < states.size());
        cf_assert(states[stateID]->isParUpdatable());
        cf_assert(stateID < nnz.size());
        nnz[stateID] += nbStates;
      }
    }
    else
    {
      for (CFuint iState = 0; iState < nbStates; ++iState)
      {
        const CFuint stateID = (*cellStates)(iCell,iState);
        cf_assert(stateID < states.size());
        cf_assert(!states[stateID]->isParUpdatable());
        cf_assert(stateID < ghostNnz.size());
        ghostNnz[stateID] += nbStates;
      }
    }
  }

  // number of inner faces
  const CFuint nbInnerFaces = faceToCells->nbRows();
  cf_assert(innerFaces->getLocalNbGeoEnts() == nbInnerFaces);

  // neighbouring cell IDs and number of states in cells
  vector< CFuint > cellIDs(2);
  vector< CFuint > nbStatesCells(2);

  // loop on faces and add contributions
  for (CFuint iFace = 0; iFace < nbInnerFaces; ++iFace)
  {
    for (CFuint iSide = 0; iSide < 2; ++iSide)
    {
      cellIDs[iSide] = (*faceToCells)(iFace,iSide);
      cf_assert(cellIDs[iSide] < nbCells);
      nbStatesCells[iSide] = cellStates->nbCols(cellIDs[iSide]);
    }

    // add neigbouring cell contributions
    for (CFuint iSide = 0; iSide < 2; ++iSide)
    {
      // variable for other side
      const CFuint iOtherSide = iSide == LEFT ? RIGHT : LEFT;

      const CFuint firstStateID = (*cellStates)(cellIDs[iOtherSide],0);
      if (states[firstStateID]->isParUpdatable())
      // if the first cell-state is parUpdatable, than all should be!!
      {
        for (CFuint iState = 0 ; iState < nbStatesCells[iSide]; ++iState)
        {
          const CFuint stateID = (*cellStates)(cellIDs[iSide],iState);
          cf_assert(stateID < states.size());
          cf_assert(stateID < nnz.size());
          nnz[stateID] += nbStatesCells[iOtherSide];
        }
      }
      else
      // if the first cell-state is ghost, then all should be!!
      {
        for (CFuint iState = 0 ; iState < nbStatesCells[iSide]; ++iState)
        {
          const CFuint stateID = (*cellStates)(cellIDs[iSide],iState);
          cf_assert(stateID < states.size());
          cf_assert(stateID < ghostNnz.size());
          ghostNnz[stateID] += nbStatesCells[iOtherSide];
        }
      }
    }
  } // loop inner faces

  ofstream fout ("nnz.dat");

  fout << " NNZ " << "\n";
  for (CFuint i = 0 ; i < nnz.size(); ++i)
  {
    fout << nnz[i] << "\n";
  }
  fout << "Ghost NNZ " << "\n";
  for (CFuint i = 0 ; i < ghostNnz.size(); ++i)
  {
    fout << ghostNnz[i] << "\n";
  }
}

//////////////////////////////////////////////////////////////////////////////

void CellCenteredSparsity::computeMatrixPattern(
  std::valarray<CFint>& nnz,
  std::valarray<CFint>& ghostNnz,
  vector< vector<CFuint> >& matrixPattern)
{
  CFAUTOTRACE;

  throw Common::NotImplementedException (FromHere(),"GlobalJacobianSparsity::computeMatrixPattern() was not implemented yet.");
}

//////////////////////////////////////////////////////////////////////////////

void CellCenteredSparsity::computeMatrixPattern
(DataSocketSink<Framework::State*, Framework::GLOBAL> statesSocket,
 Common::ConnectivityTable<CFuint>& matrixPattern)
{
  CFAUTOTRACE;

  // get states datahandle
  DataHandle < Framework::State*, Framework::GLOBAL > states = statesSocket.getDataHandle();

  SafePtr<TopologicalRegionSet> innerFaces =
    MeshDataStack::getActive()->getTrs("InnerFaces");

  // get connectivities
  Common::SafePtr< Common::ConnectivityTable<CFuint> > faceToCells =
    MeshDataStack::getActive()->getConnectivity("InnerFaces-Faces2Cells");

  // number of cells
  // number of inner faces
  const CFuint nbInnerFaces = faceToCells->nbRows();
  const CFuint nbCells = MeshDataStack::getActive()->getTrs("InnerCells")->getLocalNbGeoEnts();
  valarray<CFuint> nbRows(nbCells);
  nbRows = 0;

  // count the number ofneighbors per row
  for (CFuint iFace = 0; iFace < nbInnerFaces; ++iFace)
  {
    nbRows[(*faceToCells)(iFace,0)] += 1;
    nbRows[(*faceToCells)(iFace,1)] += 1;
  }

  cf_assert(nbRows.sum() == nbInnerFaces*2);

  matrixPattern.resize(nbRows);
  nbRows = 0;

  // loop on faces and add neighbor IDs
  for (CFuint iFace = 0; iFace < nbInnerFaces; ++iFace)
  {
    const CFuint lCellID = (*faceToCells)(iFace,0);
    const CFuint rCellID = (*faceToCells)(iFace,1);
    matrixPattern(lCellID, nbRows[lCellID]) = rCellID;
    matrixPattern(rCellID, nbRows[rCellID]) = lCellID;
    nbRows[lCellID]++;
    nbRows[rCellID]++;
  }

  cf_assert(nbRows.sum() == nbInnerFaces*2);
}
//////////////////////////////////////////////////////////////////////////////
/*
void CellCenteredSparsity::computeNNzNodeBased(std::valarray<CFint>& nnz, 
					std::valarray<CFint>& ghostNnz)
{
  CFAUTOTRACE;

  cf_assert(socket_states.isConnected());
  cf_assert(socket_nodes.isConnected());
  std::cout<<" I am inside FR CCCC *****************************************$******************$ "<<endl;
  DataHandle < Framework::Node*, Framework::GLOBAL > nodes =
    socket_nodes.getDataHandle();
  
  std::vector<std::string> tags;
  tags.push_back ( "inner" );
  tags.push_back ( "cell" );
  
  std::vector< Common::SafePtr<TopologicalRegionSet> > innerCellsList = MeshDataStack::getActive()->getFilteredTrsList(tags);
  
  const CFuint nbNodes = nodes.size();
  vector< vector< CFuint > > neighborCells(nbNodes);
  vector< vector< Common::SafePtr<TopologicalRegionSet> > > neighborCellsTrs(nbNodes);
  
  // detect all the neighbor Cells for each Node
  const CFuint nbGroups = innerCellsList.size();
  for(CFuint iGroup=0; iGroup < nbGroups; ++iGroup){
    Common::SafePtr<TopologicalRegionSet> currTrs = innerCellsList[iGroup];
    SafePtr< ConnectivityTable< CFuint > > cellNodes =
      MeshDataStack::getActive()->getConnectivity("cellNodes_" + currTrs->getName());
    
    const CFuint nbCells = cellNodes->nbRows();
    for (CFuint iCell = 0; iCell < nbCells; ++iCell) {
      const CFuint nbNodesInCell = cellNodes->nbCols(iCell);
      for (CFuint iNode = 0; iNode < nbNodesInCell; ++iNode) {
        const CFuint stateID = (*cellNodes)(iCell,iNode);
        cf_assert(stateID < nbNodes);
        neighborCells[stateID].push_back(iCell);
        neighborCellsTrs[stateID].push_back(currTrs);
      }
    }
  }

  typedef set< Node*, less< Node* > > SetOfNeighbors;

  for (CFuint iNode = 0; iNode < nbNodes; ++iNode)
  {
    const CFuint nbNeighborCells = neighborCells[iNode].size();
    SetOfNeighbors* neighborNodes = new SetOfNeighbors();


    for (CFuint iCell = 0; iCell < nbNeighborCells; ++iCell) {
      const CFuint cellTrsID = neighborCells[iNode][iCell];
      Common::SafePtr<TopologicalRegionSet> cellTrs = neighborCellsTrs[iNode][iCell];
      SafePtr< ConnectivityTable< CFuint > > cellNodes =
        MeshDataStack::getActive()->getConnectivity("cellNodes_" + cellTrs->getName());

      const CFuint nbNodesInCell = cellNodes->nbCols(cellTrsID);
      for (CFuint jNode = 0; jNode < nbNodesInCell; ++jNode) {
        // store all the neighbor nodes + the considered node itself
        Node *const currNode = nodes[(*cellNodes)(cellTrsID,jNode)];
        neighborNodes->insert(SetOfNeighbors::value_type(currNode));
      }
    }
    
    // number of neighbor nodes == number of non zero entries in the matrix
    const CFuint nbNeighNodes = neighborNodes->size();
    nnz[iNode] = nbNeighNodes;
    deletePtr( neighborNodes );

  }
}
*/
void CellCenteredSparsity::computeNNzNodeBased(std::valarray<CFint>& nnz, 
					std::valarray<CFint>& ghostNnz)
{
  CFAUTOTRACE;

  cf_assert(socket_states.isConnected());
  cf_assert(socket_nodes.isConnected());
  cf_assert(nnz.size() == ghostNnz.size());
  
  DataHandle < Framework::Node*, Framework::GLOBAL > nodes =
    socket_nodes.getDataHandle();
  
  std::vector<std::string> tags;
  tags.push_back ( "inner" );
  tags.push_back ( "cell" );
  
  std::vector< Common::SafePtr<TopologicalRegionSet> > innerCellsList = MeshDataStack::getActive()->getFilteredTrsList(tags);
  
  const CFuint nbNodes = nodes.size();
  vector< vector< CFuint > > neighborCells(nbNodes);
  vector< vector< Common::SafePtr<TopologicalRegionSet> > > neighborCellsTrs(nbNodes);
  
  // detect all the neighbor Cells for each Node
  const CFuint nbGroups = innerCellsList.size();
  for(CFuint iGroup=0; iGroup < nbGroups; ++iGroup){
    Common::SafePtr<TopologicalRegionSet> currTrs = innerCellsList[iGroup];
    SafePtr< ConnectivityTable< CFuint > > cellNodes =
      MeshDataStack::getActive()->getConnectivity("cellNodes_" + currTrs->getName());
    
    const CFuint nbCells = cellNodes->nbRows();
    for (CFuint iCell = 0; iCell < nbCells; ++iCell) {
      const CFuint nbNodesInCell = cellNodes->nbCols(iCell);
      for (CFuint iNode = 0; iNode < nbNodesInCell; ++iNode) {
        const CFuint stateID = (*cellNodes)(iCell,iNode);
        cf_assert(stateID < nbNodes);
        neighborCells[stateID].push_back(iCell);
        neighborCellsTrs[stateID].push_back(currTrs);
      }
    }
  }

  typedef set< Node*, less< Node* > > SetOfNeighbors;

  for (CFuint iNode = 0; iNode < nbNodes; ++iNode)
  {
    const CFuint nbNeighborCells = neighborCells[iNode].size();
    SetOfNeighbors* neighborNodes = new SetOfNeighbors();
    SetOfNeighbors* ghostNeighNodes = new SetOfNeighbors();

    for (CFuint iCell = 0; iCell < nbNeighborCells; ++iCell) {
      const CFuint cellTrsID = neighborCells[iNode][iCell];
      Common::SafePtr<TopologicalRegionSet> cellTrs = neighborCellsTrs[iNode][iCell];
      SafePtr< ConnectivityTable< CFuint > > cellNodes =
        MeshDataStack::getActive()->getConnectivity("cellNodes_" + cellTrs->getName());

      const CFuint nbNodesInCell = cellNodes->nbCols(cellTrsID);
      for (CFuint jNode = 0; jNode < nbNodesInCell; ++jNode) {
        // store all the neighbor nodes + the considered node itself
        Node *const currNode = nodes[(*cellNodes)(cellTrsID,jNode)];
        neighborNodes->insert(SetOfNeighbors::value_type(currNode));
        if (!currNode->isParUpdatable()) {
          // if the state is not updatable store it in ghostNeighNodes
          ghostNeighNodes->insert(SetOfNeighbors::value_type(currNode));
        }
      }
    }
    
    // number of neighbor nodes == number of non zero entries in the matrix
    const CFuint nbNeighNodes = neighborNodes->size();
    nnz[iNode] = nbNeighNodes;
    ghostNnz[iNode] = ghostNeighNodes->size();
    deletePtr( neighborNodes );
    deletePtr( ghostNeighNodes );
  }
}
//////////////////////////////////////////////////////////////////////////////

    } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
