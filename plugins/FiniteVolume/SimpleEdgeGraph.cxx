#include "FiniteVolume/SimpleEdgeGraph.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace Numerics {
    namespace FiniteVolume {
      
//////////////////////////////////////////////////////////////////////////////

void SimpleEdgeGraph::setNodeDataSocket( Framework::DataSocketSink < Framework::Node* , Framework::GLOBAL > nodes )
{
  socket_nodes = nodes;
}
 
void SimpleEdgeGraph::computeConnectivity(){

  Framework::DataHandle < Framework::Node*, Framework::GLOBAL > nodes = socket_nodes.getDataHandle();
  cf_assert(socket_nodes.isConnected());
  
  const CFuint nbNodes = nodes.size();
  m_connectivityTable.clear();
  m_connectivityTable.resize(nbNodes);


  std::vector<std::string> tags;
  tags.push_back ( "inner" );
  tags.push_back ( "cell" );
  
  std::vector< Common::SafePtr<Framework::TopologicalRegionSet> > innerCellsList = 
    Framework::MeshDataStack::getActive()->getFilteredTrsList(tags);
  
  
  std::vector< std::vector< CFuint > > neighborCells(nbNodes);
  std::vector< std::vector< Common::SafePtr<Framework::TopologicalRegionSet> > > neighborCellsTrs(nbNodes);
  
  // detect all the neighbor Cells for each Node
  const CFuint nbGroups = innerCellsList.size();
  for(CFuint iGroup=0; iGroup < nbGroups; ++iGroup){
    Common::SafePtr<Framework::TopologicalRegionSet> currTrs = innerCellsList[iGroup];
    Common::SafePtr<Common::ConnectivityTable< CFuint > > cellNodes =
      Framework::MeshDataStack::getActive()->getConnectivity("cellNodes_" + currTrs->getName());
    
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
  
  for (CFuint iNode = 0; iNode < nbNodes; ++iNode)
  {
    m_connectivityTable[iNode].clear();
    const CFuint nbNeighborCells = neighborCells[iNode].size();

    for (CFuint iCell = 0; iCell < nbNeighborCells; ++iCell)
    {
      const CFuint cellTrsID = neighborCells[iNode][iCell];
      Common::SafePtr<Framework::TopologicalRegionSet> cellTrs = neighborCellsTrs[iNode][iCell];
      Common::SafePtr<Common::ConnectivityTable< CFuint > > cellNodes =
        Framework::MeshDataStack::getActive()->getConnectivity("cellNodes_" + cellTrs->getName());

      const CFuint nbNodesInCell = cellNodes->nbCols(cellTrsID);
      for (CFuint jNode = 0; jNode < nbNodesInCell; ++jNode) 
      {
        // store all the neighbor nodes without the considered node itself
        Framework::Node *const currNode = nodes[(*cellNodes)(cellTrsID,jNode)];
        if( currNode->getLocalID() != iNode 
            && std::find(m_connectivityTable[iNode].begin(), 
                         m_connectivityTable[iNode].end(), 
                         currNode->getLocalID()) ==  m_connectivityTable[iNode].end() )
        {
          m_connectivityTable[iNode].push_back(currNode->getLocalID());
        }
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

}
}
}
