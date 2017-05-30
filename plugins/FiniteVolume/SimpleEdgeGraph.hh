#ifndef COOLFluiD_Numerics_FiniteVolume_SimpleEdge_hh
#define COOLFluiD_Numerics_FiniteVolume_SimpleEdge_hh

//////////////////////////////////////////////////////////////////////////////

#include <iostream>

#include "Framework/CellTrsGeoBuilder.hh"
#include "Framework/MeshData.hh"
#include "Framework/TopologicalRegionSet.hh"
#include "Common/ConnectivityTable.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace Numerics {
    namespace FiniteVolume {
      
//////////////////////////////////////////////////////////////////////////////

//A simple Edge with only information of the edges
class SimpleEdge{
public:
  SimpleEdge( Framework::Node* node1 , Framework::Node* node2):
    firstNode(node1),
    secondNode(node2)
  { 
  }
  
  SimpleEdge(){}
  
public:
  Framework::Node* firstNode;
  Framework::Node* secondNode;
  
};
      
class iterator;

//A simple edge graph including iterator
class SimpleEdgeGraph
{

public: 
  
  SimpleEdgeGraph():socket_nodes("socket_nodes"){;}
  
  void setNodeDataSocket( Framework::DataSocketSink < Framework::Node* , 
			  Framework::GLOBAL > nodes );
  
  void computeConnectivity();
  
  class iterator{
  public:
    
    iterator(const SimpleEdgeGraph &edgeGraph) :
      m_edgeGraph(&edgeGraph),
      m_row(0),
      m_col(0)
    {}
    
    iterator(const iterator& it):
      m_edgeGraph(it.m_edgeGraph),
      m_row(it.m_row ),
      m_col(it.m_col)
    {}
    
    iterator& operator=( const iterator& it){  
      m_edgeGraph = it.m_edgeGraph;
      m_row = it.m_row;
      m_col = it.m_col;
      return *this;
    }
    
    iterator begin(){
      iterator temp(*this);
      temp.m_row = 0;
      temp.m_col = 0; 
      return temp;
    }
    
    iterator end(){ 
      iterator temp(*this);
      temp.m_row = m_edgeGraph->m_connectivityTable.size();
      temp.m_col = 0;
      return temp; 
    }
    
    iterator& operator++(){
      ++m_col;
      while( (m_col == m_edgeGraph->m_connectivityTable[m_row].size()) 
	     && m_row < m_edgeGraph->m_connectivityTable.size() ) {
	m_col = 0;
	++m_row;
      }
      return *this;
    }
    
    iterator operator++(int){
      iterator tmp(*this);
      operator++();
      return tmp;
    }
    
    bool operator==(const iterator &rhs){
      return (m_row == rhs.m_row && m_col == rhs.m_col);
    }
    
    bool operator!=(const iterator &rhs){
      return !(m_row == rhs.m_row && m_col == rhs.m_col);
    }
    
    SimpleEdge operator*(){
      Framework::DataHandle < Framework::Node*, Framework::GLOBAL > nodes = 
	m_edgeGraph->socket_nodes.getDataHandle();
      const CFuint firstNodeID = m_row; 
      const CFuint secondNodeID = m_edgeGraph->m_connectivityTable[firstNodeID][m_col];
      SimpleEdge tmp( nodes[firstNodeID], nodes[secondNodeID]);
      return tmp;
    }
    
    SimpleEdge operator->(){
      iterator temp(*this);
      return( *temp  );
      
    }
  private:
    SimpleEdgeGraph const *m_edgeGraph;
    CFuint m_row;
    CFuint m_col;
  };
  
  iterator begin(){
    iterator temp(*this);
    return temp.begin();
  }
  
  iterator end(){ 
    iterator temp(*this);
    return temp.end();
  }
  
  
  const std::vector<Framework::Node*>& getNeighborNodesOfNode(const Framework::Node* firstNode) {
    Framework::DataHandle < Framework::Node*, Framework::GLOBAL > nodes = socket_nodes.getDataHandle();
    CFuint nodeID = firstNode->getLocalID();
    m_connectedNodes.resize(m_connectivityTable[nodeID].size());
    for(CFuint i=0; i<m_connectivityTable[nodeID].size(); ++i){
      m_connectedNodes[i] = nodes[m_connectivityTable[nodeID][i]];
    }
    return m_connectedNodes;
  }
  
private:
  /// the socket to the data handle of the state's
  Framework::DataSocketSink < Framework::Node* , Framework::GLOBAL > socket_nodes;
  
  //connectivity table
  std::vector< std::vector<CFuint> > m_connectivityTable;
  
  std::vector<Framework::Node*> m_connectedNodes;
};

//////////////////////////////////////////////////////////////////////////////

}
}
}

#endif
