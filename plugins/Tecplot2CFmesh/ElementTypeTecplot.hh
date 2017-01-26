// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_IO_Tecplot2CFmesh_ElementTypeTecplot_hh
#define COOLFluiD_IO_Tecplot2CFmesh_ElementTypeTecplot_hh

//////////////////////////////////////////////////////////////////////////////

#include <vector>
#include "Common/COOLFluiD.hh"
#include "Tecplot2CFmesh/NodeDim.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace Tecplot2CFmesh {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class stores information about the types of elements in Tecplot format
 *
 * @author Andrea Lani
 *
 */
class ElementTypeTecplot {
public:
  
  ElementTypeTecplot(CFuint nbElems, 
		     CFuint nbNodesPerElem,
		     CFuint nbNodes, 
		     CFuint nbVars) :
    m_nbElems(nbElems),
    m_nbNodesPerElem(nbNodesPerElem),
    m_nbUniqueNodesPerElem(nbNodesPerElem),
    m_nbNodes(nbNodes),
    m_nbVars(nbVars),
    m_cellNode(nbElems*nbNodesPerElem), 
    m_cellNodeBkp(nbElems*nbNodesPerElem),
    m_neighbor(nbElems, -1),
    m_variables(nbNodes*nbVars),
    m_localElemNodeIdx()
  {
  }
  
  ~ElementTypeTecplot()
  {
  }
  
  /**
   * Apply an offset to the numbering
   */
  void offsetNumbering(CFint offset)
  {
    for (CFuint i = 0; i < m_cellNode.size(); ++i) {
      m_cellNode[i] += offset;
    }
  }

  /**
   * Operator overloading to access individual element node IDs
   */
  CFuint& operator()(CFuint iElem , CFuint iNode) 
  {
    const CFuint idx = iElem*m_nbNodesPerElem + iNode;
    cf_assert(idx < m_cellNode.size()); 
    return m_cellNode[idx];
  }
  
  /**
   * Operator overloading to access individual element node IDs
   */
  CFuint operator()(CFuint iElem , CFuint iNode) const
  {
    const CFuint idx = iElem*m_nbNodesPerElem + iNode;
    cf_assert(idx < m_cellNode.size());  
    return m_cellNode[idx];
  }
  
  /**
   * Set the backup values for individual element node IDs
   */
  void setBkpElementNode(CFuint iElem , CFuint iNode, CFuint nodeID)
  {
    const CFuint idx = iElem*m_nbNodesPerElem + iNode;
    cf_assert(idx < m_cellNode.size());  
    m_cellNodeBkp[idx] = nodeID;
  }
  
  /**
   * Set the backup values for individual element node IDs
   */
  CFuint getBkpElementNode(CFuint iElem , CFuint iNode) const
  {
    const CFuint idx = iElem*m_nbNodesPerElem + iNode;
    cf_assert(idx < m_cellNode.size());  
    return m_cellNodeBkp[idx];
  }
  
  /**
   * Set variable in the nodes
   */
  void setVar(CFuint iNode, CFuint iVar, CFreal value) 
  {
    const CFuint idx = iNode*m_nbVars+iVar;
    cf_assert(idx < m_variables.size());  
    m_variables[idx] = value;
  }
   
  /**
   * Get variable in the nodes
   */
  CFreal getVar(CFuint iNode, CFuint iVar) const
  {
    const CFuint idx = iNode*m_nbVars+iVar;
    cf_assert(idx < m_variables.size());
    return m_variables[idx];
  } 
  
  /**
   * Get the pointer to the first variable for the given node
   */
  CFreal* getNodalVarPtr(CFuint iNode)
  {
    const CFuint idx = iNode*m_nbVars;
    cf_assert(idx < m_variables.size());
    return &m_variables[idx];
  } 
  
  /**
   * Set the coordinates in NodeDimGet pointer to the variables in the nodes
   */
  void setNodeDim(CFuint iNode, NodeDim& nodeDim)
  {
    const CFuint idx = iNode*m_nbVars;
    cf_assert(idx < m_variables.size());
    nodeDim.reset(&m_variables[idx]);
  }
  
  /**
   * Set neighbor ID 
   */
  void setNeighborID(CFuint iElem, CFint neighborID)
  {
    cf_assert(neighborID >= 0);
    m_neighbor[iElem] = neighborID;
  } 
  
  /**
   * Get neighbor ID
   */
  CFint getNeighborID(CFuint iElem) const
  {
    return m_neighbor[iElem];
  }
  
  /**
   * Gets the number of Nodes per cell
   */
  CFuint getNbElems() const {return m_nbElems;}

  /**
   * Gets the number of Nodes per cell
   */
  CFuint getNbNodesPerElem() const {return m_nbNodesPerElem;}
  
  /**
   * Gets the number of unique Nodes per cell
   */
  CFuint getNbUniqueNodesPerElem() const {return m_nbUniqueNodesPerElem;}
  
  /**
   * Sets the number of unique Nodes per cell
   */
  void setNbUniqueNodesPerElem() 
  {
    // consider the first element (all the other in this type are the same)
    // keep track of the local (i.e. inside the element itself) node IDs which have
    // to be written avoiding to consider duplicates
    m_localElemNodeIdx.reserve(m_nbNodesPerElem);
    std::set<CFuint> nodeIDList;
    for (CFuint j = 0; j < m_nbNodesPerElem; ++j) {
      const CFuint enodeID = (*this)(0,j);
      if (nodeIDList.count(enodeID) == 0) {
	nodeIDList.insert(enodeID);
	m_localElemNodeIdx.push_back(j);
      }
    }
    
    m_nbUniqueNodesPerElem = m_localElemNodeIdx.size();
    cf_assert(m_nbUniqueNodesPerElem <= m_nbNodesPerElem);
  }
  
  /**
   * Sets the number of unique Nodes per cell
   */
  void getUniqueNodesInSingleElem(CFuint iElem, 
				  std::set<CFuint>& nodeIDList, 
				  std::vector<CFuint>& localElemNodeIdx) 
  {
    nodeIDList.clear();
    localElemNodeIdx.clear();
    for (CFuint j = 0; j < m_nbNodesPerElem; ++j) {
      const CFuint enodeID = (*this)(iElem,j);
      if (nodeIDList.count(enodeID) == 0) {
	nodeIDList.insert(enodeID);
	localElemNodeIdx.push_back(j);
      }
    }
  }
  
  /**
   * Gets the list of local element non-duplicated node indices
   */
  const std::vector<CFuint>& getLocalElemNodeIdx() const 
  {
    return m_localElemNodeIdx;
  }
  
  /**
   * Gets the total number of Nodes
   */
  CFuint getNbNodes() const {return m_nbNodes;}
 
  /**
   * Gets the number of variables
   */
  CFuint getNbVars() const {return m_nbVars;}
  
private:
  /// nb cells with this type
  CFuint      m_nbElems;
  
  /// nb nodes per cell
  CFuint       m_nbNodesPerElem;
  
  /// unique nb nodes per cell (this is for degenerated bricks like prism, pyramids)
  CFuint       m_nbUniqueNodesPerElem;
  
  /// total nb nodes
  CFuint      m_nbNodes;
  
  /// nb variables
  CFuint       m_nbVars;
  
  /// cell to node connectivity table
  std::vector<CFuint>  m_cellNode;
  
  /// cell to node connectivity table (backup of read values, before any renumbering)
  std::vector<CFuint>  m_cellNodeBkp;
  
  /// neighbor cell (only for faces in FV type mesh)
  std::vector<CFint>  m_neighbor;
  
  /// cell to node connectivity table
  std::vector<CFreal>  m_variables;
  
  /// list of local element non-duplicated node indices
  std::vector<CFuint> m_localElemNodeIdx;
  
}; // end class ElementTypeTecplot

//////////////////////////////////////////////////////////////////////////////

    } // namespace Tecplot2CFmesh

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_IO_Tecplot2CFmesh_ElementTypeTecplot_hh
