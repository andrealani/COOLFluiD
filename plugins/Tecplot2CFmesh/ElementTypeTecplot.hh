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
  
  ElementTypeTecplot(CFuint nbElems, CFuint nbNodesPerElem, CFuint nbNodes, CFuint nbVars) :
    m_nbElems(nbElems),
    m_nbNodesPerElem(nbNodesPerElem),
    m_nbNodes(nbNodes),
    m_nbVars(nbVars),
    m_cellNode(nbElems*nbNodesPerElem),
    m_neighbor(nbElems, -1),
    m_variables(nbNodes*nbVars)
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
    assert(idx < m_cellNode.size()); 
    return m_cellNode[idx];
  }
  
  /**
   * Operator overloading to access individual element node IDs
   */
  CFuint operator()(CFuint iElem , CFuint iNode) const
  {
    const CFuint idx = iElem*m_nbNodesPerElem + iNode;
    assert(idx < m_cellNode.size());  
    return m_cellNode[idx];
  }
  
  /**
   * Set variable in the nodes
   */
  void setVar(CFuint iNode, CFuint iVar, CFreal value) 
  {
    const CFuint idx = iNode*m_nbVars+iVar;
    assert(idx < m_variables.size());  
    m_variables[idx] = value;
  }
   
  /**
   * Get variable in the nodes
   */
  CFreal getVar(CFuint iNode, CFuint iVar) const
  {
    const CFuint idx = iNode*m_nbVars+iVar;
    assert(idx < m_variables.size());
    return m_variables[idx];
  } 
  
  /**
   * Set the coordinates in NodeDimGet pointer to the variables in the nodes
   */
  template <CFuint DIM>
  void setNodeDim(CFuint iNode, NodeDim<DIM>& nodeDim)
  {
    const CFuint idx = iNode*m_nbVars;
    assert(idx < m_variables.size());
    nodeDim.reset(&m_variables[idx]);
  }
  
  /**
   * Set neighbor ID 
   */
  void setNeighborID(CFuint iElem, CFint neighborID)
  {
    assert(neighborID >= 0);
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
  
  /// total nb nodes
  CFuint      m_nbNodes;
  
  /// nb variables
  CFuint       m_nbVars;
  
  /// cell to node connectivity table
  std::vector<CFuint>  m_cellNode;
  
  /// neighbor cell (only for faces in FV type mesh)
  std::vector<CFint>  m_neighbor;
  
  /// cell to node connectivity table
  std::vector<CFreal>  m_variables;
  
}; // end class ElementTypeTecplot

//////////////////////////////////////////////////////////////////////////////

    } // namespace Tecplot2CFmesh

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_IO_Tecplot2CFmesh_ElementTypeTecplot_hh
