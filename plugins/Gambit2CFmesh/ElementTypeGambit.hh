// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_IO_Gambit2CFmesh_ElementTypeGambit_hh
#define COOLFluiD_IO_Gambit2CFmesh_ElementTypeGambit_hh

//////////////////////////////////////////////////////////////////////////////

#include "Common/ConnectivityTable.hh"


//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace IO {

    namespace Gambit2CFmesh {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class stores information about the types of elements in Gambit format
 *
 * @author Thomas Wuilbaut
 *
 */
class ElementTypeGambit {
public:

  /**
   * Typedef of a table for indexes
   */
  typedef Common::ConnectivityTable<CFuint> IndexTable;

  ElementTypeGambit() :
    _nbCellsPerType(0),
    _nbNodesPerCell(0),
    _cellNode(CFNULL),
    _typeID(0),
    _globalElemIDs()
  {
  }

  ~ElementTypeGambit()
  {
  }

  /**
   * Gets the number of Nodes per cell
   */
  CFuint getNbCellsPerType() const
  {
    return _nbCellsPerType;
  }

  /**
   * Gets the number of Nodes per cell
   */
  CFuint getNbNodesPerCell() const
  {
    return _nbNodesPerCell;
  }

  /**
   * Sets the number of Cells per type
   */
  void setNbCellsPerType(const CFuint& nbCellsPerType)
  {
    cf_assert(nbCellsPerType > 0);
    _nbCellsPerType = nbCellsPerType;
  }

  /**
   * Sets the number of nodes per cell
   */
  void setNbNodesPerCell(const CFuint& nbNodesPerCell)
  {
    cf_assert(nbNodesPerCell > 0);
    _nbNodesPerCell = nbNodesPerCell;
  }
  
  /**
   * Gets the table of connectivity
   */
  CFuint getNodeID(CFuint localTypeElemID, CFuint iNode) const
  {
    cf_assert(_cellNode != CFNULL);
    return (*_cellNode)(_globalElemIDs[localTypeElemID], iNode);
  } 
  
  /**
   * Sets the global element ID corresponding to this type and this local type element ID 
   */
  void setElementGlobalID(CFuint globalElemID, CFuint localTypeElemID)
  {
    _globalElemIDs[localTypeElemID] = globalElemID;
  }
  
  /**
   * Creates the cell to node connectivity
   */
  void createCellNodeConnectivity(IndexTable *const table)
  {
    cf_assert(_nbCellsPerType > 0);
    cf_assert(_nbNodesPerCell > 0);
    cf_assert(table != CFNULL);
    _cellNode = table;
    _globalElemIDs.resize(_nbCellsPerType);
  }

  void setTypeID(CFuint id){
    _typeID = id;
  }

  CFuint getTypeID(){
    return _typeID;
  }

 protected:

  /// nb cells with this type
  CFuint       _nbCellsPerType;

  /// nb nodes per cell
  CFuint       _nbNodesPerCell;

  /// cell to node connectivity table
  IndexTable*  _cellNode;

  /// typeID
  CFuint _typeID;
  
  /// array of row IDs corresponding to the element global IDs
  /// in the global table of connectivity
  std::vector<CFuint> _globalElemIDs;
  
 }; // end class ElementTypeGambit


//////////////////////////////////////////////////////////////////////////////

    } // namespace Gambit2CFmesh

  } // namespace IO

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_IO_Gambit2CFmesh_ElementTypeGambit_hh
