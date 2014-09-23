// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_IO_Dpl2CFmesh_ElementTypeDpl_hh
#define COOLFluiD_IO_Dpl2CFmesh_ElementTypeDpl_hh

//////////////////////////////////////////////////////////////////////////////

#include "Common/Table.hh"


//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace IO {

    namespace Dpl2CFmesh {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class stores information about the types of elements in Dpl format
 *
 * @author Thomas Wuilbaut
 *
 */
class ElementTypeDpl {
public:

  /**
   * Typedef of a table for indexes
   */
  typedef Common::Table<CFuint> IndexTable;

  ElementTypeDpl() :
    _nbCellsPerType(0),
    _nbNodesPerCell(0),
    _cellNode(CFNULL)
  {
  }

  ~ElementTypeDpl()
  {
    deletePtr(_cellNode);
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
   * Sets the current index
   */
  void setCurrentIndex(const CFuint& index)
  {
    _index = index;
  }

  /**
   * Sets the index
   */
  CFuint getCurrentIndex()
  {
    return _index;
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
  IndexTable& getTableConnectivity()
  {
    cf_assert(_cellNode != CFNULL);
    return *_cellNode;
  }

  /**
   * Creates the cell to node connectivity
   */
  void createCellNodeConnectivity()
  {
    cf_assert(_nbCellsPerType > 0);
    cf_assert(_nbNodesPerCell > 0);
    _cellNode = new IndexTable(_nbCellsPerType,
                               _nbNodesPerCell);
  }

  void setTypeID(CFuint id){
    _typeID = id;
  }

  CFuint getTypeID(){
    return _typeID;
  }

 protected:

  /// nb cells with this type
  CFuint      _nbCellsPerType;

  /// nb nodes per cell
  CFuint       _nbNodesPerCell;

  /// cell to node connectivity table
  IndexTable* _cellNode;

  /// typeID
  CFuint _typeID;

  /// local index of the element
  CFuint _index;

 }; // end class ElementTypeDpl


//////////////////////////////////////////////////////////////////////////////

    } // namespace Dpl2CFmesh

  } // namespace IO

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_IO_Dpl2CFmesh_ElementTypeDpl_hh
