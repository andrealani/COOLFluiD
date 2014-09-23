// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_IO_THOR2CFmesh_ElementTypeTHOR_hh
#define COOLFluiD_IO_THOR2CFmesh_ElementTypeTHOR_hh

//////////////////////////////////////////////////////////////////////////////

#include "Common/Table.hh"


//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace THOR2CFmesh {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class stores information about the types of elements in THOR format
 *
 * @author Andrea Lani
 * @author Tiago Quintino
 *
 */
class ElementTypeTHOR {
public:

  /**
   * Typedef of a table for indexes
   */
  typedef Common::Table<CFuint> IndexTable;

  ElementTypeTHOR() :
    m_nbCellsPerType(0),
    m_nbNodesPerCell(0),
    m_cellNode(CFNULL)
  {
  }

  ~ElementTypeTHOR()
  {
    deletePtr(m_cellNode);
  }

  /**
   * Gets the number of Nodes per cell
   */
  CFuint getNbCellsPerType() const
  {
    return m_nbCellsPerType;
  }

  /**
   * Gets the number of Nodes per cell
   */
  CFuint getNbNodesPerCell() const
  {
    return m_nbNodesPerCell;
  }

  /**
   * Sets the number of Cells per type
   */
  void setNbCellsPerType(const CFuint& nbCellsPerType)
  {
    cf_assert(nbCellsPerType > 0);
    m_nbCellsPerType = nbCellsPerType;
  }

  /**
   * Sets the number of nodes per cell
   */
  void setNbNodesPerCell(const CFuint& nbNodesPerCell)
  {
    cf_assert(nbNodesPerCell > 0);
    m_nbNodesPerCell = nbNodesPerCell;
  }

  /**
   * Gets the table of connectivity
   */
  IndexTable& getTableConnectivity()
  {
    cf_assert(m_cellNode != CFNULL);
    return *m_cellNode;
  }

  /**
   * Creates the cell to node connectivity
   */
  void createCellNodeConnectivity()
  {
    cf_assert(m_nbCellsPerType > 0);
    cf_assert(m_nbNodesPerCell > 0);
    m_cellNode = new IndexTable(m_nbCellsPerType, m_nbNodesPerCell);
  }

 protected:

  /// nb cells with this type
  CFuint      m_nbCellsPerType;

  /// nb nodes per cell
  CFuint       m_nbNodesPerCell;

  /// cell to node connectivity table
  IndexTable* m_cellNode;

 }; // end class ElementTypeTHOR

//////////////////////////////////////////////////////////////////////////////

    } // namespace THOR2CFmesh

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_IO_THOR2CFmesh_ElementTypeTHOR_hh
