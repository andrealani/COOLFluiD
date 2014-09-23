// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_CGNS2CFmesh_ElementTypeCGNS_hh
#define COOLFluiD_CGNS2CFmesh_ElementTypeCGNS_hh

//////////////////////////////////////////////////////////////////////////////

#include "Common/Table.hh"


//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace CGNS2CFmesh {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class stores information about the types of elements in CGNS format
 *
 * @author Andrea Lani
 * @author Tiago Quintino
 *
 */
class ElementTypeCGNS {
public:

  ElementTypeCGNS() :
    m_nbCellsPerType(0),
    m_nbNodesPerCell(0)
  {
  }

  ~ElementTypeCGNS()
  {
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

 protected:

  /// nb cells with this type
  CFuint      m_nbCellsPerType;

  /// nb nodes per cell
  CFuint       m_nbNodesPerCell;

 }; // end class ElementTypeCGNS

//////////////////////////////////////////////////////////////////////////////

    } // namespace CGNS2CFmesh

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_CGNS2CFmesh_ElementTypeCGNS_hh
