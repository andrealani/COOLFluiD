// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_IO_THOR2CFmesh_CheckNodeNumbering_hh
#define COOLFluiD_IO_THOR2CFmesh_CheckNodeNumbering_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/VolumeCalculator.hh"
#include "Common/Table.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace THOR2CFmesh {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a checker for node numberings coming
 * from THOR files.
 * This is important because THOR format does not specify a
 * node numbering, hence some numberings can produce negative volumes
 * in some cells.
 *
 * @author Andrea Lani
 * @author Tiago Quintino
 *
 */
class CheckNodeNumbering {
public:

  /**
   * Constructor
   */
  CheckNodeNumbering(const Common::Table<CFreal> *const nodes) :
    _nodes(nodes),
    _elementNode(CFNULL),
    _calculator(),
    _nodalCoord()
  {
    // table can only be (nbNodes x dimendion)
    // so it cannot be hybrib
    cf_assert(nodes->isNotHybrid());
  }

  /**
   * Default destructor
   */
  virtual ~CheckNodeNumbering()
  {
  }

  /**
   * Checks the numbering of the THOR read nodes and makes it stick
   * to the COOLFluiD convention.
   */
  void checkElementNodalNumbering(Common::Table<CFuint>& elementNode)
  {
    setElemNodeTable(elementNode);
    checkElementNodalNumberingImpl();
  }

protected:

  /**
   * Checks the numbering of the THOR read nodes and makes it stick
   * to the COOLFluiD convention.
   */
  virtual void checkElementNodalNumberingImpl() = 0;

protected: // functions

  /**
   * Helper function to copy the coordinates into the
   * temporary matrix _nodalCoord
   */
  template <unsigned int NBNODES,
            unsigned int NBDIM>
  void copyCoord(const CFuint& iElem)
  {
    CFuint nodeID;
    // for each local node iNode get the global ID
    for (CFuint iNode = 0; iNode < NBNODES; ++iNode) {

     nodeID = static_cast<CFuint>((*_elementNode)(iElem,iNode));

     // copy coordinates into temporary matrix for volume computation
     for (CFuint d = 0; d < NBDIM; ++d) {
       _nodalCoord(iNode,d) = (*_nodes)(nodeID,d);
     }
    }
  }

  /**
   * Helper function to
   * set the temporary element-node table
   */
  void setElemNodeTable(Common::Table<CFuint>& elementNode)
  {
    _elementNode = &elementNode;
  }

protected: // data

  // accquaintance of the nodes to check
  const Common::Table<CFreal> *const _nodes;

  /// a temporary accquaintance of the element to node table
  /// passed in checkElementNodalNumbering
  Common::Table<CFuint> * _elementNode;

  // the volume calculator to use
  Framework::VolumeCalculator _calculator;

  // the temporary matris to store element coordinates
  RealMatrix _nodalCoord;

}; // end of class CheckNodeNumbering

//////////////////////////////////////////////////////////////////////////////

    } // namespace THOR2CFmesh

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_IO_THOR2CFmesh_CheckNodeNumbering_hh
