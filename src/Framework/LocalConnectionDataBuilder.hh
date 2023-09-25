// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Framework_LocalConnectionDataBuilder_hh
#define COOLFluiD_Framework_LocalConnectionDataBuilder_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/LocalConnectionData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

/// This class represents the builder of the local connectivity data
/// for 2D and 3D meshes.
/// LocalConnectionData constructs and uses this class and should
/// be the only client of this class.
/// The local connectivity requested to be constructed is stored in
/// {@see Table} of unsigned integers (CFuint) specifying for instance
/// the local ID of a Node or a State inside a Face.
/// For the moment only connectivities face to DOFs are computed, but it
/// could be useful to compute also element-face for grid
/// adaptation
/// @author Andrea Lani
/// @author Tiago Quintino
/// @author Thomas Wuilbaut
class Framework_API LocalConnectionDataBuilder {

public:

  /// @return the table of connectivity face to DOF for a line of order 1
  static Common::Table<CFuint>* faceDofLineOrder1();

  /// @return the table of connectivity face to DOF for a triangle of order 1
  static Common::Table<CFuint>* faceDofTriagOrder1();

  /// @return the table of connectivity face to DOF for a triangle of order 2
  static Common::Table<CFuint>* faceDofTriagOrder2();

  /// @return the table of connectivity face to DOF for a triangle of order 3
  static Common::Table<CFuint>* faceDofTriagOrder3();

  /// @return the table of connectivity face to DOF for a quadrilateral of order 1
  static Common::Table<CFuint>* faceDofQuadOrder1();

  /// @return the table of connectivity face to DOF for a quadrilateral of order 2
  static Common::Table<CFuint>* faceDofQuadOrder2();

  /// @return the table of connectivity face to DOF for a quadrilateral of order 3
  static Common::Table<CFuint>* faceDofQuadOrder3();

  /// @return the table of connectivity face to DOF for a quadrilateral of order 4
  static Common::Table<CFuint>* faceDofQuadOrder4();

  /// @return the table of connectivity face to DOF for a quadrilateral of order 5
  static Common::Table<CFuint>* faceDofQuadOrder5();

  /// @return the table of connectivity face to DOF for a Tetrahedron of order 1
  /// @pre the numbering is anticlockwise on each face if you look
  ///      at the face from inside the tetra.
  static Common::Table<CFuint>* faceDofTetraOrder1();

  /// @return the table of connectivity face to DOF for a Tetrahedron of order 2
  /// @pre the numbering is anticlockwise on each face if you look
  ///      at the face from inside the tetra.
  static Common::Table<CFuint>* faceDofTetraOrder2();

  /// @return the table of connectivity face to DOF for a Pyramid of order 1
  /// @pre the numbering is anticlockwise on each face if you look
  ///      at the face from outside the pyramid.
  static Common::Table<CFuint>* faceDofPyramOrder1();

  /// @return the table of connectivity face to DOF for a Prism of order 1
  /// @pre the numbering is anticlockwise on each face if you look
  ///      at the face from outside the Prism.
  static Common::Table<CFuint>* faceDofPrismOrder1();
  
  /// @return the table of connectivity face to DOF for a Prism of order 2
  /// @pre the numbering is anticlockwise on each face if you look
  ///      at the face from outside the Prism.
  static Common::Table<CFuint>* faceDofPrismOrder2();

  /// @return the table of connectivity face to DOF for a Hexahedron of order 1
  /// @pre the numbering is anticlockwise on each face if you look
  ///      at the face from outside the hexahedron.
  static Common::Table<CFuint>* faceDofHexaOrder1();

  /// @return the table of connectivity face to DOF for a Hexahedron of order 2
  /// @pre the numbering is anticlockwise on each face if you look
  ///      at the face from outside the hexahedron.
  ///      and starts with all the corner nodes before the central nodes
  static Common::Table<CFuint>* faceDofHexaOrder2();

  /// @return the table of connectivity face to DOF for a Hexahedron of order 3
  static Common::Table<CFuint>* faceDofHexaOrder3();

  /// @return the table of connectivity face to DOF for an incomplete (20 nodes) Hexahedron of order 2
  /// @pre the numbering is anticlockwise on each face if you look
  ///      at the face from outside the hexahedron.
  ///      and starts with all the corner nodes before the central nodes
  static Common::Table<CFuint>* faceDofHexa20NodesOrder2();

  /// @return the table of connectivity edge-dof for a Tetrahedron of order 1
  static Common::Table<CFuint>* edgeDofTetraOrder1();

  /// @return the table of connectivity edge-dof for a Tetrahedron of order 2
  static Common::Table<CFuint>* edgeDofTetraOrder2();

  /// @return the table of connectivity edge-dof for a Pyramid of order 1
  static Common::Table<CFuint>* edgeDofPyramOrder1();

  /// @return the table of connectivity edge-dof for a Prism of order 1
  static Common::Table<CFuint>* edgeDofPrismOrder1();
  
  /// @return the table of connectivity edge-dof for a Prism of order 2
  static Common::Table<CFuint>* edgeDofPrismOrder2();

  /// @return the table of connectivity edge-dof for a Hexahedron of order 1
  static Common::Table<CFuint>* edgeDofHexaOrder1();

  /// @return the table of connectivity edge-dof for a Hexahedron of order 2
  static Common::Table<CFuint>* edgeDofHexaOrder2();

  /// @return the table of connectivity edge-dof for a Hexahedron of order 3
  static Common::Table<CFuint>* edgeDofHexaOrder3();

  /// @return the table of connectivity edge-dof for an incomplete Hexahedron of order 2
  static Common::Table<CFuint>* edgeDofHexa20NodesOrder2();

private:

  /// Private Constructor
  LocalConnectionDataBuilder();

  /// Private destructor
  ~LocalConnectionDataBuilder();

}; // end of class LocalConnectionDataBuilder

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Framework_LocalConnectionDataBuilder_hh
