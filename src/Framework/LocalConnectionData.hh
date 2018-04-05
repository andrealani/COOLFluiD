// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Framework_LocalConnectionData_hh
#define COOLFluiD_Framework_LocalConnectionData_hh

//////////////////////////////////////////////////////////////////////////////

#include "Common/CFMap.hh"
#include "Common/Table.hh"

#include "Framework/CFGeoShape.hh"
#include "Framework/CFPolyOrder.hh"
#include "Framework/CFPolyForm.hh"
#include "Framework/CFDofType.hh"

#include "Framework/Framework.hh"
#include "Framework/LocalConnectionDataBuilder.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

/// This class holds and handles the data concerning the local connectivity.
/// @author Andrea Lani
/// @author Tiago Quintino
/// @author Martin Vymazal
class Framework_API LocalConnectionData : public Common::NonCopyable<LocalConnectionData> {
public:

  /// Definition of a Map from a certain code to a table of local indexes
  typedef Common::CFMap<CFuint, Common::Table<CFuint>*> MapCodeTable;

  /// Definition of a Map from a certain code to a shape of a Face
  typedef Common::CFMap<CFuint, CFGeoShape::Type> MapCodeFaceShape;

  /// Definition of a Map from a certain code to a shape of a Edge
  typedef Common::CFMap<CFuint, CFGeoShape::Type> MapCodeEdgeShape;

  /// Function to access the unique instance of this singleton
  static LocalConnectionData& getInstance();

  /// Get the localConnectivity table face-dof
  /// A code is obtained by a weighted sum of order, shape and the enumerator
  /// indicating the type of connectivity requested.
  /// @parameter shape      shape of the element
  /// @parameter geomOrder  order of the geometrical representation of the element
  /// @return the table corresponding to the requested local connectivity for
  ///         an element having CFPolyOrder::Type = geomOrder and CFGeoShape::Type = shape
  /// @pre (shape < 100) && (dofType < 99) to avoid having
  ///      two different codes corresponding to the same combination
  Common::Table<CFuint>* getFaceDofLocal(const CFGeoShape::Type& shape,
                 const CFPolyOrder::Type& geomOrder,
                 const CFDofType& dofType,
                 const CFPolyForm::Type& geomPolyType);

  /// Get the localConnectivity table edge-dof
  /// A code is obtained by a weighted sum of order, shape and the enumerator
  /// indicating the type of connectivity requested.
  /// @parameter shape      shape of the element
  /// @parameter geomOrder  order of the geometrical representation of the element
  /// @return the table corresponding to the requested local connectivity for
  ///         an element having CFPolyOrder::Type = geomOrder and CFGeoShape::Type = shape
  /// @pre (shape < 100) && (dofType < 99) to avoid having
  ///      two different codes corresponding to the same combination
  Common::Table<CFuint>* getEdgeDofLocal(const CFGeoShape::Type& shape,
                 const CFPolyOrder::Type& geomOrder,
                 const CFDofType& dofType,
                 const CFPolyForm::Type& geomPolyType);

  ///  Get the code
  CFuint getCode(const CFGeoShape::Type& shape,
                 const CFPolyOrder::Type& polyOrder,
                 const CFDofType& dofType,
                 const CFPolyForm::Type& polyType);

  /// Get the shape of the given face in the given element shape
  CFGeoShape::Type getFaceShape(const CFGeoShape::Type& shape, const CFuint& iFace);

  /// Get the shape of the given edge in the given element shape
  CFGeoShape::Type getEdgeShape(const CFGeoShape::Type& shape, const CFuint& iEdge);

  /// Get the number of faces in the given element shape
  CFuint getNbFacesInShape(const CFGeoShape::Type& shape);

  /// Get the number of edges in the given element shape
  CFuint getNbEdgesInShape(const CFGeoShape::Type& shape);

  /// Set the number of faces per element
  /// @param nbFacesPerElem  list of the number of faces per element
  void setNbFacesPerElement(std::valarray<CFuint>& nbFacesPerElem);

  /// Set the number of edges per element
  /// @param nbFacesPerElem  list of the number of edges per element
  void setNbEdgesPerElement(std::valarray<CFuint>& nbEdgesPerElem);

  /// Set the shape of the face shapes per element type
  void setFaceShapesPerElemType(std::vector< std::vector<CFGeoShape::Type> >& faceShapesPerElemType);

  /// Set the shape of the edge shapes per element type
  void setEdgeShapesPerElemType(std::vector< std::vector<CFGeoShape::Type> >& edgeShapesPerElemType);

  /// Print all member maps
  void print();

private: // helper function

  /// Get the code to look up in _code2FaceShape
  static CFuint getFaceCode(const CFGeoShape::Type& shape,
                            const CFuint& iFace);

  ///  Get the code to look up in _code2EdgeShape
  static CFuint getEdgeCode(const CFGeoShape::Type& shape,
                            const CFuint& iEdge);

  /// Prepare the map _code2FaceShape
  void prepareMapCode2FaceShape();

  /// Prepare the map _code2EdgeShape
  void prepareMapCode2EdgeShape();

private: // private constructors and destructor for singleton object

  /// Private constructor
  LocalConnectionData();

  /// Private destructor
  ~LocalConnectionData();

private: // data

  /// map code to the corresponding connectivity table
  MapCodeTable _code2TableFace;

  /// map code to the corresponding connectivity table
  MapCodeTable _code2TableEdge;

  /// map code to the face shape
  MapCodeFaceShape _code2FaceShape;

  /// map code to the edge shape
  MapCodeEdgeShape _code2EdgeShape;

}; // end of class LocalConnectionData

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Framework_LocalConnectionData_hh
