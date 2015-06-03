// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Framework_ElementTypeData_hh
#define COOLFluiD_Framework_ElementTypeData_hh

//////////////////////////////////////////////////////////////////////////////

#include "Common/StringOps.hh"
#include "Framework/CFGeoShape.hh"
#include "Framework/LocalConnectionData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

/// This class represents the the data related to an ElementType
/// @author Tiago Quintino
class Framework_API ElementTypeData {

public: // functions

  /// Default constructor without arguments
  ElementTypeData();

  /// Default destructor
  ~ElementTypeData();

public: // accessors

  /// @return _nameShape
  std::string getShape() const { return _nameShape; }

  /// @return _geoShape
  CFGeoShape::Type getGeoShape() const  {  return _geoShape; }

  /// @return _nbElemsTot
  CFuint getNbTotalElems() const  {  return _nbElemsTot; }
  
  /// @return _nbElems
  CFuint getNbElems() const  {  return _nbElems; }

  /// @return _startIdx
  CFuint getStartIdx() const { return _startIdx; }

  /// @return _endIdx
  CFuint getEndIdx() const  { return _startIdx + _nbElems; }

  /// @return number of faces
  CFuint getNbFaces() const  {  return m_nbfaces;  }

  /// @return _nbNodes
  CFuint getNbNodes() const  { return _nbNodes; }

  /// @return _nbStates
  CFuint getNbStates() const {  return _nbStates; }

  /// @return _geoOrder
  CFuint getGeoOrder() const { return _geoOrder; }

  /// @return _solOrder
  CFuint getSolOrder() const  { return _solOrder; }

public: // mutators

  /// Sets _nameShape
  void setShape(const std::string& nameShape) { _nameShape = nameShape; }

  /// Sets _geoShape
  void setGeoShape(const CFGeoShape::Type& geoShape)
  {
    _geoShape = geoShape;
    m_nbfaces = LocalConnectionData::getInstance().getNbFacesInShape(_geoShape);
  }
  
  /// Sets total number of elements
  void setNbTotalElems(const CFuint nbElemsTot) {  _nbElemsTot = nbElemsTot; }
  
  /// Sets _nbElems
  void setNbElems(const CFuint nbElems) {  _nbElems = nbElems; }

  /// Sets _startIdx
  void setStartIdx(const CFuint startIdx) { _startIdx = startIdx; }

  /// Sets _nbNodes
  void setNbNodes(const CFuint nbNodes) { _nbNodes = nbNodes; }

  /// Sets _nbStates
  void setNbStates(const CFuint nbStates) { _nbStates = nbStates; }

  /// Sets _geoOrder
  void setGeoOrder(const CFuint geoOrder) { _geoOrder = geoOrder; }

  /// Sets _solOrder
  void setSolOrder(const CFuint solOrder) { _solOrder = solOrder; }

private: // data

  /// the string identifying the shape of this type
  std::string   _nameShape;
  /// the CFGeoShape::Type corresponding to the shape
  CFGeoShape::Type _geoShape;
  /// the total (across processes) number of elements with this type
  CFuint  _nbElemsTot;
  /// the local (in current process) number of elements with this type
  CFuint  _nbElems;
  /// the start index on the ordered vector of elements for this type
  CFuint _startIdx;
  /// the number of nodes in this element type
  CFuint _nbNodes;
  /// the number of states in this element type
  CFuint _nbStates;
  /// the  geometric order of this element
  CFuint _geoOrder;
  /// the  solution order of this element
  CFuint _solOrder;
  /// number of faces
  CFuint m_nbfaces;

}; // end of class ElementTypeData

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Framework_ElementTypeData_hh
