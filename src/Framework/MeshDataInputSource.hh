// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef MESHDATAINPUTSOURCE_HH
#define MESHDATAINPUTSOURCE_HH

//////////////////////////////////////////////////////////////////////////////

#include "Common/COOLFluiD.hh"
#include "Common/SafePtr.hh"

#include "Framework/CFPolyOrder.hh"
#include "Framework/ElementTypeData.hh"
#include "Framework/CFGeoEnt.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

    class TopologicalRegionSet;

//////////////////////////////////////////////////////////////////////////////

/// This class defines an abstract interface to read a mesh
/// @todo maybe have a second interface deriving from this one,
///       adding reference functions (direct access to the data without copy)
class Framework_API MeshDataInputSource {
public: // typedefs

  /// Helper type to return a node
  typedef std::valarray<CFreal> NodeDataType;

  /// Helper type to return a state
  typedef std::valarray<CFreal> StateDataType;

  /// Helper type to return an element
  typedef std::pair<std::valarray<CFuint>,
                     std::valarray<CFuint> > ElementDataType;

public: // general functions

  /// @return the space dimension
  virtual CFuint getDimension () const = 0;

  /// @return the state dimension
  virtual CFuint getNbEquations () const = 0;

  /// @return the number of nodes in the mesh
  virtual CFuint getNbNodes () const = 0;

  /// @return the number of states in the mesh
  virtual CFuint getNbStates () const = 0;

  /// @return the number of elements
  virtual CFuint getNbElements () const = 0;

  /// @return the number of element types
  virtual CFuint getNbElementTypes () const = 0;

  /// @return the geometric poly order
  virtual CFPolyOrder::Type getGeometricPolyOrder () const = 0;

  /// @return the TRS geom type
  virtual CFGeoEnt::Type getGeomType (CFuint TRS) const = 0;

  /// @return the solution poly order
  virtual CFPolyOrder::Type getSolutionPolyOrder () const = 0;

  /// determine if a solution is included
  virtual bool isWithSolution () const = 0;

  /// @return the number of TRS
  virtual CFuint getNbTRSs () const = 0;

  /// @return the number of TRs in a TRS
  virtual CFuint getNbTRs (CFuint TRS) const = 0;

  /// @return the name of a TRS
  virtual std::string getNameTRS (CFuint TRS) const = 0;

  /// @return the number of geometric entities for a TRS and TR
  virtual CFuint getNbGeomEnts(CFuint TRS, CFuint TR) const = 0;

  /// @return the element type data
  virtual void getElementType (CFuint NbElementType, ElementTypeData & Out) const =0;

  /// @return the geometric entity (ElementDataType) in OUT
  virtual void getGeoElement (Common::SafePtr<TopologicalRegionSet> trs,
                               CFuint TR,
                               CFuint Geom,
                               ElementDataType & Out) const = 0;

  /// @return the  node data (NodeDataType)
  virtual void getNode (CFuint NbNode, NodeDataType & Out) const = 0;

  /// @return the global ID of the node
  virtual CFuint getGlobalNodeID (CFuint NbNode) const = 0;

  /// @return the state data (StateDataType)
  virtual void getState (CFuint NbState, StateDataType & Out) const = 0;

  /// @return the global state ID
  virtual CFuint getGlobalStateID (CFuint State) const = 0;

  /// @return the element
  virtual void getElement (CFuint Ele, ElementDataType & Out) const = 0;

  /// @return the global element ID
  virtual CFuint getGlobalElementID (CFuint Element) const = 0;

  virtual std::vector< Common::SafePtr<TopologicalRegionSet> > getTrsList() const = 0;

  /// @return the array of global geom counts for TRS
  virtual const std::vector<std::vector<CFuint> > & getTotalTRSInfo () const = 0;

  /// @return array of names for the TRSs in TotalTRSInfo
  virtual const std::vector<std::string> & getTotalTRSNames () const = 0;

  virtual ~MeshDataInputSource ();
};

//////////////////////////////////////////////////////////////////////////////

    }
}

//////////////////////////////////////////////////////////////////////////////

#endif
