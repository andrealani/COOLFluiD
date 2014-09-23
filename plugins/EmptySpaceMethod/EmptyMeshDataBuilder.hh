// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_EmptySpaceMethod_EmptyMeshDataBuilder_hh
#define COOLFluiD_EmptySpaceMethod_EmptyMeshDataBuilder_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/MeshDataBuilder.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace EmptySpaceMethod {

//////////////////////////////////////////////////////////////////////////////

/// This class builds Empty solver data inside MeshData
/// @see MeshDataBuilder
/// @author Tiago Quintino
/// @author Pedro Maciel
class EmptyMeshDataBuilder : public Framework::MeshDataBuilder {

public: // functions

  /// Constructor
  /// @param name name of the builder used for configuration
  EmptyMeshDataBuilder(const std::string& name);

  /// Destructor
  ~EmptyMeshDataBuilder();

protected: // functions

  /// Set the max number of states in cell
  virtual void setMaxNbStatesInCell();

  /// Set the max number of nodes in cell
  virtual void setMaxNbNodesInCell();

  /// Set the max number of faces in cell
  virtual void setMaxNbFacesInCell();

  /// Create and set the mapping between faces and TRSs
  virtual void setMapGeoToTrs();

private: // data

};  // end of class EmptyMeshDataBuilder

//////////////////////////////////////////////////////////////////////////////

  }  // namespace EmptySpaceMethod

}  // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif  // COOLFluiD_EmptySpaceMethod_EmptyMeshDataBuilder_hh

