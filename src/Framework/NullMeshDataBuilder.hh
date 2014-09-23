// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Framework_NullMeshDataBuilder_hh
#define COOLFluiD_Framework_NullMeshDataBuilder_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/MeshDataBuilder.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

/// This class offers an abstract interface for NullMeshDataBuilder's.
/// ReadCFmesh aggregates the NullMeshDataBuilder.
/// @author Andrea Lani
class Framework_API NullMeshDataBuilder : public Framework::MeshDataBuilder {

public:

  /// Defines the Config Option's of this class
  /// @param options a OptionList where to add the Option's
  static void defineConfigOptions(Config::OptionList& options);
  
  /// Constructor
  /// @param name name of the builder used for configuration
  NullMeshDataBuilder(const std::string& name);
  
  /// Virtual destructor
  virtual ~NullMeshDataBuilder();
  
  /// Creates all the TopologicalRegionSet's from CFmeshData
  virtual void createTopologicalRegionSets();
  
  /// Releases temporary memory used in building the mesh
  /// Should be called by child classes
  virtual void releaseMemory();
  
protected: // methods

  /// Set the max number of states in cell
  virtual void setMaxNbStatesInCell();

  /// Set the max number of nodes in cell
  virtual void setMaxNbNodesInCell();
  
  /// Set the max number of faces in cell
  virtual void setMaxNbFacesInCell();
  
private: // methods
  
  /// Create and set the mapping between faces and TRSs.
  /// Allows to get Face info from TRS by geometric entity ID,
  /// local to the processor.
  virtual void setMapGeoToTrs();
  
}; // end of class NullMeshDataBuilder

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

Framework_Factory(NullMeshDataBuilder) // define the factoty instance

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Framework_NullMeshDataBuilder_hh
