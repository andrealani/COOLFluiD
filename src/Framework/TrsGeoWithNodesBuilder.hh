// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Framework_TrsGeoWithNodesBuilder_hh
#define COOLFluiD_Framework_TrsGeoWithNodesBuilder_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/GeometricEntity.hh"
#include "Storage.hh"
#include "Common/NonCopyable.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

    class TopologicalRegionSet;

//////////////////////////////////////////////////////////////////////////////

/// This class represents a builder of GeometricEntity's with only Node's,
/// starting from some TRS-based data. It is not supposed to be
/// overridden.
/// This builder allows to create, use and release Only one GeometricEntity
/// at a time.
/// @see GeometricEntity
/// @see Node
/// @see State
/// @see ConnectivityTable
/// @see TopologicalRegionSet
/// @see TopologicalRegion
/// @see MeshData
/// @author Andrea Lani
/// @author Tiago Quintino
class Framework_API TrsGeoWithNodesBuilder {
public:

  /// This nested struct represents groups the data needed by this builder.
  /// It must be non copyable to force the client code to use by reference
  /// the data aggregated by this TrsGeoWithNodesBuilder
  /// @see TopologicalRegionSet
  /// @author Andrea Lani
  /// @author Tiago Quintino
  struct GeoData : public Common::NonCopyable<GeoData> {
    /// Default constructor
    GeoData() {}

    /// pointer to TRS
    Common::SafePtr<TopologicalRegionSet> trs;

    /// geo index in TRS
    CFuint idx;
  };

  /// Constructor
  TrsGeoWithNodesBuilder();

  /// Destructor
  ~TrsGeoWithNodesBuilder();

  /// Set up the pool
  void setup();
  
  /// Unset up the pool
  void unsetup();
  
  /// Set up the pool to use on another namespace
  void setupInNamespace(const std::string& namespaceName);

  /// Get the data of the GeometricEntity builder.
  /// This allows the client code to set the data and then
  /// let the builder work on its own updated data.
  TrsGeoWithNodesBuilder::GeoData& getDataGE() { return m_data; }

  /// Build the GeometricEntity corresponding to the given local ID
  /// in the correspondng TopologicalRegionSet
  GeometricEntity* buildGE();

  /// Release the buildd GeometricEntity's to make them again
  /// available for creation of new ones
  void releaseGE() { m_builtGeo = CFNULL; }

private:

  /// data of this builder
  TrsGeoWithNodesBuilder::GeoData m_data;

  /// handle to the Node's storage
  DataHandle<Node*,GLOBAL> m_nodes;

  /// GeometricEntity's pool ordered by GeoType
  std::vector<GeometricEntity*> m_poolData;

  /// built GeometricEntity's
  GeometricEntity* m_builtGeo;

  ///flag to known if the builder has already been setup
  bool m_isSetup;

}; // end of class TrsGeoWithNodesBuilder

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Framework_TrsGeoWithNodesBuilder_hh
