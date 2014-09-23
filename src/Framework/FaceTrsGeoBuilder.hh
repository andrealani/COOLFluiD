// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Framework_FaceTrsGeoBuilder_hh
#define COOLFluiD_Framework_FaceTrsGeoBuilder_hh

//////////////////////////////////////////////////////////////////////////////

#include "Common/NonCopyable.hh"
#include "Framework/GeometricEntity.hh"
#include "Framework/Storage.hh"
#include "Framework/DataSocketSink.hh"
#include "Framework/State.hh"
#include "Framework/Node.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework { class TopologicalRegionSet; }

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

/// This class represents a builder of a standard GeometricEntity,
/// with Node's and State's starting from some TRS-based data. It is
/// not supposed to be overridden.
/// This builder allows to create, use and release only one GeometricEntity
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
class Framework_API FaceTrsGeoBuilder {
public:

  /// This nested struct represents groups the data needed by this builder.
  /// It must be non copyable to force the client code to use by reference
  /// the data aggregated by this FaceTrsGeoBuilder
  /// @see TopologicalRegionSet
  /// @author Andrea Lani
  struct GeoData : public Common::NonCopyable<GeoData> {

    /// Default constructor
    GeoData() {}

    /// pointer to TRS
    Common::SafePtr<Framework::TopologicalRegionSet> trs;

    /// flag telling if the face is on the boundary
    bool isBFace;

    /// geo index in TRS
    CFuint idx;
  };

  /// Constructor
  FaceTrsGeoBuilder();

  /// Destructor
  ~FaceTrsGeoBuilder();

  /// Sets the DataSockets
  void setDataSockets
  (Framework::DataSocketSink<Framework::State*, Framework::GLOBAL> statesSocket,
    Framework::DataSocketSink<Framework::State*> gstatesSocket,
  Framework::DataSocketSink<Framework::Node*, Framework::GLOBAL> nodesSocket);

  /// Set up the pool
  void setup();

  /// Unset up the pool
  void unsetup();

  /// Get the data of the GeometricEntity builder.
  /// This allows the client code to set the data and then
  /// let the builder work on its own updated data.
  FaceTrsGeoBuilder::GeoData& getDataGE()
  {
    return m_data;
  }

  /// Build the GeometricEntity corresponding to the given local ID
  /// in the correspondng TopologicalRegionSet
  Framework::GeometricEntity* buildGE();

  /// Release the buildd GeometricEntity's to make them again
  /// available for creation of new ones
  void releaseGE()
  {
    m_builtGeo = CFNULL;
  }

private:

  /// data of this builder
  FaceTrsGeoBuilder::GeoData  m_data;

  /// socket for ghost States
  Framework::DataSocketSink<Framework::State*> socket_gstates;

  /// socket for nodes
  Framework::DataSocketSink < Framework::Node* , Framework::GLOBAL > socket_nodes;

  /// socket for States
  Framework::DataSocketSink < Framework::State* , Framework::GLOBAL > socket_states;

  /// GeometricEntity's pool ordered by GeoType
  std::vector<Framework::GeometricEntity*> m_poolData;

  /// built GeometricEntity
  Framework::GeometricEntity* m_builtGeo;

  ///flag to known if the builder has already been setup
  bool m_isSetup;

  ///flag to known if the socket have already been set
  bool m_isSocketsSet;

}; // end of class FaceTrsGeoBuilder

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Framework_FaceTrsGeoBuilder_hh
