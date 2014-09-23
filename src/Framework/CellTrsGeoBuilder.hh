// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Framework_CellTrsGeoBuilder_hh
#define COOLFluiD_Framework_CellTrsGeoBuilder_hh

//////////////////////////////////////////////////////////////////////////////

#include "Common/NonCopyable.hh"
#include "Framework/GeometricEntity.hh"
#include "Framework/Storage.hh"
#include "Framework/DataSocketSink.hh"
#include "Framework/State.hh"
#include "Framework/Node.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {
    class TopologicalRegionSet;
    class MapGeoToTrsAndIdx;
  }

  namespace Common {
    template <class T> class ConnectivityTable;
  }

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

/// This class represents a builder of a standard GeometricEntity,
/// with Node's and State's starting from some TRS-based data. It is
/// not supposed to be overridden.
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
class Framework_API CellTrsGeoBuilder {
public:

  /// This nested struct represents groups the data needed by this builder.
  /// It must be non copyable to force the client code to use by reference
  /// the data aggregated by this CellTrsGeoBuilder
  /// @see TopologicalRegionSet
  /// @author Andrea Lani
  struct GeoData : public Common::NonCopyable<GeoData> {

    /// Default constructor
    GeoData() {}

    /// pointer to TRS
    Common::SafePtr<Framework::TopologicalRegionSet> trs;

    /// geo index in TRS
    CFuint idx;
  };

  /// Constructor
  CellTrsGeoBuilder();

  /// Destructor
  virtual ~CellTrsGeoBuilder();

  /// Sets the DataSockets
  void setDataSockets
  (Framework::DataSocketSink<Framework::State*,Framework::GLOBAL> statesSocket,
    Framework::DataSocketSink<Framework::State*> gstatesSocket,
  Framework::DataSocketSink<Framework::Node*,Framework::GLOBAL> nodesSocket);

  /// Set up the pool
  void setup();
  
  /// Unsetup this GeoBuilder, by deallocating the pool
  void unsetup();
  
  /// Set up the pool
  void setupInNamespace(const std::string& namespaceName);

  /// Get the data of the GeometricEntity builder.
  /// This allows the client code to set the data and then
  /// let the builder work on its own updated data.
  CellTrsGeoBuilder::GeoData& getDataGE()
  {
    return _data;
  }

  /// Build the GeometricEntity corresponding to the given local ID
  /// in the correspondng TopologicalRegionSet
  Framework::GeometricEntity* buildGE();

  /// Release the buildd GeometricEntity's to make them again
  /// available for creation of new ones
  void releaseGE()
  {
    // free all the geometric entities
    for (CFuint i = 0; i < _builtGeos.size(); ++i) {
      _builtGeos[i]->release();
    }
    _builtGeos.clear();
    _countGeo = 0;
  }
  
protected:
  
  /// data of this builder
  CellTrsGeoBuilder::GeoData _data;

  /// socket for ghost States
  Framework::DataSocketSink<Framework::State*> socket_gstates;

  /// socket for nodes
  Framework::DataSocketSink < Framework::Node* , Framework::GLOBAL > socket_nodes;

  /// socket for States
  Framework::DataSocketSink < Framework::State* , Framework::GLOBAL > socket_states;

  /// pointer to the MapGeoToTrsAndIdx for the faces
  Common::SafePtr<Framework::MapGeoToTrsAndIdx> _mapGeoToTrs;

  /// pointer to the connectivity cell-faces
  Common::SafePtr<Common::ConnectivityTable<CFuint> > _cellFaces;

  /// GeometricEntity's pool ordered by GeoType
  std::vector<std::vector<Framework::GeometricEntity*> > _poolData;

  /// count the GeometricEntity's built per type
  std::valarray<CFuint> _countGeo;

  /// list of the built GeometricEntity's
  std::vector<Framework::GeometricEntity*> _builtGeos;

  ///flag to known if the builder has already been setup
  bool _isSetup;

  ///flag to known if the sockets of the builder have already been setup
  bool _isSocketsSet;

}; // end of class CellTrsGeoBuilder

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Framework_CellTrsGeoBuilder_hh
