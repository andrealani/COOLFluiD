// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Framework_StdTrsMultiGeoBuilder_hh
#define COOLFluiD_Framework_StdTrsMultiGeoBuilder_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/GeometricEntity.hh"
#include "Storage.hh"
#include "Common/NonCopyable.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

    class TopologicalRegionSet;

//////////////////////////////////////////////////////////////////////////////

/// This class represents a builder of standard geometric entities,
/// starting from some TRS-based data. It is not supposed to be
/// overridden
/// @see GeometricEntity
/// @see Node
/// @see State
/// @see ConnectivityTable
/// @see TopologicalRegionSet
/// @see TopologicalRegion
/// @see MeshData
/// @author Andrea Lani
/// @author Tiago Quintino
class Framework_API StdTrsMultiGeoBuilder {
public:

  /// This nested struct represents groups the data needed by this builder.
  /// It must be non copyable to force the client code to use by reference
  /// the data aggregated by this StdTrsMultiGeoBuilder
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
  StdTrsMultiGeoBuilder();

  /// Destructor
  ~StdTrsMultiGeoBuilder();

  /// Set up the pool
  void setup();

  /// Get the data of the GeometricEntity builder.
  /// This allows the client code to set the data and then
  /// let the builder work on its own updated data.
  StdTrsMultiGeoBuilder::GeoData& getDataGE() {  return _data; }

  /// Build the GeometricEntity corresponding to the given local ID
  /// in the correspondng TopologicalRegionSet
  GeometricEntity* buildGE();

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

private:

  /// data of this builder
  StdTrsMultiGeoBuilder::GeoData  _data;

  /// handle to the State's storage
  DataHandle<State*,GLOBAL>  _states;

  /// handle to the Node's storage
  DataHandle<Node*,GLOBAL>   _nodes;

  /// GeometricEntity's pool ordered by GeoType
  std::vector<std::vector<GeometricEntity*> > _poolData;

  /// count the GeometricEntity's built per type
  std::valarray<CFuint> _countGeo;

  /// list of the built GeometricEntity's
  std::vector<GeometricEntity*> _builtGeos;

  ///flag to known if the builder has already been setup
  bool _isSetup;

}; // end of class StdTrsMultiGeoBuilder

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Framework_StdTrsMultiGeoBuilder_hh
