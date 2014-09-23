// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Framework_StdTrsGeoBuilder_hh
#define COOLFluiD_Framework_StdTrsGeoBuilder_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/GeometricEntity.hh"
#include "Storage.hh"
#include "Common/NonCopyable.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

    class TopologicalRegionSet;

//////////////////////////////////////////////////////////////////////////////

///  * This class represents a builder of a standard GeometricEntity,
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
/// @author Tiago Quintino
class Framework_API StdTrsGeoBuilder {
public:

  /// This nested struct represents groups the data needed by this builder.
  /// It must be non copyable to force the client code to use by reference
  /// the data aggregated by this StdTrsGeoBuilder
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
  StdTrsGeoBuilder();

  /// Destructor
  ~StdTrsGeoBuilder();

  /// Setup this GeoBuilder, by allocating the pool
  void setup();

  /// Unsetup this GeoBuilder, by deallocating the pool
  void unsetup();

  /// Set up the pool to use on another namespace
  void setupInNamespace(const std::string& namespaceName);

  /// Get the data of the GeometricEntity builder.
  /// This allows the client code to set the data and then
  /// let the builder work on its own updated data.
  StdTrsGeoBuilder::GeoData& getDataGE() {  return _data;  }

  /// Build the GeometricEntity corresponding to the given local ID
  /// in the correspondng TopologicalRegionSet
  GeometricEntity* buildGE();

  /// Release the buildd GeometricEntity's to make them again
  /// available for creation of new ones
  void releaseGE() { _builtGeo = CFNULL; }

private:

  /// data of this builder
  StdTrsGeoBuilder::GeoData  _data;

  /// handle to the State's storage
  DataHandle<State*,GLOBAL>  _states;

  /// handle to the Node's storage
  DataHandle<Node*,GLOBAL>   _nodes;

  /// GeometricEntity's pool ordered by GeoType
  std::vector<GeometricEntity*> _poolData;

  /// built GeometricEntity
  GeometricEntity* _builtGeo;

  ///flag to known if the builder has already been setup
  bool _isSetup;
  
  /// passes in setup
  CFuint pass;

  /// passes in unsetup
  CFuint unpass;

}; // end of class StdTrsGeoBuilder

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Framework_StdTrsGeoBuilder_hh
