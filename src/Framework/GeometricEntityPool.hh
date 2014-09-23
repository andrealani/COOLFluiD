// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Framework_GeometricEntityPool_hh
#define COOLFluiD_Framework_GeometricEntityPool_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/GeometricEntity.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

/// This class represents a base geometric entity pool that stores
/// GeometricEntity's (Cell, Face, Edge, ...) and return them on
/// demand. The actual job is delegated to an aggregated policy class,
/// corresponding to the template parameter GE_BUILDER.
/// GeometricEntity's are not allocated on the fly, but they are just
/// set up using connectivity info stored for instance in ConnectivityTable's
/// acquainted by TopologicalRegionSet's and aggregated by the MeshData.
/// GeometricEntityPool is meant NOT to be overridable, but to be used in
/// combination with appropriate statically bound policies, chosen by the
/// client code (Method, NumericalCommand or other class) that can select
/// therefore how to build ad-hoc GeometricEntity's in the most efficient
/// way.
/// The building steps for a GeometricEntity are the following:
/// 1. getDataGE() to get and set the GE_BUILDER::GeoData used by the builder
/// 2. buildGE()   to build the GeometricEntity based on current GE_BUILDER::GeoData
/// 3. releaseGE() to release data and make the GeometricEntity's pool fully
///                available again
/// @see GeometricEntity
/// @see Node
/// @see State
/// @see ConnectivityTable
/// @see TopologicalRegionSet
/// @see MeshData
/// @author Andrea Lani
/// @author Tiago Quintino
template <class GEOBUILDER>
class GeometricEntityPool {
public:

  /// Constructor
  GeometricEntityPool() : m_geBuilder() {}

  /// Default destructor (this class should not be overriden)
  ~GeometricEntityPool() {}

  /// Set up the pool
  void setup()   {  m_geBuilder.setup(); }

  /// Set up the pool
  void unsetup() {  m_geBuilder.unsetup(); }

  /// Set up the pool in a different namespace
  void setupInNamespace(const std::string& namespaceName)
  {
    m_geBuilder.setupInNamespace(namespaceName);
  }

  /// Get the data of the GeometricEntity builder.
  /// This allows the client code to set the data and then
  /// let the builder work on its own updated data.
  typename GEOBUILDER::GeoData& getDataGE() { return m_geBuilder.getDataGE(); }

  /// Build the GeometricEntity corresponding to the given local ID
  /// in the correspondng TopologicalRegionSet
  GeometricEntity* buildGE() {  return m_geBuilder.buildGE(); }

  /// Release the build GeometricEntity's to make them again
  /// available for creation of new ones
  void releaseGE() { m_geBuilder.releaseGE(); }

  /// Get a SafePtr to the GeoBuilder
  Common::SafePtr<GEOBUILDER> getGeoBuilder() { return &m_geBuilder; }

private: //data

  /// Geometric entity builder
  GEOBUILDER m_geBuilder;

}; // end of class GeometricEntityPool

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Framework_GeometricEntityPool_hh
