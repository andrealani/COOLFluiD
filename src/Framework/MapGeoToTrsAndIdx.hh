// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Framework_MapGeoToTrsAndIdx_hh
#define COOLFluiD_Framework_MapGeoToTrsAndIdx_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/TopologicalRegionSet.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

/// This class provides a mapping from a GeometricEntity to the corresponding
/// TopologicalRegionSet and local index in the TopologicalRegionSet.
/// It also gives info about the geo being a bounday geo or not.
/// @author Andrea Lani
/// @author Tiago Quintino
class Framework_API MapGeoToTrsAndIdx {

public: // functions

  /// Constructor
  MapGeoToTrsAndIdx();

  /// Destructor
  ~MapGeoToTrsAndIdx();

  /// Resize the mapper
  void resize(const CFuint totalNbGeos);

  /// Get the a pointer to the TopologicalRegionSet where
  /// the GeometricEntity having geoID is
  Common::SafePtr<TopologicalRegionSet> getTrs(CFuint geoID) const
  {
    cf_assert(geoID < m_trs.size());
    return m_trs[geoID];
  }

  /// Get the local GeometricEntity index in the TopologicalRegionSet
  CFuint getIdxInTrs(CFuint geoID) const
  {
    cf_assert(geoID < m_idxInTrs.size());
    return m_idxInTrs[geoID];
  }

  /// Tells if the GeometricEntity is on the boundary
  bool isBGeo(CFuint geoID) const
  {
    cf_assert(geoID < m_isBGeo.size());
    return m_isBGeo[geoID];
  }

  /// Set the TopologicalRegionSet pointer, the local index in the
  /// TopologicalRegionSet and the flag to tell if the GeometricEntity
  /// geo is on the boundary
  void setMappingData(CFuint geoID,
          Common::SafePtr<TopologicalRegionSet> trs,
          CFuint idxInTrs,
          bool isBGeo);

  /// Deallocate
  void deallocate();

private: // data

  /// TopologicalRegionSet pointer for each GeometricEntity
  std::vector< Common::SafePtr<TopologicalRegionSet> > m_trs;

  /// local GeometricEntity index in the TopologicalRegionSet
  std::vector<CFuint> m_idxInTrs;

  /// flag telling if the GeometricEntity is on the boundary
  std::vector<bool> m_isBGeo;

}; // end of class MapGeoToTrsAndIdx

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Framework_MapGeoToTrsAndIdx_hh
