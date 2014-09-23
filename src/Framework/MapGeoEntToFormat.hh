// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Framework_MapGeoEntToFormat_hh
#define COOLFluiD_Framework_MapGeoEntToFormat_hh

//////////////////////////////////////////////////////////////////////////////

#include <valarray>

#include "Common/COOLFluiD.hh"
#include "Common/NonCopyable.hh"
#include "Common/StringOps.hh"

#include "Framework/CFPolyOrder.hh"
#include "Framework/Framework.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

/// This strucutre aglomerates the information needed to uniquely
/// define an entity in a foreign format
struct Framework_API GeoEntityInfo
{
  CFuint nbNodes;
  CFuint nbStates;
  CFuint geoOrder;
  CFuint solOrder;
  CFuint dimension;
};

//////////////////////////////////////////////////////////////////////////////

/// This class kowns how to map the internal COOLFluiD representaton of
/// a GeometricEntity's to a certain file format.
/// @author Tiago Quintino
class Framework_API MapGeoEntToFormat : public Common::NonCopyable<MapGeoEntToFormat>
{
public: // interface functions

  /// Default destructor.
  virtual
  ~MapGeoEntToFormat();

  /// Identfies the Geomeric string entity by its parameters.
  virtual
  std::string identifyGeoEnt(const GeoEntityInfo& geoinfo) = 0;

  /// Writes the element connectivity into the concrete format.
  /// Useful to circunvent the fact that different formats
  /// have different element numberings.
  virtual
  void writeGeoEntConn(std::ofstream& file,
    		               std::valarray<CFuint>& stateIDs,
    		               const GeoEntityInfo& geoinfo) = 0;

  /// Computes in how many sub entities the given GeometricEntity
  /// should be partitioned for representation in the foreign format.
  /// If the format supports high-order entities, then the result will probably be 1,
  /// meaning that no sub entities are needed.
  virtual
  CFuint computeNbSubEntities(const GeoEntityInfo& geoinfo) = 0;

}; // end of class MapGeoEntToFormat

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Framework_MapGeoEntToFormat_hh
