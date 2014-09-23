// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_TecplotWriter_MapGeoEntToTecplot_hh
#define COOLFluiD_TecplotWriter_MapGeoEntToTecplot_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/MapGeoEntToFormat.hh"

#include "TecplotWriter/TecplotWriterAPI.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace TecplotWriter {

//////////////////////////////////////////////////////////////////////////////

/// This class kowns how to map the internal COOLFluiD representaton of
/// a GeometricEntity's to the Tecplot file format.
/// @author Tiago Quintino
class TecplotWriter_API MapGeoEntToTecplot : public Framework::MapGeoEntToFormat {
public:

  /// Default destructor.
  virtual
  ~MapGeoEntToTecplot();

  /// Identfies the Geomeric string entity by its parameters.
  virtual
  std::string identifyGeoEnt(const Framework::GeoEntityInfo& geoinfo);

  /// Writes the element connectivity into the Tecplot format.
  /// Useful to circunvent the fact that different formats
  /// have different element numberings.
  virtual
  void writeGeoEntConn(std::ofstream& file,
  			               std::valarray<CFuint>& stateIDs,
  			               const Framework::GeoEntityInfo& geoinfo);

  /// Computes in how many sub entities the given GeometricEntity
  /// should be partitioned for representation in the Tecplot format.
  /// If the format supports high-order entities, then the result will probably be 1,
  /// meaning that no sub entities are needed.
  virtual
  CFuint computeNbSubEntities(const Framework::GeoEntityInfo& geoinfo);
  
  /** 
   * Gets the Class name
   */
  static std::string getClassName()
  {
    return "MapGeoEntToTecplot";
  }

}; // end of class MapGeoEntToTecplot

//////////////////////////////////////////////////////////////////////////////

  } // namespace TecplotWriter

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_TecplotWriter_MapGeoEntToTecplot_hh
