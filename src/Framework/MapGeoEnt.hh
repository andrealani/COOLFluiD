// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Framework_MapGeoEnt_hh
#define COOLFluiD_Framework_MapGeoEnt_hh

//////////////////////////////////////////////////////////////////////////////

#include <valarray>

#include "Common/COOLFluiD.hh"
#include "Common/StringOps.hh"
#include "Framework/Framework.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

/// This class represents an entity that kowns how to identify Geometric
/// Entities from their.
/// In the future this might work as a factory of know types for
/// Geometric Entities.
/// @todo this class could be better implemented with as abstract class
///       and with concrete derived classes directly instanciated by the clients
///       because inthis case the clients do know the concrete type.
///       e.g. TecplotMapGeoEnt, CFmeshMapGeoEnt, etc...
/// @todo there is missing documentation in this class
/// @author Tiago Quintino
class Framework_API MapGeoEnt {
public:

  /// Identfies the Geomeric string entity by its parameters.
  /// This should be used by the CFmesh formaters
  static std::string identifyGeoEnt(const CFuint nbNodes,
    const CFuint geoOrder,
    const CFuint dim);

  /// Identfies the Geomeric string entity by its parameters.
  /// This should be used by the Tecplot formatters
  static std::string identifyGeoEntTecplot(const CFuint nbNodes,
    const CFuint geoOrder,
    const CFuint dim);

  /// Writes the element connectivity into Tecplot format.
  /// This is used to circunvent the fact that Tecplot
  /// only knows about Tetrahedra and Bricks (Hexahedra) in 3D
  static void writeTecplotGeoEntConn(std::ofstream& file,
    		     std::valarray<CFuint>& nodeIDs,
    		     const CFuint geoOrder,
    		     const CFuint dim);

  /// Writes the element connectivity into OpenDX format.
  /// This is used to circunvent the fact that OpenDX
  /// only knows about Cubes in 3D
  static void writeOpenDXGeoEntConn(std::ofstream& file,
    		    std::valarray<CFuint>& nodeIDs,
                                    const CFuint geoOrder,
    		    const CFuint dim);

private:

  /// Default constructor without arguments.
  MapGeoEnt();

  /// Default destructor.
  ~MapGeoEnt();

  /// Default constructor without arguments.
  MapGeoEnt(const MapGeoEnt&);

private:

}; // end of class MapGeoEnt

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_MapGeoEnt_hh
