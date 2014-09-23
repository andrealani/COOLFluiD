// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Framework_NullIntegrableEntity_hh
#define COOLFluiD_Framework_NullIntegrableEntity_hh

//////////////////////////////////////////////////////////////////////////////

#include "IntegrableEntity.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

    class GeometricEntity;

//////////////////////////////////////////////////////////////////////////////

/// This class represents the flux corresponding to the Euler physical
/// model 2D and assuming conservative variables
/// @author Andrea Lani
class Framework_API NullIntegrableEntity : public IntegrableEntity {
public:

  /// Constructor
  NullIntegrableEntity();

  /// Default destructor
  ~NullIntegrableEntity();

  /// Compute a contour integration of a physical quantity in the given cell
  virtual void integrateContourInGeo(GeometricEntity* const geo, RealVector& result);

  /// Compute a volume integration of a physical quantity in the given cell
  virtual void integrateVolumeInGeo(GeometricEntity* const geo, RealVector& result);

  /// Compute a contour integration of a physical quantity in the given cell
  virtual void integrateContourInGeo(GeometricEntity* const geo, RealMatrix& result);

  /// Compute a volume integration of a physical quantity   /// in the given cell
  virtual void integrateVolumeInGeo(GeometricEntity* const geo, RealMatrix& result);

}; // end of class NullIntegrableEntity

//////////////////////////////////////////////////////////////////////////////

    } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Framework_NullIntegrableEntity_hh
