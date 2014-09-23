// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Framework_IntegrableEntity_hh
#define COOLFluiD_Framework_IntegrableEntity_hh

//////////////////////////////////////////////////////////////////////////////

#include "Common/StringOps.hh"
#include "Common/OwnedObject.hh"
#include "Common/NonCopyable.hh"
#include "MathTools/RealVector.hh"
#include "MathTools/RealMatrix.hh"
#include "Environment/ConcreteProvider.hh"
#include "Framework/Framework.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

    class GeometricEntity;

//////////////////////////////////////////////////////////////////////////////

/// This class represents the basic interface for a physics-dependent
/// integrable entity
/// @author Andrea Lani
/// @author Tiago Quintino
class Framework_API IntegrableEntity : public Common::OwnedObject,
                         public Common::NonCopyable<IntegrableEntity> {
public:

  /// Type for the provider of this abstract class
  typedef Environment::ConcreteProvider<IntegrableEntity> PROVIDER;

  /// Constructor
  IntegrableEntity();

  /// Default destructor
  virtual ~IntegrableEntity();

  /// Compute a contour integration of a physical quantity
  /// in the given cell
  virtual void integrateContourInGeo(GeometricEntity* const geo, RealVector& result);

  /// Compute a volume integration of a physical quantity
  /// in the given cell
  virtual void integrateVolumeInGeo(GeometricEntity* const geo, RealVector& result);

  /// Compute a contour integration of a physical quantity
  /// in the given cell
  virtual void integrateContourInGeo(GeometricEntity* const geo, RealMatrix& result);

  /// Compute a volume integration of a physical quantity
  /// in the given cell
  virtual void integrateVolumeInGeo(GeometricEntity* const geo, RealMatrix& result);

  /// Gets the Class name
  static std::string getClassName() { return "IntegrableEntity"; }

}; // end of class IntegrableEntity

//////////////////////////////////////////////////////////////////////////////

    } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

Framework_Factory(IntegrableEntity) // define the factoty instance

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Framework_IntegrableEntity_hh
