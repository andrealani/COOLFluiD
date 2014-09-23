// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Framework_VolumeIntegrator_hh
#define COOLFluiD_Framework_VolumeIntegrator_hh

//////////////////////////////////////////////////////////////////////////////

#include "Integrator.hh"
#include "VolumeIntegratorImpl.hh"
#include "VectorialFunction.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

/// This class provides the functionality of a
/// GeometricEntity volume integrator.
/// @author Andrea Lani
/// @author Tiago Quintino
class Framework_API VolumeIntegrator : public Integrator<VolumeIntegratorImpl> {
public:

  /// Constructor
  VolumeIntegrator();

  /// Destructor
  virtual ~VolumeIntegrator();

  /// Set up the integrator for the simulation
  void setup();

  /// Compute the interpolated solution at the quadrature
  /// points on the supplied GeometricEntity
  void computeSolutionAtQuadraturePointsOnGeoEnt
  (GeometricEntity* const geo,
  std::vector<State*>& values)
  {
    getSolutionIntegrator(geo)
      ->computeSolutionAtQuadraturePoints(*geo->getStates(),values);
  }

  /// Compute the interpolated solution at the quadrature
  /// points on the supplied GeometricEntity
  void computeCoordinatesAtQuadraturePointsOnGeoEnt
  (GeometricEntity* const geo,
  std::vector<Node*>& coord)
  {
    getGeometryIntegrator(geo)
      ->computeCoordinatesAtQuadraturePoints(*geo->getNodes(),coord);
  }

  /// Compute a volume integration for on a GeometricEntity
  /// where the functor has no input paramenters
  template <class FUNCTOR, class RETURN>
  void integrateConstantFunctorOnGeoEnt(GeometricEntity* const geo,
                                     FUNCTOR& functor,
                                     RETURN& result);

  /// Compute a volume integration for on a GeometricEntity
  /// where the functor has the coordinates as paramenter
  template <class FUNCTOR, class RETURN>
  void integrateGeometricFunctorOnGeoEnt(GeometricEntity* const geo,
                                     FUNCTOR& functor,
                                     RETURN& result);

  /// Compute a volume integration for on a GeometricEntity
  /// where the functor has the solution as parameter
  template <class FUNCTOR, class RETURN>
  void integrateSolutionFunctorOnGeoEnt(GeometricEntity* const geo,
                                     FUNCTOR& functor,
                                     RETURN& result);

  /// Compute a volume integration for on a GeometricEntity
  /// where the functor has both the coordinates and
  /// the solution as paramenters
  template <class FUNCTOR, class RETURN>
  void integrateFullFunctorOnGeoEnt(GeometricEntity* const geo,
                                 FUNCTOR& functor,
                                 RETURN& result);

  /// Compute a volume integration for on a GeometricEntity
  /// where the functor takes the solution and its gradient as
  /// parameters
  template <class FUNCTOR, class RETURN>
  void integrateGradFunctorOnGeoEnt(GeometricEntity* const geo,
                                 FUNCTOR& functor,
                                 RETURN& result);

  /// Compute a volume integration for on a GeometricEntity
  /// where the functor takes the solution and its gradient
  /// and the geometric entity as parameters
  template <class FUNCTOR, class RETURN>
  void integrateGeneralFunctorOnGeoEnt(GeometricEntity* const geo,
                                 FUNCTOR& functor,
                                 RETURN& result);

  /// Compute a volume integration for on a GeometricEntity
  /// where the functor takes the states and shapefuntions as
  /// parameters
  template <class FUNCTOR, class RETURN>
  void integrateStateFunctorOnGeoEnt(GeometricEntity* const geo,
                                  FUNCTOR& functor,
                                  RETURN& result);

  /// Compute a volume integration for on a GeometricEntity
  /// where the functor takes the solution, the value of the shape function
  /// and a VectorialFunction as parameters
  template <class FUNCTOR, class RETURN>
  void integrateVectorialFunctorOnGeoEnt(GeometricEntity* const geo,
                                      FUNCTOR& functor,
                                      VectorialFunction& vf,
                                      RETURN& result);

  /// Compute a volume integration for on a GeometricEntity
  /// where the functor takes the value of the shape function
  /// as parameter
  template <class FUNCTOR, class RETURN>
  void integrateShapeFFunctorOnGeoEnt(GeometricEntity* const geo,
      FUNCTOR& functor,
      RETURN& result);

  /// Compute a volume integration for on a GeometricEntity
  /// where the functor takes the value of the shape function
  /// and the GeometricEntity as parameter
  template <class FUNCTOR, class RETURN>
  void integrateShapeFGeoFunctorOnGeoEnt(GeometricEntity* const geo,
  		     	      FUNCTOR& functor,
  			      RETURN& result);

  /// Sets the number of the solution quadrature points per each
  /// GeometricEntity
  void setNbSolQuadraturePoints(GeometricEntity *const geo,
  			CFuint& nbPoints);

private:

  /// Sets the coordinates in the States
  void buildCoordinatesAndStates(const IntegratorPattern& pattern);

  /// Copy constructor
  VolumeIntegrator(const VolumeIntegrator&);

  /// assignment operator
  VolumeIntegrator& operator= (const VolumeIntegrator&);

protected: // data

  /// VolumeIntegrator for non iso parametric geometric entities
  /// for the geometric shape function
  Common::SafePtr<VolumeIntegratorImpl> _intgGeo;

  /// VolumeIntegrator for non iso parametric geometric entities
  /// for the solution shape function
  Common::SafePtr<VolumeIntegratorImpl> _intgSol;

  /// storage for temporary interpolated solution values
  std::vector<State*> _solValues;

  /// storage for temporary interpolated coordinate values
  std::vector<Node*> _coord;

  /// storage for temporary interpolated gradient values
  /// size of vector is maximun size of quadrature points
  /// each matris is sized nb shape functions * nb of dimensions
  std::vector<RealMatrix> _gradValues;

  /// storage for temporary jacobians at quadrature points
  std::vector<RealMatrix> _jacob;

  /// storage for temporary interpolated jacobian determinant
  std::valarray<CFreal> _detJacobian;

}; // end class VolumeIntegrator

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#include "VolumeIntegrator.ci"

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Framework_VolumeIntegrator_hh
