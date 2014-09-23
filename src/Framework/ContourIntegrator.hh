// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Framework_ContourIntegrator_hh
#define COOLFluiD_Framework_ContourIntegrator_hh

//////////////////////////////////////////////////////////////////////////////

#include "Integrator.hh"
#include "ContourIntegratorImpl.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

/// This class provides the functionality of a
/// GeometricEntity contour integrator.
/// @author Andrea Lani
/// @author Tiago Quintino
class Framework_API ContourIntegrator : public Integrator<ContourIntegratorImpl> {
public:

  /// Constructor
  ContourIntegrator();

  /// Destructor
  ~ContourIntegrator();

  /// Set up the integrator for the simulation
  void setup();

  /// Compute all (interpolated values and coordinates) in all the
  /// points on the supplied GeometricEntity
  /// @pre this must be called before calling integrateSolutionFunctorOnGeo()
  void computeAllAtQuadraturePointsOnGeo(GeometricEntity* const geo)
  {
    // gets the correct implementor corresponding to the
    // shape of the element where you want to integrate
    ContourIntegratorImpl *const impl = getSolutionIntegrator(geo);
    cf_assert(impl != CFNULL);
    impl->computeAllAtQuadraturePoints(*geo->getNodes(), _coord,
                                       *geo->getStates(), _values);
  }

  /// Compute a contour integration for on a GeometricEntity
  template <class FUNCTOR>
  void integrateConstantFunctorOnGeo(GeometricEntity* const geo,
                                     FUNCTOR& functor,
                                     RealVector& result);

  /// Compute a contour integration for on a GeometricEntity
  template <class FUNCTOR>
  void integrateGeometricFunctorOnGeo(GeometricEntity* const geo,
                                      FUNCTOR& functor,
                                      RealVector& result);

  /// Compute a contour integration for on a GeometricEntity
  template <class FUNCTOR>
  void integrateSolutionFunctorOnGeo(GeometricEntity* const geo,
                                     FUNCTOR& functor,
				     const std::vector<RealVector>& pdata,
                                     RealVector& result);
  
  /// Compute a contour integration for on a GeometricEntity
  template <class FUNCTOR>
  void integrateFullFunctorOnGeo(GeometricEntity* const geo,
                                 FUNCTOR& functor,
                                 RealVector& result);
  
  /// Set a pointer to the unit face normals for the considered cell
  void setFaceNormals(std::vector<RealVector> *const faceNormals)
  {
    _unitFaceNormals = faceNormals;
  }
  
  /// Get the values
  void getValues(std::vector<State*>& v)
  {
    v.resize(_values.size());
    for (CFuint i = 0; i < v.size(); ++i) {
      v[i] = _values[i];
    }
  }
  
  /// Sets the number of the solution quadrature points per each
  /// GeometricEntity
  void setNbSolQuadraturePoints(const GeometricEntity* geo,
                                CFuint& nbPoints);

private:

  /// Sets the coordinates in the States
  void buildCoordinatesAndStates(const IntegratorPattern& pattern);

  /// Copy constructor
  ContourIntegrator(const ContourIntegrator&);

  /// assignment operator
  ContourIntegrator& operator= (const ContourIntegrator&);

protected: // data

  /// face normals
  std::vector<RealVector>* _unitFaceNormals;

  /// storage for temporary face jacobians determinants
  std::vector<RealVector> _faceJacobian;

  /// storage for temporary interpolated values
  std::vector<State*> _values;

  /// storage for temporary interpolated coordinates
  std::vector<Node*> _coord;

}; // end class ContourIntegrator

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#include "ContourIntegrator.ci"

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Framework_ContourIntegrator_hh
