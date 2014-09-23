// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Framework_ContourIntegratorImpl_hh
#define COOLFluiD_Framework_ContourIntegratorImpl_hh

//////////////////////////////////////////////////////////////////////////////



#include "MathTools/RealVector.hh"
#include "Common/NullableObject.hh"
#include "IntegratorPattern.hh"
#include "IntegratorProperties.hh"
#include "Environment/ConcreteProvider.hh"
#include "Common/OwnedObject.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

    class Node;
    class State;

//////////////////////////////////////////////////////////////////////////////

/// This class provides the shape dependent implementation
/// of an integrator
/// @author Andrea Lani
/// @author Tiago Quintino
class Framework_API ContourIntegratorImpl : public Common::NullableObject,
                              public Common::OwnedObject {
public:

  typedef Environment::ConcreteProvider<ContourIntegratorImpl> PROVIDER;

  /// Constructor
  ContourIntegratorImpl();

  /// Default destructor
  virtual ~ContourIntegratorImpl();

  /// Template class specialization for VolumeIntegratorImpl
  static void assertIntegratorType(const CFIntegration::Type& type)
  {
    cf_assert(type == CFIntegration::CONTOUR);
  }

  /// Set up the private data to prepare the simulation
  virtual void setup() = 0;

  /// Compute the coordinates of the quadrature points
  virtual const std::vector<RealVector>& getQuadraturePointsCoordinates() const = 0;

  /// Gets the ShapeFunctions values
  virtual const std::vector<RealVector>& getShapeFunctionsAtQuadraturePoints() const = 0;

  /// Computes the ShapeFunctions values
  virtual const std::vector<RealVector>& computeShapeFunctionsAtQuadraturePoints() = 0;

  /// Compute the interpolated solution at the quadrature points
  virtual void computeSolutionAtQuadraturePoints(
          const std::vector<State*>& states,
                std::vector<State*>& values) = 0;

  /// Compute the interpolated coordinates at the quadrature points
  virtual void computeCoordinatesAtQuadraturePoints(
          const std::vector<Node*>& nodes,
                std::vector<Node*>& coord) = 0;

  /// Compute the interpolated coordinates and solution at the quadrature points
  virtual void computeAllAtQuadraturePoints(
          const std::vector<Node*>& nodes,
                std::vector<Node*>& coord,
          const std::vector<State*>& states,
                std::vector<State*>& values) = 0;

  /// Compute the determinant of the face jacobian at the quadrature points
  virtual void computeFaceJacobianDetAtQuadraturePoints(
          const std::vector<Node*>& nodes,
                std::vector<RealVector>& faceJacobian) = 0;

  /// Get the quadrature points pattern.
  /// For example a integrator that loops on 4 faces, summing 2
  /// points in the first and 3 on the others has a pattern: [2,3,3,3]
  const IntegratorPattern& getIntegratorPattern() const
  {
    return _pattern;
  }

  /// Get the constant coefficients
  const std::vector<RealVector>& getCoeff() const
  {
    return _coeff;
  }

  /// Gets the Class name
  static std::string getClassName()
  {
    return "ContourIntegratorImpl";
  }

  /// Gets the IntegratorProperties of this concrete integrator
  virtual const IntegratorProperties getIntegratorProperties() const = 0;

protected: // methods

  /// Sets up the IntegratorPattern.
  const IntegratorPattern& setupIntegratorPattern()
  {
    _pattern.resize(_coeff.size());
    for(CFuint i = 0; i < _pattern.nbSteps(); ++i) {
      _pattern[i] = _coeff[i].size();
    }
    return _pattern;
  }

protected: // data

  /// constant quadrature coefficients
  std::vector<RealVector> _coeff;

  /// integrator pattern
  IntegratorPattern _pattern;

}; // end class ContourIntegratorImpl

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

Framework_Factory(ContourIntegratorImpl) // define the factoty instance

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Framework_ContourIntegratorImpl_hh
