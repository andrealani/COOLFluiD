// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Framework_VolumeIntegratorImpl_hh
#define COOLFluiD_Framework_VolumeIntegratorImpl_hh

//////////////////////////////////////////////////////////////////////////////



#include "MathTools/RealVector.hh"
#include "MathTools/RealMatrix.hh"
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
/// of a volume integrator
/// @author Andrea Lani
/// @author Tiago Quintino
class Framework_API VolumeIntegratorImpl : public Common::OwnedObject,
                              public Common::NullableObject {
public:

  typedef Environment::ConcreteProvider<VolumeIntegratorImpl> PROVIDER;

  /// Constructor
  VolumeIntegratorImpl();

  /// Default destructor
  virtual ~VolumeIntegratorImpl();

  /// Template class specialization for ContourIntegratorImpl
  static void assertIntegratorType(const CFIntegration::Type& type)
  {
    cf_assert(type == CFIntegration::VOLUME);
  }

  /// Set up the private data to prepare the simulation
  virtual void setup() = 0;

  /// Compute the interpolated solution at the quadrature points
  virtual void computeSolutionAtQuadraturePoints
  (const std::vector<State*>& states,
  std::vector<State*>& values) = 0;

    /// Compute the interpolated coordinates at the quadrature points
  virtual void computeCoordinatesAtQuadraturePoints
  (const std::vector<Node*>& nodes,
    std::vector<Node*>& coord) = 0;

  /// Compute the interpolated coordinates and solution at the quadrature points
  /// Note: This is only to be used for isoparametric elements
  virtual void computeAllAtQuadraturePoints
  (const std::vector<Node*>& nodes,
  std::vector<Node*>& coord,
  const std::vector<State*>& states,
  std::vector<State*>& values) = 0;

  /// Compute the interpolated gradients at the quadrature points
  virtual const std::vector<RealVector>& getQuadraturePointsCoordinates() const = 0;

  /// Gets the ShapeFunctions values
  virtual const std::vector<RealVector>& getShapeFunctionsAtQuadraturePoints() const = 0;

  /// Computes the ShapeFunctions values
  virtual const std::vector<RealVector>& computeShapeFunctionsAtQuadraturePoints() = 0;

  /// Compute the interpolated gradients at the quadrature points
  /// @param nodes
  /// @param mappedCoord
  /// @param jacob
  virtual void computeJacobianAtQuadraturePoints(
          const std::vector<Node*>& nodes,
          const std::vector<RealVector>& mappedCoord,
                std::vector<RealMatrix>& jacob) = 0;

  /// Compute the interpolated gradients at the quadrature points
  virtual void computeGradSolutionShapeFAtQuadraturePoints(
          const std::vector<RealMatrix>& jacob,
          const std::vector<RealVector>& mappedCoord,
                std::vector<RealMatrix>& grad) = 0;

  /// Compute the determinant of the jacobian at the quadrature points
  /// @todo this function will become pure virtual
  virtual void computeJacobianDetAtQuadraturePoints(
          const std::vector<Node*>& nodes,
                std::valarray<CFreal>& detJacobian) = 0;

  /// Get the quadrature points pattern.
  const IntegratorPattern& getIntegratorPattern() const
  {
    cf_assert(_pattern.totalNbPts() == _coeff.size());
    return _pattern;
  }

  /// Get the constant coefficients
  const std::valarray<CFreal>& getCoeff() const
  {
    return _coeff;
  }

  /// Gets the Class name
  static std::string getClassName()
  {
    return "VolumeIntegratorImpl";
  }

  /// Gets the IntegratorProperties of this concrete integrator
  virtual const IntegratorProperties getIntegratorProperties() const = 0;

protected: // methods

  /// Get the quadrature points pattern.
  const IntegratorPattern& setupIntegratorPattern()
  {
    _pattern.resize(1);
    _pattern[0] = _coeff.size();
    return _pattern;
  }

protected: // data

  /// constant quadrature coefficients
  std::valarray<CFreal> _coeff;

  /// integration pattern
  IntegratorPattern _pattern;

}; // end class VolumeIntegratorImpl

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

Framework_Factory(VolumeIntegratorImpl)

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Framework_VolumeIntegratorImpl_hh
