// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Framework_GaussLegendreVolumeIntegratorImpl_hh
#define COOLFluiD_Framework_GaussLegendreVolumeIntegratorImpl_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/CFPolyForm.hh"
#include "Framework/CFPolyOrder.hh"
#include "Framework/CFGeoShape.hh"
#include "Framework/CFQuadrature.hh"
#include "Framework/CFIntegration.hh"
#include "Framework/VolumeIntegratorImpl.hh"
#include "Framework/State.hh"
#include "Framework/PhysicalModel.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

/// This class is an interface to a Gauss Legendre point integrator
/// with shape functions supplied as template parameters.
/// @author Andrea Lani
/// @author Tiago Quintino
template <typename INTEGRATOR, typename SHAPE>
class GaussLegendreVolumeIntegratorImpl : public VolumeIntegratorImpl {
public:

  /// Constructor
  GaussLegendreVolumeIntegratorImpl() :  VolumeIntegratorImpl(), _mappedCoord() {}

  /// Default destructor
  virtual ~GaussLegendreVolumeIntegratorImpl() {}

  /// Set up the private data to prepare the simulation
  void setup()
  {
    CFLogDebugMin("Setting up " << INTEGRATOR::getName() << "\n");

    _pattern.setNbShapeFunctions(SHAPE::getNbNodes());

    setMappedCoordinates();
    setWeights();
    setupIntegratorPattern();
    setupNodalSF();
  }

  /// Set up the nodal shape function data
  void setupNodalSF()
  {
    IntegratorPattern pattern = getIntegratorPattern();
    cf_assert(pattern.nbSteps() == 1);
    const CFuint nbQPts = pattern.nbPts(0);
    _shapeFunc.resize(nbQPts);
    for (CFuint ip = 0; ip < nbQPts; ++ip) {
      _shapeFunc[ip].resize(SHAPE::getNbNodes());
      _shapeFunc[ip] = 0.0;
    }
  }

  /// Compute the interpolated values in the quadrature points
  void computeSolutionAtQuadraturePoints(const std::vector<State*>& states,
      std::vector<State*>& values)
  {
    SHAPE::computeShapeFunctions(_mappedCoord,_shapeFunc);
    SHAPE::interpolate(states, _shapeFunc, values);
  }

  /// Compute the interpolated coordinates in the quadrature points
  void computeCoordinatesAtQuadraturePoints(const std::vector<Node*>& nodes,
    			    std::vector<Node*>& coord)
  {
    SHAPE::computeShapeFunctions(_mappedCoord,_shapeFunc);
    SHAPE::interpolate(nodes, _shapeFunc, coord);
  }

  /// Compute the interpolated coordinates in the quadrature points
  void computeAllAtQuadraturePoints
  (const std::vector<Node*>& nodes,
  std::vector<Node*>& coord,
  const std::vector<State*>& states,
  std::vector<State*>& values)
  {
    SHAPE::computeShapeFunctions(_mappedCoord,_shapeFunc);
    SHAPE::interpolate(nodes, _shapeFunc, coord);
    SHAPE::interpolate(states, _shapeFunc, values);
  }

  /// Gets the quadrature points coordinates
  const std::vector<RealVector>& getQuadraturePointsCoordinates() const
  {
    return _mappedCoord;
  }

  /// Gets the ShapeFunctions values
  const std::vector<RealVector>& getShapeFunctionsAtQuadraturePoints() const
  {
    return _shapeFunc;
  }

  /// Compute the ShapeFunctions values
  const std::vector<RealVector>& computeShapeFunctionsAtQuadraturePoints()
  {
    SHAPE::computeShapeFunctions(_mappedCoord,_shapeFunc);
    return _shapeFunc;
  }

  /// Compute the interpolated gradients at the quadrature points
  /// @param nodes
  /// @param mappedCoord   * @param jacob
  void computeJacobianAtQuadraturePoints(
          const std::vector<Node*>& nodes,
          const std::vector<RealVector>& mappedCoord,
                std::vector<RealMatrix>& jacob)
  {
  const CFuint dim = PhysicalModelStack::getActive()->getDim();
  const CFuint dimensionality = SHAPE::getDimensionality();
  if(dim == dimensionality) {
    SHAPE::computeJacobian(nodes, mappedCoord, jacob);
  }
  else {
    if(dim - DIM_1D == dimensionality) {
      SHAPE::computeJacobianPlus1D(nodes, mappedCoord, jacob);
    }
    if(dim - DIM_2D == dimensionality) {
      SHAPE::computeJacobianPlus2D(nodes, mappedCoord, jacob);
    }
  }
  }

  /// Compute the interpolated gradients at the quadrature points
  void computeGradSolutionShapeFAtQuadraturePoints(
          const std::vector<RealMatrix>& jacob,
          const std::vector<RealVector>& mappedCoord,
                std::vector<RealMatrix>& grad)
  {
    SHAPE::computeGradientStates(jacob, mappedCoord, grad);
  }

  /// Compute the jacobian of transformation at the quadrature points
  /// @todo This can be optimized with a single call to compute the jacobian
  ///       including the interpolation inside.
  void computeJacobianDetAtQuadraturePoints(const std::vector<Node*>& nodes,
                                                  std::valarray<CFreal>& detJacobian)
  {
    const CFuint dim = PhysicalModelStack::getActive()->getDim();
    const CFuint dimensionality = SHAPE::getDimensionality();
    if(dim == dimensionality) {
      SHAPE::computeJacobianDeterminant(_mappedCoord,nodes,detJacobian);
    }
    else {
    if(dim - DIM_1D == dimensionality) {
      SHAPE::computeJacobianDeterminantPlus1D(_mappedCoord,nodes,detJacobian);
      }
      if(dim - DIM_2D == dimensionality) {
        SHAPE::computeJacobianDeterminantPlus2D(_mappedCoord,nodes,detJacobian);
      }
    }
  }

  /// Gets the IntegratorProperties of this concrete integrator
  const IntegratorProperties getIntegratorProperties() const
  {
    return IntegratorProperties(INTEGRATOR::getName(),
                                INTEGRATOR::getIntegrationType(),
                                INTEGRATOR::getQuadratureType(),
                                INTEGRATOR::getIntegrationOrder(),
                                INTEGRATOR::getShape(),
                                INTEGRATOR::getInterpolatorType(),
                                INTEGRATOR::getIntegrationOrder());
  }

protected: // methods

  /// Set the mapped coordinates
  virtual void setMappedCoordinates() = 0;

  /// Set the weights
  virtual void setWeights() = 0;

protected: // data

  /// mapped coordinates of the quadrature points
  std::vector<RealVector> _mappedCoord;

  /// nodal shape function values at quadrature points
  std::vector<RealVector> _shapeFunc;

}; // end class GaussLegendreVolumeIntegratorImpl

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Framework_GaussLegendreVolumeIntegratorImpl_hh
