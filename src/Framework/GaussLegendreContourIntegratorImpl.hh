// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Framework_GaussLegendreContourIntegratorImpl_hh
#define COOLFluiD_Framework_GaussLegendreContourIntegratorImpl_hh

//////////////////////////////////////////////////////////////////////////////

#include "Common/CFLog.hh"

#include "Framework/CFPolyForm.hh"
#include "Framework/CFPolyOrder.hh"
#include "Framework/CFGeoShape.hh"
#include "Framework/CFQuadrature.hh"
#include "Framework/CFIntegration.hh"
#include "Framework/ContourIntegratorImpl.hh"
#include "Framework/State.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

/// This class is an interface to a Gauss Legendre point integrator
/// with shape functions supplied as template parameters.
/// @author Andrea Lani
/// @author Tiago Quintino
template <typename INTEGRATOR, typename SHAPE>
class GaussLegendreContourIntegratorImpl : public ContourIntegratorImpl {
public:

  /// Constructor
  GaussLegendreContourIntegratorImpl() :
    ContourIntegratorImpl(),
    _mappedCoord()
  {
  }

  /// Default destructor
  virtual ~GaussLegendreContourIntegratorImpl()
  {
  }

  /// Set up the private data to prepare the simulation
  void setup()
  {
    CFLogDebugMed(
          "Setting up GaussContourIntegratorImpl for "
          << SHAPE::getName()
          << "\n");

    _pattern.setNbShapeFunctions(SHAPE::getNbNodes());

    setMappedCoordinates();
    setWeights();
    setupIntegratorPattern();
    setupNodalSF();
  }

  /// Set up the nodal shape function data
  void setupNodalSF()
  {
    /// @todo use of sum() is inefficient !!!
    const CFuint nbQPts = _pattern.sum();
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
  void computeAllAtQuadraturePoints(
          const std::vector<Node*>& nodes,
                std::vector<Node*>& coord,
          const std::vector<State*>& states,
                std::vector<State*>& values)
  {
    SHAPE::computeShapeFunctions(_mappedCoord,_shapeFunc);
    SHAPE::interpolate(nodes, _shapeFunc, coord);
    SHAPE::interpolate(states, _shapeFunc, values);
  }

  /// Compute the determinant of the face jacobian at the quadrature points
  void computeFaceJacobianDetAtQuadraturePoints(
          const std::vector<Node*>& nodes,
                std::vector<RealVector>& faceJacobian)
  {
    SHAPE::computeFaceJacobianDeterminant(_mappedCoord,
                                          nodes,
                                          _pattern,
                                          faceJacobian);
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

  /// Computes the ShapeFunctions values
  const std::vector<RealVector>& computeShapeFunctionsAtQuadraturePoints()
  {
    SHAPE::computeShapeFunctions(_mappedCoord,_shapeFunc);
    return _shapeFunc;
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

}; // end class GaussLegendreContourIntegratorImpl

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Framework_GaussLegendreContourIntegratorImpl_hh
