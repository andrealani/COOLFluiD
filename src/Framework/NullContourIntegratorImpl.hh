// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Framework_NullContourIntegratorImpl_hh
#define COOLFluiD_Framework_NullContourIntegratorImpl_hh

//////////////////////////////////////////////////////////////////////////////

#include "Common/CFLog.hh"
#include "ContourIntegratorImpl.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

/// This class provides the shape dependent implementation
/// of a Null volume integrator
/// @author Tiago Quintino
class Framework_API NullContourIntegratorImpl : public ContourIntegratorImpl {
public:

  /// Constructor
  NullContourIntegratorImpl()
  {
  }

  /// Default destructor
  ~NullContourIntegratorImpl()
  {
  }

  /// This is a null object so isNull is set to true
  bool isNull() const
  {
    return true;
  }

  /// Gets the name of this concrete integrator
  static const std::string getName()
  {
    return "Null";
  }

  /// Get the number of quadrature points
  static CFuint getNbQuadraturePoints()
  {
    return 0;
  }

  /// Get the number of faces
  static CFuint getNbFaces()
  {
    return 0;
  }

  /// Gets integration type
  static CFIntegration::Type getIntegrationType()
  {
    return CFIntegration::CONTOUR;
  }

  /// Gets quadrature type of this integrator
  static CFQuadrature::Type getQuadratureType()
  {
    return CFQuadrature::INVALID;
  }

  /// Gets the shape on which this integrator operates
  static CFPolyOrder::Type getIntegrationOrder()
  {
    return CFPolyOrder::MAXORDER;
  }

  /// Gets the shape on which this integrator operates
  static CFGeoShape::Type getShape()
  {
    return CFGeoShape::INVALID;
  }

  /// Gets the interpolator order with whom this integrator operates
  static CFPolyOrder::Type getInterpolatorOrder()
  {
    return CFPolyOrder::ORDER0;
  }

  /// Gets the interpolator type with whom this integrator operates
  static CFPolyForm::Type getInterpolatorType()
  {
    return CFPolyForm::INVALID;
  }

  /// Set up the private data to prepare the simulation
  void setup()
  {
    CFLog(VERBOSE,"Calling setup() on a NULL integrator\n");
  }

  /// Compute the interpolated gradients at the quadrature points
  const std::vector<RealVector>& getQuadraturePointsCoordinates() const
  {
    CFLog(VERBOSE,"Calling getQuadraturePointsCoordinates() on a NULL integrator\n");
    return _dummy;
  }

  /// Gets the ShapeFunctions values
  const std::vector<RealVector>& getShapeFunctionsAtQuadraturePoints() const
  {
    CFLog(VERBOSE,"Calling getShapeFunctionsAtQuadraturePoints() on a NULL integrator\n");
    return _dummy;
  }

  /// Computes the ShapeFunctions values
  const std::vector<RealVector>& computeShapeFunctionsAtQuadraturePoints()
  {
    CFLog(VERBOSE,"Calling computeShapeFunctionsAtQuadraturePoints() on a NULL integrator\n");
    return _dummy;
  }

  /// Compute the interpolated solution at the quadrature points
  void computeSolutionAtQuadraturePoints(
          const std::vector<State*>& states,
                std::vector<State*>& values)
  {
    CFLog(VERBOSE,"Calling computeSolutionAtQuadraturePoints() on a NULL integrator\n");
  }

  /// Compute the interpolated coordinates at the quadrature points
  void computeCoordinatesAtQuadraturePoints(
          const std::vector<Node*>& nodes,
                std::vector<Node*>& coord)
  {
    CFLog(VERBOSE,"Calling computeCoordinatesAtQuadraturePoints() on a NULL integrator\n");
  }

  /// Compute the interpolated coordinates and solution at the quadrature points
  void computeAllAtQuadraturePoints(
          const std::vector<Node*>& nodes,
                std::vector<Node*>& coord,
          const std::vector<State*>& states,
                std::vector<State*>& values)
  {
    CFLog(VERBOSE,"Calling computeAllAtQuadraturePoints() on a NULL integrator\n");
  }

  /// Compute the determinant of the face jacobian at the quadrature points
  void computeFaceJacobianDetAtQuadraturePoints(
          const std::vector<Node*>& nodes,
                std::vector<RealVector>& faceJacobian)
  {
    CFLog(VERBOSE,"Calling computeFaceJacobianDetAtQuadraturePoints() on a NULL integrator\n");
  }

  /// Gets the Class name
  static std::string getClassName()
  {
    return "NullContourIntegratorImpl";
  }

  /// Gets the IntegratorProperties of this concrete integrator
  const IntegratorProperties getIntegratorProperties() const
  {
    return IntegratorProperties(getName(),
                                getIntegrationType(),
                                getQuadratureType(),
                                getIntegrationOrder(),
                                getShape(),
                                getInterpolatorType(),
                                getIntegrationOrder());
  }

  private:

    /// Dummy data
    std::vector<RealVector> _dummy;

}; // end class NullContourIntegratorImpl

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Framework_NullContourIntegratorImpl_hh
