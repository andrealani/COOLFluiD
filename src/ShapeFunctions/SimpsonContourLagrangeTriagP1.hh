// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_ShapeFunctions_SimpsonContourLagrangeTriagP1_hh
#define COOLFluiD_ShapeFunctions_SimpsonContourLagrangeTriagP1_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/CFPolyForm.hh"
#include "Framework/CFGeoShape.hh"
#include "Framework/CFPolyOrder.hh"
#include "Framework/CFQuadrature.hh"
#include "Framework/CFIntegration.hh"
#include "Framework/ContourIntegratorImpl.hh"
#include "Common/NotImplementedException.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



  namespace ShapeFunctions {

//////////////////////////////////////////////////////////////////////////////

/// This class
/// @author Andrea Lani
/// @author Tiago Quintino
class SimpsonContourLagrangeTriagP1 : public Framework::ContourIntegratorImpl {
public:

  /// Gets the Framework::IntegratorProperties of this concrete integrator
  const Framework::IntegratorProperties getIntegratorProperties() const
  {
    Framework::IntegratorProperties prop(getName(),
                              getIntegrationType(),
                              getQuadratureType(),
                              getIntegrationOrder(),
                              getShape(),
                              getInterpolatorType(),
                              getInterpolatorOrder());
    return prop;
  }

  /// Gets the name of this concrete integrator
  static const std::string getName()
  {
    return "SimpsonContourLagrangeTriagP1";
  }

  /// Get the number of quadrature points
  static CFuint getNbQuadraturePoints()
  {
    return 9;
  }

  /// Gets integration type
  static CFIntegration::Type getIntegrationType()
  {
    return CFIntegration::CONTOUR;
  }

  /// Gets quadrature type of this integrator
  static CFQuadrature::Type getQuadratureType()
  {
    return CFQuadrature::SIMPSON;
  }

  /// Gets the shape on which this integrator operates
  static CFPolyOrder::Type getIntegrationOrder()
  {
    // @TODO AL: here it should be CFPolyOrder::ORDER2 but then the current integrator is not selected
    return CFPolyOrder::ORDER1;
  }

  /// Gets the shape on which this integrator operates
  static CFGeoShape::Type getShape()
  {
    return CFGeoShape::TRIAG;
  }

  /// Gets the interpolator order with whom this integrator operates
  static CFPolyOrder::Type getInterpolatorOrder()
  {
    return CFPolyOrder::ORDER1;
  }

  /// Gets the interpolator type with whom this integrator operates
  static CFPolyForm::Type getInterpolatorType()
  {
    return CFPolyForm::LAGRANGE;
  }

  /// Set up the private data to prepare the simulation
  void setup();

  /// Compute the interpolated gradients at the quadrature points
  const std::vector<RealVector>& getQuadraturePointsCoordinates() const
  {
    throw Common::NotImplementedException (FromHere(),"SimpsonContourLagrangeTriagP1::getQuadraturePointsCoordinates()");
  }

  /// Gets the ShapeFunctions values
  const std::vector<RealVector>& getShapeFunctionsAtQuadraturePoints() const
  {
    throw Common::NotImplementedException (FromHere(),"SimpsonContourLagrangeTriagP1::getShapeFunctionsAtQuadraturePoints()");
  }

  /// Computes the ShapeFunctions values
  const std::vector<RealVector>& computeShapeFunctionsAtQuadraturePoints()
  {
    throw Common::NotImplementedException (FromHere(),"SimpsonContourLagrangeTriagP1::computeShapeFunctionsAtQuadraturePoints()");
  }

  /// Compute the interpolated values in the quadrature points
  void computeSolutionAtQuadraturePoints(
      const std::vector<Framework::State*>& states,
             std::vector<Framework::State*>& values);

  /// Compute the interpolated coordinates in the quadrature points
  void computeCoordinatesAtQuadraturePoints(
      const std::vector<Framework::Node*>& nodes,
            std::vector<Framework::Node*>& coord)
  {
    throw Common::NotImplementedException (FromHere(),"SimpsonContourLagrangeTriagP1::computeCoordinatesAtQuadraturePoints()");
  }

  /// Compute the interpolated coordinates in the quadrature points
  void computeAllAtQuadraturePoints(
          const std::vector<Framework::Node*>& nodes,
                std::vector<Framework::Node*>& coord,
          const std::vector<Framework::State*>& states,
                std::vector<Framework::State*>& values)
  {
    computeSolutionAtQuadraturePoints(states,values);
  }

  /// Compute the determinant of the face jacobian at the quadrature points
  void computeFaceJacobianDetAtQuadraturePoints(
          const std::vector<Framework::Node*>& nodes,
  std::vector<RealVector>& faceJacobian);

}; // end class SimpsonContourLagrangeTriagP1

//////////////////////////////////////////////////////////////////////////////

    } // namespace ShapeFunctions



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_ShapeFunctions_SimpsonContourLagrangeTriagP1_hh
