// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_ShapeFunctions_ContourGaussLegendre3LagrangePyram_hh
#define COOLFluiD_ShapeFunctions_ContourGaussLegendre3LagrangePyram_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/GaussLegendreContourIntegratorImpl.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



  namespace ShapeFunctions {

//////////////////////////////////////////////////////////////////////////////

/// This class performs a Gauss3 integration on an P1 Pyramid element
/// with lagrange shape function
/// @author Tiago Quintino
template <typename INTERPOLATOR>
class ContourGaussLegendre3LagrangePyram : public Framework::GaussLegendreContourIntegratorImpl<ContourGaussLegendre3LagrangePyram<INTERPOLATOR>,INTERPOLATOR> {
public:

  /// Gets the name of this concrete integrator
  static const std::string getName()
  {
    return "ContourGaussLegendre3" + INTERPOLATOR::getName();
  }

  /// Get the number of quadrature points
  static CFuint getNbQuadraturePoints()
  {
    return 16;
  }

  /// Get the number of faces
  static CFuint getNbFaces()
  {
    return 5;
  }

  /// Gets integration type
  static CFIntegration::Type getIntegrationType()
  {
    return CFIntegration::CONTOUR;
  }

  /// Gets quadrature type of this integrator
  static CFQuadrature::Type getQuadratureType()
  {
    return CFQuadrature::GAUSSLEGENDRE;
  }

  /// Gets the shape on which this integrator operates
  static CFPolyOrder::Type getIntegrationOrder()
  {
    return CFPolyOrder::ORDER3;
  }

  /// Gets the shape on which this integrator operates
  static CFGeoShape::Type getShape()
  {
    return INTERPOLATOR::getShape();
  }

  /// Gets the interpolator order with whom this integrator operates
  static CFPolyOrder::Type getInterpolatorOrder()
  {
    return INTERPOLATOR::getInterpolatorOrder();
  }

  /// Gets the interpolator type with whom this integrator operates
  static CFPolyForm::Type getInterpolatorType()
  {
    return INTERPOLATOR::getInterpolatorType();
  }

protected: // methods

using Framework::GaussLegendreContourIntegratorImpl<ContourGaussLegendre3LagrangePyram<INTERPOLATOR>,INTERPOLATOR>::_coeff;
using Framework::GaussLegendreContourIntegratorImpl<ContourGaussLegendre3LagrangePyram<INTERPOLATOR>,INTERPOLATOR>::_mappedCoord;

  /// Set the mapped coordinates
  void setMappedCoordinates();

  /// Set the weights
  void setWeights();

}; // end ContourGaussLegendre3LagrangePyram

//////////////////////////////////////////////////////////////////////////////

    } // namespace ShapeFunctions



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_ShapeFunctions_ContourGaussLegendre3LagrangePyram_hh
