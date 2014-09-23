// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_ShapeFunctions_ContourGaussLegendre1LagrangeTetra_hh
#define COOLFluiD_ShapeFunctions_ContourGaussLegendre1LagrangeTetra_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/GaussLegendreContourIntegratorImpl.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
namespace ShapeFunctions {

//////////////////////////////////////////////////////////////////////////////

/// This class performs a Gauss1 integration on an P1 tetrahedron element
/// with lagrange shape function
/// @author Andrea Lani
/// @author Tiago Quintino
template <typename INTERPOLATOR>
class ContourGaussLegendre1LagrangeTetra : public Framework::GaussLegendreContourIntegratorImpl<ContourGaussLegendre1LagrangeTetra<INTERPOLATOR>,INTERPOLATOR> {
public:

  /// Gets the name of this concrete integrator
  static const std::string getName()
  {
    return "ContourGaussLegendre1" + INTERPOLATOR::getName();
  }

  /// Get the number of quadrature points
  CFuint getNbQuadraturePoints() const
  {
    return 4;
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
    return CFPolyOrder::ORDER1;
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

  using Framework::GaussLegendreContourIntegratorImpl<ContourGaussLegendre1LagrangeTetra<INTERPOLATOR>,INTERPOLATOR>::_coeff;
  using Framework::GaussLegendreContourIntegratorImpl<ContourGaussLegendre1LagrangeTetra<INTERPOLATOR>,INTERPOLATOR>::_mappedCoord;

  /// Set the mapped coordinates
  void setMappedCoordinates();

  /// Set the weights
  void setWeights();

}; // end class ContourGaussLegendre1LagrangeTetra

//////////////////////////////////////////////////////////////////////////////

    } // namespace ShapeFunctions



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_ShapeFunctions_ContourGaussLegendre1LagrangeTetra_hh
