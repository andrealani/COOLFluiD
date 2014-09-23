// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_ShapeFunctions_VolumeGaussLegendre2LagrangeTetra_hh
#define COOLFluiD_ShapeFunctions_VolumeGaussLegendre2LagrangeTetra_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/GaussLegendreVolumeIntegratorImpl.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



  namespace ShapeFunctions {

//////////////////////////////////////////////////////////////////////////////

/// This class performs a VolumeGaussLegendre2 integration on an
/// tetrahedron element
/// with lagrange shape function
/// @author Thomas Wuilbaut
template <typename INTERPOLATOR>
class VolumeGaussLegendre2LagrangeTetra :
public Framework::GaussLegendreVolumeIntegratorImpl<VolumeGaussLegendre2LagrangeTetra<INTERPOLATOR>,INTERPOLATOR> {
public:

  /// Gets the name of this concrete integrator
  static const std::string getName()
  {
    return "VolumeGaussLegendre2" + INTERPOLATOR::getName();
  }

  /// Gets the number of quadrature points
  static CFuint getNbQuadraturePoints()
  {
    return 4;
  }

  /// Gets integration type
  static CFIntegration::Type getIntegrationType()
  {
    return CFIntegration::VOLUME;
  }

  /// Gets quadrature type of this integrator
  static CFQuadrature::Type getQuadratureType()
  {
    return CFQuadrature::GAUSSLEGENDRE;
  }

  /// Gets the shape on which this integrator operates
  static CFPolyOrder::Type getIntegrationOrder()
  {
    return CFPolyOrder::ORDER2;
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
  using Framework::GaussLegendreVolumeIntegratorImpl<VolumeGaussLegendre2LagrangeTetra<INTERPOLATOR>,INTERPOLATOR>::_mappedCoord;
  using Framework::GaussLegendreVolumeIntegratorImpl<VolumeGaussLegendre2LagrangeTetra<INTERPOLATOR>,INTERPOLATOR>::_coeff;

  /// Set the mapped coordinates
  void setMappedCoordinates();

  /// Set the weights
  void setWeights();

}; // end VolumeGaussLegendre2LagrangeTetra

//////////////////////////////////////////////////////////////////////////////

    } // namespace ShapeFunctions



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_ShapeFunctions_VolumeGaussLegendre2LagrangeTetra_hh
