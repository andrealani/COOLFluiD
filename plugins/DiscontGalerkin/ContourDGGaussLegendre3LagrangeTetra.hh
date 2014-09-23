#ifndef COOLFluiD_Numerics_ShapeFunctions_ContourDGGaussLegendre3LagrangeTetra_hh
#define COOLFluiD_Numerics_ShapeFunctions_ContourDGGaussLegendre3LagrangeTetra_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/GaussLegendreContourIntegratorImpl.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



   namespace ShapeFunctions {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class performs a Gauss3 integration on an P1 triangular element
 * with lagrange shape function
 *
 * @author Andrea Lani
 * @author Tiago Quintino
 *
 */
template <typename INTERPOLATOR>
class ContourDGGaussLegendre3LagrangeTetra :
  public Framework::GaussLegendreContourIntegratorImpl<ContourDGGaussLegendre3LagrangeTetra<INTERPOLATOR>,INTERPOLATOR> {
public:

  /**
   * Gets the name of this concrete integrator
   */
  static const std::string getName()
  {
    return "ContourDGGaussLegendre3" + INTERPOLATOR::getName();
  }

  /**
   * Get the number of quadrature points
   */
  static const CFuint getNbQuadraturePoints()
  {
    return 16;
  }

  /**
   * Get the number of faces
   */
  static const CFuint getNbFaces()
  {
    return 4;
  }

  /**
   * Gets integration type
   */
  static const CFIntegration::Type getIntegrationType()
  {
    return CFIntegration::CONTOUR;
  }

  /**
   * Gets quadrature type of this integrator
   */
  static const CFQuadrature::Type getQuadratureType()
  {
    return CFQuadrature::DGGAUSSLEGENDRE;
  }

  /**
   * Gets the shape on which this integrator operates
   */
  static const CFPolyOrder::Type getIntegrationOrder()
  {
    return CFPolyOrder::ORDER3;
  }

  /**
   * Gets the shape on which this integrator operates
   */
  static const CFGeoShape::Type getShape()
  {
    return INTERPOLATOR::getShape();
  }

  /**
   * Gets the interpolator order with whom this integrator operates
   */
  static const CFPolyOrder::Type getInterpolatorOrder()
  {
    return INTERPOLATOR::getInterpolatorOrder();
  }

  /**
   * Gets the interpolator type with whom this integrator operates
   */
  static const CFPolyForm::Type getInterpolatorType()
  {
    return INTERPOLATOR::getInterpolatorType();
  }

protected: // methods

  using Framework::GaussLegendreContourIntegratorImpl<ContourDGGaussLegendre3LagrangeTetra<INTERPOLATOR>,INTERPOLATOR>::_coeff;
  using Framework::GaussLegendreContourIntegratorImpl<ContourDGGaussLegendre3LagrangeTetra<INTERPOLATOR>,INTERPOLATOR>::_mappedCoord;

  /**
   * Set the mapped coordinates
   */
  void setMappedCoordinates();

  /**
   * Set the weights
   */
  void setWeights();

}; // end ContourDGGaussLegendre3LagrangeTetra

//////////////////////////////////////////////////////////////////////////////

    } // namespace ShapeFunctions



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_ShapeFunctions_ContourDGGaussLegendre3LagrangeTetra_hh
