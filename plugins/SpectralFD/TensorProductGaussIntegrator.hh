#ifndef COOLFluiD_Numerics_SpectralFD_TensorProductGaussIntegrator_hh
#define COOLFluiD_Numerics_SpectralFD_TensorProductGaussIntegrator_hh

//////////////////////////////////////////////////////////////////////////////

#include "Common/COOLFluiD.hh"

#include "Common/CFLog.hh"

#include "Framework/CFPolyOrder.hh"

#include "MathTools/RealVector.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace SpectralFD {

//////////////////////////////////////////////////////////////////////////////

/// This class implements an integrator over an arbitrary tensor product cell/face (line, quadrangle, box)
class TensorProductGaussIntegrator {

public: // functions

  /// Constructor
  TensorProductGaussIntegrator();

  /// Constructor
  TensorProductGaussIntegrator(CFDim dimensionality, CFPolyOrder::Type integratorOrder);

  /// Destructor
  ~TensorProductGaussIntegrator();

  /// @return m_dimensionality
  CFDim getDimensionality()
  {
    return m_dimensionality;
  }

  /// @return m_integratorOrder
  CFPolyOrder::Type getIntegratorOrder()
  {
    return m_integratorOrder;
  }

  /// @return m_nbrQuadPnts
  CFuint getNbrQuadPnts()
  {
    return m_nbrQuadPnts;
  }

  /// @return m_quadPntsWheights
  std::vector< CFreal > getQuadPntsWheights()
  {
    return m_quadPntsWheights;
  }

  /// @return m_quadPntsMappedCoor
  std::vector< RealVector > getQuadPntsMappedCoords()
  {
    return m_quadPntsMappedCoor;
  }

  /**
   * Returns a vector with the (global) coordinates of the quadrature points,
   * which can than be used to evaluate a general (non-polynomial) function.
   * Primarily intended for solution initialization (where there can be non-polynomial functions).
   * @par nodes the nodes of the line/quadrangle/box over which to integrate.
   * @return a vector of quadrature point coordinates
   */
  std::vector< RealVector > getQuadPntsCoords(const std::vector< RealVector >& nodeCoord);

  /**
   * Returns the vector with the quadrature wheights multiplied with the jacobian determinant
   * @return vector of quadrature wheights multiplied with the jacobian determinant
   */
  std::vector< CFreal > getQuadPntsWheights(const std::vector< RealVector >& nodeCoord);

  /**
   * Returns a vector with the (global) coordinates of the quadrature points,
   * which can than be used to evaluate a general (non-polynomial) function.
   * @return a vector of quadrature point coordinates
   */
  std::vector< RealVector > getQuadPntsCoordsPlus1D(const std::vector< RealVector >& nodeCoord);

  /**
   * Returns a vector with the local face unit normals at the quadrature points.
   * @return a vector of quadrature point coordinates
   */
  std::vector< RealVector > getQuadPntsUnitNormalsPlus1D(const std::vector< RealVector >& nodeCoord);

  /**
   * Returns the vector with the quadrature wheights multiplied with the jacobian determinant
   * @return vector of quadrature wheights multiplied with the jacobian determinant
   */
  std::vector< CFreal > getQuadPntsWheightsPlus1D(const std::vector< RealVector >& nodeCoord);

  /**
   * Function that sets the variable m_dimensionality and
   * resets the variables m_quadPntsMappedCoor and m_quadPntsWheights.
   * @par dimensionality: the new dimensionality
   */
  void setDimensionality(const CFDim dimensionality);

  /**
   * Function that sets the variable m_integratorOrder and
   * resets the variables m_quadPntsMappedCoor and m_quadPntsWheights.
   * @par integratorOrder: the new order
   */
  void setIntegratorOrder(const CFPolyOrder::Type integratorOrder);

private: // functions

  /**
   * Function that resets m_quadPntsMappedCoor and m_quadPntsWheights.
   * @par integratorOrder: the new order
   * @pre resetQuadPntsMappedCoorAndWheights1D
   */
  void resetQuadPntsMappedCoorAndWheights(const CFDim dimensionality);

  /**
   * Function that resets m_quadPntsMappedCoor1D and m_quadPntsWheights1D.
   * @par integratorOrder: the new order
   */
  void resetQuadPntsMappedCoorAndWheights1D(const CFPolyOrder::Type integratorOrder);

  /**
   * Function that calculates the number of terms in the highest order polynomial
   * that is integrated exactly, based on m_dimensionality and m_integratorOrder,
   * and puts it in m_nbrPolyTerms4ExIntegr.
   */
  void calcNbrPolyTerms4ExIntegr();

  /**
   * @return the spatial basisfunctions evaluated at the given local coordinates
   */
  std::vector< CFreal > evaluateBasisFunction(const CFDim dimensionality,
                                              const CFPolyOrder::Type spatialOrder,
                                              const RealVector& localCoord);

  /**
   * @return the derivatives of the spatial basisfunctions evaluated at the given local coordinates
   */
  std::vector< RealVector > evaluateDerivBasisFunction(const CFDim dimensionality,
                                                       const CFPolyOrder::Type spatialOrder,
                                                       const RealVector& localCoord);

private: // data

  /// Dimensionality of the simplex
  CFDim m_dimensionality;

  /// Maximum polynomial order which is integrated exactly
  CFPolyOrder::Type m_integratorOrder;

  /// The number of terms in the highest degree polynomial that is exactly integrated
  CFuint m_nbrPolyTerms4ExIntegr;

  /// The number of quadrature points
  CFuint m_nbrQuadPnts;

  /// Vector for the mapped coordinates of the quadrature points
  std::vector< RealVector > m_quadPntsMappedCoor;

  /// Vector for the quadrature wheights
  std::vector< CFreal > m_quadPntsWheights;

  /// The number of quadrature points in 1D
  CFuint m_nbrQuadPnts1D;

  /// mapped coordinates of the quadrature points in 1D
  std::vector< CFreal > m_quadPntsMappedCoor1D;

  /// Vector for the quadrature wheights in 1D
  std::vector< CFreal > m_quadPntsWheights1D;

};

//////////////////////////////////////////////////////////////////////////////

  }  // namespace SpectralFD
}  // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_SpectralFD_TensorProductGaussIntegrator_hh
