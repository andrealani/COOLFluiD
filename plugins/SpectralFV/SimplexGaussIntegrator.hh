#ifndef COOLFluiD_Numerics_SpectralFV_SimplexGaussIntegrator_hh
#define COOLFluiD_Numerics_SpectralFV_SimplexGaussIntegrator_hh

//////////////////////////////////////////////////////////////////////////////



#include "Common/COOLFluiD.hh"

#include "Common/CFLog.hh"

#include "Framework/CFPolyOrder.hh"

#include "MathTools/RealVector.hh"

#include "SpectralFV/SimplexGaussIntegrator.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace SpectralFV {

//////////////////////////////////////////////////////////////////////////////

/// This class implements an integrator over an arbitrary simplex
class SimplexGaussIntegrator {

public: // functions

  /// Constructor
  SimplexGaussIntegrator();

  /// Constructor
  SimplexGaussIntegrator(CFDim dimensionality, CFPolyOrder::Type integratorOrder);

  /// Destructor
  ~SimplexGaussIntegrator();

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

  /// @return m_quadPntsNormals
  std::vector< RealVector > getQuadPntsNormals()
  {
    return m_quadPntsNormals;
  }

  /**
   * Returns a vector with the (global) coordinates of the quadrature points,
   * which can than be used to evaluate a general (non-polynomial) function.
   * Primarily intended for solution initialization (where there can be non-polynomial functions).
   * Only to be used for P1 simplices.
   * @par simplexNodes the nodes of the simplex over which to integrate.
   * @return a vector of quadrature point coordinates
   */
  std::vector< RealVector > getQuadPntsCoords(const std::vector< RealVector >& simplexNodeCoord);

  /**
   * Returns the vector with the quadrature wheights multiplied with the jacobian determinant
   * @return vector of quadrature wheights multiplied with the jacobian determinant
   */
  std::vector< CFreal > getQuadPntsWheights(const std::vector< RealVector >& simplexNodeCoord);

  /**
   * Returns a vector with the (global) coordinates of the quadrature points,
   * which can than be used to evaluate a general (non-polynomial) function.
   * Only to be used for P1 simplices in a space with dimensionality of one degree higher.
   * @par simplexNodes the nodes of the simplex over which to integrate.
   * @return a vector of quadrature point coordinates
   */
  std::vector< RealVector > getQuadPntsCoordsPlus1D(const std::vector< RealVector >& simplexNodeCoord);

  /**
   * Returns the vector with the quadrature wheights multiplied with the jacobian determinant
   * @return vector of quadrature wheights multiplied with the jacobian determinant
   */
  std::vector< CFreal > getQuadPntsWheightsPlus1D(const std::vector< RealVector >& simplexNodeCoord);

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
   */
  void resetQuadPntsMappedCoorAndWheights(const CFDim dimensionality, const CFPolyOrder::Type integratorOrder);

  /**
   * Function that calculates the number of terms in the highest order polynomial
   * that is integrated exactly, based on m_dimensionality and m_integratorOrder,
   * and puts it in m_nbrPolyTerms4ExIntegr.
   */
  void calcNbrPolyTerms4ExIntegr();

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

  /// Vector for the normals in the quadrature points
  std::vector< RealVector > m_quadPntsNormals;

};

//////////////////////////////////////////////////////////////////////////////

  }  // namespace SpectralFV
}  // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_SpectralFV_SimplexGaussIntegrator_hh
