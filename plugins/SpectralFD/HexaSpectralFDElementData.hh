#ifndef COOLFluiD_SpectralFD_HexaSpectralFDElementData_hh
#define COOLFluiD_SpectralFD_HexaSpectralFDElementData_hh

//////////////////////////////////////////////////////////////////////////////



#include "Common/COOLFluiD.hh"
#include "Framework/CFPolyOrder.hh"
#include "Framework/CFGeoShape.hh"
#include "Common/SafePtr.hh"
#include "MathTools/RealVector.hh"
#include "SpectralFD/SpectralFDElementData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace SpectralFD {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a hexahedral spectral finite difference element
 *
 * @author Kris Van den Abeele
 *
 */
class HexaSpectralFDElementData : public SpectralFDElementData {
public:

  /**
   * Default constructor without arguments.
   */
  HexaSpectralFDElementData();

  /**
   * Constructor initializing polyOrder.
   */
  HexaSpectralFDElementData(CFPolyOrder::Type polyOrder);

  /**
   * Default destructor.
   */
  ~HexaSpectralFDElementData();

protected: // functions
  /**
   * create vector with flux points local coordinates
   * @pre createFlxPntsLocalCoord1D()
   */
  void createFlxPntsLocalCoords();

  /**
   * create vector with solution points local coordinates
   * @pre createSolPntsLocalCoord1D()
   */
  void createSolPntsLocalCoords();

  /**
   * create vector with face flux points local coordinates (coordinate system local to face)
   * @pre createSolPntsLocalCoord1D()
   */
  void createFaceFlxPntsFaceLocalCoords();

  /**
   * Creates a vector containing the exponents of the terms in the solution polynomials.
   */
  void createSolPolyExponents();

  /**
   * Creates the derivation direction of the flux points
   */
  void createFlxPntDerivDir();

  /**
   * Creates list of internal flux points
   */
  void createIntFlxPntIdxs();

  /**
   * Creates the connectivity between faces and flux points
   */
  void createFaceFluxPntsConn();

  /**
   * Creates the connectivity between faces and flux points,
   * arranged per face connectivity orientation.
   * @pre createFaceFluxPntsConn()
   */
  void createFaceFluxPntsConnPerOrient();

  /**
   * Creates a vector containing the exponents of the terms in the flux polynomials.
   */
  void createFlxPolyExponents();

  /**
   * Creates the flux point index (the row index in m_recCoefsFlxPnts1D)
   * for the solution reconstruction in each flux point
   */
  void createFlxPntMatrixIdxForReconstruction();

  /**
   * Creates the solution point indices for the solution reconstruction in each flux point
   */
  void createSolPntIdxsForReconstruction();

   /**
   * Creates the solution point index (the row index in m_derivCoefsSolPnts1D)
   * for the flux derivation in each solution point
   */
  void createSolPntMatrixIdxForDerivation();

  /**
   * Creates the flux point index (the column index in m_derivCoefsSolPnts1D)
   * for the flux derivation in each solution point
   */
  void createFlxPntMatrixIdxForDerivation();

  /**
   * Creates the solution point index for the flux derivation in each flux point
   */
  void createSolPntIdxsForDerivation();

  /**
   * Creates the flux point index for the flux derivation in each solution point
   */
  void createFlxPntIdxsForDerivation();

  /**
   * Creates the local coorinates of the cell nodes
   */
  void createCellNodeCoords();

  /**
   * Creates face-node connectivity
   */
  void createFaceNodeConnectivity();

  /**
   * Creates face-mapped coordinate direction
   */
  void createFaceMappedCoordDir();

  /**
   * Creates a list with the different possible connectivities of faces
   * @pre createFaceNodeConnectivity()
   */
  void createFaceNodeConnectivityPerOrient();

  /**
   * Create the coefficients for the integration over a face
   */
  void createFaceIntegrationCoefs();

  /**
   * Create the coefficients for the cell average solution
   */
  void createCellAvgSolCoefs();

  /**
   * Create the coefficients for the cell center derivatives
   */
  void createCellCenterDerivCoefs();

  /**
   * sets the convective/diffusive cfl ratio
   */
  void setCFLConvDiffRatio();

  /**
   * Creates a set of nodes for interpolation with the requested polynomial degree
   */
  void setInterpolationNodeSet(const CFPolyOrder::Type order,std::vector< RealVector >& nodalSet);

  /**
   * create the cell mapped coordinates of a uniform distribution of points on the cell faces (for output)
   */
  void createFaceOutputPntCellMappedCoords();

  /**
   * create the connectivity in a uniform distribution of points on the cell faces (for output)
   */
  void createFaceOutputPntConn();

}; // end of class HexaSpectralFDElementData

//////////////////////////////////////////////////////////////////////////////

  } // namespace SpectralFD

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_SpectralFD_LineSpectralFDElementData_hh
