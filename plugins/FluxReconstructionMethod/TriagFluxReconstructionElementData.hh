#ifndef COOLFluiD_FluxReconstructionMethod_TriagFluxReconstructionElementData_hh
#define COOLFluiD_FluxReconstructionMethod_TriagFluxReconstructionElementData_hh

//////////////////////////////////////////////////////////////////////////////



#include "Common/COOLFluiD.hh"
#include "Framework/CFPolyOrder.hh"
#include "Framework/CFGeoShape.hh"
#include "Common/SafePtr.hh"
#include "MathTools/RealVector.hh"
#include "FluxReconstructionMethod/FluxReconstructionElementData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a triangular Flux Reconstruction element
 *
 * @author Ray Vandenhoeck
 * @author Alexander Papen
 *
 */
class TriagFluxReconstructionElementData : public FluxReconstructionElementData {
public:

  /**
   * Default constructor without arguments.
   */
  TriagFluxReconstructionElementData();

  /**
   * Constructor initializing polyOrder.
   */
  TriagFluxReconstructionElementData(CFPolyOrder::Type polyOrder);
  
  /**
   * Constructor initializing polyOrder, solution and flux point distribution
   */
  TriagFluxReconstructionElementData(CFPolyOrder::Type polyOrder, Common::SafePtr< BasePointDistribution > solPntDist, 
				    Common::SafePtr< BasePointDistribution > flxPntDist);

  /**
   * Default destructor.
   */
  ~TriagFluxReconstructionElementData();

protected: // functions

  /**
   * Creates the wheight coordinates of the flux points in a SV face
   * @pre createFluxPolyNodeCoord() and createFaceFluxPntsConn()
   */
  virtual void createFaceFluxPolyNodeWheightCoord();

  /**
   * Creates a list with the different possible connectivities of faces,
   * not taking into account possible symmetries
   * (--> more possible orientations than with function above)
   * @pre createSVFaceNodeConnectivity()
   */
  virtual void createSVFaceNodeConnectivityPerOrientNoSymm();

  /**
   * Computes the polynomial coefficients of the spectral FV basis functions
   * @pre createLocalCVNodeConn() and createLocalNodeCoord() must be called first
   */
  virtual void computePolyCoefs();

  /**
   * Creates the connectivity between SV faces and flux points
   */
  virtual void createFaceFluxPntsConn();

  /**
   * Creates a vector containing the exponents of the terms in the flux basis polynomials.
   */
  virtual void createFluxPolyExponents();

  /**
   * Creates a set of nodes for interpolation with the requested polynomial degree
   */
  void setInterpolationNodeSet(const CFPolyOrder::Type order,std::vector< RealVector >& nodalSet);

  /**
   * sets the convective/diffusive cfl ratio
   */
  void setCFLConvDiffRatio();

  /**
   * create the cell mapped coordinates of a uniform distribution of points on the cell faces (for output)
   */
  void createFaceOutputPntCellMappedCoords();

  /**
   * create the connectivity in a uniform distribution of points on the cell faces (for output)
   */
  void createFaceOutputPntConn();
  
  //Added from the FRElementData:
  
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
   * @todo check if still necessary
   */
  void createFaceFlxPntsFaceLocalCoords();

  /**
   * Creates a vector containing the exponents of the terms in the solution polynomials.
   */
  void createSolPolyExponents();

  /**
   * Computes the polynomial coefficients of the solution polynomial basis functions.
   */
  void computeSolPolyCoefs();

  /**
   * Computes the coordinates used for the initialization of the states.
   */
  void computeInitPntsCoords();

  /**
   * Computes the transformation matrix used for the initialization of the states.
   */
  void computeInitTransfMatrix();

  /**
   * Creates the derivation direction of the flux points
   * @todo Should be deleted intirely
   */
  void createFlxPntDerivDir();

  /**
   * Creates list of internal flux points
   */
  void createIntFlxPntIdxs();
  
  /**
   * Creates face normals
   */
  void createFaceNormals();

  /**
   * Creates the connectivity between faces and flux points,
   * arranged per face connectivity orientation.
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
   * create the dimensions on which the flux must be projected in the flux points
   */
  void createFluxPntFluxDim() {}


}; // end of class TriagFluxReconstructionElementData

//////////////////////////////////////////////////////////////////////////////

  } // namespace FluxReconstructionMethod

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_FluxReconstructionMethod_LineFluxReconstructionElementData_hh
