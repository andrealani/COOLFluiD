#ifndef COOLFluiD_FluxReconstructionMethod_TriagFluxReconstructionElementData_hh
#define COOLFluiD_FluxReconstructionMethod_TriagFluxReconstructionElementData_hh

//////////////////////////////////////////////////////////////////////////////



#include "Common/COOLFluiD.hh"
#include "Framework/CFPolyOrder.hh"
#include "Framework/CFGeoShape.hh"
#include "Common/SafePtr.hh"
#include "MathTools/RealVector.hh"
#include "FluxReconstructionMethod/FluxReconstructionElementData.hh"
#include "FluxReconstructionMethod/BasePointDistribution.hh"


//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a triangular Flux Reconstruction element
 *
 * @author Ray Vandenhoeck
 * @author Alexander Papen
 * @author Firas Ben Ameur
 * @author Rayan Dhib
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
   * create vector with face flux points local coordinates (coordinate system local to face) per face type
   */
  void createFaceFlxPntsLocalCoordsPerType();

  /**
   * Creates a vector containing the exponents of the terms in the solution polynomials.
   */
  void createSolPolyExponents();
  
  /**
   * Creates a vector containing the exponents of the terms in the node associated base polynomials.
   */
  void createNodePolyExponents();

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
   * Creates face normals
   */
  void createFaceNormals();

  /**
   * Creates the connectivity between faces and flux points,
   * arranged per face connectivity orientation.
   */
  void createFaceFluxPntsConnPerOrient();

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
   * Create the coefficients for the integration over a face per face type
   */
  void createFaceIntegrationCoefsPerType();

  /**
   * Create the coefficients for the cell average solution
   */
  void createCellAvgSolCoefs();

  /**
   * Create the coefficients for the cell center derivatives
   */
  void createCellCenterDerivCoefs();

   void createFlxSolDependencies();
  
  /**
   * create the dimensions on which the flux must be projected in the flux points
   */
  void createFluxPntFluxDim();

  void createVandermondeMatrix() ;

  /**
   * create the connectivity of flux points connectivty to faces
   */
  void createFluxPntsFaceConn();

private:

  std::vector<std::vector < std::vector < CFreal> > > getLocalCoords2D(CFPolyOrder::Type solOrder);
  std::vector<CFreal> getPercentage(CFPolyOrder::Type solOrder);

  /// Compute the orthonormal normalized jacobi poly (for triag)
  CFreal ComputeJacobi(CFuint N,CFuint alpha, CFuint beta,CFreal x);
  
  /// Computes the factorial
  CFreal factorial(CFreal n);

private:
  
  /// Distribution of flux points
  Common::SafePtr<BasePointDistribution> flxPntDist;

  /// solution point local coordinate in 2D
  std::vector<std::vector < CFreal > > solPntsLocalCoord2D;

  std::vector<CFreal> coordsFluxPointsInit;


  /// flux point local coordinate in 2D
  std::vector<std::vector < std::vector < CFreal> > > flxPntsLocalCoord2D;

  /// Distribution of flux points
  //Common::SafePtr<BasePointDistribution> m_flxPntDistribution;

}; // end of class TriagFluxReconstructionElementData

//////////////////////////////////////////////////////////////////////////////

  } // namespace FluxReconstructionMethod

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_FluxReconstructionMethod_LineFluxReconstructionElementData_hh
