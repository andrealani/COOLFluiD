#ifndef COOLFluiD_FluxReconstructionMethod_PrismFluxReconstructionElementData_hh
#define COOLFluiD_FluxReconstructionMethod_PrismFluxReconstructionElementData_hh

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
 * This class represents a Prism flux reconstruction element
 * 
 * @author Rayan Dhib
 *
 */
class PrismFluxReconstructionElementData : public FluxReconstructionElementData {
public:

  /**
   * Default constructor without arguments.
   */
  PrismFluxReconstructionElementData();

  /**
   * Constructor initializing polyOrder.
   */
  PrismFluxReconstructionElementData(CFPolyOrder::Type polyOrder);
  
  /**
   * Constructor initializing polyOrder, solution and flux point distribution
   */
  PrismFluxReconstructionElementData(CFPolyOrder::Type polyOrder, Common::SafePtr< BasePointDistribution > solPntDist, 
				    Common::SafePtr< BasePointDistribution > flxPntDist);

  /**
   * Default destructor.
   */
  ~PrismFluxReconstructionElementData();

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
  
  /**
   * Creates face normals
   */
  void createFaceNormals();

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
  
  /**
   * create the dimensions on which the flux must be projected in the flux points
   */
  void createFluxPntFluxDim();
  
  /**
   * create the vandermonde matrix of the transformation to modal basis
   */
  virtual void createVandermondeMatrix();
  
  /**
   * create the sol/flx and sol/sol dependencies
   */
  void createFlxSolDependencies();

  
  std::vector<CFreal> getPercentage(CFPolyOrder::Type solOrder);

private:
  

  /// solution point local coordinate in 2D
  std::vector<std::vector < CFreal > > solPntsLocalCoord2D;

  /// Compute the orthonormal normalized jacobi poly (for triag)
  CFreal ComputeJacobi(CFuint N,CFuint alpha, CFuint beta,CFreal x);
  
  /// Computes the factorial
  CFreal factorial(CFreal n);


}; // end of class PrismFluxReconstructionElementData

//////////////////////////////////////////////////////////////////////////////

  } // namespace FluxReconstructionMethod

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_FluxReconstructionMethod_LineFluxReconstructionElementData_hh
