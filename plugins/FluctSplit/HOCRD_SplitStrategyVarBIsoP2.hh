#ifndef COOLFluiD_Numerics_FluctSplit_HOCRD_SplitStrategyVarBIsoP2_hh
#define COOLFluiD_Numerics_FluctSplit_HOCRD_SplitStrategyVarBIsoP2_hh

//////////////////////////////////////////////////////////////////////////////

#include "MathTools/CFMat.hh"
#include "Framework/StdTrsGeoBuilder.hh"
#include "FluctSplit/FluctuationSplitStrategy.hh"
#include "FluctSplit/P2Normal.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

/// This class represent a fluctuation splitting strategy of higher order
/// elements which residual is computed with a contour integral.
/// @author Martin Vymazal
/// @author Nabil Abderrahaman
/// @author Tiago Quintino
class HOCRD_SplitStrategyVarBIsoP2 : public FluctuationSplitStrategy {

public: // methods

  /// Constructor.
  HOCRD_SplitStrategyVarBIsoP2(const std::string& name);

  /// Destructor.
  virtual ~HOCRD_SplitStrategyVarBIsoP2();

  /// Configure the object
  virtual void configure ( Config::ConfigArgs& args )
  {
    FluctuationSplitStrategy::configure(args);
  }

  /// Set up private data
  virtual void setup();

  /// Unsetup private data
  virtual void unsetup();

  /// Returns the DataSocket's that this numerical strategy needs as sinks
  /// @return a vector of SafePtr with the DataSockets
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();

  /// Compute the fluctuation
  /// @param residual the residual for each variable to distribute in each state
  virtual void computeFluctuation(std::vector<RealVector>& residual);

protected: // methods

  /// Compute the integral of the fluxes
  void computeHOCurvedFluctuation();
  
 /** The function getVelQuadPoints saves in 'vel' the advection velocity 
 * at quadrature points 'coord'
 * @param vel Vector where the velocity will be saved
 * @param coord Coordinates of the quadrature points
 * @note both are matrix of dimension [spacedimension, number of quad points]
 **/
  void getVelQuadPoints(std::vector<RealVector>& vel, const std::vector<Framework::Node*>& coord);

 /** The function getJacQuadPoints saves in 'jac' the jacobian
  * at quadrature points 'coord'
  * @param jac RealVector where the jacobian will be saved
  * @param coord Coordinates of the quadrature points
  * @note jac is vector of dimension [number of quadrature points]
  * @note coord is a matrix of dimension [spacedimension, number of quad points]
  **/
  void computeJacQuadPoints(RealVector& jac, const std::vector<Framework::Node*>& coord);

  
  /** The function computesShapeFunctionQuadPoints saves in 'shapeFun' the 
  * value of the shape function at quadrature points 'coord'
  * @param shapeFun Matrix storing the values at the cuadrature points [nodes, quadrature points]
  * @param coord Coordinates of the quadrature points [spacedimension, number of quad points]
  **/
  void computesShapeFunctionQuadPoints(std::vector<RealVector>& shapeFun, std::vector<RealVector>& coord);

  
  /** The function computesGradShapeFunctionQuadPoints saves in 'gradShapeFun' the 
  * value of the gradient at the quadrature points 'coord'
  * @param gradShapeFun 3DMatrix storing the values at the cuadrature points [nodes, quadrature points, spacial direction(x,y)]
  * @param coord Coordinates of the quadrature points [spacedimension, number of quad points]
  **/
  void computesGradShapeFunctionQuadPoints(const std::vector<Framework::Node*>& coord, RealMatrix& dPhidX, RealMatrix& dPhidY);
 
  /** The function computesBetaQuadPoints saves in 'beta' the 
  * value of beta at the quadrature points.
  * @param beta matrix storing the value of beta_i at the quadrature points.[node, quadrature point]
  * @param vel advection velocity at the quadrature points. [spacedimension, number of quad points]
  * @param gradShapeFun 3D matrix containg the gradient of the shape functions at the quadrature points. [node, quadrature points, spacial dimension]
  * */
  void computesBetaQuadPoints(const std::vector<RealVector>& vel,
                              const RealMatrix& dPhidX, const RealMatrix& dPhidY,
                              std::vector<RealVector>& beta);

  /** The function computesNablaF stores in 'nablaF' the 
  * value of the advection flux at the quadrature points.
  * @param nablaF vector storing the value of the flux at the quadrature points.[quadrature point]
  * @param vel advection velocity at the quadrature points. [spacedimension, number of quad points]
  * @param gradShapeFun 3D matrix containg the gradient of the shape functions at the quadrature points. [node, quadrature points, spacial dimension]
  * */
  //void computesNablaF(RealVector& nablaF, std::vector<RealVector> vel, std::vector<std::vector<RealVector > > gradShapeFun);
  void computesNablaF(const std::vector<RealVector>& vel, const RealMatrix& dPhidX, const RealMatrix& dPhidY, RealVector& nablaF);


  /** The function computeResiduals stores in 'resid' the residuals 
  * @param resid matrix storing the value of the residuals.
  * @param jac RealVector where the jacobian at the quadrature points. [number of quadrature points]
  * @param beta matrix containing the value of beta_i at the quadrature points. [node, quadrature point]
  * @param nablaF vector containg the value of the flux at the quadrature points.[quadrature point]
  * @param weight vector containg nodes' weight 
  * */
  //void computeRes(std::vector<RealVector>& resid, RealVector jac, std::vector<RealVector> beta, RealVector nablaF, RealVector weight);
  void computeRes(const RealVector& jac, const RealVector& weight,
                  const std::vector<RealVector>& beta, const RealVector& nablaF,
                  std::vector<RealVector>& resid);

  
  
  
protected: // data

  /// The sockets for the update coefficients
  Framework::DataSocketSink<CFreal> socket_updateCoeff;

  /// solution variable set
  Common::SafePtr<Framework::ConvectiveVarSet> m_solutionVar;

  /// update variable set
  Common::SafePtr<Framework::ConvectiveVarSet> m_updateVar;

  /// matrix storing the unit sized face normals
  std::vector<RealVector> m_unitFaceNormals;

  /// matrix storing the scaled face normals generated on the fly
  std::vector<RealVector> m_scaledFaceNormals;

  /// fluctuation on each sub element
  std::vector<RealVector*> m_phisubT;

  /// direction of the faces
  MathTools::CFMat<CFreal> subelemfacedir;

  /// faces that compose each sub element
  MathTools::CFMat<CFuint> subelemtable;

  /// states that compose each sub face
  MathTools::CFMat<CFuint> subfacetable;

  /// states that compose each sub element
  std::vector<Framework::State*> substates;

  /// residuals for each sub element
  std::vector<RealVector> subresidual;

  /// fluxes on each sub face
  std::vector<RealVector> faceflux;

  /// quadrature points per face
  RealVector qd0; /// Quadrature points: Coordinate x
  RealVector qd1; /// Quadrature points: Coordinate y
  RealVector wqd; /// Quadrature points: Points weight

  /// Linear shape functions
  CFreal L1, L2, L3;

  /// Coordinates of reference triangle
  RealVector xi_ref;
  RealVector eta_ref;

  /// Shape function values in quadrature points
  /// This is a matrix of size nbQdPts x nbShapeFunctions
  /// Each row contains all shape functions evaluated at
  /// one quadrature points
  RealMatrix sfValues;

  /// Derivatives of shape functions in quadrature points
  /// with respect to xi
  /// One row contains dphi_j/dxi for all shape functions phi_j
  /// evaluated at one quadrature point
  RealMatrix sfDerivativesXi;

  /// Derivatives of shape functions in quadrature points
  /// with respect to eta
  RealMatrix sfDerivativesEta;

  /// Derivatives of shape functions in quadrature points
  /// with respect to X and Y (i.e. in physical space)
  RealMatrix sfDerivativesX;

  RealMatrix sfDerivativesY;

  /// states at quadrature points
  std::vector<Framework::State*> qdstates;

  /// extra values at quadrature points
  std::vector<RealVector*> m_qdExtraVars;

  /// coordinates of these states
  std::vector<Framework::Node*> qdnodes;

  /// temporary face normal
  RealVector facenormal;

  /// temporary subcell normals
  InwardNormalsData * m_subcell_normals;

  RealMatrix matrix_face_norms;
  RealMatrix matrix_node_norms;

  RealVector vector_face_areas;
  RealVector vector_node_areas;

  ///Object to compute normals of P2P2 triangle on the fly
  P2Normal m_CP2N;

private: // data

  ///Number of Gauss points used in the quadrature
  enum { nbQdPts = 3 };

  /// the single splitter
  Common::SafePtr<Splitter> m_splitter;

}; // class HOCRD_SplitStrategyVarBIsoP2

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FluctSplit_HOCRD_SplitStrategyVarBIsoP2_hh
