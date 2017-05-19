#ifndef COOLFluiD_FluxReconstructionMethod_FluxReconstructionElementData_hh
#define COOLFluiD_FluxReconstructionMethod_FluxReconstructionElementData_hh

//////////////////////////////////////////////////////////////////////////////



#include "Common/COOLFluiD.hh"
#include "Framework/CFPolyOrder.hh"
#include "Framework/CFGeoShape.hh"
#include "Framework/BadFormatException.hh"
#include "Framework/CFSide.hh"
#include "MathTools/RealVector.hh"
#include "MathTools/RealMatrix.hh"
#include "MathTools/ZeroDeterminantException.hh"
#include "Common/SafePtr.hh"
#include "FluxReconstructionMethod/BasePointDistribution.hh"

#include "Framework/DataSocketSink.hh"

#include "Framework/BaseMethodStrategyProvider.hh"
#include "FluxReconstructionMethod/FluxReconstructionSolverData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a general Flux Reconstruction element
 *
 * @author Ray Vandenhoeck
 * @author Alexander Papen
 *
 */
class FluxReconstructionElementData {
public:

  /**
   * Default constructor without arguments.
   */
  FluxReconstructionElementData();

  /**
   * Default destructor.
   */
  virtual ~FluxReconstructionElementData();
  
//   /**
//    * Returns the DataSocket's that this command needs as sinks
//    * @return a vector of SafePtr with the DataSockets
//    */
//   std::vector< Common::SafePtr< Framework::BaseDataSocketSink > >
//       needsSockets();

  /**
   * @return m_shape
   */
  CFGeoShape::Type getShape()
  {
    return m_shape;
  }

  /**
   * @return m_dimensionality
   */
  CFDim getDimensionality()
  {
    return m_dimensionality;
  }

  /**
   * @return m_polyOrder
   */
  CFPolyOrder::Type getPolyOrder()
  {
    return m_polyOrder;
  }

  /**
   * Set m_polyOrder
   * @par polyOrder is the new polynomial order
   */
  void setPolyOrder(CFPolyOrder::Type polyOrder);
  
  /**
   * Set flux point distribution
   * @par flxPntDist is the new 1D flux point distribution
   */
  void setFlxPntDistribution(Common::SafePtr< BasePointDistribution > flxPntDist);
  
  /**
   * Set solution point distribution
   * @par solPntDist is the new 1D solution point distribution
   */
  void setSolPntDistribution(Common::SafePtr< BasePointDistribution > solPntDist);

  /**
   * @return number of solution points
   */
  CFuint getNbrOfSolPnts()
  {
    cf_assert(m_solPntsLocalCoords.size() > 0);
    return m_solPntsLocalCoords.size();
  }

  /**
   * @return number of flux points
   */
  CFuint getNbrOfFlxPnts()
  {
    cf_assert(m_flxPntsLocalCoords.size() > 0);
    return m_flxPntsLocalCoords.size();
  }

  /**
   * @return number of internal flux points
   */
  CFuint getNbrOfIntFlxPnts()
  {
//     cf_assert(m_solPntsLocalCoord1D.size() > 0);
//     cf_assert(m_flxPntsLocalCoord1D.size() > 0);
//     CFuint nbrIntFlx = m_dimensionality*(m_flxPntsLocalCoord1D.size()-2);
//     const CFuint dim = static_cast<CFuint>(m_dimensionality);
//     for (CFuint iDim = 1; iDim < dim; ++iDim)
//     {
//       nbrIntFlx *= m_solPntsLocalCoord1D.size();
//     }
    return 0;
  }

  /**
   * @return number of flux points on a face
   */
  CFuint getNbrOfFaceFlxPnts()
  {
    return (getNbrOfFlxPnts() - getNbrOfIntFlxPnts())/getNbrCellFaces();
  }

  /**
   * @return the number of faces
   */
  CFuint getNbrCellFaces()
  {
    cf_assert(m_faceNodeConn.size() > 0);
    return m_faceNodeConn.size();
  }

  /**
   * @return the number of corner nodes
   */
  CFuint getNbrCornerNodes()
  {
    cf_assert(m_cellNodeCoords.size() > 0);
    return m_cellNodeCoords.size();
  }

  Common::SafePtr< std::vector< CFreal > > getSolPntsLocalCoord1D()
  {
    return &m_solPntsLocalCoord1D;
  }
  
  Common::SafePtr< std::vector< CFreal > > getFlxPntsLocalCoord1D()
  {
    return &m_flxPntsLocalCoord1D;
  }

  /**
   * @return m_recCoefsFlxPnts1D
   */
  Common::SafePtr< RealMatrix > getRecCoefsFlxPnts1D()
  {
    return &m_recCoefsFlxPnts1D;
  }

  /**
   * @return m_recCoefsFlxPnts1DOptim
   */
  Common::SafePtr< RealMatrix > getRecCoefsFlxPnts1DOptim()
  {
    return &m_recCoefsFlxPnts1DOptim;
  }

  /**
   * @return m_derivCoefsSolPnts1D
   */
  Common::SafePtr< RealMatrix > getDerivCoefsSolPnts1D()
  {
    return &m_derivCoefsSolPnts1D;
  }

  /**
   * @return m_solPolyDerivCoefsFlxPnts1D
   */
  Common::SafePtr< RealMatrix > getSolPolyDerivCoefsFlxPnts1D()
  {
    return &m_solPolyDerivCoefsFlxPnts1D;
  }

  /**
   * @return m_solPolyDerivCoefsSolPnts1D
   */
  Common::SafePtr< RealMatrix > getSolPolyDerivCoefsSolPnts1D()
  {
    return &m_solPolyDerivCoefsSolPnts1D;
  }

  /**
   * @return m_solPntsLocalCoords
   */
  Common::SafePtr< std::vector< RealVector > > getSolPntsLocalCoords()
  {
    return &m_solPntsLocalCoords;
  }

  /**
   * @return m_flxPntsLocalCoords
   */
  Common::SafePtr< std::vector< RealVector > > getFlxPntsLocalCoords()
  {
    return &m_flxPntsLocalCoords;
  }

  /**
   * @return m_faceFlxPntsFaceLocalCoords
   */
  Common::SafePtr< std::vector< RealVector > > getFaceFlxPntsFaceLocalCoords()
  {
    return &m_faceFlxPntsFaceLocalCoords;
  }

  /**
   * @return m_flxPntMatrixIdxForReconstruction
   */
  Common::SafePtr< std::vector< CFuint > > getFlxPntMatrixIdxForReconstruction()
  {
    return &m_flxPntMatrixIdxForReconstruction;
  }

  /**
   * @return m_solPntIdxsForReconstruction
   */
  Common::SafePtr< std::vector< std::vector< CFuint > > > getSolPntIdxsForReconstruction()
  {
    return &m_solPntIdxsForReconstruction;
  }

  /**
   * @return m_solPntIdxsForRecOptim
   */
  Common::SafePtr< std::vector< std::vector< CFuint > > > getSolPntIdxsForRecOptim()
  {
    return &m_solPntIdxsForRecOptim;
  }

  /**
   * @return m_solPntMatrixIdxForDerivation
   */
  Common::SafePtr< std::vector< std::vector< CFuint > > > getSolPntMatrixIdxForDerivation()
  {
    return &m_solPntMatrixIdxForDerivation;
  }

  /**
   * @return m_flxPntMatrixIdxForDerivation
   */
  Common::SafePtr< std::vector< CFuint > > getFlxPntMatrixIdxForDerivation()
  {
    return &m_flxPntMatrixIdxForDerivation;
  }

  /**
   * @return m_solPntIdxsForDerivation
   */
  Common::SafePtr< std::vector< std::vector< CFuint > > > getSolPntIdxsForDerivation()
  {
    return &m_solPntIdxsForDerivation;
  }

  /**
   * @return m_flxPntIdxsForDerivation
   */
  Common::SafePtr< std::vector< std::vector< std::vector< CFuint > > > > getFlxPntIdxsForDerivation()
  {
    return &m_flxPntIdxsForDerivation;
  }

  /**
   * @return m_flxPntDerivDir
   */
  Common::SafePtr< std::vector< CFuint > > getFlxPntDerivDir()
  {
    return &m_flxPntDerivDir;
  }

  /**
   * @return m_allSolPntIdxs
   */
  Common::SafePtr< std::vector< CFuint > > getAllSolPntIdxs()
  {
    return &m_allSolPntIdxs;
  }

  /**
   * @return m_allFlxPntIdxs
   */
  Common::SafePtr< std::vector< CFuint > > getAllFlxPntIdxs()
  {
    return &m_allFlxPntIdxs;
  }

  /**
   * @return m_intFlxPntIdxs
   */
  Common::SafePtr< std::vector< CFuint > > getIntFlxPntIdxs()
  {
    return &m_intFlxPntIdxs;
  }

  /**
   * @return m_faceFlxPntConn
   */
  Common::SafePtr< std::vector< std::vector< CFuint > > > getFaceFlxPntConn()
  {
    return &m_faceFlxPntConn;
  }

  /**
   * @return m_faceFlxPntConnPerOrient
   */
  Common::SafePtr< std::vector< std::vector< std::vector< CFuint > > > > getFaceFlxPntConnPerOrient()
  {
    return &m_faceFlxPntConnPerOrient;
  }

  /**
   * @return m_faceFlxPntCellMappedCoords
   */
  Common::SafePtr< std::vector< std::vector< RealVector > > > getFaceFlxPntCellMappedCoords()
  {
    return &m_faceFlxPntCellMappedCoords;
  }

  /**
   * @return m_faceFlxPntCellMappedCoordsPerOrient
   */
  Common::SafePtr< std::vector< std::vector< std::vector< RealVector > > > > getFaceFlxPntCellMappedCoordsPerOrient()
  {
    return &m_faceFlxPntCellMappedCoordsPerOrient;
  }

  /**
   * @return m_cellNodeCoords
   */
  Common::SafePtr< std::vector< RealVector > > getCellNodeCoords()
  {
    return &m_cellNodeCoords;
  }

  /**
   * @return m_faceNodeCoords
   */
  Common::SafePtr< std::vector< std::vector< RealVector > > > getFaceNodeCoords()
  {
    return &m_faceNodeCoords;
  }

  /**
   * @return m_faceNodeConn
   */
  Common::SafePtr< std::vector< std::vector< CFuint > > > getFaceNodeConn()
  {
    return &m_faceNodeConn;
  }

  /**
   * @return m_faceMappedCoordDir
   */
  Common::SafePtr< std::vector< CFint > > getFaceMappedCoordDir()
  {
    return &m_faceMappedCoordDir;
  }
  
  /**
   * @return m_faceNormals
   */
  Common::SafePtr< std::vector<RealVector> > getFaceNormals()
  {
    return &m_faceNormals;
  }

  /**
   * @return m_faceNodeConnPerOrient
   */
  Common::SafePtr< std::vector< std::vector< std::vector< CFuint > > > > getFaceNodeConnPerOrient()
  {
    return &m_faceNodeConnPerOrient;
  }

  /**
   * @return m_faceNodeConnPerOrientNoSymm
   */
//   Common::SafePtr< std::vector< std::vector< std::vector< CFuint > > > > getFaceNodeConnPerOrientNoSymm()
//   {
//     return &m_faceNodeConnPerOrientNoSymm;
//   }

  /**
   * @return m_faceNodeCoordsPerOrient
   */
  Common::SafePtr< std::vector< std::vector< std::vector< RealVector > > > > getFaceNodeCoordsPerOrient()
  {
    return &m_faceNodeCoordsPerOrient;
  }

  /**
   * @return m_faceConnPerOrient
   */
  Common::SafePtr< std::vector< std::vector< CFuint > > > getFaceConnPerOrient()
  {
    return &m_faceConnPerOrient;
  }

  /**
   * @return m_faceMappedCoordDirPerOrient
   */
  Common::SafePtr< std::vector< std::vector< CFint > > > getFaceMappedCoordDirPerOrient()
  {
    return &m_faceMappedCoordDirPerOrient;
  }

  /**
   * @return m_faceNodeCoordsPerOrientNoSymm
   */
//   Common::SafePtr< std::vector< std::vector< std::vector< RealVector > > > > getFaceNodeCoordsPerOrientNoSymm()
//   {
//     return &m_faceNodeCoordsPerOrientNoSymm;
//   }

  /**
   * @return m_solPolyExponents
   */
  Common::SafePtr< std::vector< std::vector< CFint > > > getSolPolyExponents()
  {
    return &m_solPolyExponents;
  }

  /**
   * @return m_solPolyCoefs
   */
  Common::SafePtr< std::vector< std::vector< CFreal > > > getSolPolyCoefs()
  {
    return &m_solPolyCoefs;
  }

  /// @return m_initPntsCoords
  Common::SafePtr< std::vector< RealVector > > getInitPntsCoords()
  {
    return &m_initPntsCoords;
  }

  /// @return m_initTransfMatrix
  Common::SafePtr< RealMatrix > getInitTransfMatrix()
  {
    return &m_initTransfMatrix;
  }

  /**
   * @return m_faceIntegrationCoefs
   */
  Common::SafePtr< RealVector > getFaceIntegrationCoefs()
  {
    return &m_faceIntegrationCoefs;
  }

  /**
   * @return m_cellAvgSolCoefs
   */
  Common::SafePtr< RealVector > getCellAvgSolCoefs()
  {
    return &m_cellAvgSolCoefs;
  }

  /**
   * @return m_cellCenterDerivCoefs
   */
  Common::SafePtr< std::vector< RealVector > > getCellCenterDerivCoefs()
  {
    return &m_cellCenterDerivCoefs;
  }

  /**
   * @return m_flxPolyExponents
   */
  Common::SafePtr< std::vector< std::vector< std::vector< CFint > > > > getFlxPolyExponents()
  {
    return &m_flxPolyExponents;
  }

  /**
   * @return m_flxPolyCoefs
   */
  Common::SafePtr< std::vector< std::vector< std::vector< CFreal > > > > getFlxPolyCoefs()
  {
    return &m_flxPolyCoefs;
  }

  /**
   * @return m_coefSolPolyDerivInSolPnts
   */
  Common::SafePtr< std::vector< std::vector< std::vector< CFreal > > > > getCoefSolPolyDerivInSolPnts()
  {
    return &m_coefSolPolyDerivInSolPnts;
  }
  
  /**
   * @return m_coefSolPolyInFlxPnts
   */
  Common::SafePtr< std::vector< std::vector< CFreal > > > getCoefSolPolyInFlxPnts()
  {
    return &m_coefSolPolyInFlxPnts;
  }

  /**
   * @return m_coefSolPolyDerivInFlxPntsOptim
   */
  Common::SafePtr< std::vector< std::vector< std::vector< CFreal > > > > getCoefSolPolyDerivInFlxPntsOptim()
  {
    return &m_coefSolPolyDerivInFlxPntsOptim;
  }

  /**
   * @return m_solPntIdxsSolPolyDerivInFlxPntsOptim
   */
  Common::SafePtr< std::vector< std::vector< std::vector< CFuint > > > > getSolPntIdxsSolPolyDerivInFlxPntsOptim()
  {
    return &m_solPntIdxsSolPolyDerivInFlxPntsOptim;
  }

  /**
   * @return m_cflConvDiffRatio
   */
  CFreal getCFLConvDiffRatio()
  {
    return m_cflConvDiffRatio;
  }

  /**
   * @return m_faceOutputPntFaceMappedCoords
   */
  Common::SafePtr< std::vector< RealVector > > getFaceOutputPntFaceMappedCoords()
  {
    return &m_faceOutputPntFaceMappedCoords;
  }

  /**
   * @return m_faceOutputPntCellMappedCoords
   */
  Common::SafePtr< std::vector< std::vector< RealVector > > > getFaceOutputPntCellMappedCoords()
  {
    return &m_faceOutputPntCellMappedCoords;
  }

  /**
   * @return m_faceOutputPntSolPolyCoef
   */
  Common::SafePtr< std::vector< std::vector< std::vector< CFreal > > > > getFaceOutputPntSolPolyCoef()
  {
    return &m_faceOutputPntSolPolyCoef;
  }

  /**
   * @return m_faceOutputPntSolDerivCoef
   */
  Common::SafePtr< std::vector< std::vector< std::vector< std::vector< CFreal > > > > > getFaceOutputPntSolDerivCoef()
  {
    return &m_faceOutputPntSolDerivCoef;
  }

  /**
   * @return m_faceOutputPntConn
   */
  Common::SafePtr< std::vector< std::vector< CFuint > > > getFaceOutputPntConn()
  {
    return &m_faceOutputPntConn;
  }
  
  /**
   * @return m_flxPntFlxDim
   */
  Common::SafePtr< std::vector< CFuint > > getFluxPntFluxDim()
  {
    return &m_flxPntFlxDim;
  }

  /**
   * This function evaluates the solution polynomials at the given nodes
   * @param nodeLocalCoords is a vector with the local coordinates of the nodes
   * @return the solution polynomial values at the given nodes
   * @pre createSolPolyExponents, createSolPolyExponents
   */
  std::vector< std::vector< CFreal > > getSolPolyValsAtNode(std::vector< RealVector > nodeLocalCoords);

  /**
   * This function evaluates the solution polynomial derivatives at the given nodes
   * @param nodeLocalCoords is a vector with the local coordinates of the nodes
   * @return the solution polynomial values at the given nodes
   * @pre createSolPolyExponents, createSolPolyExponents
   */
  std::vector< std::vector< std::vector< CFreal > > >
      getSolPolyDerivsAtNode(std::vector< RealVector > nodeLocalCoords);

protected: // functions

  /**
   * Calls all the functions below to recreate the entire datastructure.
   */
  void resetFluxReconstructionElementData();

//   /**
//    * create vector with flux points local coordinate in 1D
//    */
//   void createFlxPntsLocalCoord1D();
// 
//   /**
//    * create vector with solution points local coordinate in 1D
//    */
//   void createSolPntsLocalCoord1D();

  /**
   * compute coefficients for reconstruction of solution in flux points in 1D
   * @pre createSolPntsLocalCoord1D(), createFlxPntsLocalCoord1D()
   */
  void computeRecCoefsFlxPnts1D();

  /**
   * compute coefficients for derivative of flux in solution points in 1D
   * @pre createSolPntsLocalCoord1D(), createFlxPntsLocalCoord1D()
   */
  void computeDerivCoefsSolPnts1D();

  /**
   * compute coefficients for derivative of solution polynomial in solution points in 1D
   * @pre createSolPntsLocalCoord1D(), createFlxPntsLocalCoord1D()
   */
  void computeSolPolyDerivCoefsSolPnts1D();

  /**
   * compute coefficients for derivative of solution polynomial in flux points in 1D
   * @pre createSolPntsLocalCoord1D(), createFlxPntsLocalCoord1D()
   */
  void computeSolPolyDerivCoefsFlxPnts1D();

  /**
   * create vector with flux points local coordinates
   * @pre createFlxPntsLocalCoord1D()
   */
  virtual void createFlxPntsLocalCoords() = 0;

  /**
   * create vector with solution points local coordinates
   * @pre createSolPntsLocalCoord1D()
   */
  virtual void createSolPntsLocalCoords() = 0;

  /**
   * create vector with face flux points local coordinates (coordinate system local to face)
   * @pre createSolPntsLocalCoord1D()
   */
  virtual void createFaceFlxPntsFaceLocalCoords() = 0;

  /**
   * Creates a vector containing the exponents of the terms in the solution polynomials.
   */
  virtual void createSolPolyExponents() = 0;

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
   */
  virtual void createFlxPntDerivDir() = 0;

  /**
   * Creates list of all solution points
   */
  void createAllSolPntIdxs();

  /**
   * Creates list of all flux points
   */
  void createAllFlxPntIdxs();

  /**
   * Creates list of internal flux points
   */
  virtual void createIntFlxPntIdxs() = 0;

  /**
   * Creates the connectivity between faces and flux points
   */
  virtual void createFaceFluxPntsConn() = 0;

  /**
   * Creates the connectivity between faces and flux points,
   * arranged per face connectivity orientation.
   */
  virtual void createFaceFluxPntsConnPerOrient() = 0;

  /**
   * create vector with face flux points cell mapped coordinates (coordinate system local to cell)
   * @pre createFaceFluxPntsConn, createFaceFluxPntsConnPerOrient()
   */
  void createFaceFlxPntsCellLocalCoords();

  /**
   * Creates a vector containing the exponents of the terms in the flux polynomials.
   */
  virtual void createFlxPolyExponents() = 0;

  /**
   * Computes the polynomial coefficients of the flux polynomial basis functions.
   */
  void computeFlxPolyCoefs();

  /**
   * Creates the flux point index (the row index in m_recCoefsFlxPnts1D)
   * for the solution reconstruction in each flux point
   */
  virtual void createFlxPntMatrixIdxForReconstruction() = 0;

  /**
   * Creates the solution point indices for the solution reconstruction in each flux point
   */
  virtual void createSolPntIdxsForReconstruction() = 0;

  /**
   * Creates the solution point index (the row index in m_derivCoefsSolPnts1D)
   * for the flux derivation in each solution point
   */
  virtual void createSolPntMatrixIdxForDerivation() = 0;

  /**
   * Creates the flux point index (the column index in m_derivCoefsSolPnts1D)
   * for the flux derivation in each solution point
   */
  virtual void createFlxPntMatrixIdxForDerivation() = 0;

  /**
   * Creates the solution point index for the flux derivation in each flux point
   */
  virtual void createSolPntIdxsForDerivation() = 0;

  /**
   * Creates the flux point index for the flux derivation in each solution point
   */
  virtual void createFlxPntIdxsForDerivation() = 0;

  /**
   * Creates the local coorinates of the cell nodes
   */
  virtual void createCellNodeCoords() = 0;

  /**
   * Creates face-node connectivity
   */
  virtual void createFaceNodeConnectivity() = 0;

  /**
   * Creates face-mapped coordinate direction
   */
  virtual void createFaceMappedCoordDir() = 0;

  /**
   * Creates a list with the different possible connectivities of faces
   * @pre createFaceNodeConnectivity()
   */
  virtual void createFaceNodeConnectivityPerOrient() = 0;

  /**
   * Creates a list with the different possible connectivities of faces,
   * not taking into account possible symmetries
   * (--> more possible orientations than with function above)
   * @pre createFaceNodeConnectivity()
   */
//   virtual void createFaceNodeConnectivityPerOrientNoSymm() = 0;

  /**
   * Creates a vector containing the face node coordinates for each face.
   * @pre computeNodeLocalCoords() and createFaceNodeConnectivity()
   */
  virtual void createFaceNodeCoords();

  /**
   * Creates a vector containing the face node coordinates for each face,
   * arranged per orientation of the face connecitivity.
   * @pre computeNodeLocalCoords() and createFaceNodeConnectivityPerOrient()
   */
  virtual void createFaceNodeCoordsPerOrient();

  /**
   * Creates a vector containing the face node coordinates for each face,
   * arranged per orientation of the face connecitivity (not taking into account any symmetries).
   * @pre computeNodeLocalCoords() and createFaceNodeConnectivityPerOrient()
   */
//   virtual void createFaceNodeCoordsPerOrientNoSymm();

  /**
   * Create the coefficients for the integration over a face
   */
  virtual void createFaceIntegrationCoefs() = 0;

  /**
   * Create the coefficients for the cell average solution
   */
  virtual void createCellAvgSolCoefs() = 0;

  /**
   * Create the coefficients for the cell center derivatives
   */
  virtual void createCellCenterDerivCoefs() = 0;

  /**
   * sets the convective/diffusive cfl ratio
   */
  virtual void setCFLConvDiffRatio() = 0;

  /**
   * Creates a set of nodes for interpolation with the requested polynomial degree
   */
  virtual void setInterpolationNodeSet(const CFPolyOrder::Type order,std::vector< RealVector >& nodalSet) = 0;

  /**
   * create coefficients for computation of solution polynomial derivatives
   */
  void createCoefSolPolyDerivInSolPnts();
  
  /**
   * create coefficients for computation of solution polynomials in the flx pnts
   */
  void createCoefSolPolyInFlxPnts();

  /**
   * create coefficients for optimized computation of solution polynomial derivatives
   */
  void createCoefSolPolyDerivInFlxPntsOptim();

  /**
   * create the cell mapped coordinates of a uniform distribution of points on the cell faces (for output)
   */
  virtual void createFaceOutputPntCellMappedCoords() = 0;

  /**
   * create the solution polynomial values and derivation coefficients
   * in a uniform distribution of points on the cell faces (for output).
   */
  void createFaceOutputPntSolPolyAndDerivCoef();

  /**
   * create the connectivity in a uniform distribution of points on the cell faces (for output)
   */
  virtual void createFaceOutputPntConn() = 0;
  
  /**
   * create the normals for each face
   */
  virtual void createFaceNormals() = 0;

  /**
   * A function that computes the inverse of the matrix A and puts it in the patrix AI.
   * @par A is the matrix to be inverted.
   * @par AI is the matrix where the inverse of A will be put.
   * @note this function can probably be replaced by a linear system solver
   */
  void InvertMatrix(RealMatrix A, RealMatrix& AI);

  /**
   * A function that swaps two given rows in a matrix.
   * @par A is the matrix in which to swap the rows.
   * @par row1 and row2 are the numbers of the rows to be swapped.
   */
  void SwapRows(RealMatrix& A, CFuint row1, CFuint row2);
  
  /**
   * create the dimensions on which the flux must be projected in the flux points
   */
  virtual void createFluxPntFluxDim() = 0;

protected: // protected data

  /// dimensionality
  CFDim m_dimensionality;

  /// shape
  CFGeoShape::Type m_shape;

  /// polynomial order
  CFPolyOrder::Type m_polyOrder;

  /// solution point local coordinate in 1D
  std::vector< CFreal > m_solPntsLocalCoord1D;

  /// flux point local coordinate in 1D
  std::vector< CFreal > m_flxPntsLocalCoord1D;
  
  /// Distribution of solution points
  Common::SafePtr<BasePointDistribution> m_solPntDistribution;
  
  /// Distribution of flux points
  Common::SafePtr<BasePointDistribution> m_flxPntDistribution;

  /// coefficients for solution reconstruction in the flux points
  RealMatrix m_recCoefsFlxPnts1D;

  /// coefficients for solution reconstruction in the flux points in an optimized way
  RealMatrix m_recCoefsFlxPnts1DOptim;

  /// indexes of solution points corresponding to the coefficients
  /// for solution reconstruction in the flux points in an optimized way
  std::vector< std::vector< CFuint > > m_solPntIdxsForRecFlxPnts1DOptim;

  /// coefficients for derivative computation in the solution points
  RealMatrix m_derivCoefsSolPnts1D;

  /// coefficients for computation of solution polynomial derivaties in the flux points
  RealMatrix m_solPolyDerivCoefsFlxPnts1D;

  /// coefficients for computation of solution polynomial derivaties in the solution points
  RealMatrix m_solPolyDerivCoefsSolPnts1D;

  /// solution point local coordinates
  std::vector< RealVector > m_solPntsLocalCoords;

  /// flux point local coordinates
  std::vector< RealVector > m_flxPntsLocalCoords;

  /// face flux point local coordinates (face coordinates)
  std::vector< RealVector > m_faceFlxPntsFaceLocalCoords;

  /// flux point index (row index in m_recCoefsFlxPnts1D) for solution reconstruction
  std::vector< CFuint > m_flxPntMatrixIdxForReconstruction;

  /// solution point indices for solution reconstruction
  std::vector< std::vector< CFuint > > m_solPntIdxsForReconstruction;

  /// solution point indices for solution reconstruction in an optimized way
  std::vector< std::vector< CFuint > > m_solPntIdxsForRecOptim;

  /// solution point index (row index in m_derivCoefsSolPnts1D) for flux derivation
  std::vector< std::vector< CFuint > > m_solPntMatrixIdxForDerivation;

  /// flux point index (column index in m_derivCoefsSolPnts1D) for flux derivation
  std::vector< CFuint > m_flxPntMatrixIdxForDerivation;

  /// solution point indices for flux derivation
  std::vector< std::vector< CFuint > > m_solPntIdxsForDerivation;

  /// flux point indices for flux derivation
  std::vector< std::vector< std::vector< CFuint > > > m_flxPntIdxsForDerivation;

  /// flux point local coordinate derivative direction
  std::vector< CFuint > m_flxPntDerivDir;

  /// all solution points indexes
  std::vector< CFuint > m_allSolPntIdxs;

  /// all flux points indexes
  std::vector< CFuint > m_allFlxPntIdxs;

  /// internal flux points indexes
  std::vector< CFuint > m_intFlxPntIdxs;

  /// local cell face - flux point connectivity
  std::vector< std::vector< CFuint > > m_faceFlxPntConn;

  /// local cell face - flux point connectivity per face connection orientation
  std::vector< std::vector< std::vector< CFuint > > > m_faceFlxPntConnPerOrient;

  /// local cell face - cell flux point mapped coordinates
  std::vector< std::vector< RealVector > > m_faceFlxPntCellMappedCoords;

  /// local cell face - cell flux point mapped coordinates per face connection orientation
  std::vector< std::vector< std::vector< RealVector > > > m_faceFlxPntCellMappedCoordsPerOrient;

  /// cell node coordinates
  std::vector< RealVector > m_cellNodeCoords;

  /// cell face - node coordinates
  std::vector< std::vector< RealVector > > m_faceNodeCoords;

  /// cell face - node connectivity
  std::vector< std::vector< CFuint > > m_faceNodeConn;

  /// cell face - mapped coordinate direction
  std::vector< CFint > m_faceMappedCoordDir;

  /// cell face - node connectivity per orientation
  std::vector< std::vector< std::vector< CFuint > > > m_faceNodeConnPerOrient;

  /// cell face connectivity per orientation
  std::vector< std::vector< CFuint > > m_faceConnPerOrient;

  /// cell face - mapped coordinate direction per orientation
  std::vector< std::vector< CFint > > m_faceMappedCoordDirPerOrient;
  
  /// cell face - mapped coordinate direction per orientation
  std::vector< RealVector > m_faceNormals;

  /// cell face - node connectivity per orientation, not taking into account symmetries
//   std::vector< std::vector< std::vector< CFuint > > > m_faceNodeConnPerOrientNoSymm;

  /// cell face - node coordinates per orientation
  std::vector< std::vector< std::vector< RealVector > > > m_faceNodeCoordsPerOrient;

  /// cell face - node coordinates per orientation
//   std::vector< std::vector< std::vector< RealVector > > > m_faceNodeCoordsPerOrientNoSymm;

  /// polynomial exponents
  std::vector< std::vector< CFint > > m_solPolyExponents;

  /// spectral finite difference polynomial coefficients
  std::vector< std::vector< CFreal > > m_solPolyCoefs;

  /// local coordinates of points used for initialization of the solution
  std::vector< RealVector > m_initPntsCoords;

  /// transformation matrix for initialization of the states
  RealMatrix m_initTransfMatrix;

  /// coefficients for integration over a face
  RealVector m_faceIntegrationCoefs;

  /// coefficients for cell average solution
  RealVector m_cellAvgSolCoefs;

  /// coefficients for cell center gradient (mapped coordinates)
  std::vector< RealVector > m_cellCenterDerivCoefs;

  /// flux polynomial exponents
  std::vector< std::vector< std::vector< CFint > > > m_flxPolyExponents;

  /// flux polynomial coefficients
  std::vector< std::vector< std::vector< CFreal > > > m_flxPolyCoefs;

  /// coefficients for derivatives of solution polynomials in the soltution points
  std::vector< std::vector< std::vector< CFreal > > > m_coefSolPolyDerivInSolPnts;
  
  /// coefficients for derivatives of solution polynomials in the flux points
  std::vector< std::vector< std::vector< CFreal > > > m_coefSolPolyDerivInFlxPnts;
  
  /// coefficients for solution polynomials in the flux points
  std::vector< std::vector < CFreal > > m_coefSolPolyInFlxPnts;

  /// coefficients for derivatives of solution polynomials in the flux points (optimized computation)
  std::vector< std::vector< std::vector< CFreal > > > m_coefSolPolyDerivInFlxPntsOptim;

  /// solution point indexes for optimized computation of derivatives of solution polynomials in the flux points
  std::vector< std::vector< std::vector< CFuint > > > m_solPntIdxsSolPolyDerivInFlxPntsOptim;

  /// ratio between convective and diffusive CFL limit (results from a trial and error procedure...)
  CFreal m_cflConvDiffRatio;

  /// local cell face - output point face mapped coordinates
  std::vector< RealVector > m_faceOutputPntFaceMappedCoords;

  /// local cell face - output point cell mapped coordinates
  std::vector< std::vector< RealVector > > m_faceOutputPntCellMappedCoords;

  /// local cell face - output point solution polynomial values
  std::vector< std::vector< std::vector< CFreal > > > m_faceOutputPntSolPolyCoef;

  /// local cell face - output point solution derivation coefficients
  std::vector< std::vector< std::vector< std::vector< CFreal > > > > m_faceOutputPntSolDerivCoef;

  /// local cell face - output point connectivity
  std::vector< std::vector< CFuint > > m_faceOutputPntConn;
  
  /// dimension on which the flux must be projected in the flux points
  std::vector<CFuint> m_flxPntFlxDim;
  
  
//   /// socket for solution coordinates in 1D
//   Framework::DataSocketSink< std::vector< CFreal > > socket_solCoords1D;
//   
//   /// socket for flux coordinates in 1D
//   Framework::DataSocketSink< std::vector< CFreal > > socket_flxCoords1D;

}; // end of class FluxReconstructionElementData

//////////////////////////////////////////////////////////////////////////////

  } // namespace FluxReconstructionMethod

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_FluxReconstructionMethod_FluxReconstructionElementData_hh
