#ifndef COOLFluiD_SpectralFV_SpectralFVElementData_hh
#define COOLFluiD_SpectralFV_SpectralFVElementData_hh

//////////////////////////////////////////////////////////////////////////////



#include "Common/COOLFluiD.hh"
#include "Framework/CFPolyOrder.hh"
#include "Framework/CFGeoShape.hh"
#include "Framework/BadFormatException.hh"
#include "Framework/CFSide.hh"
#include "Common/NotImplementedException.hh"
#include "Common/SafePtr.hh"
#include "MathTools/RealVector.hh"
#include "MathTools/RealMatrix.hh"
#include "MathTools/ZeroDeterminantException.hh"
#include "SpectralFV/SimplexGaussIntegrator.hh"
#include "SpectralFV/TensorProductGaussIntegrator.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace SpectralFV {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a general spectral finite volume element
 *
 * @author Kris Van den Abeele
 *
 */
class SpectralFVElementData {
public:

  /**
   * Default constructor without arguments.
   */
  SpectralFVElementData();

  /**
   * Default destructor.
   */
  virtual ~SpectralFVElementData();

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
   * @return number of local nodes
   */
  CFuint getNbrOfLocalNodes()
  {
    return m_localNodeCoord.size();
  }

  /**
   * @return number of local faces
   */
  CFuint getNbrOfLocalFaces()
  {
    return m_localFaceNodeConn.size();
  }

  /**
   * @return number of local internal faces
   */
  CFuint getNbrOfLocalIntFaces()
  {
    return m_localIntFaceCVConn.size();
  }

  /**
   * @return the number of SV faces
   */
  CFuint getNbrSVFaces()
  {
    return m_faceFracCV.size();
  }

  /**
   * @return number of control volumes
   */
  CFuint getNbrOfCVs()
  {
    return m_localCVNodeConn.size();
  }

  /**
   * @return m_localNodeCoord
   */
  Common::SafePtr< std::vector< RealVector > > getLocalNodeCoord()
  {
    return &m_localNodeCoord;
  }

  /**
   * @return m_localFaceNodeConn
   */
  Common::SafePtr< std::vector< std::vector< CFuint > > > getLocalFaceNodeConn()
  {
    return &m_localFaceNodeConn;
  }

  /**
   * @return m_localCVNodeConn
   */
  Common::SafePtr< std::vector< std::vector< CFuint > > > getLocalCVNodeConn()
  {
    return &m_localCVNodeConn;
  }

  /**
   * @return m_localIntFaceCVConn
   */
  Common::SafePtr< std::vector< std::vector< CFuint > > > getLocalIntFaceCVConn()
  {
    return &m_localIntFaceCVConn;
  }

  /**
   * @return m_intFaceQuadPntNorm
   */
  Common::SafePtr< std::vector< std::vector< RealVector > > > getLocalIntFaceNorm()
  {
    return &m_intFaceQuadPntNorm;
  }

  /**
   * @return m_extFaceLocalNorm
   */
  Common::SafePtr< std::vector< RealVector > > getExtFaceLocalNorm()
  {
    return &m_extFaceLocalNorm;
  }

  /**
   * @return m_volFracCV
   */
  Common::SafePtr< std::vector< CFreal > > getVolFracCV()
  {
    return &m_volFracCV;
  }

  /**
   * @return m_invVolFracCV
   */
  Common::SafePtr< std::vector< CFreal > > getInvVolFracCV()
  {
    return &m_invVolFracCV;
  }

  /**
   * @return m_localExtFaceCVConn
   */
  Common::SafePtr< std::vector< CFuint > > getLocalExtFaceCVConn()
  {
    return &m_localExtFaceCVConn;
  }

  /**
   * @return m_svFaceLocalExtFaceConn
   */
  Common::SafePtr< std::vector< std::vector< CFuint > > > getSVFaceLocalExtFaceConn()
  {
    return &m_svFaceLocalExtFaceConn;
  }

  /// @return m_extSVFaceCVConn
  Common::SafePtr< std::vector< std::vector< CFuint > > > getExtSVFaceCVConn()
  {
    return &m_extSVFaceCVConn;
  }

  /// @return m_extSVFaceCVConnPerOrient
  Common::SafePtr< std::vector< std::vector< std::vector< CFuint > > > > getExtSVFaceCVConnPerOrient()
  {
    return &m_extSVFaceCVConnPerOrient;
  }

  /// @return m_eSVFaceLocESVFaceConnPerOrient
  Common::SafePtr< std::vector< std::vector< CFuint > > > getESVFaceLocESVFaceConnPerOrient()
  {
    return &m_eSVFaceLocESVFaceConnPerOrient;
  }

  /// @return m_svNodeCoords
  Common::SafePtr< std::vector< RealVector > >  getSVNodeCoords()
  {
    return &m_svNodeCoords;
  }

  /// @return m_svFaceNodeCoords
  Common::SafePtr< std::vector< std::vector< RealVector > > > getSVFaceNodeCoords()
  {
    return &m_svFaceNodeCoords;
  }

  /// @return m_svFaceNodeConn
  Common::SafePtr< std::vector< std::vector< CFuint > > > getSVFaceNodeConn()
  {
    return &m_svFaceNodeConn;
  }

  /// @return m_svFaceNodeConnPerOrientation
  Common::SafePtr< std::vector< std::vector< std::vector< CFuint > > > > getSVFaceNodeConnPerOrientation()
  {
    return &m_svFaceNodeConnPerOrientation;
  }

  /// @return m_svFaceNodeConnPerOrientationNoSymm
  Common::SafePtr< std::vector< std::vector< std::vector< CFuint > > > > getSVFaceNodeConnPerOrientationNoSymm()
  {
    return &m_svFaceNodeConnPerOrientationNoSymm;
  }

  /// @return m_svFaceNodeCoordsPerOrient
  Common::SafePtr< std::vector< std::vector< std::vector< RealVector > > > > getSVFaceNodeCoordsPerOrient()
  {
    return &m_svFaceNodeCoordsPerOrient;
  }

  /// @return m_svFaceNodeCoordsPerOrientNoSymm
  Common::SafePtr< std::vector< std::vector< std::vector< RealVector > > > > getSVFaceNodeCoordsPerOrientNoSymm()
  {
    return &m_svFaceNodeCoordsPerOrientNoSymm;
  }

  /**
   * @return m_faceFracCV
   */
  Common::SafePtr< std::vector< std::vector< CFreal > > > getFaceFracCV()
  {
    return &m_faceFracCV;
  }

  /**
   * @return m_invVolFracFaceCV
   */
  Common::SafePtr< std::vector< std::vector< CFreal > > > getInvVolFracFaceCV()
  {
    return &m_invVolFracFaceCV;
  }

  /**
   * @return m_invVolFracFaceCVPerOrient
   */
  Common::SafePtr< std::vector< std::vector< std::vector< CFreal > > > > getInvVolFracFaceCVPerOrient()
  {
    return &m_invVolFracFaceCVPerOrient;
  }

  /**
   * @return m_polyCoefSFV
   */
  Common::SafePtr< std::vector< std::vector< CFreal > > > getPolyCoefSFV()
  {
    return &m_polyCoefSFV;
  }

  /**
   * @return m_polyExponents
   */
  Common::SafePtr< std::vector< std::vector< CFint > > > getPolyExponents()
  {
    return &m_polyExponents;
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
   * @return m_intFaceQuadPntCoords
   */
  Common::SafePtr< std::vector< std::vector< RealVector > > > getIntFaceQuadPntCoords()
  {
    return &m_intFaceQuadPntCoords;
  }

  /**
   * @return m_intFaceQuadWheights
   */
  Common::SafePtr< std::vector< std::vector< CFreal > > > getIntFaceQuadWheights()
  {
    return &m_intFaceQuadWheights;
  }

  /**
   * @return m_intFaceQuadPntPolyVals
   */
  Common::SafePtr< std::vector< std::vector< std::vector< CFreal > > > > getIntFaceQuadPntPolyVals()
  {
    return &m_intFaceQuadPntPolyVals;
  }

  /**
   * @return m_intQuadPntCoords
   */
  Common::SafePtr< std::vector< RealVector > > getIntQuadPntCoords()
  {
    return &m_intQuadPntCoords;
  }

  /**
   * @return m_intQuadWheights
   */
  Common::SafePtr< std::vector< CFreal > > getIntQuadWheights()
  {
    return &m_intQuadWheights;
  }

  /**
   * @return m_intQuadPntPolyVals
   */
  Common::SafePtr< std::vector< std::vector< CFreal > > > getIntQuadPntPolyVals()
  {
    return &m_intQuadPntPolyVals;
  }

  /**
   * @return m_extFaceQPntWheightCoords
   */
  Common::SafePtr< std::vector< std::vector< RealVector > > >  getExtFaceQPntWheightCoords()
  {
    return &m_extFaceQPntWheightCoords;
  }

  /**
   * @return m_extQPntWheightCoords
   */
  Common::SafePtr< std::vector< RealVector > >  getExtQPntWheightCoords()
  {
    return &m_extQPntWheightCoords;
  }

  /**
   * @return m_extFaceQuadPntCoords
   */
  Common::SafePtr< std::vector< std::vector< std::vector< RealVector > > > > getExtFaceQuadPntCoords()
  {
    return &m_extFaceQuadPntCoords;
  }

  /**
   * @return m_extFaceQuadWheights
   */
  Common::SafePtr< std::vector< std::vector< std::vector< CFreal > > > > getExtFaceQuadWheights()
  {
    return &m_extFaceQuadWheights;
  }

  /**
   * @return m_extFaceQuadPntPolyVals
   */
  Common::SafePtr< std::vector< std::vector< std::vector< std::vector< CFreal > > > > >
    getExtFaceQuadPntPolyVals()
  {
    return &m_extFaceQuadPntPolyVals;
  }

  /**
   * @return m_extQPntCoords
   */
  Common::SafePtr< std::vector< std::vector< RealVector > > > getExtQPntCoords()
  {
    return &m_extQPntCoords;
  }

  /**
   * @return m_extQWheights
   */
  Common::SafePtr< std::vector< std::vector< CFreal > > > getExtQWheights()
  {
    return &m_extQWheights;
  }

  /**
   * @return m_extQPntPolyVals
   */
  Common::SafePtr< std::vector< std::vector< std::vector< CFreal > > > > getExtQPntPolyVals()
  {
    return &m_extQPntPolyVals;
  }

  /**
   * @return m_extFaceQPntCoordsPerOrient
   */
  Common::SafePtr< std::vector< std::vector< std::vector< std::vector< RealVector > > > > >
                                                                      getExtFaceQPntCoordsPerOrient()
  {
    return &m_extFaceQPntCoordsPerOrient;
  }

  /**
   * @return m_extFaceQWheightsPerOrient
   */
  Common::SafePtr< std::vector< std::vector< std::vector< CFreal > > > > getExtFaceQWheightsPerOrient()
  {
    return &m_extFaceQWheightsPerOrient;
  }

  /**
   * @return m_extFaceQPntPolyValsPerOrient
   */
  Common::SafePtr< std::vector< std::vector< std::vector< std::vector< std::vector< CFreal > > > > > >
                                                                        getExtFaceQPntPolyValsPerOrient()
  {
    return &m_extFaceQPntPolyValsPerOrient;
  }

  /**
   * @return m_extFaceQPntPolyValsPerOrientNoSymm
   */
  Common::SafePtr< std::vector< std::vector< std::vector< std::vector< std::vector< CFreal > > > > > >
      getExtFaceQPntPolyValsPerOrientNoSymm()
      {
        return &m_extFaceQPntPolyValsPerOrientNoSymm;
      }

  /**
   * @return m_extQPntCoordsPerOrient
   */
  Common::SafePtr< std::vector< std::vector< std::vector< RealVector > > > > getExtQPntCoordsPerOrient()
  {
    return &m_extQPntCoordsPerOrient;
  }

  /**
   * @return m_extQWheightsPerOrient
   */
  Common::SafePtr< std::vector< std::vector< CFreal > > > getExtQWheightsPerOrient()
  {
    return &m_extQWheightsPerOrient;
  }

  /**
   * @return m_extQPntPolyValsPerOrient
   */
  Common::SafePtr< std::vector< std::vector< std::vector< std::vector< CFreal > > > > >
                                                                        getExtQPntPolyValsPerOrient()
  {
    return &m_extQPntPolyValsPerOrient;
  }

  /**
   * @return m_extQPntPolyValsPerOrientNoSymm
   */
  Common::SafePtr< std::vector< std::vector< std::vector< std::vector< CFreal > > > > >
      getExtQPntPolyValsPerOrientNoSymm()
      {
        return &m_extQPntPolyValsPerOrientNoSymm;
      }

  /**
   * @return m_svFaceAvgPolyVals
   */
  Common::SafePtr< std::vector< std::vector< CFreal > > > getSVFaceAvgPolyVals()
  {
    return &m_svFaceAvgPolyVals;
  }

  /**
   * @return m_svFaceAvgPolyValsPerOrient
   */
  Common::SafePtr< std::vector< std::vector< std::vector< CFreal > > > > getSVFaceAvgPolyValsPerOrient()
  {
    return &m_svFaceAvgPolyValsPerOrient;
  }

  /**
   * @return m_fluxPolyNodeCoord
   */
  Common::SafePtr< std::vector< RealVector > > getFluxPolyNodeCoord()
  {
    return &m_fluxPolyNodeCoord;
  }

  /**
   * @return m_faceFluxPolyNodeWheightCoord
   */
  Common::SafePtr< std::vector< RealVector > > getFaceFluxPolyNodeWheightCoord()
  {
    return &m_faceFluxPolyNodeWheightCoord;
  }

  /**
   * @return m_fluxPolyExponents
   */
  Common::SafePtr< std::vector< std::vector< CFint > > > getFluxPolyExponents()
  {
    return &m_fluxPolyExponents;
  }

  /**
   * @return m_fluxPolyCoefs
   */
  Common::SafePtr< std::vector< std::vector< CFreal > > > getFluxPolyCoefs()
  {
    return &m_fluxPolyCoefs;
  }

  /**
   * @return m_solInFluxPntCoef
   */
  Common::SafePtr< std::vector< std::vector< CFreal > > > getSolInFluxCoef()
  {
    return &m_solInFluxPntCoef;
  }

  /**
   * @return m_volTermTensor
   */
  Common::SafePtr< std::vector< std::vector< std::vector< CFreal > > > > getVolTermTensor()
  {
    return &m_volTermTensor;
  }

  /**
   * @return number of flux points
   */
  CFuint getNbrOfFlxPnts()
  {
    return m_fluxPolyNodeCoord.size();
  }

  /**
   * @return m_faceFlxPntsConn
   */
  Common::SafePtr< std::vector< std::vector< CFuint > > > getFaceFlxPntsConn()
  {
    return &m_faceFlxPntsConn;
  }

  /**
   * @return m_solInFaceFluxPntCoef
   */
  Common::SafePtr< std::vector< std::vector< std::vector< CFreal > > > > getSolInFaceFluxPntCoef()
  {
    return &m_solInFaceFluxPntCoef;
  }

  /**
   * @return m_solInFaceFluxPntCoefPerOrient

   */
  Common::SafePtr< std::vector< std::vector< std::vector< std::vector< CFreal > > > > >
    getSolInFaceFluxPntCoefPerOrient()
  {
    return &m_solInFaceFluxPntCoefPerOrient;
  }

  /**
   * @return m_solInFaceFluxPntCoefPerOrientNoSymm

   */
  Common::SafePtr< std::vector< std::vector< std::vector< std::vector< CFreal > > > > >
      getSolInFaceFluxPntCoefPerOrientNoSymm()
      {
        return &m_solInFaceFluxPntCoefPerOrientNoSymm;
      }

  /**
   * @return m_cvExtFaceFluxCoef
   */
  Common::SafePtr< std::vector< std::vector< CFreal > > > getCVExtFaceFluxCoef()
  {
    return &m_cvExtFaceFluxCoef;
  }

  /**
   * @return m_avgSolInSVFaceCoef
   */
  Common::SafePtr< std::vector< CFreal > > getAvgSolInSVFaceCoef()
  {
    return &m_avgSolInSVFaceCoef;
  }

  /**
   * @return number of face flux points
   */
  CFuint getNbrOfFaceFlxPnts()
  {
    cf_assert(m_faceFlxPntsConn.size() > 0);
    return m_faceFlxPntsConn[0].size();
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
   * This function evaluates the SV-polynomials at the given nodes
   * @param nodeLocalCoords is a vector with the local coordinates of the nodes
   * @return the SV-polynomial values at the given nodes
   */
  std::vector< std::vector< CFreal > > getSVPolyValsAtNode(std::vector< RealVector > nodeLocalCoords);

protected: // functions
  /**
   * vector with the local coordinates of the SV nodes
   */
  virtual void computeSVNodeLocalCoords() = 0;

  /**
   * Creates the vector containing the local node coordinates of a SV.
   */
  virtual void createLocalNodeCoord() = 0;

  /**
   * Creates the connectivity between a local face and its local nodes in a SV.
   * External faces are put first.
   */
  virtual void createLocalFaceNodeConn() = 0;

  /**
   * Creates local CV - node connectivity.
   */
  virtual void createLocalCVNodeConn() = 0;

  /**
   * Computes the local internal face normals.
   * @pre createLocalNodeCoord() and createLocalFaceNodeConn()
   */
  virtual void computeLocalFaceNormals() = 0;

  /**
   * Creates the connectivity between a local internal face and its local neighbouring CVs.
   * @pre createLocalNodeCoord(), createLocalFaceNodeConn(),
   *      createLocalCVNodeConn() and computeLocalFaceNormals()
   */
  virtual void createLocalIntFaceCVConn() = 0;

  /**
   * Computes the external face local normals
   */
  virtual void computeExtFaceLocalNormals() = 0;

  /**
   * Computes the volume fraction that each CV occupates inside the SV.
   * @pre createLocalCVNodeConn() and createLocalNodeCoord()
   */
  virtual void computeVolumeFractionsOfCVs() = 0;

  /**
   * Computes the inverse volume fraction of each CV
   * @pre computeVolumeFractionsOfCVs()
   */
  void computeInverseVolFracsOfCVs();

  /**
   * Creates the connectivity between a local external face and its local neighbouring CV.
   */
  virtual void createLocalExtFaceCVConn() = 0;

  /**
   * Creates the connectivity between SV faces and local external faces.
   */
  virtual void createSVFaceLocalExtFaceConn() = 0;

  /**
   * Creates the external face CV connectivity
   * @pre createLocalExtFaceCVConn() and createSVFaceLocalExtFaceConn()
   */
  void createExtFaceCVConn();

  /**
   * Creates the external face CV connectivity per orientation
   * @pre createLocalExtFaceCVConn() and createSVFaceLocalExtFaceConn()
   * @warning should be overwritten for elements other than LINEs and TRIAGs
   */
  virtual void createExtFaceCVConnPerOrient();

  /**
   * Creates SV face-node connectivity
   */
  virtual void createSVFaceNodeConnectivity() = 0;

  /**
   * Creates a list with the different possible connectivities of faces
   * @pre createSVFaceNodeConnectivity()
   */
  virtual void createSVFaceNodeConnectivityPerOrient() = 0;

  /**
   * Creates a list with the different possible connectivities of faces,
   * not taking into account possible symmetries
   * (--> more possible orientations than with function above)
   * @pre createSVFaceNodeConnectivity()
   */
  virtual void createSVFaceNodeConnectivityPerOrientNoSymm() = 0;

  /**
   * Creates a vector containing the face node coordinates for each SV face.
   * @pre computeSVNodeLocalCoords() and createSVFaceNodeConnectivity()
   */
  virtual void createSVFaceNodeCoords();

  /**
   * Creates a vector containing the face node coordinates for each SV face, arranged per orientation of the
   * face connecitivity.
   * @pre computeSVNodeLocalCoords() and createSVFaceNodeConnectivityPerOrient()
   */
  virtual void createSVFaceNodeCoordsPerOrient();

  /**
   * Creates a vector containing the face node coordinates for each SV face, arranged per orientation of the
   * face connecitivity (not taking into account any symmetries).
   * @pre computeSVNodeLocalCoords() and createSVFaceNodeConnectivityPerOrient()
   */
  virtual void createSVFaceNodeCoordsPerOrientNoSymm();

  /**
   * Computes the fraction of the face that each CV occupates at a SV boundary.
   * @pre createLocalFaceNodeConn()
   */
  virtual void computeFaceFractionsOfCVs() = 0;

  /**
   * Computes the volume fraction of that each CV at a SV boundary occupates.
   * @pre computeInverseVolFracsOfCVs() and createLocalExtFaceCVConn()
   */
  void computeInvVolFracsOfFaceCVs();

  /**
   * Computes the volume fraction of that each CV at a SV boundary occupates, per orientation.
   * @pre computeInverseVolFracsOfCVs() and createLocalExtFaceCVConnPerOrient()
   */
  void computeInvVolFracsOfFaceCVsPerOrient();

  /**
   * Creates a vector containing the exponents of the terms in the Spectral FV polynomials.
   */
  virtual void createPolyExponents() = 0;

  /**
   * Computes the polynomial coefficients of the spectral FV basis functions.
   */
  virtual void computePolyCoefs() = 0;

  /**
   * Computes the coordinates used for the initialization of the states.
   */
  void computeInitPntsCoords();

  /**
   * Computes the transformation matrix used for the initialization of the states.
   */
  void computeInitTransfMatrix();

  /**
   * Computes coordinates, wheights and polynomial values in the quadrature points on internal faces
   */
  virtual void computeIntFaceQuadPntsData();

  /**
   * Computes coordinates, wheights and polynomial values in the quadrature points on internal faces
   * but ordered in one list (instead of a separate index for faces and quadrature points)
   * @pre computeIntFaceQuadPntsData()
   */
  void computeIntQuadPntsData();

  /**
   * Computes the the external face node coordinates local to the SV face
   */
  virtual void createExtFaceNodeLocalCoords() = 0;

  /**
   * Computes coordinates, wheights and polynomial values in the quadrature points on external faces
   * @pre all data for external faces and SV polynomials has to be computed first
   */
  virtual void computeExtFaceQuadPntsData();

  /**
   * Computes coordinates, wheights and polynomial values in the external quadrature points
   * @pre computeExtFaceQuadPntsData()
   */
  virtual void computeExtQuadPntsData();

  /**
   * Computes coordinates, wheights and polynomial values in the quadrature points on external faces
   * ordered per orientation of the face connectivity
   * @pre all data for external faces and SV polynomials has to be computed first
   */
  virtual void computeExtFaceQuadPntsDataPerOrientation();

  /**
   * Computes polynomial values in the quadrature points on external faces
   * ordered per orientation of the face connectivity, not taking into account any symmetries
   * @pre all data for external faces and SV polynomials has to be computed first
   */
  virtual void computeExtFaceQuadPntsDataPerOrientationNoSymm();

  /**
   * Computes coordinates, wheights and polynomial values in the external quadrature points
   * ordered per orientation of the face connectivity
   * @pre computeExtFaceQuadPntsDataPerOrientation()
   */
  virtual void computeExtQuadPntsDataPerOrientation();

  /**
   * Computes polynomial values in the external quadrature points
   * ordered per orientation of the face connectivity, not taking into account any symmetries
   * @pre computeExtFaceQuadPntsDataPerOrientation()
   */
  virtual void computeExtQuadPntsDataPerOrientationNoSymm();

  /**
   * Computes the polynomial values in the face centers
   * @pre computePolyCoefs() and createSVFaceNodeCoords()
   */
  void computeFaceCenPolyVals();

  /**
   * Computes the polynomial values in the face centers, ordered per orientation of the face connectivity
   * @pre computePolyCoefs() and createSVFaceNodeCoordsPerOrient()
   */
  void computeFaceCenPolyValsPerOrient();

  /**
   * Computes the polynomial values averaged over the faces
   * @pre computePolyCoefs() and createSVFaceNodeCoords()
   */
  void computeFaceAvgPolyVals();

  /**
   * Computes the polynomial values averaged over the faces, ordered per orientation of the face connectivity
   * @pre computePolyCoefs() and createSVFaceNodeCoordsPerOrient()
   */
  void computeFaceAvgPolyValsPerOrient();

  /**
   * Creates the vector containing the local flux polynomial node coordinates of a SV.
   */
  virtual void createFluxPolyNodeCoord();

  /**
   * Creates the connectivity between SV faces and flux points
   */
  virtual void createFaceFluxPntsConn() = 0;

  /**
   * Creates the wheight coordinates of the flux points in a SV face
   * @pre createFluxPolyNodeCoord() and createFaceFluxPntsConn()
   */
  virtual void createFaceFluxPolyNodeWheightCoord() = 0;

  /**
   * Creates a vector containing the exponents of the terms in the flux basis polynomials.
   */
  virtual void createFluxPolyExponents() = 0;

  /**
   * Computes the flux polynomial coefficients of the spectral FV basis functions.
   * @pre createFluxPolyNodeCoord() and createFluxPolyExponents()
   */
  virtual void createFluxPolyCoef();

  /**
   * Computes the coefficients for the computation of the numerical solution in the flux points.
   * @pre createFluxPolyNodeCoord() and computePolyCoefs()
   */
  virtual void createSolInFlxPntsCoefs();

  /**
   * Computes the tensor for the evaluation of the volume terms
   * @pre createFluxPolyCoef must be called first, and basic SV data must be built as well.
   */
  virtual void createVolumeTermsTensor();

  /**
   * Computes the coefficients for the computation of the numerical solution in the flux points that lie on a
   * SV face.
   * @pre createSolInFlxPntsCoefs() and createFaceFluxPntsConn()
   */
  virtual void createSolInFaceFlxPntsCoefs();

  /**
   * Computes the coefficients for the computation of the numerical solution in the flux points that lie on a
   * SV face, arranged per face connectivity orientation.
   * @pre createFluxPolyCoef() and createFaceFluxPolyNodeWheightCoord().
   */
  virtual void createSolInFaceFlxPntsCoefsPerOrient();

  /**
   * Computes the coefficients for the computation of the numerical solution in the flux points that lie on a
   * SV face, arranged per face connectivity orientation, not taking into account any symmetries.
   * @pre createFluxPolyCoef() and createFaceFluxPolyNodeWheightCoord().
   */
  virtual void createSolInFaceFlxPntsCoefsPerOrientNoSymm();

  /**
   * Computes the coefficients for the computation of the averaged solution over CV faces,
   * starting from solution values in the SV face flux points.
   * @pre createFluxPolyCoef, createFluxPolyNodeCoord(),
   * createFaceFluxPntsConn() and computeInvVolFractionsOfFaceCVs()
   * and basic SV data must be built as well.
   */
  virtual void createCVExtFaceFluxCoef();

  /**
   * Computes the coefficients for the computation of the averaged solution over SV faces,
   * starting from solution values in the SV face flux points.
   * @pre createAvgSolInCVFaceCoef() and computeFaceFractionsOfCVs()
   */
  virtual void createAvgSolInSVFaceCoef();

  /**
   * sets the convective/diffusive cfl ratio
   */
  virtual void setCFLConvDiffRatio() = 0;

  /**
   * Creates a set of nodes for interpolation with the requested polynomial degree
   */
  virtual void setInterpolationNodeSet(const CFPolyOrder::Type order,std::vector< RealVector >& nodalSet) = 0;

  /**
   * Calls all the functions above to recreate the entire datastructure.
   */
  void resetSpectralFVElementData();

  /**
   * A function that computes the surface of a triangle from the node positions.
   */
  inline CFreal computeTriangleSurface(const RealVector& node1Coord,
                                       const RealVector& node2Coord,
                                       const RealVector& node3Coord);

  /**
   * A function that computes the surface of a polygon from the node positions.
   */
  inline CFreal computePolygonSurface(const std::vector< RealVector >& nodesCoord);

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

protected: // protected data

  /// dimensionality
  CFDim m_dimensionality;

  /// shape
  CFGeoShape::Type m_shape;

  /// polynomial order
  CFPolyOrder::Type m_polyOrder;

  /// simplex integrator
  SimplexGaussIntegrator m_sIntegrator;

  /// tensor product integrator
  TensorProductGaussIntegrator m_tpIntegrator;

  /// local node coordinates
  std::vector< RealVector > m_localNodeCoord;

  /// local face - local node connectivity
  std::vector< std::vector< CFuint > > m_localFaceNodeConn;

  /// local CV - local node connectivity
  std::vector< std::vector< CFuint > > m_localCVNodeConn;

  /// local internal face - CV connectivity
  std::vector< std::vector< CFuint > > m_localIntFaceCVConn;

  /// local external face normals (whole SV face)
  std::vector< RealVector > m_extFaceLocalNorm;

  /// CV volume fraction
  std::vector< CFreal > m_volFracCV;

  /// inverse CV volume fraction
  std::vector< CFreal > m_invVolFracCV;

  /// SV face - local external face connectivity
  std::vector< std::vector< CFuint > > m_svFaceLocalExtFaceConn;

  /// local external face - CV connectivity
  std::vector< CFuint > m_localExtFaceCVConn;

  /// local external SV face - CV connectivity
  std::vector< std::vector< CFuint > > m_extSVFaceCVConn;

  /// external face - CV connectivity per face connection orientation
  std::vector< std::vector< std::vector< CFuint > > > m_extSVFaceCVConnPerOrient;

  /// external face - local external SV face
  std::vector< std::vector< CFuint > > m_eSVFaceLocESVFaceConnPerOrient;

  /// SV node coordinates
  std::vector< RealVector > m_svNodeCoords;

  /// SV face - node coordinates
  std::vector< std::vector< RealVector > > m_svFaceNodeCoords;

  /// SV face - node connectivity
  std::vector< std::vector< CFuint > > m_svFaceNodeConn;

  /// SV face - node connectivity per orientation
  std::vector< std::vector< std::vector< CFuint > > > m_svFaceNodeConnPerOrientation;

  /// SV face - node connectivity per orientation, not taking into account symmetries
  std::vector< std::vector< std::vector< CFuint > > > m_svFaceNodeConnPerOrientationNoSymm;

  /// SV face - node coordinates per orientation
  std::vector< std::vector< std::vector< RealVector > > > m_svFaceNodeCoordsPerOrient;

  /// SV face - node coordinates per orientation
  std::vector< std::vector< std::vector< RealVector > > > m_svFaceNodeCoordsPerOrientNoSymm;

  /// ext face node coordinates local to the SV face
  std::vector< std::vector< RealVector > > m_extFaceNodeLocalCoords;

  /// CV face fraction
  std::vector< std::vector< CFreal > > m_faceFracCV;

  /// face-CVs inverse volume fraction
  std::vector< std::vector< CFreal > > m_invVolFracFaceCV;

  /// face-CVs inverse volume fraction per orientation
  std::vector< std::vector< std::vector< CFreal > > > m_invVolFracFaceCVPerOrient;

  /// polynomial exponents
  std::vector< std::vector< CFint > > m_polyExponents;

  /// spectral finite volume polynomial coefficients
  std::vector< std::vector< CFreal > > m_polyCoefSFV;

  /// local coordinates of points used for initialization of the solution
  std::vector< RealVector > m_initPntsCoords;

  /// transformation matrix for initialization of the states
  RealMatrix m_initTransfMatrix;

  /// local internal face normals in each quadrature point
  std::vector< std::vector< RealVector > > m_intFaceQuadPntNorm;

  /// local internal face quadrature point coordinates
  std::vector< std::vector< RealVector > > m_intFaceQuadPntCoords;

  /// local internal face quadrature wheights
  std::vector< std::vector< CFreal > > m_intFaceQuadWheights;

  /// local internal face quadrature point polynomial values
  std::vector< std::vector< std::vector< CFreal > > > m_intFaceQuadPntPolyVals;

  /// local internal quadrature point coordinates (ordered in one list)
  std::vector< RealVector > m_intQuadPntCoords;

  /// local internal quadrature wheights (ordered in one list)
  std::vector< CFreal > m_intQuadWheights;

  /// local internal quadrature point polynomial values (ordered in one list)
  std::vector< std::vector< CFreal > > m_intQuadPntPolyVals;

  /// wheight coordinates of face quadrature points
  std::vector< RealVector >  m_extQPntWheightCoords;

  /// local external quadrature point coordinates, per SV face
  std::vector< std::vector< RealVector > > m_extQPntCoords;

  /// local external quadrature wheights, per SV face
  std::vector< std::vector< CFreal > > m_extQWheights;

  /// local external quadrature quadrature point polynomial values, per SV face
  std::vector< std::vector< std::vector< CFreal > > > m_extQPntPolyVals;

  /// wheight coordinates of face quadrature points
  std::vector< std::vector< RealVector > >  m_extFaceQPntWheightCoords;

  /// local external face quadrature point coordinates, per SV face
  std::vector< std::vector< std::vector< RealVector > > > m_extFaceQuadPntCoords;

  /// local external face quadrature wheights, per SV face
  std::vector< std::vector< std::vector< CFreal > > > m_extFaceQuadWheights;

  /// local external face quadrature quadrature point polynomial values, per SV face
  std::vector< std::vector< std::vector< std::vector< CFreal > > > > m_extFaceQuadPntPolyVals;

  /// local external face quadrature point coordinates, per orientation of the face connectivity
  std::vector< std::vector< std::vector< std::vector< RealVector > > > > m_extFaceQPntCoordsPerOrient;

  /// local external face quadrature wheights, per orientation of the face connectivity
  std::vector< std::vector< std::vector< CFreal > > > m_extFaceQWheightsPerOrient;

  /// local external face quadrature quadrature point polynomial values,
  /// per orientation of the face connectivity
  std::vector< std::vector< std::vector< std::vector< std::vector< CFreal > > > > >
                                                                        m_extFaceQPntPolyValsPerOrient;

  /// local external face quadrature quadrature point polynomial values,
  /// per orientation of the face connectivity, not taking into account any possible symmetries
  std::vector< std::vector< std::vector< std::vector< std::vector< CFreal > > > > >
      m_extFaceQPntPolyValsPerOrientNoSymm;

  /// local external quadrature point coordinates, per orientation of the face connectivity
  std::vector< std::vector< std::vector< RealVector > > > m_extQPntCoordsPerOrient;

  /// local external quadrature wheights, per orientation of the face connectivity
  std::vector< std::vector< CFreal > > m_extQWheightsPerOrient;

  /// local external quadrature quadrature point polynomial values,
  /// per orientation of the face connectivity
  std::vector< std::vector< std::vector< std::vector< CFreal > > > > m_extQPntPolyValsPerOrient;

  /// local external quadrature quadrature point polynomial values,
  /// per orientation of the face connectivity,
  /// not taking into account any symmetries
  std::vector< std::vector< std::vector< std::vector< CFreal > > > > m_extQPntPolyValsPerOrientNoSymm;

  /// polynomial values averaged over SV face (for computation of wave speed)
  std::vector< std::vector< CFreal > > m_svFaceAvgPolyVals;

  /// polynomial values averaged over SV face, ordered per orientation of the face connectivity
  std::vector< std::vector< std::vector< CFreal > > > m_svFaceAvgPolyValsPerOrient;


  // Data for quadrature free implementation
  /// local coordinates of nodes for flux polynomial
  std::vector< RealVector > m_fluxPolyNodeCoord;

  /// face - flux points connectivity
  std::vector< std::vector< CFuint > > m_faceFlxPntsConn;

  /// wheight coordinates of flux polynomial node coordinates in a face
  std::vector< RealVector > m_faceFluxPolyNodeWheightCoord;

  /// flux polynomial exponents
  std::vector< std::vector< CFint > > m_fluxPolyExponents;

  /// flux polynomial coefficients
  std::vector< std::vector< CFreal > > m_fluxPolyCoefs;

  /// coefficients for solution reconstruction in flux points
  std::vector< std::vector< CFreal > > m_solInFluxPntCoef;

  /// tensor for volume terms of quadrature free approach
  std::vector< std::vector< std::vector< CFreal > > > m_volTermTensor;

  /// coefficients for solution reconstruction in flux points on SV face
  std::vector< std::vector< std::vector< CFreal > > > m_solInFaceFluxPntCoef;

  /// coefficients for solution reconstruction in flux points on SV face,
  /// arranged per face connectivity orientation
  std::vector< std::vector< std::vector< std::vector< CFreal > > > > m_solInFaceFluxPntCoefPerOrient;

  /// coefficients for solution reconstruction in flux points on SV face,
  /// arranged per face connectivity orientation, not taking into account any symmetries
  std::vector< std::vector< std::vector< std::vector< CFreal > > > > m_solInFaceFluxPntCoefPerOrientNoSymm;

  /// coefficients for computation of integralover CV face on SV face, starting from SV face flux points
  std::vector< std::vector< CFreal > > m_cvExtFaceFluxCoef;

  /// coefficients for computation of averaged solution over SV face, starting from SV face flux points
  /// Actually, there are more flux points than necessary to compute the solution to the desired order of
  /// accuracy, but the solution in the flux points is available, and it's cheaper to use these values than to
  /// use the CV averaged solutions.
  std::vector< CFreal > m_avgSolInSVFaceCoef;

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

}; // end of class SpectralFVElementData

//////////////////////////////////////////////////////////////////////////////

// Inline functions

//////////////////////////////////////////////////////////////////////////////

inline CFreal SpectralFVElementData::computeTriangleSurface(const RealVector& node1Coord,
                                                            const RealVector& node2Coord,
                                                            const RealVector& node3Coord)
{
  return 0.5*( (node2Coord[0]-node1Coord[0])*(node3Coord[1]-node1Coord[1]) -
               (node3Coord[0]-node1Coord[0])*(node2Coord[1]-node1Coord[1]) );
}

//////////////////////////////////////////////////////////////////////////////

inline CFreal SpectralFVElementData::computePolygonSurface(const std::vector< RealVector >& nodesCoord)
{
  CFreal surface = 0.0;

  const CFuint nbrTriangles = nodesCoord.size()-2;

  for (CFuint iTriangle = 0; iTriangle < nbrTriangles; ++iTriangle)
  {
    surface += computeTriangleSurface(nodesCoord[0],nodesCoord[iTriangle+1],nodesCoord[iTriangle+2]);
  }

  if (surface < 0) throw Framework::BadFormatException (FromHere(),"Polygon with negative surface. Check if the node-order is correct");

  return surface;
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace SpectralFV

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_SpectralFV_SpectralFVElementData_hh
