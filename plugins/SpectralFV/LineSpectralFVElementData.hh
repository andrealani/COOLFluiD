#ifndef COOLFluiD_SpectralFV_LineSpectralFVElementData_hh
#define COOLFluiD_SpectralFV_LineSpectralFVElementData_hh

//////////////////////////////////////////////////////////////////////////////



#include "Common/COOLFluiD.hh"
#include "Framework/CFPolyOrder.hh"
#include "Framework/CFGeoShape.hh"
#include "Common/SafePtr.hh"
#include "MathTools/RealVector.hh"
#include "SpectralFV/SpectralFVElementData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace SpectralFV {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a linear spectral finite volume element
 *
 * @author Kris Van den Abeele
 *
 */
class LineSpectralFVElementData : public SpectralFVElementData {
public:

  /**
   * Default constructor without arguments.
   */
  LineSpectralFVElementData();

  /**
   * Constructor initializing polyOrder.
   */
  LineSpectralFVElementData(CFPolyOrder::Type polyOrder);

  /**
   * Default destructor.
   */
  ~LineSpectralFVElementData();

protected: // functions
  /**
   * vector with the local coordinates of the SV nodes
   */
  virtual void computeSVNodeLocalCoords();

  /**
   * Creates the vector containing the local node coordinates of a SV.
   */
  virtual void createLocalNodeCoord();

  /**
   * Creates the connectivity between a local face and its local nodes in a SV.
   * External faces are put first.
   */
  virtual void createLocalFaceNodeConn();

  /**
   * Creates local CV - node connectivity.
   */
  virtual void createLocalCVNodeConn();

  /**
   * Creates the connectivity between a local internal face and its local neighbouring CVs.
   */
  virtual void createLocalIntFaceCVConn();

  /**
   * Computes the local internal face normals.
   * @pre createLocalNodeCoord() and createLocalFaceNodeConn() have to be executed first.
   */
  virtual void computeLocalFaceNormals();

  /**
   * Computes the external face local normals
   */
  virtual void computeExtFaceLocalNormals();

  /**
   * Computes the volume fraction that each CV occupates inside the SV.
   * @pre createLocalCVNodeConn() and createLocalNodeCoord() must be called first.
   */
  virtual void computeVolumeFractionsOfCVs();

  /**
   * Creates the connectivity between a local external face and its local neighbouring CV.
   */
  virtual void createLocalExtFaceCVConn();

  /**
   * Creates the connectivity between SV faces and local external faces.
   */
  virtual void createSVFaceLocalExtFaceConn();

  /**
   * Creates SV face-node connectivity
   */
  virtual void createSVFaceNodeConnectivity();

  /**
   * Creates a list with the different possible connectivities of faces
   * @pre createSVFaceNodeConnectivity()
   */
  virtual void createSVFaceNodeConnectivityPerOrient();

  /**
   * Creates a list with the different possible connectivities of faces,
   * not taking into account possible symmetries
   * (--> more possible orientations than with function above)
   * @pre createSVFaceNodeConnectivity()
   */
  virtual void createSVFaceNodeConnectivityPerOrientNoSymm();

  /**
   * Computes the fraction of the face that each CV occupates at a SV boundary.
   * @pre createLocalFaceNodeConn() and createLocalNodeCoord() must be called first.
   */
  virtual void computeFaceFractionsOfCVs();

  /**
   * Computes the the external face node coordinates local to the SV face
   */
  virtual void createExtFaceNodeLocalCoords();

  /**
   * Creates a vector containing the exponents of the terms in the Spectral FV polynomials
   */
  virtual void createPolyExponents();

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
   * Creates the wheight coordinates of the flux points in a SV face
   * @pre createFluxPolyNodeCoord() and createFaceFluxPntsConn()
   */
  virtual void createFaceFluxPolyNodeWheightCoord() {};

  /**
   * Creates a vector containing the exponents of the terms in the flux basis polynomials.
   */
  virtual void createFluxPolyExponents();

  /**
   * sets the convective/diffusive cfl ratio
   */
  virtual void setCFLConvDiffRatio();

  /**
   * Creates a set of nodes for interpolation with the requested polynomial degree
   */
  virtual void setInterpolationNodeSet(const CFPolyOrder::Type order,std::vector< RealVector >& nodalSet);

  /**
   * create the cell mapped coordinates of a uniform distribution of points on the cell faces (for output)
   */
  void createFaceOutputPntCellMappedCoords();

  /**
   * create the connectivity in a uniform distribution of points on the cell faces (for output)
   */
  void createFaceOutputPntConn();

}; // end of class LineSpectralFVElementData

//////////////////////////////////////////////////////////////////////////////

  } // namespace SpectralFV

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_SpectralFV_LineSpectralFVElementData_hh
