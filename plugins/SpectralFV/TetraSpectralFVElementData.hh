#ifndef COOLFluiD_SpectralFV_TetraSpectralFVElementData_hh
#define COOLFluiD_SpectralFV_TetraSpectralFVElementData_hh

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
 * This class represents a tetrahedral spectral finite volume element
 *
 * @author Kris Van den Abeele
 *
 */
class TetraSpectralFVElementData : public SpectralFVElementData {
public:

  /**
   * Default constructor without arguments.
   */
  TetraSpectralFVElementData();

  /**
   * Constructor initializing polyOrder.
   */
  TetraSpectralFVElementData(CFPolyOrder::Type polyOrder);

  /**
   * Default destructor.
   */
  ~TetraSpectralFVElementData();

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
   * Creates the external face CV connectivity per orientation
   * @pre createLocalExtFaceCVConn() and createSVFaceLocalExtFaceConn()
   * @warning should be overwritten for elements other than LINEs and TRIAGs
   */
  virtual void createExtFaceCVConnPerOrient();

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
   * Computes the the external face node coordinates local to the SV face
   */
  virtual void createExtFaceNodeLocalCoords();

  /**
   * Computes the fraction of the face that each CV occupates at a SV boundary.
   * @pre createLocalFaceNodeConn() and createLocalNodeCoord() must be called first.
   */
  virtual void computeFaceFractionsOfCVs();

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
   * Computes coordinates, wheights and polynomial values in the quadrature points on internal faces
   */
  virtual void computeIntFaceQuadPntsData();

  /**
   * Computes the tensor for the evaluation of the volume terms
   * @pre createFluxPolyCoef must be called first, and basic SV data must be built as well.
   */
  virtual void createVolumeTermsTensor();

  /**
   * Creates the wheight coordinates of the flux points in a SV face
   * @pre createFluxPolyNodeCoord() and createFaceFluxPntsConn()
   */
  virtual void createFaceFluxPolyNodeWheightCoord();

  /**
   * Creates a vector containing the exponents of the terms in the flux basis polynomials.
   */
  virtual void createFluxPolyExponents();

  /**
   * Creates a set of nodes for interpolation with the requested polynomial degree
   */
  virtual void setInterpolationNodeSet(const CFPolyOrder::Type order,std::vector< RealVector >& nodalSet);

  /**
   * sets the convective/diffusive cfl ratio
   */
  virtual void setCFLConvDiffRatio();

private: // functions
  /**
   * Computes the normal to a triangular face
   */
  inline RealVector computeTriangleNormal(const RealVector& node0,
                                          const RealVector& node1,
                                          const RealVector& node2);

  /**
   * Computes the normal to a polygonal face
   */
  inline RealVector computePolygonNormal(const std::vector< RealVector >& nodes);

  /**
   * Computes the volume of a tetrahedron
   */
  inline CFreal computeTetraVolume(const RealVector& node0,
                                   const RealVector& node1,
                                   const RealVector& node2,
                                   const RealVector& node3);

  /**
   * Computes the volume of a polyhedron
   */
  inline CFreal computePolyHedronVolume(const std::vector< RealVector >& nodes);

  /**
   * sets the decomposition of a polyhedron into tetrahedrons
   * (assuming a certain ordering of the polyhedron nodes!!!)
   */
  void setPolyHedronTetraDecomposition(const CFuint nbrNodes,std::vector< std::vector< CFuint > >& tetraDecomp);

  /**
   * sets the edge node connectivity for all edges in a certain polyhedron
   * (assuming a certain ordering of the polyhedron nodes!!!)
   */
  void setPolyHedronEdgeNodeConn(const CFuint nbrNodes, std::vector< std::vector< CFuint > >& edgeToNodes);

  /**
   * sets the decomposition of a polyhedron into boxes
   * (assuming a certain ordering of the polyhedron nodes!!!)
   */
  void setPolyHedronBoxDecomposition(const std::vector< RealVector > polyHedronNodes,
                                     std::vector< std::vector< RealVector > >& boxDecomp);

  /**
   * sets the decomposition of a polygon into quadrilaterals
   * (assuming a certain ordering of the polygon nodes!!!)
   */
  void setPolyGonQuadrilateralDecomposition(const std::vector< RealVector > polygonNodes,
                                            std::vector< std::vector< RealVector > >& quadDecomp);

  /**
   * create the cell mapped coordinates of a uniform distribution of points on the cell faces (for output)
   */
  void createFaceOutputPntCellMappedCoords();

  /**
   * create the connectivity in a uniform distribution of points on the cell faces (for output)
   */
  void createFaceOutputPntConn();

}; // end of class TetraSpectralFVElementData

//////////////////////////////////////////////////////////////////////////////

// Inline functions

//////////////////////////////////////////////////////////////////////////////

inline RealVector TetraSpectralFVElementData::computeTriangleNormal(const RealVector& node0,
                                                                    const RealVector& node1,
                                                                    const RealVector& node2)
{
  cf_assert(node0.size() == 3);
  cf_assert(node1.size() == 3);
  cf_assert(node2.size() == 3);

  // triangle vectors
  RealVector vec1 = node1 - node0;
  RealVector vec2 = node2 - node0;

  // return variable
  RealVector normal(0.0,m_dimensionality);
  normal[XX] = 0.5*(vec1[YY]*vec2[ZZ] - vec1[ZZ]*vec2[YY]);
  normal[YY] = 0.5*(vec1[ZZ]*vec2[XX] - vec1[XX]*vec2[ZZ]);
  normal[ZZ] = 0.5*(vec1[XX]*vec2[YY] - vec1[YY]*vec2[XX]);

  return normal;
}

//////////////////////////////////////////////////////////////////////////////

inline RealVector TetraSpectralFVElementData::computePolygonNormal(const std::vector< RealVector >& nodes)
{
  // return variable
  RealVector normal(0.0,m_dimensionality);

  // number of triangles in polygon
  const CFuint nbrTriangles = nodes.size()-2;

  // compute polygon normal
  for (CFuint iTriangle = 0; iTriangle < nbrTriangles; ++iTriangle)
  {
    normal += computeTriangleNormal(nodes[0],nodes[iTriangle+1],nodes[iTriangle+2]);
  }

  return normal;
}

//////////////////////////////////////////////////////////////////////////////

inline CFreal TetraSpectralFVElementData::computeTetraVolume(const RealVector& node0,
                                                             const RealVector& node1,
                                                             const RealVector& node2,
                                                             const RealVector& node3)
{
  return (+ (node1[XX]-node0[XX])*(node2[YY]-node0[YY])*(node3[ZZ]-node0[ZZ])
          + (node2[XX]-node0[XX])*(node3[YY]-node0[YY])*(node1[ZZ]-node0[ZZ])
          + (node3[XX]-node0[XX])*(node1[YY]-node0[YY])*(node2[ZZ]-node0[ZZ])
          - (node3[XX]-node0[XX])*(node2[YY]-node0[YY])*(node1[ZZ]-node0[ZZ])
          - (node2[XX]-node0[XX])*(node1[YY]-node0[YY])*(node3[ZZ]-node0[ZZ])
          - (node1[XX]-node0[XX])*(node3[YY]-node0[YY])*(node2[ZZ]-node0[ZZ]))/6.0;
}

//////////////////////////////////////////////////////////////////////////////

inline CFreal TetraSpectralFVElementData::computePolyHedronVolume(const std::vector< RealVector >& nodes)
{
  // return variable
  CFreal volume = 0.0;

  // number of nodes
  const CFuint nbrNodes = nodes.size();

  // set decomposition into tetrahedrons
  std::vector< std::vector< CFuint > > tetraDecomp;
  setPolyHedronTetraDecomposition(nbrNodes,tetraDecomp);

  // loop over tetrahedrons to compute volume
  const CFuint nbrTetras = tetraDecomp.size();
  for (CFuint iTetra = 0; iTetra < nbrTetras; ++iTetra)
  {
    // get node IDs
    const CFuint node0ID = tetraDecomp[iTetra][0];
    const CFuint node1ID = tetraDecomp[iTetra][1];
    const CFuint node2ID = tetraDecomp[iTetra][2];
    const CFuint node3ID = tetraDecomp[iTetra][3];

    // compute volume
    volume += computeTetraVolume(nodes[node0ID],nodes[node1ID],nodes[node2ID],nodes[node3ID]);
  }

  return volume;
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace SpectralFV

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_SpectralFV_LineSpectralFVElementData_hh
