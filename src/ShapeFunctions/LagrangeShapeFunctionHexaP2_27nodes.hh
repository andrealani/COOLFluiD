// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_ShapeFunctions_LagrangeShapeFunctionHexaP2_27Nodes_hh
#define COOLFluiD_ShapeFunctions_LagrangeShapeFunctionHexaP2_27Nodes_hh

//////////////////////////////////////////////////////////////////////////////

#include "Common/StringOps.hh"
#include "ShapeFunctions/LagrangeShapeFunction.hh"
#include "ShapeFunctions/ShapeFunctionsAPI.hh"
#include "Framework/VolumeCalculator.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace ShapeFunctions {

//////////////////////////////////////////////////////////////////////////////

/// This class provides the lagrangian shape function
/// of the second order for a hexahedron.
/// @author Tiago Quintino
class ShapeFunctions_API LagrangeShapeFunctionHexaP2_27Nodes : public ShapeFunctions::LagrangeShapeFunction {
public:

  /// Get the inheritant dimensionality of the ShapeFunction
  static CFuint getDimensionality()
  {
    return DIM_3D;
  }

  /// Get the number of nodal shape functions in the element
  static CFuint getNbNodes()
  {
    return 27;
  }

  /// Get the number of faces
  static CFuint getNbFaces()
  {
    return 6;
  }

  /// Gets the type of CFGeoShape::Type
  static CFGeoShape::Type getShape()
  {
    return CFGeoShape::HEXA;
  }

  /// Gets the type of Interpolator
  static CFPolyForm::Type getInterpolatorType()
  {
    return CFPolyForm::LAGRANGE;
  }

  /// Gets the Interpolator order
  static CFPolyOrder::Type getInterpolatorOrder()
  {
    return CFPolyOrder::ORDER2;
  }

  static void getStatesMappedCoordinates(std::vector<RealVector>& mappedCoords)
  {
    cf_assert (mappedCoords.size() == 27);

    mappedCoords[0][KSI] = -1. ;
    mappedCoords[0][ETA] = -1. ;
    mappedCoords[0][ZTA] = -1. ;

    mappedCoords[1][KSI] = 1. ;
    mappedCoords[1][ETA] = -1. ;
    mappedCoords[1][ZTA] = -1. ;

    mappedCoords[2][KSI] = 1. ;
    mappedCoords[2][ETA] = 1. ;
    mappedCoords[2][ZTA] = -1. ;

    mappedCoords[3][KSI] = -1. ;
    mappedCoords[3][ETA] = 1. ;
    mappedCoords[3][ZTA] = -1. ;

    mappedCoords[4][KSI] = -1. ;
    mappedCoords[4][ETA] = -1. ;
    mappedCoords[4][ZTA] = 1. ;

    mappedCoords[5][KSI] = 1. ;
    mappedCoords[5][ETA] = -1. ;
    mappedCoords[5][ZTA] = 1. ;

    mappedCoords[6][KSI] = 1. ;
    mappedCoords[6][ETA] = 1. ;
    mappedCoords[6][ZTA] = 1. ;

    mappedCoords[7][KSI] = -1. ;
    mappedCoords[7][ETA] = 1. ;
    mappedCoords[7][ZTA] = 1. ;

    mappedCoords[8][KSI] = 0. ;
    mappedCoords[8][ETA] = -1. ;
    mappedCoords[8][ZTA] = -1. ;

    mappedCoords[9][KSI] = 1. ;
    mappedCoords[9][ETA] = 0. ;
    mappedCoords[9][ZTA] = -1. ;

    mappedCoords[10][KSI] = 0. ;
    mappedCoords[10][ETA] = 1. ;
    mappedCoords[10][ZTA] = -1. ;

    mappedCoords[11][KSI] = -1. ;
    mappedCoords[11][ETA] = 0. ;
    mappedCoords[11][ZTA] = -1. ;

    mappedCoords[12][KSI] = 0. ;
    mappedCoords[12][ETA] = 0. ;
    mappedCoords[12][ZTA] = -1. ;

    mappedCoords[13][KSI] = -1. ;
    mappedCoords[13][ETA] = -1. ;
    mappedCoords[13][ZTA] = 0. ;

    mappedCoords[14][KSI] = 0. ;
    mappedCoords[14][ETA] = -1. ;
    mappedCoords[14][ZTA] = 0. ;

    mappedCoords[15][KSI] = 1. ;
    mappedCoords[15][ETA] = -1. ;
    mappedCoords[15][ZTA] = 0. ;

    mappedCoords[16][KSI] = 1. ;
    mappedCoords[16][ETA] = 0. ;
    mappedCoords[16][ZTA] = 0. ;

    mappedCoords[17][KSI] = 1. ;
    mappedCoords[17][ETA] = 1. ;
    mappedCoords[17][ZTA] = 0. ;

    mappedCoords[18][KSI] = 0. ;
    mappedCoords[18][ETA] = 1. ;
    mappedCoords[18][ZTA] = 0. ;

    mappedCoords[19][KSI] = -1. ;
    mappedCoords[19][ETA] = 1. ;
    mappedCoords[19][ZTA] = 0. ;

    mappedCoords[20][KSI] = -1. ;
    mappedCoords[20][ETA] = 0. ;
    mappedCoords[20][ZTA] = 0. ;

    mappedCoords[21][KSI] = 0. ;
    mappedCoords[21][ETA] = 0. ;
    mappedCoords[21][ZTA] = 0. ;

    mappedCoords[22][KSI] = 0. ;
    mappedCoords[22][ETA] = -1. ;
    mappedCoords[22][ZTA] = 1. ;

    mappedCoords[23][KSI] = 1. ;
    mappedCoords[23][ETA] = 0. ;
    mappedCoords[23][ZTA] = 1. ;

    mappedCoords[24][KSI] = 0. ;
    mappedCoords[24][ETA] = 1. ;
    mappedCoords[24][ZTA] = 1. ;

    mappedCoords[25][KSI] = -1. ;
    mappedCoords[25][ETA] = 0. ;
    mappedCoords[25][ZTA] = 1. ;

    mappedCoords[26][KSI] = 0. ;
    mappedCoords[26][ETA] = 0. ;
    mappedCoords[26][ZTA] = 1. ;
  }

  /// Compute the shape functions corresponding to the given
  /// mapped coordinates
  static void computeShapeFunctions(
          const std::vector<RealVector>& mappedCoord, std::vector<RealVector>& shapeFunc)
  {
    for (CFuint ip = 0; ip < mappedCoord.size(); ++ip) {
      computeShapeFunction(mappedCoord[ip],shapeFunc[ip]);
    }
  }

  /// Compute the shape functions corresponding to the given
  /// mapped coordinates
  static void computeShapeFunction(const RealVector& mappedCoord,
                                         RealVector& shapeFunc)
  {
    const CFreal xi   = mappedCoord[0];
    const CFreal eta  = mappedCoord[1];
    const CFreal zeta = mappedCoord[2];
    const CFreal xi2  = xi*xi;
    const CFreal eta2 = eta*eta;
    const CFreal zeta2= zeta*zeta;

    const CFreal Nxi_m1 = -xi * (1. - xi);
    const CFreal Nxi_0  =  2. * (1. - xi2);
    const CFreal Nxi_p1 =  xi * (1. + xi);

    const CFreal Neta_m1 = -eta * (1. - eta);
    const CFreal Neta_0  =  2.  * (1. - eta2);
    const CFreal Neta_p1 =  eta * (1. + eta);

    const CFreal Nzeta_m1 = -zeta * (1. - zeta);
    const CFreal Nzeta_0  =   2.  * (1. - zeta2);
    const CFreal Nzeta_p1 =  zeta * (1. + zeta);

    shapeFunc[ 0] = Nxi_m1 * Neta_m1 * Nzeta_m1 ;
    shapeFunc[ 1] = Nxi_p1 * Neta_m1 * Nzeta_m1 ;
    shapeFunc[ 2] = Nxi_p1 * Neta_p1 * Nzeta_m1 ;
    shapeFunc[ 3] = Nxi_m1 * Neta_p1 * Nzeta_m1 ;
    shapeFunc[ 4] = Nxi_m1 * Neta_m1 * Nzeta_p1 ;
    shapeFunc[ 5] = Nxi_p1 * Neta_m1 * Nzeta_p1 ;
    shapeFunc[ 6] = Nxi_p1 * Neta_p1 * Nzeta_p1 ;
    shapeFunc[ 7] = Nxi_m1 * Neta_p1 * Nzeta_p1 ;
    shapeFunc[ 8] = Nxi_0  * Neta_m1 * Nzeta_m1 ;
    shapeFunc[ 9] = Nxi_p1 * Neta_0  * Nzeta_m1 ;
    shapeFunc[10] = Nxi_0  * Neta_p1 * Nzeta_m1 ;
    shapeFunc[11] = Nxi_m1 * Neta_0  * Nzeta_m1 ;
    shapeFunc[12] = Nxi_0  * Neta_0  * Nzeta_m1 ;
    shapeFunc[13] = Nxi_m1 * Neta_m1 * Nzeta_0  ;
    shapeFunc[14] = Nxi_0  * Neta_m1 * Nzeta_0  ;
    shapeFunc[15] = Nxi_p1 * Neta_m1 * Nzeta_0  ;
    shapeFunc[16] = Nxi_p1 * Neta_0  * Nzeta_0  ;
    shapeFunc[17] = Nxi_p1 * Neta_p1 * Nzeta_0  ;
    shapeFunc[18] = Nxi_0  * Neta_p1 * Nzeta_0  ;
    shapeFunc[19] = Nxi_m1 * Neta_p1 * Nzeta_0  ;
    shapeFunc[20] = Nxi_m1 * Neta_0  * Nzeta_0  ;
    shapeFunc[21] = Nxi_0  * Neta_0  * Nzeta_0  ;
    shapeFunc[22] = Nxi_0  * Neta_m1 * Nzeta_p1 ;
    shapeFunc[23] = Nxi_p1 * Neta_0  * Nzeta_p1 ;
    shapeFunc[24] = Nxi_0  * Neta_p1 * Nzeta_p1 ;
    shapeFunc[25] = Nxi_m1 * Neta_0  * Nzeta_p1 ;
    shapeFunc[26] = Nxi_0  * Neta_0  * Nzeta_p1 ;

    shapeFunc *= 0.125;
  }

  /// Compute the Gradient of the Shape Function
  static void computeGradientStates(
         const std::vector<RealMatrix>& jacob,
         const std::vector<RealVector>& mappedCoord,
               std::vector<RealMatrix>& grad);

  /// Computes the normal to a face at the given mapped coordinates,
  /// scaled with the 'face Jacobian determinant'.
  /// (Normal has the dimensionality of the Face + 1)
  static void computeFaceJacobDetVectorAtMappedCoords(const std::vector<RealVector>& mappedCoord,
      const std::vector<Framework::Node*>& nodes,
      std::vector<RealVector>& normal)
  {
    throw Common::NotImplementedException (FromHere(),getName() + "::computeFaceJacobDetVectorAtMappedCoords()");
  }

  /// Computes the normal to a given mapped coordinate plane, at the given mapped coordinates
  static void computeMappedCoordPlaneNormal(const std::vector<CFuint>& planeIdx,
                                            const std::vector<RealVector>& mappedCoord,
                                            const std::vector<Framework::Node*>& nodes,
                                            std::vector<RealVector>& normal);

  /// Compute the Jacobian
  static void computeJacobian(
         const std::vector<Framework::Node*>& nodes,
         const std::vector<RealVector>& mappedCoord,
               std::vector<RealMatrix>& jacob);

  /// Compute the Jacobian
  static void computeJacobianPlus1D(
         const std::vector<Framework::Node*>& nodes,
         const std::vector<RealVector>& mappedCoord,
               std::vector<RealMatrix>& jacob)
  {
    throw Common::ShouldNotBeHereException (FromHere(),getName() + "::computeJacobianPlus1D()");
  }

  /// Compute the Jacobian
  static void computeJacobianPlus2D(
         const std::vector<Framework::Node*>& nodes,
         const std::vector<RealVector>& mappedCoord,
               std::vector<RealMatrix>& jacob)
  {
    throw Common::ShouldNotBeHereException (FromHere(),getName() + "::computeJacobianPlus2D()");
  }

  /// Compute the jacobian determinant at the given
  /// mapped coordinates
  static void computeJacobianDeterminant(
         const std::vector<RealVector>& mappedCoord,
         const std::vector<Framework::Node*>& nodes,
               std::valarray<CFreal>& detJacobian)
  {
    /// @note implemented in a not very efficient way
    cf_assert(nodes.size() == getNbNodes());

    const CFuint nbrPnts = mappedCoord.size();
    cf_assert(detJacobian.size() == nbrPnts);

    std::vector< RealMatrix > pntJacob(nbrPnts,RealMatrix(3,3));
    computeJacobian(nodes,mappedCoord,pntJacob);
    for (CFuint ip = 0; ip < nbrPnts; ++ip)
    {
      detJacobian[ip] = pntJacob[ip].determ3();
    }
  }

  /// Compute the jacobian determinant at the given
  /// mapped coordinates
  static void computeJacobianDeterminantPlus1D(
         const std::vector<RealVector>& mappedCoord,
         const std::vector<Framework::Node*>& nodes,
               std::valarray<CFreal>& detJacobian)
  {
    throw Common::NotImplementedException (FromHere(),getName() + "::computeJacobianDeterminantPlus1D()");
  }

  /// Compute the jacobian determinant at the given
  /// mapped coordinates
  static void computeJacobianDeterminantPlus2D(
         const std::vector<RealVector>& mappedCoord,
         const std::vector<Framework::Node*>& nodes,
               std::valarray<CFreal>& detJacobian)
  {
    throw Common::ShouldNotBeHereException (FromHere(),getName() + "::computeJacobianDeterminantPlus2D()");
  }

  /// Get the name of this shape function
  static const std::string getName()
  {
    static const std::string name = "LagrangeHexaP2_27nodes";
    return name;
  }

  /// Compute the determinant of the face jacobian at the quadrature points
  static void computeFaceJacobianDeterminant(
          const std::vector<RealVector>& mappedCoord,
          const std::vector<Framework::Node*>& nodes,
          const Framework::IntegratorPattern& pattern,
                std::vector<RealVector>& faceJacobian)
  {
    throw Common::NotImplementedException (FromHere(),getName() + "::computeFaceJacobianDeterminant()");
  }


  /// Get the ID for the solution integration
  static Framework::InterpolatorID getInterpolatorID()
  {
    return _interpolatorID;
  }

  /// Get the ID for the solution integration
  static void setInterpolatorID(const Framework::InterpolatorID& id)
  {
    _interpolatorID = id;
  }

  /// Get the volume
  static CFreal computeVolume(const std::vector<Framework::Node*>& nodes)
  {
    /// @todo this is only an approximate volume, of a P1 hexahedron
    static Framework::VolumeCalculator volumeCalc;
    std::vector<Framework::Node*> nodesP1(8);
    for (CFuint in = 0; in < 8; ++in)
    {
      nodesP1[in] = nodes[in];
    }
    return volumeCalc.calculateHexaVolume(nodesP1);
  }

  /// Get the centroid
  /// If Ai are the nodes, and A0 arbitrary node, then
  /// vect(A0G) = 1/n sum(vect(A0Ai))
  static RealVector computeCentroid(const std::vector<Framework::Node*>& nodes)
  {
    return *nodes[21];
  }

  /// Compute Mapped Coordinates
  /// @param coord contains the coordinates to be mapped
  /// @param nodes contains the nodes
  /// @return RealVector containing the Mapped Coordinates
  static RealVector computeMappedCoordinates(const RealVector& coord, const std::vector<Framework::Node*>& nodes)
  {
    throw Common::NotImplementedException
      (FromHere(), getName()+"::computeMappedCoordinates()");
  }

  /// Compute Mapped Coordinates
  /// @param coord contains the coordinates to be mapped
  /// @param nodes contains the nodes
  /// @return RealVector containing the Mapped Coordinates
  static RealVector computeMappedCoordinatesPlus1D(const RealVector& coord, const std::vector<Framework::Node*>& nodes)
  {
    throw Common::NotImplementedException
      (FromHere(), getName()+"::computeMappedCoordinatesPlus1D()");
  }

  /// Compute Mapped Coordinates
  /// @param coord contains the coordinates to be mapped
  /// @param nodes contains the nodes
  /// @return RealVector containing the Mapped Coordinates
  static RealVector computeMappedCoordinatesPlus2D(const RealVector& coord, const std::vector<Framework::Node*>& nodes)
  {
    throw Common::NotImplementedException
      (FromHere(), getName()+"::computeMappedCoordinatesPlus2D()");
  }

  /// Compute Face Average Normals
  /// @param nodes contains the nodes
  /// @return vector<RealVector> containing the average face normals
  static std::vector<RealVector> computeAvgFaceNormals(const std::vector<Framework::Node*>& nodes)
  {
    cf_assert(nodes.size() == 6);

    // Face 0321
    computeFaceQuadNormal(0,nodes,0,3,2,1);

    // Face 4567
    computeFaceQuadNormal(1,nodes,4,5,6,7);

    // Face 0154
    computeFaceQuadNormal(2,nodes,0,1,5,4);

    // Face 1265
    computeFaceQuadNormal(3,nodes,1,2,6,5);

    // Face 3762
    computeFaceQuadNormal(4,nodes,3,7,6,2);

    // Face 0473
    computeFaceQuadNormal(5,nodes,0,4,7,3);

    return _normals;
  }

  /// Compute Face Normals
  /// @param mappedCoord contains the coordinates of the location at which the face normal should be computed
  /// @param nodes contains the nodes
  /// @return vector<RealVector> containing the average face normals
  static std::vector<RealVector> computeFaceNormals(const RealVector mappedCoord, const std::vector<Framework::Node*>& nodes)
  {
    throw Common::NotImplementedException
      (FromHere(), getName()+"::computeMappedCoordinates()");
  }

  /// Compute Cell Average Normal
  /// @param nodes contains the nodes
  /// @return RealVector containing the average cell normal
  static RealVector computeAvgCellNormal(const std::vector<Framework::Node*>& nodes)
  {
    throw Common::ShouldNotBeHereException (FromHere(),getName()+"::computeAvgCellNormal()");
  }

  /// Compute Cell Normal at a given mapped coord
  /// @param nodes contains the nodes
  /// @return RealVector containing the cell normal
  static RealVector computeCellNormal(const RealVector& mappedCoord, const std::vector<Framework::Node*>& nodes)
  {
    throw Common::ShouldNotBeHereException (FromHere(),getName()+"::computeCellNormal()");
  }

  /// Check if a point (defined with mapped coordinates) is inside an element
  static bool isInMappedElement(const RealVector& mappedCoord)
  {
    cf_assert(mappedCoord.size() == 3);
    if((mappedCoord[0] > 1.) ||
      (mappedCoord[1] > 1.) ||
      (mappedCoord[2] > 1.) ||
      (mappedCoord[0] < -1.) ||
      (mappedCoord[1] < -1.) ||
      (mappedCoord[2] < -1.))
    {
      return false;
    }
    else
    {
      return true;
    }
  }

  /// Check if a point is inside an element
  static bool isInElement(const std::vector<Framework::Node*>& nodes, const RealVector& coord)
  {
    throw Common::NotImplementedException
      (FromHere(), getName()+"::isInElement()");
  }

private:

  /// Constructor without arguments
  /// @see LagrangeShapeFunction()
  LagrangeShapeFunctionHexaP2_27Nodes()
  {
    for(CFuint i=0;i<_normals.size();i++)
    {
      _normals[i].resize(DIM_3D);
    }
  }

  /// Default destructor
  ~LagrangeShapeFunctionHexaP2_27Nodes();

  /// Computes the normal to the quadrilateral face.
  /// Note that the order of indexes is important
  /// Cyclic permutations are allowed:
  /// Fijkl = Fjkli = Fklij = Flijk
  static void computeFaceQuadNormal(const CFuint& n,
                          const std::vector<Framework::Node*>& nodes,
                         const CFuint& i,
                         const CFuint& j,
                         const CFuint& k,
                         const CFuint& l)
  {
    _vec1[XX] = (*nodes[k])[XX] - (*nodes[i])[XX];
    _vec1[YY] = (*nodes[k])[YY] - (*nodes[i])[YY];
    _vec1[ZZ] = (*nodes[k])[ZZ] - (*nodes[i])[ZZ];

    _vec2[XX] = (*nodes[l])[XX] - (*nodes[j])[XX];
    _vec2[YY] = (*nodes[l])[YY] - (*nodes[j])[YY];
    _vec2[ZZ] = (*nodes[l])[ZZ] - (*nodes[j])[ZZ];

    MathTools::MathFunctions::crossProd(_vec1,_vec2,_vec3);

    _vec3 *= 0.5;

    _normals[n][XX] = _vec3[XX];
    _normals[n][YY] = _vec3[YY];
    _normals[n][ZZ] = _vec3[ZZ];
  }

private:

  /// solution integrator ID
  static CFuint _interpolatorID;

  /// Temporary matrix for the holding the 3D jacobian inverse
  static RealMatrix m_invJ;

  /// Temporary Vector for the computation of normals
  static RealVector _vec1;
  static RealVector _vec2;
  static RealVector _vec3;

  /// gradients of shape functions
  static std::vector< RealVector > _gradShapFunc;

  /// Vector of normals
  static std::vector<RealVector> _normals;

}; // end of class LagrangeShapeFunctionHexaP2_27Nodes

//////////////////////////////////////////////////////////////////////////////

    } // namespace ShapeFunctions

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

// #include "ShapeFunctions/LagrangeShapeFunctionHexaP2_27Nodes.ci"

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_ShapeFunctions_LagrangeShapeFunctionHexaP2_27Nodes_hh
