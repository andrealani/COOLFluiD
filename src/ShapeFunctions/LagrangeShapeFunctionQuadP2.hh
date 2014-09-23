// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_ShapeFunctions_LagrangeShapeFunctionQuadP2_hh
#define COOLFluiD_ShapeFunctions_LagrangeShapeFunctionQuadP2_hh

//////////////////////////////////////////////////////////////////////////////

#include "ShapeFunctions/LagrangeShapeFunction.hh"
#include "ShapeFunctions/ShapeFunctionsAPI.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {
    class Node;
  }



  namespace ShapeFunctions {

//////////////////////////////////////////////////////////////////////////////

/// This class provides the lagrangian shape function describing the
/// representation of the solution and/or the geometry in a quadrilateral
/// P2 element.
/// @author Thomas Wuilbaut
class ShapeFunctions_API LagrangeShapeFunctionQuadP2 : public ShapeFunctions::LagrangeShapeFunction {
public:

  /// Get the inheritant dimensionality of the ShapeFunction
  static CFuint getDimensionality()
  {
    return DIM_2D;
  }

  /// Get the number of nodal shape functions in the element
  static CFuint getNbNodes()
  {
    return 9;
  }

  /// Get the number of faces
  static CFuint getNbFaces()
  {
    return 4;
  }

  /// Gets the type of CFGeoShape::Type
  static CFGeoShape::Type getShape()
  {
    return CFGeoShape::QUAD;
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

    cf_assert (mappedCoords.size() == 9);

    mappedCoords[0][KSI] = -1. ;
    mappedCoords[0][ETA] = -1. ;

    mappedCoords[1][KSI] =  1. ;
    mappedCoords[1][ETA] = -1. ;

    mappedCoords[2][KSI] =  1. ;
    mappedCoords[2][ETA] =  1. ;

    mappedCoords[3][KSI] = -1. ;
    mappedCoords[3][ETA] =  1. ;

    mappedCoords[4][KSI] =  0. ;
    mappedCoords[4][ETA] = -1. ;

    mappedCoords[5][KSI] =  1. ;
    mappedCoords[5][ETA] =  0. ;

    mappedCoords[6][KSI] =  0. ;
    mappedCoords[6][ETA] =  1. ;

    mappedCoords[7][KSI] = -1. ;
    mappedCoords[7][ETA] =  0. ;

    mappedCoords[8][KSI] =  0. ;
    mappedCoords[8][ETA] =  0. ;

  }

  /// Compute the shape functions corresponding to the given
  /// mapped coordinates
  static void computeShapeFunctions
  (const std::vector<RealVector>& mappedCoord,
  std::vector<RealVector>& shapeFunc)
  {
    for (CFuint ip = 0; ip < mappedCoord.size(); ++ip) {
      computeShapeFunction(mappedCoord[ip],shapeFunc[ip]);
    }
  }


  /// Compute the shape functions corresponding to the given
  /// mapped coordinates
  static void computeShapeFunction
  (const RealVector& mappedCoord,
    RealVector& shapeFunc)
  {
      const CFreal xi  = mappedCoord[0];
      const CFreal eta = mappedCoord[1];
      const CFreal xi2 = xi*xi;
      const CFreal eta2 = eta*eta;
      const CFreal xiEta = xi*eta;

      shapeFunc[0] = 0.25 * (1.0 - xi)  * (1.0 - eta) * xiEta;
      shapeFunc[1] = -0.25 * (1.0 + xi)  * (1.0 - eta) * xiEta;
      shapeFunc[2] = 0.25 * (1.0 + xi)  * (1.0 + eta) * xiEta;
      shapeFunc[3] = -0.25 * (1.0 - xi)  * (1.0 + eta) * xiEta;
      shapeFunc[4] = -0.5  * (1.0 - xi2) * (1.0 - eta) * eta;
      shapeFunc[5] = 0.5  * (1.0 + xi)  * (1.0 - eta2) * xi;
      shapeFunc[6] = 0.5  * (1.0 - xi2) * (1.0 + eta) * eta;
      shapeFunc[7] = -0.5  * (1.0 - xi)  * (1.0 - eta2) * xi;
      shapeFunc[8] =        (1.0 - xi2) * (1.0 - eta2);
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
                                                      std::vector<RealVector>& normal);

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
               std::vector<RealMatrix>& jacob);

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

    std::vector< RealMatrix > pntJacob(nbrPnts,RealMatrix(2,2));
    computeJacobian(nodes,mappedCoord,pntJacob);
    for (CFuint ip = 0; ip < nbrPnts; ++ip)
    {
      detJacobian[ip] = pntJacob[ip].determ2();
    }
  }

  /// Compute the jacobian determinant at the given
  /// mapped coordinates
  /// @note Implemented in very inefficient way. Constant allocation and deallocation of memory.
  static void computeJacobianDeterminantPlus1D(
         const std::vector<RealVector>& mappedCoord,
         const std::vector<Framework::Node*>& nodes,
               std::valarray<CFreal>& detJacobian)
  {
    throw Common::NotImplementedException (FromHere(),getName() + "::computeJacobianDeterminantPlus1D()");
  }

  /// Compute the jacobian determinant at the given
  /// mapped coordinatesnodes
  static void computeJacobianDeterminantPlus2D(
         const std::vector<RealVector>& mappedCoord,
         const std::vector<Framework::Node*>& nodes,
               std::valarray<CFreal>& detJacobian)
  {
    throw Common::ShouldNotBeHereException (FromHere(),getName() + "::computeJacobianDeterminantPlus2D()");
  }

  /// Compute the determinant of the face jacobian at the quadrature points
  static void computeFaceJacobianDeterminant(
          const std::vector<RealVector>& mappedCoord,
          const std::vector<Framework::Node*>& nodes,
          const Framework::IntegratorPattern& pattern,
                std::vector<RealVector>& faceJacobian);

  /// Get the name of this shape function
  static const std::string getName()
  {
    return "LagrangeQuadP2";
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
    const CFreal x0 = (*nodes[0])[XX];
    const CFreal x1 = (*nodes[1])[XX];
    const CFreal x2 = (*nodes[2])[XX];
    const CFreal x3 = (*nodes[3])[XX];
    const CFreal x4 = (*nodes[4])[XX];
    const CFreal x5 = (*nodes[5])[XX];
    const CFreal x6 = (*nodes[6])[XX];
    const CFreal x7 = (*nodes[7])[XX];
    const CFreal y0 = (*nodes[0])[YY];
    const CFreal y1 = (*nodes[1])[YY];
    const CFreal y2 = (*nodes[2])[YY];
    const CFreal y3 = (*nodes[3])[YY];
    const CFreal y4 = (*nodes[4])[YY];
    const CFreal y5 = (*nodes[5])[YY];
    const CFreal y6 = (*nodes[6])[YY];
    const CFreal y7 = (*nodes[7])[YY];

    return (4*((x7 - x4)*y0 + (x4 - x5)*y1 + (x5 - x6)*y2 + (x6 - x7)*y3) +
            x1*(y0 - y2 - 4*y4 + 4*y5) + x2*(y1 - y3 - 4*y5 + 4*y6) +
            x0*(y3 - y1 + 4*y4 - 4*y7) + x3*(y2 - y0 - 4*y6 + 4*y7)
           )/6.;
  }

  /// Get the centroid
  static RealVector computeCentroid(const std::vector<Framework::Node*>& nodes)
  {
    return *nodes[8];
  }

  /// Compute Mapped Coordinates
  /// @param coord contains the coordinates to be mapped
  /// @param nodes contains the nodes
  /// @return RealVector containing the Mapped Coordinates
  static RealVector computeMappedCoordinates(const RealVector& coord, const std::vector<Framework::Node*>& nodes);

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
    cf_assert(nodes.size() == 4);

    // Face 01
    computeFaceLineNormal(0,nodes,0,1);

    // Face 12
    computeFaceLineNormal(1,nodes,1,2);

    // Face 23
    computeFaceLineNormal(2,nodes,2,3);

    // Face 30
    computeFaceLineNormal(3,nodes,3,0);

    return _normals;
  }

  /// Compute Face Normals
  /// @param mappedCoord contains the coordinates of the location at which the face normal should be computed
  /// @param nodes contains the nodes
  /// @return vector<RealVector> containing the average face normals
  static std::vector<RealVector> computeFaceNormals(const RealVector mappedCoord, const std::vector<Framework::Node*>& nodes)
  {
    throw Common::NotImplementedException
      (FromHere(), getName()+"::computeFaceNormals()");
  }

  /// Compute Cell Average Normal
  /// @param nodes contains the nodes
  /// @return RealVector containing the average cell normal
  static RealVector computeAvgCellNormal(const std::vector<Framework::Node*>& nodes)
  {
    _3Dvec1[XX] = (*nodes[2])[XX] - (*nodes[0])[XX];
    _3Dvec1[YY] = (*nodes[2])[YY] - (*nodes[0])[YY];
    _3Dvec1[ZZ] = (*nodes[2])[ZZ] - (*nodes[0])[ZZ];

    _3Dvec2[XX] = (*nodes[3])[XX] - (*nodes[1])[XX];
    _3Dvec2[YY] = (*nodes[3])[YY] - (*nodes[1])[YY];
    _3Dvec2[ZZ] = (*nodes[3])[ZZ] - (*nodes[1])[ZZ];

    MathTools::MathFunctions::crossProd(_3Dvec1,_3Dvec2,_tempVector3D);

    _tempVector3D *= 0.5;

    return _tempVector3D;
  }

  /// Compute Cell Normal at a given mappedCoord
  /// @param nodes contains the nodes
  /// @return RealVector containing the cell normal
  static RealVector computeCellNormal(const RealVector& mappedCoord, const std::vector<Framework::Node*>& nodes)
  {
    throw Common::NotImplementedException
      (FromHere(), getName()+"::computeCellNormal()");
  }

  /// Check if a point (defined with mapped coordinates) is inside an element
  static bool isInMappedElement(const RealVector& mappedCoord)
  {
    cf_assert(mappedCoord.size() == 2);
    if((mappedCoord[0] > 1.) ||
      (mappedCoord[1] > 1.) ||
      (mappedCoord[0] < -1.) ||
      (mappedCoord[1] < -1.))
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

  /// Default constructor without arguments
  LagrangeShapeFunctionQuadP2()
  {
    for(CFuint i=0;i<_normals.size();i++)
    {
      _normals[i].resize(DIM_3D);
    }
  }

  /// Default destructor
  ~LagrangeShapeFunctionQuadP2();

  /// Computes the normal to the line.
  /// Note that the order of indexes is important
  static void computeFaceLineNormal(const CFuint& n,
                          const std::vector<Framework::Node*>& nodes,
                          const CFuint& i,
                          const CFuint& j)
  {
    _normals[n][XX] = (*nodes[j])[YY] - (*nodes[i])[YY];
    _normals[n][YY] = (*nodes[i])[XX] - (*nodes[j])[XX];
  }

  /// Computes the coordinates form the mapped coordinates
  static void computeCoordFromMappedCoord(const RealVector& mappedCoord, const std::vector<Framework::Node*>& nodes, RealVector& coord)
  {

    throw Common::NotImplementedException
      (FromHere(), getName()+"::computeCoordFromMappedCoord()");

    const CFuint nbNodes = nodes.size();
    const CFreal xi  = mappedCoord[0];
    const CFreal eta = mappedCoord[1];
    const CFreal xi2 = xi*xi;
    const CFreal eta2 = eta*eta;
    const CFreal xiEta = xi*eta;

    _shapeFunc[0] = 0.25 * (1.0 - xi)  * (1.0 - eta) * xiEta;
    _shapeFunc[1] = -0.25 * (1.0 + xi)  * (1.0 - eta) * xiEta;
    _shapeFunc[2] = 0.25 * (1.0 + xi)  * (1.0 + eta) * xiEta;
    _shapeFunc[3] = -0.25 * (1.0 - xi)  * (1.0 + eta) * xiEta;
    _shapeFunc[4] = -0.5  * (1.0 - xi2) * (1.0 - eta) * eta;
    _shapeFunc[5] = 0.5  * (1.0 + xi)  * (1.0 - eta2) * xi;
    _shapeFunc[6] = 0.5  * (1.0 - xi2) * (1.0 + eta) * eta;
    _shapeFunc[7] = -0.5  * (1.0 - xi)  * (1.0 - eta2) * xi;
    _shapeFunc[8] =        (1.0 - xi2) * (1.0 - eta2);

    for(CFuint i=0;i < DIM_2D; ++i)
    {
      coord[i] = _shapeFunc[0]*(*(nodes[0]))[i];
      for(CFuint j=1;j<nbNodes;++j)
      {
      coord[i] += _shapeFunc[j]*(*(nodes[j]))[i];
      }
    }

  }

private:

  /// solution integrator ID
  static CFuint _interpolatorID;

  /// Temporary Vector for the computation of normals
  static RealVector _2Dvec1;
  static RealVector _2Dvec2;
  static RealVector _3Dvec1;
  static RealVector _3Dvec2;

  /// Vector of normals
  static std::vector<RealVector> _normals;

  /// temporary vector for 2D
  static RealVector _shapeFunc;

  /// temporary vector for 2D
  static RealVector _tempVector2D;

  /// temporary vector for 3D
  static RealVector _tempVector3D;

  /// Vector for the mapped coordinates
  static RealVector _mappedCoord;

  /// Temporary Vector for the computation of mapped coord
  static RealVector _vec10;
  static RealVector _vec20;
  static RealVector _vecP0;

  /// vector for shape function derivatives
  static std::vector< RealVector > _shapeFuncDerivs;

}; // end of class LagrangeShapeFunctionQuadP2

//////////////////////////////////////////////////////////////////////////////

    } // namespace ShapeFunctions



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_ShapeFunctions_LagrangeShapeFunctionQuadP2_hh
