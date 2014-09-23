// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_ShapeFunctions_LagrangeShapeFunctionQuadP1_hh
#define COOLFluiD_ShapeFunctions_LagrangeShapeFunctionQuadP1_hh

//////////////////////////////////////////////////////////////////////////////

#include "MathTools/MathConsts.hh"
#include "ShapeFunctions/LagrangeShapeFunction.hh"
#include "ShapeFunctions/ShapeFunctionsAPI.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework { class Node; }

  namespace ShapeFunctions {

//////////////////////////////////////////////////////////////////////////////

/// This class provides the lagrangian shape function describing the
/// representation of the solution and/or the geometry in a quadrilateral
/// P1 element.
/// @author Andrea Lani
/// @author Tiago Quintino
class ShapeFunctions_API LagrangeShapeFunctionQuadP1 : public ShapeFunctions::LagrangeShapeFunction {
public:

  /// Get the inheritant dimensionality of the ShapeFunction
  static CFuint getDimensionality()
  {
    return DIM_2D;
  }

  /// Get the number of nodal shape functions in the element
  static CFuint getNbNodes()
  {
    return 4;
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
    return CFPolyOrder::ORDER1;
  }

  static void getStatesMappedCoordinates(std::vector<RealVector>& mappedCoords)
  {

    throw Common::ShouldNotBeHereException (FromHere(),getName() + "::getStatesMappedCoordinates()");

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

      shapeFunc[0] = 0.25 * (1.0 - xi) * (1.0 - eta);
      shapeFunc[1] = 0.25 * (1.0 + xi) * (1.0 - eta);
      shapeFunc[2] = 0.25 * (1.0 + xi) * (1.0 + eta);
      shapeFunc[3] = 0.25 * (1.0 - xi) * (1.0 + eta);
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

    cf_assert(nodes.size() == getNbNodes());

    const CFreal x0 = (*nodes[0])[XX];
    const CFreal y0 = (*nodes[0])[YY];

    const CFreal x1 = (*nodes[1])[XX];
    const CFreal y1 = (*nodes[1])[YY];

    const CFreal x2 = (*nodes[2])[XX];
    const CFreal y2 = (*nodes[2])[YY];

    const CFreal x3 = (*nodes[3])[XX];
    const CFreal y3 = (*nodes[3])[YY];

    const CFreal A0 = 0.125 * ((y3-y1)*(x2-x0)-(y2-y0)*(x3-x1));
    const CFreal A1 = 0.125 * ((y2-y3)*(x1-x0)-(y1-y0)*(x2-x3));
    const CFreal A2 = 0.125 * ((y3-y0)*(x2-x1)-(y2-y1)*(x3-x0));

    const CFuint nbQdPts = mappedCoord.size();

    for (CFuint ip = 0; ip < nbQdPts; ++ip) {
      cf_assert(mappedCoord[ip].size() == DIM_2D);

      const CFreal xi   = mappedCoord[ip][KSI];
      const CFreal eta  = mappedCoord[ip][ETA];

      detJacobian[ip] = A0 + A1 * xi + A2 * eta;
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
//    throw Common::NotImplementedException (FromHere(),getName() + "::computeJacobianDeterminantPlus1D()");

    static std::vector<RealMatrix> jacob;
    const CFuint nbQdPts = mappedCoord.size();
    jacob.resize(nbQdPts);
    for(CFuint i = 0; i < jacob.size(); ++i) {
      jacob[i].resize(DIM_3D,DIM_3D);
    }
    computeJacobianPlus1D(nodes,mappedCoord,jacob);
    for(CFuint i = 0; i < jacob.size(); ++i) {
      detJacobian = jacob[i].determ3();
    }
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
    return "LagrangeQuadP1";
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
    using namespace std;
    using namespace COOLFluiD::Framework;

    const CFreal diagonalsProd = ((*nodes[2])[0] - (*nodes[0])[0])*
      ((*nodes[3])[1] - (*nodes[1])[1])-
      ((*nodes[2])[1] - (*nodes[0])[1])*
      ((*nodes[3])[0] - (*nodes[1])[0]);
    return 0.5*std::abs(diagonalsProd);
  }

  /// Get the centroid
  /// If Ai are the nodes, and A0 arbitrary node, then
  /// vect(A0G) = 1/n sum(vect(A0Ai))
  static RealVector computeCentroid(const std::vector<Framework::Node*>& nodes)
  {
    using namespace std;
    using namespace COOLFluiD::Framework;

    CFreal nbNodes = getNbNodes();
    RealVector coordCentroid(PhysicalModelStack::getActive()->getDim());
    coordCentroid = 0.;
    for (CFuint i=1; i < nbNodes; ++i){
      coordCentroid += (*nodes[i])-(*nodes[0]);
    }
    coordCentroid /= nbNodes;
    coordCentroid += (*nodes[0]);

    return coordCentroid;
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
    /// Linear element -> same as computeAvgFaceNormals for all mappedCoord
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
    /// Linear element -> same as computeAvgCellNormal for all mappedCoord
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
    using namespace MathTools;
    cf_assert(nodes.size() == 4);
    cf_assert(nodes[0]->size() == DIM_2D);
    cf_assert(nodes[1]->size() == DIM_2D);
    cf_assert(nodes[2]->size() == DIM_2D);
    cf_assert(nodes[3]->size() == DIM_2D);
    CFuint scalarProduct = 0;

    // Face 01
    computeFaceLineNormal(0,nodes,0,1);
    _tempVector2D = (*nodes[0]) - coord;

    //scalar product of AP with outward pointing normal
    CFreal inside = 0.;
    for(CFuint iDim=0; iDim < DIM_2D; ++iDim){
      inside += _normals[0][iDim]*_tempVector2D[iDim];
    }
    if(inside > -MathTools::MathConsts::CFrealEps()) scalarProduct++;

    // Face 12
    computeFaceLineNormal(1,nodes,1,2);
    _tempVector2D = (*nodes[1]) - coord;

    //scalar product of AP with outward pointing normal
    inside =0.;
    for(CFuint iDim=0; iDim < DIM_2D; ++iDim){
      inside += _normals[1][iDim]*_tempVector2D[iDim];
    }
    if(inside > -MathTools::MathConsts::CFrealEps()) scalarProduct++;

    // Face 23
    computeFaceLineNormal(2,nodes,2,3);
    _tempVector2D = (*nodes[2]) - coord;

    //scalar product of AP with outward pointing normal
    inside = 0.;
    for(CFuint iDim=0; iDim < DIM_2D; ++iDim){
      inside += _normals[2][iDim]*_tempVector2D[iDim];
    }
    if(inside > -MathTools::MathConsts::CFrealEps()) scalarProduct++;


    // Face 30
    computeFaceLineNormal(3,nodes,3,0);
    _tempVector2D = (*nodes[3]) - coord;

    //scalar product of AP with outward pointing normal
    inside = 0.;
    for(CFuint iDim=0; iDim < DIM_2D; ++iDim){
      inside += _normals[3][iDim]*_tempVector2D[iDim];
    }
    if(inside > -MathTools::MathConsts::CFrealEps()) scalarProduct++;


    if(scalarProduct == 4) return true;

    return false;
  }

private:

  /// Default constructor without arguments
  LagrangeShapeFunctionQuadP1()
  {
    for(CFuint i=0;i<_normals.size();i++)
    {
      _normals[i].resize(DIM_3D);
    }
  }

  /// Default destructor
  ~LagrangeShapeFunctionQuadP1();

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
    const CFuint nbNodes = nodes.size();
    const CFreal xi  = mappedCoord[0];
    const CFreal eta = mappedCoord[1];

    _shapeFunc[0] = 0.25 * (1.0 - xi) * (1.0 - eta);
    _shapeFunc[1] = 0.25 * (1.0 + xi) * (1.0 - eta);
    _shapeFunc[2] = 0.25 * (1.0 + xi) * (1.0 + eta);
    _shapeFunc[3] = 0.25 * (1.0 - xi) * (1.0 + eta);

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

}; // end of class LagrangeShapeFunctionQuadP1

//////////////////////////////////////////////////////////////////////////////

    } // namespace ShapeFunctions



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_ShapeFunctions_LagrangeShapeFunctionQuadP1_hh
