// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_ShapeFunctions_LagrangeShapeFunctionPyramP1_hh
#define COOLFluiD_ShapeFunctions_LagrangeShapeFunctionPyramP1_hh

//////////////////////////////////////////////////////////////////////////////

#include "ShapeFunctions/LagrangeShapeFunction.hh"
#include "ShapeFunctions/ShapeFunctionsAPI.hh"
#include "Framework/VolumeCalculator.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



  namespace ShapeFunctions {

    class Node;

//////////////////////////////////////////////////////////////////////////////

/// This class provides the lagrangian shape function describing the
/// representation of the solution and/or the geometry in a triangular
/// element.
/// This is an explicit specialization for order = CFPolyOrder::ORDER1
/// @author Tiago Quintino
class ShapeFunctions_API LagrangeShapeFunctionPyramP1 : public ShapeFunctions::LagrangeShapeFunction {
public:

  /// Get the inheritant dimensionality of the ShapeFunction
  static CFuint getDimensionality()
  {
    return DIM_3D;
  }

  /// Get the number of nodal shape functions in the element
  static CFuint getNbNodes()
  {
    return 5;
  }

  /// Get the number of faces
  static CFuint getNbFaces()
  {
    return 5;
  }

  /// Gets the type of CFGeoShape::Type
  static CFGeoShape::Type getShape()
  {
    return CFGeoShape::PYRAM;
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
  static void computeShapeFunctions(
          const std::vector<RealVector>& mappedCoord, std::vector<RealVector>& shapeFunc)
  {
    for (CFuint ip = 0; ip < mappedCoord.size(); ++ip) {
      computeShapeFunction(mappedCoord[ip],shapeFunc[ip]);
    }
  }

  /// Compute the shape functions corresponding to the given
  /// mapped coordinates
  static void computeShapeFunction(
          const RealVector& mappedCoord, RealVector& shapeFunc)
  {
    const CFreal xi   = mappedCoord[0];
    const CFreal eta  = mappedCoord[1];
    const CFreal zeta = mappedCoord[2];

    const CFreal xi_eta = xi*eta;

    shapeFunc[0] = 0.25 * (0.5 + xi - eta - xi_eta - 0.5*zeta);
    shapeFunc[1] = 0.25 * (0.5 + xi + eta + xi_eta - 0.5*zeta);
    shapeFunc[2] = 0.25 * (0.5 - xi + eta - xi_eta - 0.5*zeta);
    shapeFunc[3] = 0.25 * (0.5 - xi - eta + xi_eta - 0.5*zeta);
    shapeFunc[4] = 0.50 * (1.0 + zeta);
  }

  /// Compute the Gradient of the Shape Function
  static void computeGradientStates(
         const std::vector<RealMatrix>& jacob,
         const std::vector<RealVector>& mappedCoord,
               std::vector<RealMatrix>& grad)
  {
    throw Common::NotImplementedException (FromHere(),getName() + "::computeGradientStates()");
  }

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
                                            std::vector<RealVector>& normal)
  {
    throw Common::NotImplementedException (FromHere(),getName() + "::computeMappedCoordPlaneNormal()");
  }

  /// Compute the Jacobian
  static void computeJacobian(
         const std::vector<Framework::Node*>& nodes,
         const std::vector<RealVector>& mappedCoord,
               std::vector<RealMatrix>& jacob)
  {
    throw Common::NotImplementedException (FromHere(),getName() + "::computeJacobian()");
  }

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
    cf_assert(nodes.size() == getNbNodes());

    const CFreal x0 = (*nodes[0])[XX];
    const CFreal y0 = (*nodes[0])[YY];
    const CFreal z0 = (*nodes[0])[ZZ];

    const CFreal x1 = (*nodes[1])[XX];
    const CFreal y1 = (*nodes[1])[YY];
    const CFreal z1 = (*nodes[1])[ZZ];

    const CFreal x2 = (*nodes[2])[XX];
    const CFreal y2 = (*nodes[2])[YY];
    const CFreal z2 = (*nodes[2])[ZZ];

    const CFreal x3 = (*nodes[3])[XX];
    const CFreal y3 = (*nodes[3])[YY];
    const CFreal z3 = (*nodes[3])[ZZ];

    const CFreal x4 = (*nodes[4])[XX];
    const CFreal y4 = (*nodes[4])[YY];
    const CFreal z4 = (*nodes[4])[ZZ];

  const CFreal A0 =
  2*x2*y1*z0 + 2*x3*y1*z0 - 4*x4*y1*z0 - 2*x1*y2*z0 + 2*x3*y2*z0 - 2*x1*y3*z0 - 2*x2*y3*z0 + 4*x4*y3*z0 + 4*x1*y4*z0 - 4*x3*y4*z0 -
  2*x2*y0*z1 - 2*x3*y0*z1 + 4*x4*y0*z1 + 2*x0*y2*z1 + 2*x3*y2*z1 - 4*x4*y2*z1 + 2*x0*y3*z1 - 2*x2*y3*z1 - 4*x0*y4*z1 +
  4*x2*y4*z1 + 2*x1*y0*z2 - 2*x3*y0*z2 - 2*x0*y1*z2 - 2*x3*y1*z2 + 4*x4*y1*z2 + 2*x0*y3*z2 + 2*x1*y3*z2 - 4*x4*y3*z2 -
  4*x1*y4*z2 + 4*x3*y4*z2 + 2*x1*y0*z3 + 2*x2*y0*z3 - 4*x4*y0*z3 - 2*x0*y1*z3 + 2*x2*y1*z3 - 2*x0*y2*z3 - 2*x1*y2*z3 +
  4*x4*y2*z3 + 4*x0*y4*z3 - 4*x2*y4*z3 - 4*x1*y0*z4 + 4*x3*y0*z4 + 4*x0*y1*z4 - 4*x2*y1*z4 + 4*x1*y2*z4 - 4*x3*y2*z4 -
  4*x0*y3*z4 + 4*x2*y3*z4;

  const CFreal A2 =
      2*( x2*y1*z0 - x3*y1*z0 - x1*y2*z0 - x3*y2*z0 - x2*y0*z1 + x1*y3*z0 + x2*y3*z0 - x0*y3*z1 + x0*y2*z3 - x1*y2*z3
        + x2*y1*z3 + x3*y0*z1 + x0*y2*z1 + x3*y2*z1 - x2*y3*z1 + x1*y0*z2 + x3*y0*z2 - x0*y1*z2 - x3*y1*z2 - x0*y3*z2
        + x1*y3*z2 - x1*y0*z3 - x2*y0*z3 + x0*y1*z3 )
    + 4*( x3*y4*z0 - x4*y3*z0 - x2*y4*z0 + x4*y2*z0 - x4*y1*z3 - x0*y4*z3 + x1*y4*z3 + x2*y0*z4 - x3*y0*z4 - x2*y1*z4
        + x3*y1*z4 - x4*y0*z2 + x0*y4*z2 - x1*y4*z2 - x0*y2*z4 + x1*y2*z4 + x0*y3*z4 - x1*y3*z4 - x4*y2*z1 + x4*y3*z1
        + x2*y4*z1 - x3*y4*z1 + x4*y1*z2 + x4*y0*z3 );

  const CFreal A1 =
      x2*y1*z0 + x3*y1*z0 - 2*x4*y1*z0 - x1*y2*z0 - x3*y2*z0 + 2*x4*y2*z0 - x1*y3*z0 + x2*y3*z0 + 2*x1*y4*z0 - 2*x2*y4*z0 -
      x2*y0*z1 - x3*y0*z1 + 2*x4*y0*z1 + x0*y2*z1 - x3*y2*z1 + x0*y3*z1 + x2*y3*z1 - 2*x4*y3*z1 - 2*x0*y4*z1 + 2*x3*y4*z1 +
      x1*y0*z2 + x3*y0*z2 - 2*x4*y0*z2 - x0*y1*z2 + x3*y1*z2 - x0*y3*z2 - x1*y3*z2 + 2*x4*y3*z2 + 2*x0*y4*z2 - 2*x3*y4*z2 +
      x1*y0*z3 - x2*y0*z3 - x0*y1*z3 - x2*y1*z3 + 2*x4*y1*z3 + x0*y2*z3 + x1*y2*z3 - 2*x4*y2*z3 - 2*x1*y4*z3 + 2*x2*y4*z3 -
      2*x1*y0*z4 + 2*x2*y0*z4 + 2*x0*y1*z4 - 2*x3*y1*z4 - 2*x0*y2*z4 + 2*x3*y2*z4 + 2*x1*y3*z4 - 2*x2*y3*z4;

    const CFuint nbQdPts = mappedCoord.size();

    for (CFuint ip = 0; ip < nbQdPts; ++ip) {
      cf_assert(mappedCoord[ip].size() == DIM_3D);

      const CFreal xi   = mappedCoord[ip][KSI];
      const CFreal eta  = mappedCoord[ip][ETA];
      //  not used // const CFreal zeta = mappedCoord[ip][ZTA];

      detJacobian[ip] = 0.01562*(A0 + A1 * xi + A2 * eta);
    }
  }

  /// Compute the determinant of the face jacobian at the quadrature points
  static void computeFaceJacobianDeterminant(
          const std::vector<RealVector>& mappedCoord,
          const std::vector<Framework::Node*>& nodes,
          const Framework::IntegratorPattern& pattern,
                std::vector<RealVector>& faceJacobian);

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
    static const std::string name = "LagrangePyramP1";
    return name;
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
    static Framework::VolumeCalculator volumeCalc;
    return volumeCalc.calculatePyramVolume(nodes);
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
    cf_assert(nodes.size() == 5);

    // Face 0321
    computeFaceQuadNormal(0,nodes,0,3,2,1);

    // Face 014
    computeFaceTriagNormal(1,nodes,0,1,4);

    // Face 124
    computeFaceTriagNormal(2,nodes,1,2,4);

    // Face 234
    computeFaceTriagNormal(3,nodes,2,3,4);

    // Face 043
    computeFaceTriagNormal(4,nodes,0,4,3);

    return _normals;
  }

  /// Compute Face Normals
  /// @param mappedCoord contains the coordinates of the location at which the face normal should be computed
  /// @param nodes contains the nodes
  /// @return vector<RealVector> containing the average face normals
  static std::vector<RealVector> computeFaceNormals(const RealVector mappedCoord, const std::vector<Framework::Node*>& nodes)
  {
    /// Linear element -> same as computeAvgFaceNormals for all mappedCoord
    cf_assert(nodes.size() == 5);

    // Face 0321
    computeFaceQuadNormal(0,nodes,0,3,2,1);

    // Face 014
    computeFaceTriagNormal(1,nodes,0,1,4);

    // Face 124
    computeFaceTriagNormal(2,nodes,1,2,4);

    // Face 234
    computeFaceTriagNormal(3,nodes,2,3,4);

    // Face 043
    computeFaceTriagNormal(4,nodes,0,4,3);

    return _normals;
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

    CFreal zeta = mappedCoord[2];
    //Bounds for eta and xi for a given zeta
    CFreal posBound = 0.5*(1-zeta);
    CFreal negBound = -posBound;
    if((mappedCoord[0] > posBound) ||
      (mappedCoord[1] > posBound) ||
      (mappedCoord[2] > 1.) ||
      (mappedCoord[0] < negBound) ||
      (mappedCoord[1] < negBound) ||
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
  LagrangeShapeFunctionPyramP1()
  {
    for(CFuint i=0;i<_normals.size();i++)
    {
      _normals[i].resize(DIM_3D);
    }
  }

  /// Default destructor
  ~LagrangeShapeFunctionPyramP1();

  /// Computes the normal to the triangle.
  /// Note that the order of indexes is important
  /// Cycliv permutations are allowed:
  /// Fijk = Fjki = Fkij
  static void computeFaceTriagNormal(const CFuint& n,
                          const std::vector<Framework::Node*>& nodes,
                          const CFuint& i,
                          const CFuint& j,
                          const CFuint& k)
  {
    _vec1[XX] = (*nodes[j])[XX] - (*nodes[i])[XX];
    _vec1[YY] = (*nodes[j])[YY] - (*nodes[i])[YY];
    _vec1[ZZ] = (*nodes[j])[ZZ] - (*nodes[i])[ZZ];

    _vec2[XX] = (*nodes[k])[XX] - (*nodes[j])[XX];
    _vec2[YY] = (*nodes[k])[YY] - (*nodes[j])[YY];
    _vec2[ZZ] = (*nodes[k])[ZZ] - (*nodes[j])[ZZ];

    MathTools::MathFunctions::crossProd(_vec1,_vec2,_vec3);

    _vec3 *= 0.5;

    _normals[n][XX] = _vec3[XX];
    _normals[n][YY] = _vec3[YY];
    _normals[n][ZZ] = _vec3[ZZ];
  }

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

  /// Temporary Vector for the computation of normals
  static RealVector _vec1;
  static RealVector _vec2;
  static RealVector _vec3;

  /// Vector of normals
  static std::vector<RealVector> _normals;
}; // end of class LagrangeShapeFunctionPyramP1

//////////////////////////////////////////////////////////////////////////////

    } // namespace ShapeFunctions



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

// #include "ShapeFunctions/LagrangeShapeFunctionPyramP1.ci"

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_ShapeFunctions_LagrangeShapeFunctionPyramP1_hh
