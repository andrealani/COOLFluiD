// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_ShapeFunctions_LagrangeShapeFunctionHexaP1_hh
#define COOLFluiD_ShapeFunctions_LagrangeShapeFunctionHexaP1_hh

//////////////////////////////////////////////////////////////////////////////

#include "Common/StringOps.hh"
#include "ShapeFunctions/LagrangeShapeFunction.hh"
#include "Framework/VolumeCalculator.hh"

#include "ShapeFunctions/ShapeFunctionsAPI.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



  namespace ShapeFunctions {

//////////////////////////////////////////////////////////////////////////////

/// This class provides the lagrangian shape function
/// of the first order for a hexahedron.
/// @author Tiago Quintino
class ShapeFunctions_API LagrangeShapeFunctionHexaP1 : public ShapeFunctions::LagrangeShapeFunction {
public:

  /// Get the inheritant dimensionality of the ShapeFunction
  static CFuint getDimensionality()
  {
    return DIM_3D;
  }

  /// Get the number of nodal shape functions in the element
  static CFuint getNbNodes()
  {
    return 8;
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
  static void computeShapeFunction(const RealVector& mappedCoord,
                                         RealVector& shapeFunc)
  {
    const CFreal xi   = mappedCoord[0];
    const CFreal eta  = mappedCoord[1];
    const CFreal zeta = mappedCoord[2];

    const CFreal a1 = (1 + xi);
    const CFreal a2 = (1 - xi);

    const CFreal b1 = (1 + eta);
    const CFreal b2 = (1 - eta);

    const CFreal c1 = (1 + zeta);
    const CFreal c2 = (1 - zeta);

    shapeFunc[0] = a2*b2*c2;
    shapeFunc[1] = a1*b2*c2;
    shapeFunc[2] = a1*b1*c2;
    shapeFunc[3] = a2*b1*c2;
    shapeFunc[4] = a2*b2*c1;
    shapeFunc[5] = a1*b2*c1;
    shapeFunc[6] = a1*b1*c1;
    shapeFunc[7] = a2*b1*c1;

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

    const CFreal x5 = (*nodes[5])[XX];
    const CFreal y5 = (*nodes[5])[YY];
    const CFreal z5 = (*nodes[5])[ZZ];

    const CFreal x6 = (*nodes[6])[XX];
    const CFreal y6 = (*nodes[6])[YY];
    const CFreal z6 = (*nodes[6])[ZZ];

    const CFreal x7 = (*nodes[7])[XX];
    const CFreal y7 = (*nodes[7])[YY];
    const CFreal z7 = (*nodes[7])[ZZ];

    const CFuint nbQdPts = mappedCoord.size();

    for (CFuint ip = 0; ip < nbQdPts; ++ip) {
      cf_assert(mappedCoord[ip].size() == DIM_3D);

      const CFreal xi   = mappedCoord[ip][KSI];
      const CFreal eta  = mappedCoord[ip][ETA];
      const CFreal zeta = mappedCoord[ip][ZTA];

      const CFreal etaxi = (eta*xi);
      const CFreal etazeta = (eta*zeta);
      const CFreal xizeta = (xi*zeta);

      const CFreal f1 = 0.125*(-1.0 + eta + xi - etaxi);
      const CFreal f2 = 0.125*(-1.0 + eta - xi + etaxi);
      const CFreal f3 = 0.125*(-1.0 - eta - xi - etaxi);
      const CFreal f4 = 0.125*(-1.0 - eta + xi + etaxi);
      const CFreal f5 = 0.125*(+1.0 - eta - xi + etaxi);
      const CFreal f6 = 0.125*(+1.0 - eta + xi - etaxi);
      const CFreal f7 = 0.125*(+1.0 + eta + xi + etaxi);
      const CFreal f8 = 0.125*(+1.0 + eta - xi - etaxi);

      const CFreal f9  = 0.125*(-1.0 - eta - zeta - etazeta);
      const CFreal f10 = 0.125*(+1.0 + eta - zeta - etazeta);
      const CFreal f11 = 0.125*(+1.0 - eta + zeta - etazeta);
      const CFreal f12 = 0.125*(-1.0 + eta + zeta - etazeta);
      const CFreal f13 = 0.125*(+1.0 - eta - zeta + etazeta);
      const CFreal f14 = 0.125*(-1.0 + eta - zeta + etazeta);
      const CFreal f15 = 0.125*(-1.0 - eta + zeta + etazeta);
      const CFreal f16 = 0.125*(+1.0 + eta + zeta + etazeta);

      const CFreal f17 = 0.125*(-1.0 - xi - zeta - xizeta);
      const CFreal f18 = 0.125*(+1.0 + xi - zeta - xizeta);
      const CFreal f19 = 0.125*(+1.0 - xi + zeta - xizeta);
      const CFreal f20 = 0.125*(-1.0 + xi + zeta - xizeta);
      const CFreal f21 = 0.125*(+1.0 - xi - zeta + xizeta);
      const CFreal f22 = 0.125*(-1.0 + xi - zeta + xizeta);
      const CFreal f23 = 0.125*(-1.0 - xi + zeta + xizeta);
      const CFreal f24 = 0.125*(+1.0 + xi + zeta + xizeta);

      detJacobian[ip] = (f1*z0 + f2*z1 + f3*z2 + f4*z3 + f5*z4 + f6*z5 + f7*z6 + f8*z7)*
                        (-((y7*f9 + y2*f10 + y5*f11 + y0*f12 + y1*f13 + y4*f14 + y3*f15 + y6*f16)*
                        (x5*f17 + x2*f18 + x7*f19 + x0*f20 + x3*f21 + x4*f22 + x1*f23 + x6*f24)) +
                        (x7*f9 + x2*f10 + x5*f11 + x0*f12 + x1*f13 + x4*f14 + x3*f15 + x6*f16)*
                        (y5*f17 + y2*f18 + y7*f19 + y0*f20 + y3*f21 + y4*f22 + y1*f23 + y6*f24)) -
                        (f1*y0 + f2*y1 + f3*y2 + f4*y3 + f5*y4 + f6*y5 + f7*y6 + f8*y7)*
                        (-((z7*f9 + z2*f10 + z5*f11 + z0*f12 + z1*f13 + z4*f14 + z3*f15 + z6*f16)*
                        (x5*f17 + x2*f18 + x7*f19 + x0*f20 + x3*f21 + x4*f22 + x1*f23 + x6*f24)) +
                        (x7*f9 + x2*f10 + x5*f11 + x0*f12 + x1*f13 + x4*f14 + x3*f15 + x6*f16)*
                        (z5*f17 + z2*f18 + z7*f19 + z0*f20 + z3*f21 + z4*f22 + z1*f23 + z6*f24)) +
                        (x2*f3 + x7*f8 + x5*f6 + x0*f1 + x4*f5 + x1*f2 + x3*f4 + x6*f7)*
                        (-((z7*f9 + z2*f10 + z5*f11 + z0*f12 + z1*f13 + z4*f14 + z3*f15 + z6*f16)*
                        (y5*f17 + y2*f18 + y7*f19 + y0*f20 + y3*f21 + y4*f22 + y1*f23 + y6*f24)) +
                        (y7*f9 + y2*f10 + y5*f11 + y0*f12 + y1*f13 + y4*f14 + y3*f15 + y6*f16)*
                        (z5*f17 + z2*f18 + z7*f19 + z0*f20 + z3*f21 + z4*f22 + z1*f23 + z6*f24));
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
    static const std::string name = "LagrangeHexaP1";
    return name;
  }

  /// Compute the determinant of the face jacobian at the quadrature points
  static void computeFaceJacobianDeterminant(
          const std::vector<RealVector>& mappedCoord,
          const std::vector<Framework::Node*>& nodes,
          const Framework::IntegratorPattern& pattern,
                std::vector<RealVector>& faceJacobian);

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
    return volumeCalc.calculateHexaVolume(nodes);
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
/*    RealVector test(3);
    test = 0.;
    return test;*/
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
    /// Linear element -> same as computeAvgFaceNormals for all mappedCoord
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
  LagrangeShapeFunctionHexaP1()
  {
    for(CFuint i=0;i<_normals.size();i++)
    {
      _normals[i].resize(DIM_3D);
    }
  }

  /// Default destructor
  ~LagrangeShapeFunctionHexaP1();

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

  /// gradients of shape functions
  static std::vector< RealVector > _gradShapFunc;

  /// Vector of normals
  static std::vector<RealVector> _normals;

}; // end of class LagrangeShapeFunctionHexaP1

//////////////////////////////////////////////////////////////////////////////

    } // namespace ShapeFunctions



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

// #include "ShapeFunctions/LagrangeShapeFunctionHexaP1.ci"

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_ShapeFunctions_LagrangeShapeFunctionHexaP1_hh
