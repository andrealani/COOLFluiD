// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_ShapeFunctions_LagrangeShapeFunctionPrismP1_hh
#define COOLFluiD_ShapeFunctions_LagrangeShapeFunctionPrismP1_hh

//////////////////////////////////////////////////////////////////////////////

#include "ShapeFunctions/LagrangeShapeFunction.hh"
#include "ShapeFunctions/ShapeFunctionsAPI.hh"
#include "Framework/VolumeCalculator.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



  namespace ShapeFunctions {

//////////////////////////////////////////////////////////////////////////////

/// This class provides the lagrangian shape function describing the
/// representation of the solution and/or the geometry in a prismatic
/// element.
/// This is an explicit specialization for order = CFPolyOrder::ORDER1
/// @author Tiago Quintino
class ShapeFunctions_API LagrangeShapeFunctionPrismP1 : public ShapeFunctions::LagrangeShapeFunction {
public:

  /// Get the inheritant dimensionality of the ShapeFunction
  static CFuint getDimensionality()
  {
    return DIM_3D;
  }

  /// Get the number of nodal shape functions in the element
  static CFuint getNbNodes()
  {
    return 6;
  }

  /// Get the number of faces
  static CFuint getNbFaces()
  {
    return 5;
  }

  /// Gets the type of CFGeoShape::Type
  static CFGeoShape::Type getShape()
  {
    return CFGeoShape::PRISM;
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
    const CFreal xi   = mappedCoord[KSI];
    const CFreal eta  = mappedCoord[ETA];
    const CFreal zeta = mappedCoord[ZTA];

    const CFreal lbd = 1.0 - xi - eta;
    const CFreal a   = 0.5 * (1.0 - zeta);
    const CFreal b   = 0.5 * (1.0 + zeta);

    shapeFunc[0] = lbd * a;
    shapeFunc[1] = xi  * a;
    shapeFunc[2] = eta * a;
    shapeFunc[3] = lbd * b;
    shapeFunc[4] = xi  * b;
    shapeFunc[5] = eta * b;
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
                                            std::vector<RealVector>& normal);

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
               std::valarray<CFreal>& detJacobian);
  /// Compute the determinant of the face jacobian at the quadrature points
  /// @todo fix this function
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
    static const std::string name = "LagrangePrismP1";
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
    return volumeCalc.calculatePrismVolume(nodes);
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
    throw Common::NotImplementedException (FromHere(), getName()+"::computeMappedCoordinates()");
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

    // Face 021
    computeFaceTriagNormal(0,nodes,0,2,1);

    // Face 345
    computeFaceTriagNormal(1,nodes,3,4,5);

    // Face 0143
    computeFaceQuadNormal(2,nodes,0,1,4,3);

    // Face 1254
    computeFaceQuadNormal(3,nodes,1,2,5,4);

    // Face 0352
    computeFaceQuadNormal(4,nodes,0,3,5,2);

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

    // Face 021
    computeFaceTriagNormal(0,nodes,0,2,1);

    // Face 345
    computeFaceTriagNormal(1,nodes,3,4,5);

    // Face 0143
    computeFaceQuadNormal(2,nodes,0,1,4,3);

    // Face 1254
    computeFaceQuadNormal(3,nodes,1,2,5,4);

    // Face 0352
    computeFaceQuadNormal(4,nodes,0,3,5,2);

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
    if( (mappedCoord[0] >= 0.) &&
        (mappedCoord[1] >= 0.) &&
        (mappedCoord[0]+mappedCoord[1] <= 1.) &&
        (mappedCoord[2] >= -1.) &&
        (mappedCoord[2] <= 1.))
    {
      return true;
    }
    else
    {
      return false;
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
  LagrangeShapeFunctionPrismP1()
  {
    for(CFuint i=0;i<_normals.size();i++)
    {
      _normals[i].resize(DIM_3D);
    }
  }

  /// Default destructor
  ~LagrangeShapeFunctionPrismP1();

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
  
  /// gradients of shape functions
  static std::vector< RealVector > _gradShapFunc;

  /// Vector of normals
  static std::vector<RealVector> _normals;

}; // end of class LagrangeShapeFunctionPrismP1

//////////////////////////////////////////////////////////////////////////////

    } // namespace ShapeFunctions



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

// #include "ShapeFunctions/LagrangeShapeFunctionPrismP1.ci"

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_ShapeFunctions_LagrangeShapeFunctionPrismP1_hh
