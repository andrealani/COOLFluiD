// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_ShapeFunctions_LagrangeShapeFunctionTetraP1_hh
#define COOLFluiD_ShapeFunctions_LagrangeShapeFunctionTetraP1_hh

//////////////////////////////////////////////////////////////////////////////

#include "MathTools/MatrixInverterT.hh"
#include "MathTools/MathConsts.hh"

#include "ShapeFunctions/LagrangeShapeFunction.hh"
#include "Framework/VolumeCalculator.hh"
#include "ShapeFunctions/ShapeFunctionsAPI.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace ShapeFunctions {

//////////////////////////////////////////////////////////////////////////////

/// This class provides the lagrangian shape function
/// of the first order for a tetrahedron.
/// @author Andrea Lani
/// @author Tiago Quintino
class ShapeFunctions_API LagrangeShapeFunctionTetraP1 : public ShapeFunctions::LagrangeShapeFunction {
public:

  /// Get the inheritant dimensionality of the ShapeFunction
  static CFuint getDimensionality()
  {
    return DIM_3D;
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
    return CFGeoShape::TETRA;
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

    cf_assert (mappedCoords.size() == 4);

    mappedCoords[0][KSI] = 0. ;
    mappedCoords[0][ETA] = 0. ;
    mappedCoords[0][ZTA] = 0. ;

    mappedCoords[1][KSI] = 1. ;
    mappedCoords[1][ETA] = 0. ;
    mappedCoords[1][ZTA] = 0. ;

    mappedCoords[2][KSI] = 0. ;
    mappedCoords[2][ETA] = 1. ;
    mappedCoords[2][ZTA] = 0. ;

    mappedCoords[3][KSI] = 0. ;
    mappedCoords[3][ETA] = 0. ;
    mappedCoords[3][ZTA] = 1. ;

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
         const RealVector& mappedCoord,
               RealVector& shapeFunc)
  {
    cf_assert(mappedCoord.size() == DIM_3D);
    shapeFunc[0] = 1.0 - mappedCoord.sum();
    shapeFunc[1] = mappedCoord[0];
    shapeFunc[2] = mappedCoord[1];
    shapeFunc[3] = mappedCoord[2];
  }

  /// Compute the Gradient of the Shape Function
  static void computeGradientStates(
         const std::vector<RealMatrix>& jacob,
         const std::vector<RealVector>& mappedCoord,
               std::vector<RealMatrix>& grad)
  {
    static RealMatrix invJ(3,3);
    MathTools::MatrixInverterT<3> inverter;

    // Derivatives of shape functions are constant
    // hence Gradients are independent of the mappedCoord
    for (CFuint ip = 0; ip < mappedCoord.size(); ++ip) {

      inverter.invert(jacob[ip], invJ);
      RealMatrix& lgrad = grad[ip];

      lgrad(0,XX) = -(invJ(0,0) + invJ(0,1) + invJ(0,2));
      lgrad(0,YY) = -(invJ(1,0) + invJ(1,1) + invJ(1,2));
      lgrad(0,ZZ) = -(invJ(2,0) + invJ(2,1) + invJ(2,2));

      lgrad(1,XX) = invJ(0,0);
      lgrad(1,YY) = invJ(1,0);
      lgrad(1,ZZ) = invJ(2,0);

      lgrad(2,XX) = invJ(0,1);
      lgrad(2,YY) = invJ(1,1);
      lgrad(2,ZZ) = invJ(2,1);

      lgrad(3,XX) = invJ(0,2);
      lgrad(3,YY) = invJ(1,2);
      lgrad(3,ZZ) = invJ(2,2);
    }
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

    const CFreal dxdksi = -x0 + x1;
    const CFreal dydksi = -y0 + y1;
    const CFreal dzdksi = -z0 + z1;

    const CFreal dxdeta = -x0 + x2;
    const CFreal dydeta = -y0 + y2;
    const CFreal dzdeta = -z0 + z2;

    const CFreal dxdzta = -x0 + x3;
    const CFreal dydzta = -y0 + y3;
    const CFreal dzdzta = -z0 + z3;

    // Derivatives of shape functions are constant
    // hence Jacobians are independent of the mappedCoord
    for (CFuint ip = 0; ip < mappedCoord.size(); ++ip) {
      RealMatrix& pointJacob = jacob[ip];

      pointJacob(KSI,XX) = dxdksi;
      pointJacob(KSI,YY) = dydksi;
      pointJacob(KSI,ZZ) = dzdksi;

      pointJacob(ETA,XX) = dxdeta;
      pointJacob(ETA,YY) = dydeta;
      pointJacob(ETA,ZZ) = dzdeta;

      pointJacob(ZTA,XX) = dxdzta;
      pointJacob(ZTA,YY) = dydzta;
      pointJacob(ZTA,ZZ) = dzdzta;
    }
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

    const CFuint nbQdPts = mappedCoord.size();

    const CFreal jacob =
        x2*y1*z0 - x3*y1*z0 - x1*y2*z0 + x3*y2*z0 + x1*y3*z0 -
        x2*y3*z0 - x2*y0*z1 + x3*y0*z1 + x0*y2*z1 - x3*y2*z1 -
        x0*y3*z1 + x2*y3*z1 + x1*y0*z2 - x3*y0*z2 - x0*y1*z2 +
        x3*y1*z2 + x0*y3*z2 - x1*y3*z2 - x1*y0*z3 + x2*y0*z3 +
        x0*y1*z3 - x2*y1*z3 - x0*y2*z3 + x1*y2*z3;

    for (CFuint ip = 0; ip < nbQdPts; ++ip) {
      cf_assert(mappedCoord[ip].size() == DIM_3D);

      // could be dependent on mappedCoord[ip][KSI] mappedCoord[ip][ETA] mappedCoord[ip][ZTA]
      detJacobian[ip] = jacob;
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

  /// Compute the determinant of the face jacobian at the quadrature points
  /// @todo fix this function
  static void computeFaceJacobianDeterminant(
          const std::vector<RealVector>& mappedCoord,
          const std::vector<Framework::Node*>& nodes,
          const Framework::IntegratorPattern& pattern,
                std::vector<RealVector>& faceJacobian)
  {
    cf_assert(pattern.nbSteps() == getNbFaces());
    cf_assert(faceJacobian.size() >= getNbFaces());

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

    cf_assert(pattern.nbSteps() == 4);

    const CFreal jacob0 =
      Framework::FaceJacobiansDeterminant::compute3DFaceJacobDet(x1-x0,y1-y0,z1-z0,x2-x0,y2-y0,z2-z0);

    const CFreal jacob1 =
      Framework::FaceJacobiansDeterminant::compute3DFaceJacobDet(x3-x0,y3-y0,z3-z0,x1-x0,y1-y0,z1-z0);

    /// sqrt(2) scales this face, since its parametrizations range is [0,sqrt(2)] not [0,1]
    const CFreal A = 1. / sqrt(2.0);
    const CFreal B = 1. / sqrt(3.0/2.0);
    const CFreal jacob2 =
      Framework::FaceJacobiansDeterminant::compute3DFaceJacobDet(A*(x2-x1),A*(y2-y1),A*(z2-z1),B*(x3-x1),B*(y3-y1),B*(z3-z1));

    const CFreal jacob3 =
      Framework::FaceJacobiansDeterminant::compute3DFaceJacobDet(x2-x0,y2-y0,z2-z0,x3-x0,y3-y0,z3-z0);

    const CFuint iFace0 = 0;
    for(CFuint iPoint = 0; iPoint < pattern.nbPts(iFace0); ++iPoint) {
        faceJacobian[iFace0][iPoint] = jacob0;
    }
    const CFuint iFace1 = 1;
    for(CFuint iPoint = 0; iPoint < pattern.nbPts(iFace1); ++iPoint) {
        faceJacobian[iFace1][iPoint] = jacob1;
    }
    const CFuint iFace2 = 2;
    for(CFuint iPoint = 0; iPoint < pattern.nbPts(iFace2); ++iPoint) {
        faceJacobian[iFace2][iPoint] = jacob2;
    }
    const CFuint iFace3 = 3;
    for(CFuint iPoint = 0; iPoint < pattern.nbPts(iFace3); ++iPoint) {
        faceJacobian[iFace3][iPoint] = jacob3;
    }
  }

  /// Get the name of this shape function
  static const std::string getName()
  {
    return "LagrangeTetraP1";
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

    const CFuint nbNodes = nodes.size();

    cf_assert(nodes.size() == 4);

    static RealMatrix coord(nbNodes,DIM_3D);
    static VolumeCalculator volumeCalc;

    for (CFuint iNode = 0; iNode < nbNodes; ++iNode) {
      for (CFuint iCoord = 0; iCoord < DIM_3D; ++iCoord) {
        coord(iNode, iCoord) = (*nodes[iNode])[iCoord];
      }
    }
    return volumeCalc.calculateTetraVolume(coord);
  }

  /// Get the centroid
  static RealVector computeCentroid(const std::vector<Framework::Node*>& nodes)
  {
    RealVector centroid(*(nodes[0]));

    centroid += *(nodes[1]);
    centroid += *(nodes[2]);
    centroid += *(nodes[3]);
    centroid *= 0.25;

    return centroid;
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
      (FromHere(), getName()+"::computeMappedCoordinatesPlus1D() ");
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

    // Face 021
    computeFaceTriagNormal(0,nodes,0,2,1);

    // Face 013
    computeFaceTriagNormal(1,nodes,0,1,3);

    // Face 123
    computeFaceTriagNormal(2,nodes,1,2,3);

    // Face 032
    computeFaceTriagNormal(3,nodes,0,3,2);

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

    // Face 021
    computeFaceTriagNormal(0,nodes,0,2,1);

    // Face 013
    computeFaceTriagNormal(1,nodes,0,1,3);

    // Face 123
    computeFaceTriagNormal(2,nodes,1,2,3);

    // Face 032
    computeFaceTriagNormal(3,nodes,0,3,2);

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
        (mappedCoord[2] >= 0.) &&
        (mappedCoord.sum() <= 1.))
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
    using namespace MathTools;
    cf_assert(nodes.size() == 4);
    cf_assert(nodes[0]->size() == DIM_3D);
    cf_assert(nodes[1]->size() == DIM_3D);
    cf_assert(nodes[2]->size() == DIM_3D);
    cf_assert(nodes[3]->size() == DIM_3D);
    CFuint scalarProduct = 0;

    // Face 021
    computeFaceTriagNormal(0,nodes,0,2,1);
    _vec1 = (*nodes[0]) - coord;

    //scalar product of AP with outward pointing normal
    CFreal inside = 0.;
    for(CFuint iDim=0; iDim < DIM_2D; ++iDim){
      inside += _normals[0][iDim]*_vec1[iDim];
    }
    if(inside > - MathTools::MathConsts::CFrealEps()) scalarProduct++;

    // Face 013
    computeFaceTriagNormal(1,nodes,0,1,3);
    _vec1 = (*nodes[0]) - coord;

    //scalar product of AP with outward pointing normal
    inside =0.;
    for(CFuint iDim=0; iDim < DIM_2D; ++iDim){
      inside += _normals[1][iDim]*_vec1[iDim];
    }
    if(inside > -MathTools::MathConsts::CFrealEps()) scalarProduct++;

    // Face 123
    computeFaceTriagNormal(2,nodes,1,2,3);
    _vec1 = (*nodes[1]) - coord;

    //scalar product of AP with outward pointing normal
    inside = 0.;
    for(CFuint iDim=0; iDim < DIM_2D; ++iDim){
      inside += _normals[2][iDim]*_vec1[iDim];
    }
    if(inside > -MathTools::MathConsts::CFrealEps()) scalarProduct++;


    // Face 032
    computeFaceTriagNormal(3,nodes,0,3,2);
    _vec1 = (*nodes[3]) - coord;

    //scalar product of AP with outward pointing normal
    inside = 0.;
    for(CFuint iDim=0; iDim < DIM_2D; ++iDim){
      inside += _normals[3][iDim]*_vec1[iDim];
    }
    if(inside > -MathTools::MathConsts::CFrealEps()) scalarProduct++;


    if(scalarProduct == 4) return true;

    return false;
  }

private:

  /// Constructor without arguments
  LagrangeShapeFunctionTetraP1()
  {
    for(CFuint i=0;i<_normals.size();i++)
    {
      _normals[i].resize(DIM_3D);
    }
  }

  /// Default destructor
  ~LagrangeShapeFunctionTetraP1();

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

private:

  /// solution integrator ID
  static CFuint _interpolatorID;

  /// Temporary Vector for the computation of normals
  static RealVector _vec1;
  static RealVector _vec2;
  static RealVector _vec3;
  static RealVector _vec4;

  /// gradients of shape functions
  static std::vector< RealVector > _gradShapFunc;
    
  /// Vector of normals
  static RealVector m_mappedCoord;

  /// Temporary Matrices
  static RealMatrix m_matrix1;
  static RealMatrix m_matrix2;

  /// Vector of normals
  static std::vector<RealVector> _normals;

}; // end of class LagrangeShapeFunctionTetraP1

//////////////////////////////////////////////////////////////////////////////

    } // namespace ShapeFunctions



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

// #include "ShapeFunctions/LagrangeShapeFunctionTetraP1.ci"

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_ShapeFunctions_LagrangeShapeFunctionTetraP1_hh
