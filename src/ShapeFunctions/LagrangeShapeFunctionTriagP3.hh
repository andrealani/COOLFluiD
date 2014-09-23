// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_ShapeFunctions_LagrangeShapeFunctionTriagP3_hh
#define COOLFluiD_ShapeFunctions_LagrangeShapeFunctionTriagP3_hh

//////////////////////////////////////////////////////////////////////////////

#include "Common/StringOps.hh"
#include "ShapeFunctions/LagrangeShapeFunction.hh"
#include "Framework/FaceJacobiansDeterminant.hh"
#include "Common/NotImplementedException.hh"
#include "ShapeFunctions/ShapeFunctionsAPI.hh"
#include "MathTools/MatrixInverterT.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace ShapeFunctions {

//////////////////////////////////////////////////////////////////////////////

/// This class provides the Lagrangian shape function describing the
/// representation of the solution and/or the geometry in a P3 (cubic)
/// triangular element.
/// @author Thomas Wuilbaut
/// @author Martin Vymazal

class ShapeFunctions_API LagrangeShapeFunctionTriagP3 : public ShapeFunctions::LagrangeShapeFunction {
public:

  /// Get the inheritant dimensionality of the ShapeFunction
  static CFuint getDimensionality()
  {
    return DIM_2D;
  }

  /// Get the number of nodal shape functions in the element
  static CFuint getNbNodes()
  {
    return 10;
  }

  /// Get the number of faces
  static CFuint getNbFaces()
  {
    return 3;
  }

  /// Gets the type of CFGeoShape::Type
  static CFGeoShape::Type getShape()
  {
    return CFGeoShape::TRIAG;
  }

  /// Gets the type of Interpolator
  static CFPolyForm::Type getInterpolatorType()
  {
    return CFPolyForm::LAGRANGE;
  }

  /// Gets the Interpolator order
  static CFPolyOrder::Type getInterpolatorOrder()
  {
    return CFPolyOrder::ORDER3;
  }

  static void getStatesMappedCoordinates(std::vector<RealVector>& mappedCoords)
  {
    cf_assert (mappedCoords.size() == 10);

    CFreal onethird = 1.0/3.0;
    CFreal twothird = 2.0*onethird;

    mappedCoords[0][KSI] = 0.       ;
    mappedCoords[0][ETA] = 0.       ;

    mappedCoords[1][KSI] = 1.       ;
    mappedCoords[1][ETA] = 0.       ;

    mappedCoords[2][KSI] = 0.       ;
    mappedCoords[2][ETA] = 1.       ;

    mappedCoords[3][KSI] = onethird ;
    mappedCoords[3][ETA] = 0.       ;

    mappedCoords[4][KSI] = twothird ;
    mappedCoords[4][ETA] = 0.0      ;

    mappedCoords[5][KSI] = twothird ;
    mappedCoords[5][ETA] = onethird ;

    mappedCoords[6][KSI] = onethird ;
    mappedCoords[6][ETA] = twothird ;

    mappedCoords[7][KSI] = 0.0      ;
    mappedCoords[7][ETA] = twothird ;

    mappedCoords[8][KSI] = 0.0      ;
    mappedCoords[8][ETA] = onethird ;

    mappedCoords[9][KSI] = onethird ;
    mappedCoords[9][ETA] = onethird ;

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

      CFreal ninehalf = 9./2.;

      shapeFunc[0] = 0.5*(3.*(1.0 - mappedCoord.sum()) - 1.0)*(3.0*(1.0 - mappedCoord.sum()) - 2.0)*(1.0 - mappedCoord.sum());
      shapeFunc[1] = 0.5*(3.*mappedCoord[0] - 1.)*(3. * mappedCoord[0] - 2.)*mappedCoord[0];
      shapeFunc[2] = 0.5*(3.*mappedCoord[1] - 1.)*(3. * mappedCoord[1] - 2.)*mappedCoord[1];
      shapeFunc[3] = ninehalf*(1. -  mappedCoord.sum())*mappedCoord[0]*(3.*(1. -  mappedCoord.sum()) - 1.) ;
      shapeFunc[4] = ninehalf*(1. -  mappedCoord.sum())*mappedCoord[0]*(3.*mappedCoord[0] - 1.) ;
      shapeFunc[5] = ninehalf*mappedCoord[0]*mappedCoord[1]*(3.*mappedCoord[0] - 1.);
      shapeFunc[6] = ninehalf*mappedCoord[0]*mappedCoord[1]*(3.*mappedCoord[1] - 1.);
      shapeFunc[7] = ninehalf*(1. -  mappedCoord.sum())*mappedCoord[1]*(3.*mappedCoord[1] - 1.) ;
      shapeFunc[8] = ninehalf*(1. -  mappedCoord.sum())*mappedCoord[1]*(3.*(1. -  mappedCoord.sum()) - 1.) ;
      shapeFunc[9] = 27.*mappedCoord[0]*mappedCoord[1]*(1.0 - mappedCoord.sum());
  }

  /// Compute the Gradient of the Shape Function
  static void computeGradientStates(
         const std::vector<RealMatrix>& jacob,
         const std::vector<RealVector>& mappedCoord,
               std::vector<RealMatrix>& grad)
  {
    MathTools::MatrixInverterT<2> inverter;
    for (CFuint ip = 0; ip < mappedCoord.size(); ++ip)
    {
      inverter.invert(jacob[ip], m_invJ);
      RealMatrix& lgrad = grad[ip];

      const CFreal xi =  mappedCoord[ip][KSI];
      const CFreal eta = mappedCoord[ip][ETA];

      CFreal ninehalf = 9./2.;

      const CFreal dN0dxi = 0.5*(9.0*(1.0 - xi - eta)*(-3.0*(1.0 - xi - eta) + 2.0) - 2.0);
      const CFreal dN1dxi = 0.5*(27.0*xi*xi - 18.0*xi + 2.0);
      const CFreal dN2dxi = 0.0;
      const CFreal dN3dxi = ninehalf*((1.0 - xi - eta)*(-6.0*xi + 3.0*(1.0 - xi - eta) - 1.0) + xi);
      const CFreal dN4dxi = ninehalf*(6.0*(1.0 - xi -eta)*xi - 3.0*xi*(xi - 1.0) - (1.0 - xi - eta));
      const CFreal dN5dxi = ninehalf*(6.0*xi*eta - eta);
      const CFreal dN6dxi = ninehalf*(eta*(3.0*eta - 1.0));
      const CFreal dN7dxi = ninehalf*(-eta*(3.0*eta - 1.0));
      const CFreal dN8dxi = ninehalf*(-6.0*(1.0 - xi - eta) + eta);
      const CFreal dN9dxi = -27.0*xi*eta + 27.0*(1.0 - xi - eta)*eta ;

      const CFreal dN0deta = 0.5*(9.0*(1.0 - xi - eta)*(-3.0*(1.0 - xi - eta) + 2.0) - 2.0);
      const CFreal dN1deta = 0.0;
      const CFreal dN2deta = 0.5*(27.0*eta*eta - 18.0*eta + 2.0);
      const CFreal dN3deta = ninehalf*(-6.0*(1.0 - xi - eta)*xi + xi);
      const CFreal dN4deta = -ninehalf*xi*(3.0*xi - 1.0);
      const CFreal dN5deta = ninehalf*xi*(3.0*xi - 1.0);
      const CFreal dN6deta = ninehalf*(6.0*xi*eta - xi);
      const CFreal dN7deta = ninehalf*(-eta*(3.0*eta - 1.0) + (1.0 - xi - eta)*(3.0*eta - 1.0) +
                                       3.0*(1.0 - xi - eta)*eta);
      const CFreal dN8deta = ninehalf*(-6.0*(1.0 - xi - eta)*eta + eta +
                                       3.0*(1.0 - xi - eta)*(1.0 - xi - eta) -
                                       - (1.0 - xi - eta)) ;
      const CFreal dN9deta = -27.0*xi*eta + 27.0*(1.0 - xi - eta)*xi ;

      const CFreal JXX = m_invJ(0,0);
      const CFreal JXY = m_invJ(0,1);
      const CFreal JYX = m_invJ(1,0);
      const CFreal JYY = m_invJ(1,1);

      lgrad(0,XX) = dN0dxi * JXX + dN0deta * JXY;
      lgrad(0,YY) = dN0dxi * JYX + dN0deta * JYY;

      lgrad(1,XX) = dN1dxi * JXX + dN1deta * JXY;
      lgrad(1,YY) = dN1dxi * JYX + dN1deta * JYY;

      lgrad(2,XX) = dN2dxi * JXX + dN2deta * JXY;
      lgrad(2,YY) = dN2dxi * JYX + dN2deta * JYY;

      lgrad(3,XX) = dN3dxi * JXX + dN3deta * JXY;
      lgrad(3,YY) = dN3dxi * JYX + dN3deta * JYY;

      lgrad(4,XX) = dN4dxi * JXX + dN4deta * JXY;
      lgrad(4,YY) = dN4dxi * JYX + dN4deta * JYY;

      lgrad(5,XX) = dN5dxi * JXX + dN5deta * JXY;
      lgrad(5,YY) = dN5dxi * JYX + dN5deta * JYY;

      lgrad(6,XX) = dN6dxi * JXX + dN6deta * JXY;
      lgrad(6,YY) = dN6dxi * JYX + dN6deta * JYY;

      lgrad(7,XX) = dN7dxi * JXX + dN7deta * JXY;
      lgrad(7,YY) = dN7dxi * JYX + dN7deta * JYY;

      lgrad(8,XX) = dN8dxi * JXX + dN8deta * JXY;
      lgrad(8,YY) = dN8dxi * JYX + dN8deta * JYY;

      lgrad(9,XX) = dN9dxi * JXX + dN9deta * JXY;
      lgrad(9,YY) = dN9dxi * JYX + dN9deta * JYY;

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

    throw Common::NotImplementedException
      (FromHere(), getName()+"::computeJacobian()");

  }

  /// Compute the Jacobian
  static void computeJacobianPlus1D(
         const std::vector<Framework::Node*>& nodes,
         const std::vector<RealVector>& mappedCoord,
               std::vector<RealMatrix>& jacob)
  {

    throw Common::NotImplementedException
      (FromHere(), getName()+"::computeJacobianPlus1D()");
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
    throw Common::NotImplementedException
      (FromHere(), getName()+"::computeJacobianDeterminant()");
  }

  /// Compute the jacobian determinant at the given
  /// mapped coordinates
  static void computeJacobianDeterminantPlus1D(
         const std::vector<RealVector>& mappedCoord,
         const std::vector<Framework::Node*>& nodes,
               std::valarray<CFreal>& detJacobian)
  {
    throw Common::NotImplementedException
      (FromHere(), getName()+"::computeJacobianDeterminantPlus1D()");

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
                std::vector<RealVector>& faceJacobian);

  /// Get the name of this shape function
  static const std::string getName()
  {
    return "LagrangeTriagP3";
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

    CFreal volume = 0.0;

    const Framework::Node& n0 = (*nodes[0]);
    const Framework::Node& n1 = (*nodes[1]);
    const Framework::Node& n2 = (*nodes[2]);
    const Framework::Node& n3 = (*nodes[3]);
    const Framework::Node& n4 = (*nodes[4]);
    const Framework::Node& n5 = (*nodes[5]);
    const Framework::Node& n6 = (*nodes[6]);
    const Framework::Node& n7 = (*nodes[7]);
    const Framework::Node& n8 = (*nodes[8]);
    const Framework::Node& n9 = (*nodes[9]);


    for(CFuint iqdpt=0;iqdpt<m_nbQdPts;++iqdpt) {

      const CFreal xi = m_xiQdPtsTriagP3[iqdpt];
      const CFreal eta = m_etaQdPtsTriagP3[iqdpt];

      const CFreal L0 = 1.0 - xi - eta;
      const CFreal L1 = xi;
      const CFreal L2 = eta;

      const CFreal dN0dxi = -0.5*(9.0*L0*(2.0*L0-1.0) + (3.0*L0-1.0)*(3.0*L0-2.0));
      const CFreal dN1dxi = 0.5*(9.0*L1*(2.0*L1-1.0) + (3.0*L1-1.0)*(3.0*L1-2.0));
      const CFreal dN2dxi = 0.0;
      const CFreal dN3dxi = 4.5*(3.0*L0*L0 - 6.0*L0*L1 + L1 - L0);
      const CFreal dN4dxi = 4.5*(-3.0*L1*L1 + 6.0*L0*L1 + L1 - L0);
      const CFreal dN5dxi = 4.5*(6.0*L1*L2 - L2);
      const CFreal dN6dxi = 4.5*(3.0*L2*L2 - L2);
      const CFreal dN7dxi = 4.5*(L2 - 3.0*L2*L2);
      const CFreal dN8dxi = 4.5*(L2 - 6.0*L0*L2);
      const CFreal dN9dxi = 27.0*L2*(L0-L1);

      const CFreal dN0deta = -0.5*(9.0*L0*(2.0*L0-1.0) + (3.0*L0-1.0)*(3.0*L0-2.0));
      const CFreal dN1deta = 0.0;
      const CFreal dN2deta = 0.5*(9.0*L2*(2.0*L2-1.0) + (3.0*L2-1.0)*(3.0*L2-2.0));
      const CFreal dN3deta = 4.5*(L1 - 6.0*L0*L1);
      const CFreal dN4deta = 4.5*(L1 - 3.0*L1*L1);
      const CFreal dN5deta = 4.5*(3.0*L1*L1 - L1);
      const CFreal dN6deta = 4.5*(6.0*L1*L2 - L1);
      const CFreal dN7deta = 4.5*(-3.0*L2*L2 + 6.0*L0*L2 - L0 + L2);
      const CFreal dN8deta = 4.5*(3.0*L0*L0 - 6.0*L0*L2 + L2 - L0);
      const CFreal dN9deta = 27.0*L1*(L0-L2);


      const CFreal JXX = dN0dxi*n0[XX] + dN1dxi*n1[XX] + dN2dxi*n2[XX] + dN3dxi*n3[XX] + dN4dxi*n4[XX] + \
                        dN5dxi*n5[XX] + dN6dxi*n6[XX] + dN7dxi*n7[XX] + dN8dxi*n8[XX] + dN9dxi*n9[XX];

      const CFreal JXY = dN0deta*n0[XX] + dN1deta*n1[XX] + dN2deta*n2[XX] + dN3deta*n3[XX] + dN4deta*n4[XX] + \
                        dN5deta*n5[XX] + dN6deta*n6[XX] + dN7deta*n7[XX] + dN8deta*n8[XX] + dN9deta*n9[XX];

      const CFreal JYX = dN0dxi*n0[YY] + dN1dxi*n1[YY] + dN2dxi*n2[YY] + dN3dxi*n3[YY] + dN4dxi*n4[YY] + \
                        dN5dxi*n5[YY] + dN6dxi*n6[YY] + dN7dxi*n7[YY] + dN8dxi*n8[YY] + dN9dxi*n9[YY];

      const CFreal JYY = dN0deta*n0[YY] + dN1deta*n1[YY] + dN2deta*n2[YY] + dN3deta*n3[YY] + dN4deta*n4[YY] + \
                        dN5deta*n5[YY] + dN6deta*n6[YY] + dN7deta*n7[YY] + dN8deta*n8[YY] + dN9deta*n9[YY];

      volume += std::abs(JXX*JYY-JXY*JYX)*m_qdWeightsTriagP3[iqdpt]*0.5; //0.5 is the volume of the triangle in reference space

    } //loop over quadrature points

    return volume;
  }

  /// Get the centroid
  static RealVector computeCentroid(const std::vector<Framework::Node*>& nodes)
  {
    throw Common::NotImplementedException
      (FromHere(), getName()+"::computeCentroid()");
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
    cf_assert(nodes.size() == 3);
    // Face 01
    computeFaceLineNormal(0,nodes,0,1);

    // Face 12
    computeFaceLineNormal(1,nodes,1,2);

    // Face 20
    computeFaceLineNormal(2,nodes,2,0);

    return m_normals;
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

  /// Compute Cell Average Normals
  /// @param nodes contains the nodes
  /// @return RealVector containing the average cell normal
  static RealVector computeAvgCellNormal(const std::vector<Framework::Node*>& nodes)
  {
    m_vec1[XX] = (*nodes[1])[XX] - (*nodes[0])[XX];
    m_vec1[YY] = (*nodes[1])[YY] - (*nodes[0])[YY];
    m_vec1[ZZ] = (*nodes[1])[ZZ] - (*nodes[0])[ZZ];

    m_vec2[XX] = (*nodes[2])[XX] - (*nodes[1])[XX];
    m_vec2[YY] = (*nodes[2])[YY] - (*nodes[1])[YY];
    m_vec2[ZZ] = (*nodes[2])[ZZ] - (*nodes[1])[ZZ];

    MathTools::MathFunctions::crossProd(m_vec1,m_vec2,m_vnormal);

    m_vnormal *= 0.5;

    return m_vnormal;
  }

  /// Compute Cell Normal at the supplied mapped coordinate
  /// @param mappedCoord vector with xsi, eta and zeta coordinates where to compute the normal (XSI,ETA,ZTA)
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
    if( (mappedCoord[0] >= 0.) &&
        (mappedCoord[1] >= 0.) &&
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

    throw Common::NotImplementedException
      (FromHere(), getName()+"::isInElement()");
  }

private:

  /// Default constructor without arguments
  LagrangeShapeFunctionTriagP3();

  /// Default destructor
  ~LagrangeShapeFunctionTriagP3();

  /// Computes the normal to the line.
  /// Note that the order of indexes is important
  static void computeFaceLineNormal(const CFuint& n,
                          const std::vector<Framework::Node*>& nodes,
                          const CFuint& i,
                          const CFuint& j)
  {
    m_normals[n][XX] = (*nodes[j])[YY] - (*nodes[i])[YY];
    m_normals[n][YY] = (*nodes[i])[XX] - (*nodes[j])[XX];
  }

private: // data

  /// solution integrator ID
  static CFuint _interpolatorID;

  /// Temporary Vector for the computation of jacobian
  static RealVector m_dPhidxi;

  /// Temporary Vector for the computation of the jacobian
  static RealVector m_dPhideta;

  /// Temporary Vector for the computation of normals
  static RealVector m_vec1;
  /// Temporary Vector for the computation of normals
  static RealVector m_vec2;

  /// Temporary matrix for the computation of volume
  static RealMatrix m_vmatrix;

  /// Temporary matrix for the holding the 3D jacobian
  static RealMatrix m_jac3d;

  /// Temporary matrix for the holding the 2D jacobian inverse
  static RealMatrix m_invJ;

  /// Vector of normals
  static std::vector<RealVector> m_normals;

  /// normal to the triangle in 3D
  static RealVector m_vnormal;

  /// storage of the m_cross product in 3D
  static RealVector m_cross;

  /// temporary vector for mapped coordinates in 2D
  static RealVector m_mappedCoord;

  /// temporary vector for mapped coordinates in 3D
  static RealVector m_mappedCoord3D;

  /// temporary vector for computaiton of mapped coordinates
  static RealVector m_vBA;
  static RealVector m_vCA;
  static RealVector m_vPA;

  /// temporary vector for computaiton of mapped coordinates
  static RealVector m_vBA3D;
  static RealVector m_vCA3D;
  static RealVector m_vPA3D;

  /// temporary vector for computaiton of mapped coordinates
  static RealVector m_temp3DVector;

  /// temporary vector for a 3 by 3 matrix
  static RealMatrix m_temp3DMatrix;

  ///data for numerical quadrature
  static CFreal m_qdWeightsTriagP3[7];
  static CFreal m_xiQdPtsTriagP3[7];
  static CFreal m_etaQdPtsTriagP3[7];
  static CFuint m_nbQdPts;
  static CFreal m_a1;
  static CFreal m_b1;
  static CFreal m_a2;
  static CFreal m_b2;
  static CFreal m_w1;
  static CFreal m_w2;
  static CFreal m_w3;

}; // end of class LagrangeShapeFunctionTriagP3

//////////////////////////////////////////////////////////////////////////////

    } // namespace ShapeFunctions



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

// #include "ShapeFunctions/LagrangeShapeFunctionTriagP3.ci"

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_ShapeFunctions_LagrangeShapeFunctionTriagP2_hh
