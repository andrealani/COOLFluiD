// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_ShapeFunctions_LagrangeShapeFunctionTriagP1_hh
#define COOLFluiD_ShapeFunctions_LagrangeShapeFunctionTriagP1_hh

//////////////////////////////////////////////////////////////////////////////

#include "Common/NotImplementedException.hh"

#include "MathTools/MatrixInverterT.hh"

#include "ShapeFunctions/LagrangeShapeFunction.hh"
#include "Framework/FaceJacobiansDeterminant.hh"
#include "ShapeFunctions/ShapeFunctionsAPI.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace ShapeFunctions {

//////////////////////////////////////////////////////////////////////////////

/// This class provides the lagrangian shape function describing the
/// representation of the solution and/or the geometry in a P1 (linear)
/// triangular element.
/// @author Andrea Lani
/// @author Geoffrey Deliege
/// @author Tiago Quintino
class ShapeFunctions_API LagrangeShapeFunctionTriagP1 : public ShapeFunctions::LagrangeShapeFunction {
public:

  /// Get the inheritant dimensionality of the ShapeFunction
  static CFuint getDimensionality()
  {
    return DIM_2D;
  }

  /// Get the number of nodal shape functions in the element
  static CFuint getNbNodes()
  {
    return 3;
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
    return CFPolyOrder::ORDER1;
  }

  static void getStatesMappedCoordinates(std::vector<RealVector>& mappedCoords)
  {

    cf_assert (mappedCoords.size() == 3);

    mappedCoords[0][KSI] = 0. ;
    mappedCoords[0][ETA] = 0. ;

    mappedCoords[1][KSI] = 1. ;
    mappedCoords[1][ETA] = 0. ;

    mappedCoords[2][KSI] = 0. ;
    mappedCoords[2][ETA] = 1. ;

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
      cf_assert(shapeFunc.size() == 3);
      shapeFunc[0] = 1.0 - mappedCoord.sum();
      shapeFunc[1] = mappedCoord[0];
      shapeFunc[2] = mappedCoord[1];
  }

  /// Compute the Gradient of the Shape Function
  static void computeGradientStates(
         const std::vector<RealMatrix>& jacob,
         const std::vector<RealVector>& mappedCoord,
               std::vector<RealMatrix>& grad)
  {
    MathTools::MatrixInverterT<2> inverter;
    // Derivatives of shape functions are constant
    // hence Gradients are independent of the mappedCoord
    for (CFuint ip = 0; ip < mappedCoord.size(); ++ip) {

      inverter.invert(jacob[ip], m_invJ);

      RealMatrix& lgrad = grad[ip];

      lgrad(0,XX) = -(m_invJ(0,0) + m_invJ(0,1));
      lgrad(0,YY) = -(m_invJ(1,0) + m_invJ(1,1));

      lgrad(1,XX) = m_invJ(0,0);
      lgrad(1,YY) = m_invJ(1,0);

      lgrad(2,XX) = m_invJ(0,1);
      lgrad(2,YY) = m_invJ(1,1);
    }
  }

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
                                            std::vector<RealVector>& normal)
  {
      const CFreal x0 = (*nodes[0])[XX];
      const CFreal y0 = (*nodes[0])[YY];

      const CFreal x1 = (*nodes[1])[XX];
      const CFreal y1 = (*nodes[1])[YY];

      const CFreal x2 = (*nodes[2])[XX];
      const CFreal y2 = (*nodes[2])[YY];

      for (CFuint ip = 0; ip < mappedCoord.size(); ++ip)
      {
        RealVector& pointNormal = normal[ip];

        const CFreal xi =  mappedCoord[ip][KSI];
        const CFreal eta = mappedCoord[ip][ETA];
        
          if (planeIdx[ip] == 0)  // normal to face nb 1
          {
            const CFreal dN0dxi = -1.;
            const CFreal dN1dxi = 1.;
            const CFreal dN2dxi = 0.;
            
            pointNormal[XX] = +(y0*dN0dxi  + y1*dN1dxi  + y2*dN2dxi);
            pointNormal[YY] = -(x0*dN0dxi  + x1*dN1dxi  + x2*dN2dxi);
          }
          else if (planeIdx[ip] == 1) // normal to face nb 2
          {
            const CFreal dN0dxi = 0.;
            const CFreal dN1dxi = 1.;
            const CFreal dN2dxi = -1.;
            
            pointNormal[XX] = -(y0*dN0dxi  + y1*dN1dxi  + y2*dN2dxi);
            pointNormal[YY] = +(x0*dN0dxi  + x1*dN1dxi  + x2*dN2dxi);
          }
          else if (planeIdx[ip] == 2) // normal to face nb 3
          {
            const CFreal dN0deta = -1.;
            const CFreal dN1deta = 0.;
            const CFreal dN2deta = 1.;
            
            pointNormal[XX] = -(y0*dN0deta  + y1*dN1deta  + y2*dN2deta);
            pointNormal[YY] = +(x0*dN0deta  + x1*dN1deta  + x2*dN2deta);
          }
          else if (planeIdx[ip] == 3) // vector ~ in the x direction
          {
            const CFreal dN0deta = -1.;
            const CFreal dN1deta = 0.;
            const CFreal dN2deta = 1.;

            pointNormal[XX] = +(y0*dN0deta + y1*dN1deta + y2*dN2deta);
            pointNormal[YY] = -(x0*dN0deta + x1*dN1deta + x2*dN2deta);
          }
          
        else if (planeIdx[ip] == 4) // vector ~ in the y direction
        {
            const CFreal dN0dxi = -1.;
          const CFreal dN1dxi = 1.;
          const CFreal dN2dxi = 0.;
          
          pointNormal[XX] = -(y0*dN0dxi  + y1*dN1dxi  + y2*dN2dxi);
          pointNormal[YY] = +(x0*dN0dxi  + x1*dN1dxi  + x2*dN2dxi);
        }
          
      }
  }
  /// Compute the Jacobian
  static void computeJacobian(
         const std::vector<Framework::Node*>& nodes,
         const std::vector<RealVector>& mappedCoord,
               std::vector<RealMatrix>& jacob)
  {
    const CFreal x0 = (*nodes[0])[XX];
    const CFreal y0 = (*nodes[0])[YY];

    const CFreal x1 = (*nodes[1])[XX];
    const CFreal y1 = (*nodes[1])[YY];

    const CFreal x2 = (*nodes[2])[XX];
    const CFreal y2 = (*nodes[2])[YY];

    const CFreal dxdxi  = -x0 + x1;
    const CFreal dydxi  = -y0 + y1;

    const CFreal dxdeta = -x0 + x2;
    const CFreal dydeta = -y0 + y2;

    // Derivatives of shape functions are constant
    // hence Jacobians are independent of the mappedCoord
    for (CFuint ip = 0; ip < mappedCoord.size(); ++ip) {
      RealMatrix& pointJacob = jacob[ip];
      pointJacob(KSI,XX) = dxdxi;
      pointJacob(KSI,YY) = dydxi;
      pointJacob(ETA,XX) = dxdeta;
      pointJacob(ETA,YY) = dydeta;
    }
  }

  /// Compute the Jacobian
  static void computeJacobianPlus1D(
         const std::vector<Framework::Node*>& nodes,
         const std::vector<RealVector>& mappedCoord,
               std::vector<RealMatrix>& jacob)
  {
    const CFreal x1x0 = (*nodes[1])[XX] - (*nodes[0])[XX];
    const CFreal y1y0 = (*nodes[1])[YY] - (*nodes[0])[YY];
    const CFreal z1z0 = (*nodes[1])[ZZ] - (*nodes[0])[ZZ];

    const CFreal x2x0 = (*nodes[2])[XX] - (*nodes[0])[XX];
    const CFreal y2y0 = (*nodes[2])[YY] - (*nodes[0])[YY];
    const CFreal z2z0 = (*nodes[2])[ZZ] - (*nodes[0])[ZZ];

    m_cross[XX] = y1y0 * z2z0 - z1z0 * y2y0;
    m_cross[YY] = z1z0 * x2x0 - x1x0 * z2z0;
    m_cross[ZZ] = x1x0 * y2y0 - y1y0 * x2x0;

    m_cross.normalize();

    m_jac3d(XX,KSI) = x1x0;
    m_jac3d(XX,ETA) = y1y0;
    m_jac3d(XX,ZTA) = z1z0;
    m_jac3d(YY,KSI) = x2x0;
    m_jac3d(YY,ETA) = y2y0;
    m_jac3d(YY,ZTA) = z2z0;
    m_jac3d(ZZ,KSI) = m_cross[XX];
    m_jac3d(ZZ,ETA) = m_cross[YY];
    m_jac3d(ZZ,ZTA) = m_cross[ZZ];

    for (CFuint ip = 0; ip < mappedCoord.size(); ++ip) {
      jacob[ip] = m_jac3d;
    }
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

    const CFreal x1 = (*nodes[1])[XX];
    const CFreal y1 = (*nodes[1])[YY];

    const CFreal x2 = (*nodes[2])[XX];
    const CFreal y2 = (*nodes[2])[YY];

    const CFreal jacob = (x1-x0)*(y2-y0) - (x2-x0)*(y1-y0);

    detJacobian = jacob;
  }

  /// Compute the jacobian determinant at the given
  /// mapped coordinates
  static void computeJacobianDeterminantPlus1D(
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

    detJacobian = Framework::FaceJacobiansDeterminant::compute3DFaceJacobDet(x1-x0,y1-y0,z1-z0,x2-x0,y2-y0,z2-z0);
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
    return "LagrangeTriagP1";
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

    for (CFuint i = 0; i < 3; ++i) {
      for (CFuint j = 0; j < 3; ++j) {
        if (j > 0) {
          m_vmatrix(i,j) = (*nodes[i])[j-1];
        }
        else {
          m_vmatrix(i,j) = 1.0;
        }
      }
    }
    return (0.5*m_vmatrix.determ3());
  }

  /// Get the centroid
  static RealVector computeCentroid(const std::vector<Framework::Node*>& nodes)
  {
    RealVector centroid(*(nodes[0]));

    centroid += *(nodes[1]);
    centroid += *(nodes[2]);
    centroid /= 3.;

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
  static RealVector computeMappedCoordinatesPlus1D(const RealVector& coord, const std::vector<Framework::Node*>& nodes);

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
    /// Linear element -> same as computeAvgFaceNormals for all mappedCoord
    return computeAvgFaceNormals(nodes);
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
    return computeAvgCellNormal(nodes);
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
    using namespace MathTools;
    cf_assert(nodes.size() == 3);
    cf_assert(nodes.size() == 3);
    CFuint scalarProduct = 0;

    // Face 01
    computeFaceLineNormal(0,nodes,0,1);
    m_temp2DVector = (*nodes[0]) - coord;

    //scalar product of AP with outward pointing normal
    CFreal inside = 0.;
    for(CFuint iDim=0; iDim < DIM_2D; ++iDim){
      inside += m_normals[0][iDim]*m_temp2DVector[iDim];
    }
    if(inside > -MathTools::MathConsts::CFrealEps()) scalarProduct++;

    // Face 12
    computeFaceLineNormal(1,nodes,1,2);
    m_temp2DVector = (*nodes[1]) - coord;

    //scalar product of AP with outward pointing normal
    inside =0.;
    for(CFuint iDim=0; iDim < DIM_2D; ++iDim){
      inside += m_normals[1][iDim]*m_temp2DVector[iDim];
    }
    if(inside > -MathTools::MathConsts::CFrealEps()) scalarProduct++;

    // Face 20
    computeFaceLineNormal(2,nodes,2,0);
    m_temp2DVector = (*nodes[2]) - coord;

    //scalar product of AP with outward pointing normal
    inside = 0.;
    for(CFuint iDim=0; iDim < DIM_2D; ++iDim){
      inside += m_normals[2][iDim]*m_temp2DVector[iDim];
    }
    if(inside > -MathTools::MathConsts::CFrealEps()) scalarProduct++;

    if(scalarProduct == 3) return true;

    return false;
  }

private:

  /// Default constructor without arguments
  LagrangeShapeFunctionTriagP1();

  /// Default destructor
  ~LagrangeShapeFunctionTriagP1();

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

  /// temporary vectors
  static RealVector m_temp2DVector;
  static RealVector m_temp3DVector;

  /// temporary vector for a 3 by 3 matrix
  static RealMatrix m_temp3DMatrix;
  static RealMatrix m_temp3DMatrix2;

}; // end of class LagrangeShapeFunctionTriagP1

//////////////////////////////////////////////////////////////////////////////

    } // namespace ShapeFunctions



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

// #include "ShapeFunctions/LagrangeShapeFunctionTriagP1.ci"

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_ShapeFunctions_LagrangeShapeFunctionTriagP1_hh
