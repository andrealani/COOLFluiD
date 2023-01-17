// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_ShapeFunctions_LagrangeShapeFunctionTriagP2_hh
#define COOLFluiD_ShapeFunctions_LagrangeShapeFunctionTriagP2_hh

//////////////////////////////////////////////////////////////////////////////

#include "Common/StringOps.hh"
#include "ShapeFunctions/LagrangeShapeFunction.hh"
#include "Framework/FaceJacobiansDeterminant.hh"
#include "Common/NotImplementedException.hh"
#include "ShapeFunctions/ShapeFunctionsAPI.hh"
#include "MathTools/RealMatrix.hh"
#include "MathTools/MatrixInverterT.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace ShapeFunctions {

//////////////////////////////////////////////////////////////////////////////

/// This class provides the lagrangian shape function describing the
/// representation of the solution and/or the geometry in a P2 (quadratic)
/// triangular element.
/// @author Thomas Wuilbaut

class ShapeFunctions_API LagrangeShapeFunctionTriagP2 : public ShapeFunctions::LagrangeShapeFunction {
public:

  /// Get the inheritant dimensionality of the ShapeFunction
  static CFuint getDimensionality()
  {
    return DIM_2D;
  }

  /// Get the number of nodal shape functions in the element
  static CFuint getNbNodes()
  {
    return 6;
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
    return CFPolyOrder::ORDER2;
  }

  static void getStatesMappedCoordinates(std::vector<RealVector>& mappedCoords)
  {

    cf_assert (mappedCoords.size() == 6);

    mappedCoords[0][KSI] = 0. ;
    mappedCoords[0][ETA] = 0. ;

    mappedCoords[1][KSI] = 1. ;
    mappedCoords[1][ETA] = 0. ;

    mappedCoords[2][KSI] = 0. ;
    mappedCoords[2][ETA] = 1. ;

    mappedCoords[3][KSI] = 0.5 ;
    mappedCoords[3][ETA] = 0.  ;

    mappedCoords[4][KSI] = 0.5 ;
    mappedCoords[4][ETA] = 0.5 ;

    mappedCoords[5][KSI] = 0.  ;
    mappedCoords[5][ETA] = 0.5 ;
    
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
      shapeFunc[0] = (1.0 - mappedCoord.sum())*(2.*(1.0 - mappedCoord.sum()) - 1.);
      shapeFunc[1] = mappedCoord[0] * (2. * mappedCoord[0] - 1.);
      shapeFunc[2] = mappedCoord[1] * (2. * mappedCoord[1] - 1.);
      shapeFunc[3] = 4.*mappedCoord[0]*(1.0 - mappedCoord.sum());
      shapeFunc[4] = 4.*mappedCoord[0]*mappedCoord[1];
      shapeFunc[5] = 4.*mappedCoord[1]*(1.0 - mappedCoord.sum());
  }

  /// Compute the Gradient of the Shape Function
  static void computeGradientStates(
         const std::vector<RealMatrix>& jacob,
         const std::vector<RealVector>& mappedCoord,
               std::vector<RealMatrix>& grad)
  {
    MathTools::MatrixInverterT<2> inverter;

    for (CFuint ip = 0; ip < mappedCoord.size(); ++ip) {
      inverter.invert(jacob[ip], m_invJ);
      RealMatrix& lgrad = grad[ip];

      const CFreal xi =  mappedCoord[ip][KSI];
      const CFreal eta = mappedCoord[ip][ETA];

      const CFreal dN0dxi = -3. + 4.*eta + 4.*xi;
      const CFreal dN1dxi = 4.*xi - 1.;
      const CFreal dN2dxi = 0.;
      const CFreal dN3dxi = 4. - 8.*xi - 4.*eta;
      const CFreal dN4dxi = 4.*eta;
      const CFreal dN5dxi = -4.*eta;

      const CFreal dN0deta = -3. + 4.*eta + 4.*xi;
      const CFreal dN1deta = 0.;
      const CFreal dN2deta = 4.*eta - 1.;
      const CFreal dN3deta = -4.*xi;
      const CFreal dN4deta = 4.*xi;
      const CFreal dN5deta = 4. - 4.*xi - 8.*eta;
      
      
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

    const CFreal x3 = (*nodes[3])[XX];
    const CFreal y3 = (*nodes[3])[YY];
    
    const CFreal x4 = (*nodes[4])[XX];
    const CFreal y4 = (*nodes[4])[YY];

    const CFreal x5 = (*nodes[5])[XX];
    const CFreal y5 = (*nodes[5])[YY];
    
    RealVector _shapeFunc;
    
    _shapeFunc.resize(9);
    
    for (CFuint ip = 0; ip < mappedCoord.size(); ++ip)
    {
      RealVector& pointNormal = normal[ip];
      
      const CFreal xi =  mappedCoord[ip][KSI];
      const CFreal eta = mappedCoord[ip][ETA];
      
      if (planeIdx[ip] == 0)  // normal to face nb 1
      {
          const CFreal dN0dxi = -3. + 4.*eta + 4.*xi;
      const CFreal dN1dxi = 4.*xi - 1.;
      const CFreal dN2dxi = 0.;
      const CFreal dN3dxi = 4. - 8.*xi - 4.*eta;
      const CFreal dN4dxi = 4.*eta;
      const CFreal dN5dxi = -4.*eta;
        
        pointNormal[XX] = +(y0*dN0dxi  + y1*dN1dxi +y2*dN2dxi +y3*dN3dxi + y4*dN4dxi + y5*dN5dxi);
        pointNormal[YY] = -(x0*dN0dxi  + x1*dN1dxi +x2*dN2dxi +x3*dN3dxi + x4*dN4dxi + x5*dN5dxi);
        }
      else if (planeIdx[ip] == 1) // normal to face nb 2
      {
        const CFreal dN0deta = -3. + 4.*eta + 4.*xi;
      const CFreal dN1deta = 0.;
      const CFreal dN2deta = 4.*eta - 1.;
      const CFreal dN3deta = -4.*xi;
      const CFreal dN4deta = 4.*xi;
      const CFreal dN5deta = 4. - 4.*xi - 8.*eta;

        pointNormal[XX] = +(y0*dN0deta  + y1*dN1deta +y2*dN2deta +y3*dN3deta + y4*dN4deta + y5*dN5deta);
        pointNormal[YY] = -(x0*dN0deta + x1*dN1deta + x2*dN2deta + x3*dN3deta + x4*dN4deta + x5*dN5deta);
        
        const CFreal dN0dxi = -3. + 4.*eta + 4.*xi;
      const CFreal dN1dxi = 4.*xi - 1.;
      const CFreal dN2dxi = 0.;
      const CFreal dN3dxi = 4. - 8.*xi - 4.*eta;
      const CFreal dN4dxi = 4.*eta;
      const CFreal dN5dxi = -4.*eta;
        
        pointNormal[XX] += -(y0*dN0dxi  + y1*dN1dxi +y2*dN2dxi +y3*dN3dxi + y4*dN4dxi + y5*dN5dxi);
        pointNormal[YY] += +(x0*dN0dxi  + x1*dN1dxi +x2*dN2dxi +x3*dN3dxi + x4*dN4dxi + x5*dN5dxi);
      }
      else if (planeIdx[ip] == 2) // normal to face nb 3
      {
      const CFreal dN0deta = -3. + 4.*eta + 4.*xi;
      const CFreal dN1deta = 0.;
      const CFreal dN2deta = 4.*eta - 1.;
      const CFreal dN3deta = -4.*xi;
      const CFreal dN4deta = 4.*xi;
      const CFreal dN5deta = 4. - 4.*xi - 8.*eta;

        pointNormal[XX] = -(y0*dN0deta  + y1*dN1deta +y2*dN2deta +y3*dN3deta + y4*dN4deta + y5*dN5deta);
        pointNormal[YY] = +(x0*dN0deta + x1*dN1deta + x2*dN2deta + x3*dN3deta + x4*dN4deta + x5*dN5deta);
      }
      else if (planeIdx[ip] == 3) // vector ~ in the x direction
      {
      const CFreal dN0deta = -3. + 4.*eta + 4.*xi;
      const CFreal dN1deta = 0.;
      const CFreal dN2deta = 4.*eta - 1.;
      const CFreal dN3deta = -4.*xi;
      const CFreal dN4deta = 4.*xi;
      const CFreal dN5deta = 4. - 4.*xi - 8.*eta;

        pointNormal[XX] = +(y0*dN0deta  + y1*dN1deta +y2*dN2deta +y3*dN3deta + y4*dN4deta + y5*dN5deta);
        pointNormal[YY] = -(x0*dN0deta + x1*dN1deta + x2*dN2deta + x3*dN3deta + x4*dN4deta + x5*dN5deta);
        
        }
      else if (planeIdx[ip] == 4) // vector ~ in the y direction
      {
        const CFreal dN0dxi = -3. + 4.*eta + 4.*xi;
      const CFreal dN1dxi = 4.*xi - 1.;
      const CFreal dN2dxi = 0.;
      const CFreal dN3dxi = 4. - 8.*xi - 4.*eta;
      const CFreal dN4dxi = 4.*eta;
      const CFreal dN5dxi = -4.*eta;
        
        pointNormal[XX] = -(y0*dN0dxi  + y1*dN1dxi +y2*dN2dxi +y3*dN3dxi + y4*dN4dxi + y5*dN5dxi);
        pointNormal[YY] = +(x0*dN0dxi  + x1*dN1dxi +x2*dN2dxi +x3*dN3dxi + x4*dN4dxi + x5*dN5dxi);
      }
    }
  }

  

  /// Compute the Jacobian
  static void computeJacobian(
         const std::vector<Framework::Node*>& nodes,
         const std::vector<RealVector>& mappedCoord,
               std::vector<RealMatrix>& jacob)
  {
    /// @todo must be finished
    throw Common::NotImplementedException (FromHere(),getName()+"::computeJacobian()");

    static RealVector x,y;

    x.resize(6);
    y.resize(6);

    static RealVector dNdxi, dNdeta;
    dNdxi.resize(6);
    dNdeta.resize(6);

    for(CFuint ip = 0; ip < x.size(); ++ip) {
      x[ip] = (*nodes[ip])[XX];
      y[ip] = (*nodes[ip])[YY];
}

    for (CFuint ip = 0; ip < mappedCoord.size(); ++ip) {

      const CFreal xi =  mappedCoord[ip][KSI];
      const CFreal eta = mappedCoord[ip][ETA];

      dNdxi[0] = -3. + 4.*eta + 4.*xi;
      dNdxi[1] = 4.*xi - 1.;
      dNdxi[2] = 0.;
      dNdxi[3] = 4. - 8.*xi - 4.*eta;
      dNdxi[4] = 4.*eta;
      dNdxi[5] = -4.*eta;

      dNdeta[0] = -3. + 4.*eta + 4.*xi;
      dNdeta[1] = 0.;
      dNdeta[2] = 4.*eta - 1.;
      dNdeta[3] = -4.*xi;
      dNdeta[4] = 4.*xi;
      dNdeta[5] = 4. - 4.*xi - 8.*eta;

CFreal dxdxi  = 0.0;
CFreal dydxi  = 0.0;

CFreal dxdeta = 0.0;
CFreal dydeta = 0.0;

for(CFuint intp = 0; intp < x.size(); ++intp) {
      dxdxi += dNdxi[intp]*x[intp];
      dydxi += dNdxi[intp]*y[intp];
      dxdeta += dNdeta[intp]*x[intp];
      dydeta += dNdeta[intp]*y[intp];
}
    }

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
        cf_assert(nodes.size() == getNbNodes());
        
    static RealVector x,y;

    x.resize(6);
    y.resize(6);

    static RealVector dNdxi, dNdeta;
    dNdxi.resize(6);
    dNdeta.resize(6);
    
    for(CFuint ip = 0; ip < x.size(); ++ip) {
          x[ip] = (*nodes[ip])[XX];
          y[ip] = (*nodes[ip])[YY];
      }
    
    for (CFuint ip = 0; ip < mappedCoord.size(); ++ip) {
        const CFreal xi =  mappedCoord[ip][KSI];
        const CFreal eta = mappedCoord[ip][ETA];
        
        dNdxi[0] = -3. + 4.*eta + 4.*xi;
        dNdxi[1] = 4.*xi - 1.;
        dNdxi[2] = 0.;
        dNdxi[3] = 4. - 8.*xi - 4.*eta;
        dNdxi[4] = 4.*eta;
        dNdxi[5] = -4.*eta;
        
        dNdeta[0] = -3. + 4.*eta + 4.*xi;
        dNdeta[1] = 0.;
        dNdeta[2] = 4.*eta - 1.;
        dNdeta[3] = -4.*xi;
        dNdeta[4] = 4.*xi;
        dNdeta[5] = 4. - 4.*xi - 8.*eta;
        
    CFreal dxdxi  = 0.0;
    CFreal dydxi  = 0.0;

    CFreal dxdeta = 0.0;
    CFreal dydeta = 0.0;

    for(CFuint intp = 0; intp < x.size(); ++intp) {
          dxdxi += dNdxi[intp]*x[intp];
          dydxi += dNdxi[intp]*y[intp];
          dxdeta += dNdeta[intp]*x[intp];
          dydeta += dNdeta[intp]*y[intp];
    }
    const CFreal jacob = (dxdxi*dydeta)-(dxdeta*dydxi);

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
    return "LagrangeTriagP2";
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
    const Framework::Node& n1 = (*nodes[0]);
    const Framework::Node& n2 = (*nodes[1]);
    const Framework::Node& n3 = (*nodes[2]);
    const Framework::Node& n4 = (*nodes[3]);
    const Framework::Node& n5 = (*nodes[4]);
    const Framework::Node& n6 = (*nodes[5]);

  return 1./6. * (n2[XX]*n1[YY]-n1[XX]*n2[YY]+n3[XX]*n2[YY]-n2[XX]*n3[YY]-n3[XX]*n1[YY]+n1[XX]*n3[YY]) + \
             + 2./3. * (n1[XX]*n4[YY]-n4[XX]*n1[YY]+n2[XX]*n5[YY]-n5[XX]*n2[YY]-n1[XX]*n6[YY]+n6[XX]*n1[YY] + \
  n4[XX]*n2[YY]-n2[XX]*n4[YY]+n5[XX]*n3[YY]-n3[XX]*n5[YY]-n6[XX]*n3[YY]+n3[XX]*n6[YY]);

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
  LagrangeShapeFunctionTriagP2();

  /// Default destructor
  ~LagrangeShapeFunctionTriagP2();

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


}; // end of class LagrangeShapeFunctionTriagP2

//////////////////////////////////////////////////////////////////////////////

    } // namespace ShapeFunctions



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

// #include "ShapeFunctions/LagrangeShapeFunctionTriagP2.ci"

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_ShapeFunctions_LagrangeShapeFunctionTriagP2_hh
