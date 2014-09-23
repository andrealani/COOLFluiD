// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "ShapeFunctions/LagrangeShapeFunctionQuadP2.hh"
#include "MathTools/MathFunctions.hh"
#include "MathTools/MatrixInverterT.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

   namespace ShapeFunctions {

//////////////////////////////////////////////////////////////////////////////

CFuint LagrangeShapeFunctionQuadP2::_interpolatorID = 0;
RealVector LagrangeShapeFunctionQuadP2::_vec10 = RealVector(DIM_2D);
RealVector LagrangeShapeFunctionQuadP2::_vec20 = RealVector(DIM_2D);
RealVector LagrangeShapeFunctionQuadP2::_vecP0 = RealVector(DIM_2D);
RealVector LagrangeShapeFunctionQuadP2::_shapeFunc = RealVector(9);
RealVector LagrangeShapeFunctionQuadP2::_2Dvec1 = RealVector(DIM_2D);
RealVector LagrangeShapeFunctionQuadP2::_2Dvec2 = RealVector(DIM_2D);
RealVector LagrangeShapeFunctionQuadP2::_3Dvec1 = RealVector(DIM_3D);
RealVector LagrangeShapeFunctionQuadP2::_3Dvec2 = RealVector(DIM_3D);
RealVector LagrangeShapeFunctionQuadP2::_mappedCoord  = RealVector(DIM_2D);
RealVector LagrangeShapeFunctionQuadP2::_tempVector2D = RealVector(DIM_2D);
RealVector LagrangeShapeFunctionQuadP2::_tempVector3D = RealVector(DIM_3D);
std::vector<RealVector> LagrangeShapeFunctionQuadP2::_normals =
    std::vector<RealVector>(4,LagrangeShapeFunctionQuadP2::_tempVector2D);
std::vector<RealVector> LagrangeShapeFunctionQuadP2::_shapeFuncDerivs =
    std::vector<RealVector>(9,RealVector(DIM_2D));
//////////////////////////////////////////////////////////////////////////////

void LagrangeShapeFunctionQuadP2::computeGradientStates(
         const std::vector<RealMatrix>& jacob,
         const std::vector<RealVector>& mappedCoord,
               std::vector<RealMatrix>& grad)
{
  static RealMatrix invJ(2,2);
  MathTools::MatrixInverterT<2> inverter;

  for (CFuint ip = 0; ip < mappedCoord.size(); ++ip) {

      inverter.invert(jacob[ip], invJ);

    RealMatrix& lgrad = grad[ip];

    const CFreal xi =  mappedCoord[ip][KSI];
    const CFreal eta = mappedCoord[ip][ETA];
    const CFreal xi2 = xi*xi;
    const CFreal eta2 = eta*eta;
    const CFreal xiEta = xi*eta;
    const CFreal xi2Eta = xi2*eta;
    const CFreal xiEta2 = xi*eta2;

    const CFreal dN0dxi =  0.25 * (eta - 2.*xiEta - eta2 + 2.*xiEta2);
    const CFreal dN1dxi = -0.25 * (eta + 2.*xiEta - eta2 - 2.*xiEta2);
    const CFreal dN2dxi =  0.25 * (eta + 2.*xiEta + eta2 + 2.*xiEta2);
    const CFreal dN3dxi = -0.25 * (eta - 2.*xiEta + eta2 - 2.*xiEta2);
    const CFreal dN4dxi = -0.5  * (-2.*xiEta + 2.*xiEta2);
    const CFreal dN5dxi =  0.5  * (1. - eta2 + 2.*xi - 2.*xiEta2);
    const CFreal dN6dxi =  0.5  * (-2.*xiEta - 2.*xiEta2);
    const CFreal dN7dxi = -0.5  * (1. - eta2 - 2.*xi + 2.*xiEta2);
    const CFreal dN8dxi =  2.*xiEta2 - 2.*xi;

    const CFreal dN0deta =  0.25 * (xi - xi2 - 2.*xiEta + 2.*xi2Eta);
    const CFreal dN1deta = -0.25 * (xi + xi2 - 2.*xiEta - 2.*xi2Eta);
    const CFreal dN2deta =  0.25 * (xi + xi2 + 2.*xiEta + 2.*xi2Eta);
    const CFreal dN3deta = -0.25 * (xi - xi2 + 2.*xiEta - 2.*xi2Eta);
    const CFreal dN4deta = -0.5 * (1. - xi2 - 2.*eta + 2.*xi2Eta);
    const CFreal dN5deta =  0.5 * (-2.*xiEta - 2.*xi2Eta);
    const CFreal dN6deta =  0.5 * (1. - xi2 + 2.*eta - 2.*xi2Eta);
    const CFreal dN7deta = -0.5 * (-2.*xiEta + 2.*xi2Eta);
    const CFreal dN8deta =  2.*xi2Eta - 2.*eta;

    const CFreal JXX = invJ(0,0);
    const CFreal JXY = invJ(0,1);
    const CFreal JYX = invJ(1,0);
    const CFreal JYY = invJ(1,1);

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

  }
}

//////////////////////////////////////////////////////////////////////////////

void LagrangeShapeFunctionQuadP2::computeFaceJacobDetVectorAtMappedCoords(const std::vector<RealVector>& mappedCoord,
                                                                          const std::vector<Framework::Node*>& nodes,
                                                                          std::vector<RealVector>& normal)
{
  for (CFuint ip = 0; ip < mappedCoord.size(); ++ip)
  {
    RealVector& pointNormal = normal[ip];

    // get mapped coordinates
    const CFreal xi =  mappedCoord[ip][KSI];
    const CFreal eta = mappedCoord[ip][ETA];
    const CFreal xi2 = xi*xi;
    const CFreal eta2 = eta*eta;
    const CFreal xiEta = xi*eta;
    const CFreal xi2Eta = xi2*eta;
    const CFreal xiEta2 = xi*eta2;

    // compute shape function gradients
    _shapeFuncDerivs[0][KSI] =  (eta - eta2 - 2.*(xiEta - xiEta2));
    _shapeFuncDerivs[1][KSI] = -(eta - eta2 + 2.*(xiEta - xiEta2));
    _shapeFuncDerivs[2][KSI] =  (eta + eta2 + 2.*(xiEta + xiEta2));
    _shapeFuncDerivs[3][KSI] = -(eta + eta2 - 2.*(xiEta + xiEta2));
    _shapeFuncDerivs[4][KSI] = -4. * (-xiEta + xiEta2);
    _shapeFuncDerivs[5][KSI] =  2. * (1. - eta2 + 2.*(xi - xiEta2));
    _shapeFuncDerivs[6][KSI] =  4. * (-xiEta - xiEta2);
    _shapeFuncDerivs[7][KSI] = -2. * (1. - eta2 - 2.*(xi - xiEta2));
    _shapeFuncDerivs[8][KSI] =  8. * (xiEta2 - xi);

    _shapeFuncDerivs[0][ETA] =  (xi - xi2 - 2.*(xiEta - xi2Eta));
    _shapeFuncDerivs[1][ETA] = -(xi + xi2 - 2.*(xiEta + xi2Eta));
    _shapeFuncDerivs[2][ETA] =  (xi + xi2 + 2.*(xiEta + xi2Eta));
    _shapeFuncDerivs[3][ETA] = -(xi - xi2 + 2.*(xiEta - xi2Eta));
    _shapeFuncDerivs[4][ETA] = -2. * (1. - xi2 - 2.*(eta - xi2Eta));
    _shapeFuncDerivs[5][ETA] =  4. * (-xiEta - xi2Eta);
    _shapeFuncDerivs[6][ETA] =  2. * (1. - xi2 + 2.*(eta - xi2Eta));
    _shapeFuncDerivs[7][ETA] = -4. * (-xiEta + xi2Eta);
    _shapeFuncDerivs[8][ETA] =  8. * (xi2Eta - eta);

    // compute gradients
    _3Dvec1 = (*nodes[0])*_shapeFuncDerivs[0][KSI];
    _3Dvec2 = (*nodes[0])*_shapeFuncDerivs[0][ETA];
    for (CFuint in = 1; in < 9; ++in)
    {
      _3Dvec1 += (*nodes[in])*_shapeFuncDerivs[in][KSI];
      _3Dvec2 += (*nodes[in])*_shapeFuncDerivs[in][ETA];
    }

    // compute face jacobian vector
    MathTools::MathFunctions::crossProd(_3Dvec1,_3Dvec2,pointNormal);
    pointNormal *= 0.0625;
  }
}

//////////////////////////////////////////////////////////////////////////////

void LagrangeShapeFunctionQuadP2::computeMappedCoordPlaneNormal(const std::vector<CFuint>& planeIdx,
                                                                const std::vector<RealVector>& mappedCoord,
                                                                const std::vector<Framework::Node*>& nodes,
                                                                std::vector<RealVector>& normal)
{
  for (CFuint ip = 0; ip < mappedCoord.size(); ++ip)
  {
    RealVector& pointNormal = normal[ip];

    // get mapped coordinates
    const CFreal xi =  mappedCoord[ip][KSI];
    const CFreal eta = mappedCoord[ip][ETA];
    const CFreal xi2 = xi*xi;
    const CFreal eta2 = eta*eta;
    const CFreal xiEta = xi*eta;

    cf_assert(planeIdx[ip] == 0 || planeIdx[ip] == 1);
    if (planeIdx[ip] == 0)
    {
      const CFreal xi2Eta = xi2*eta;

      /// @note below, the derivatives of shapefunctions are computed, not the shapefunctions themselves
      _shapeFunc[0] =  (xi - xi2 - 2.*(xiEta - xi2Eta));
      _shapeFunc[1] = -(xi + xi2 - 2.*(xiEta + xi2Eta));
      _shapeFunc[2] =  (xi + xi2 + 2.*(xiEta + xi2Eta));
      _shapeFunc[3] = -(xi - xi2 + 2.*(xiEta - xi2Eta));
      _shapeFunc[4] = -2. * (1. - xi2 - 2.*(eta - xi2Eta));
      _shapeFunc[5] =  4. * (-xiEta - xi2Eta);
      _shapeFunc[6] =  2. * (1. - xi2 + 2.*(eta - xi2Eta));
      _shapeFunc[7] = -4. * (-xiEta + xi2Eta);
      _shapeFunc[8] =  8. * (xi2Eta - eta);

      pointNormal[XX] = +(*nodes[0])[YY]*_shapeFunc[0];
      pointNormal[YY] = -(*nodes[0])[XX]*_shapeFunc[0];
      for (CFuint in = 1; in < 9; ++in)
      {
        pointNormal[XX] += (*nodes[in])[YY]*_shapeFunc[in];
        pointNormal[YY] -= (*nodes[in])[XX]*_shapeFunc[in];
      }
    }
    else
    {
      const CFreal xiEta2 = xi*eta2;

      /// @note below, the derivatives of shapefunctions are computed, not the shapefunctions themselves
      _shapeFunc[0] =  (eta - eta2 - 2.*(xiEta - xiEta2));
      _shapeFunc[1] = -(eta - eta2 + 2.*(xiEta - xiEta2));
      _shapeFunc[2] =  (eta + eta2 + 2.*(xiEta + xiEta2));
      _shapeFunc[3] = -(eta + eta2 - 2.*(xiEta + xiEta2));
      _shapeFunc[4] = -4. * (-xiEta + xiEta2);
      _shapeFunc[5] =  2. * (1. - eta2 + 2.*(xi - xiEta2));
      _shapeFunc[6] =  4. * (-xiEta - xiEta2);
      _shapeFunc[7] = -2. * (1. - eta2 - 2.*(xi - xiEta2));
      _shapeFunc[8] =  8. * (xiEta2 - xi);

      pointNormal[XX] = -(*nodes[0])[YY]*_shapeFunc[0];
      pointNormal[YY] = +(*nodes[0])[XX]*_shapeFunc[0];
      for (CFuint in = 1; in < 9; ++in)
      {
        pointNormal[XX] -= (*nodes[in])[YY]*_shapeFunc[in];
        pointNormal[YY] += (*nodes[in])[XX]*_shapeFunc[in];
      }
    }
    pointNormal *= 0.25;
  }
}

//////////////////////////////////////////////////////////////////////////////

void LagrangeShapeFunctionQuadP2::computeJacobian(const std::vector<Framework::Node*>& nodes,
                                                  const std::vector<RealVector>& mappedCoord,
                                                  std::vector<RealMatrix>& jacob)
{
  std::vector< RealVector > shapeFuncDerivs(9,RealVector(2));

  for (CFuint ip = 0; ip < mappedCoord.size(); ++ip) {

    // dereference jacobian
    RealMatrix& pointJacob = jacob[ip];

    // get mapped coordinates
    const CFreal xi =  mappedCoord[ip][KSI];
    const CFreal eta = mappedCoord[ip][ETA];
    const CFreal xi2 = xi*xi;
    const CFreal eta2 = eta*eta;
    const CFreal xiEta = xi*eta;
    const CFreal xi2Eta = xi2*eta;
    const CFreal xiEta2 = xi*eta2;

    // set shape function derivatives
    _shapeFuncDerivs[0][KSI] =  0.25 * (eta - 2.*xiEta - eta2 + 2.*xiEta2);
    _shapeFuncDerivs[1][KSI] = -0.25 * (eta + 2.*xiEta - eta2 - 2.*xiEta2);
    _shapeFuncDerivs[2][KSI] =  0.25 * (eta + 2.*xiEta + eta2 + 2.*xiEta2);
    _shapeFuncDerivs[3][KSI] = -0.25 * (eta - 2.*xiEta + eta2 - 2.*xiEta2);
    _shapeFuncDerivs[4][KSI] = -0.5  * (-2.*xiEta + 2.*xiEta2);
    _shapeFuncDerivs[5][KSI] =  0.5  * (1. - eta2 + 2.*xi - 2.*xiEta2);
    _shapeFuncDerivs[6][KSI] =  0.5  * (-2.*xiEta - 2.*xiEta2);
    _shapeFuncDerivs[7][KSI] = -0.5  * (1. - eta2 - 2.*xi + 2.*xiEta2);
    _shapeFuncDerivs[8][KSI] =  2.*xiEta2 - 2.*xi;

    _shapeFuncDerivs[0][ETA] =  0.25 * (xi - xi2 - 2.*xiEta + 2.*xi2Eta);
    _shapeFuncDerivs[1][ETA] = -0.25 * (xi + xi2 - 2.*xiEta - 2.*xi2Eta);
    _shapeFuncDerivs[2][ETA] =  0.25 * (xi + xi2 + 2.*xiEta + 2.*xi2Eta);
    _shapeFuncDerivs[3][ETA] = -0.25 * (xi - xi2 + 2.*xiEta - 2.*xi2Eta);
    _shapeFuncDerivs[4][ETA] = -0.5 * (1. - xi2 - 2.*eta + 2.*xi2Eta);
    _shapeFuncDerivs[5][ETA] =  0.5 * (-2.*xiEta - 2.*xi2Eta);
    _shapeFuncDerivs[6][ETA] =  0.5 * (1. - xi2 + 2.*eta - 2.*xi2Eta);
    _shapeFuncDerivs[7][ETA] = -0.5 * (-2.*xiEta + 2.*xi2Eta);
    _shapeFuncDerivs[8][ETA] =  2.*xi2Eta - 2.*eta;

    // evaluate Jacobian
    pointJacob(KSI,XX) = _shapeFuncDerivs[0][KSI]*(*nodes[0])[XX];
    pointJacob(ETA,XX) = _shapeFuncDerivs[0][ETA]*(*nodes[0])[XX];

    pointJacob(KSI,YY) = _shapeFuncDerivs[0][KSI]*(*nodes[0])[YY];
    pointJacob(ETA,YY) = _shapeFuncDerivs[0][ETA]*(*nodes[0])[YY];
    for (CFuint in = 1; in < 9; ++in)
    {
      pointJacob(KSI,XX) += _shapeFuncDerivs[in][KSI]*(*nodes[in])[XX];
      pointJacob(ETA,XX) += _shapeFuncDerivs[in][ETA]*(*nodes[in])[XX];

      pointJacob(KSI,YY) += _shapeFuncDerivs[in][KSI]*(*nodes[in])[YY];
      pointJacob(ETA,YY) += _shapeFuncDerivs[in][ETA]*(*nodes[in])[YY];
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void LagrangeShapeFunctionQuadP2::computeJacobianPlus1D(
         const std::vector<Framework::Node*>& nodes,
         const std::vector<RealVector>& mappedCoord,
               std::vector<RealMatrix>& jacob)
{
  throw Common::NotImplementedException (FromHere(),getName() + "::computeJacobianPlus1D()");
}

//////////////////////////////////////////////////////////////////////////////

void LagrangeShapeFunctionQuadP2::computeFaceJacobianDeterminant(
          const std::vector<RealVector>& mappedCoord,
          const std::vector<Framework::Node*>& nodes,
          const Framework::IntegratorPattern& pattern,
                std::vector<RealVector>& faceJacobian)
{
  throw Common::NotImplementedException (FromHere(),getName() + "::computeFaceJacobianDeterminant()");

}

//////////////////////////////////////////////////////////////////////////////

RealVector LagrangeShapeFunctionQuadP2::computeMappedCoordinates(const RealVector& coord,
                                                                 const std::vector<Framework::Node*>& nodes)
{

  throw Common::NotImplementedException (FromHere(),getName() + "::computeMappedCoordinates()");

}
//////////////////////////////////////////////////////////////////////

    } // namespace ShapeFunctions

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
