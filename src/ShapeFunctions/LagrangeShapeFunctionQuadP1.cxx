// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "ShapeFunctions/LagrangeShapeFunctionQuadP1.hh"
#include "MathTools/MathFunctions.hh"
#include "MathTools/MatrixInverterT.hh"


//////////////////////////////////////////////////////////////////////////////

using namespace COOLFluiD::MathTools;

namespace COOLFluiD {

   namespace ShapeFunctions {

//////////////////////////////////////////////////////////////////////////////

CFuint LagrangeShapeFunctionQuadP1::_interpolatorID = 0;
RealVector LagrangeShapeFunctionQuadP1::_vec10 = RealVector(DIM_2D);
RealVector LagrangeShapeFunctionQuadP1::_vec20 = RealVector(DIM_2D);
RealVector LagrangeShapeFunctionQuadP1::_vecP0 = RealVector(DIM_2D);
RealVector LagrangeShapeFunctionQuadP1::_shapeFunc = RealVector(4);
RealVector LagrangeShapeFunctionQuadP1::_2Dvec1 = RealVector(DIM_2D);
RealVector LagrangeShapeFunctionQuadP1::_2Dvec2 = RealVector(DIM_2D);
RealVector LagrangeShapeFunctionQuadP1::_3Dvec1 = RealVector(DIM_3D);
RealVector LagrangeShapeFunctionQuadP1::_3Dvec2 = RealVector(DIM_3D);
RealVector LagrangeShapeFunctionQuadP1::_mappedCoord  = RealVector(DIM_2D);
RealVector LagrangeShapeFunctionQuadP1::_tempVector2D = RealVector(DIM_2D);
RealVector LagrangeShapeFunctionQuadP1::_tempVector3D = RealVector(DIM_3D);
std::vector<RealVector> LagrangeShapeFunctionQuadP1::_normals =
 std::vector<RealVector>(4,LagrangeShapeFunctionQuadP1::_tempVector2D);

//////////////////////////////////////////////////////////////////////////////

void LagrangeShapeFunctionQuadP1::computeGradientStates(
         const std::vector<RealMatrix>& jacob,
         const std::vector<RealVector>& mappedCoord,
               std::vector<RealMatrix>& grad)
{
  static RealMatrix invJ(2,2);
  MathTools::MatrixInverterT<2> inverter;

  for (CFuint ip = 0; ip < mappedCoord.size(); ++ip) {

    inverter.invert(jacob[ip], invJ);

    RealMatrix& lgrad = grad[ip];

    const CFreal mxi   = 0.25 * (1.0 - mappedCoord[ip][KSI]);
    const CFreal meta  = 0.25 * (1.0 - mappedCoord[ip][ETA]);
    const CFreal pxi   = 0.25 * (1.0 + mappedCoord[ip][KSI]);
    const CFreal peta  = 0.25 * (1.0 + mappedCoord[ip][ETA]);

    const CFreal JXX = invJ(0,0);
    const CFreal JXY = invJ(0,1);
    const CFreal JYX = invJ(1,0);
    const CFreal JYY = invJ(1,1);

    lgrad(0,XX) = - JXX * meta - JXY * mxi;
    lgrad(0,YY) = - JYX * meta - JYY * mxi;

    lgrad(1,XX) = + JXX * meta - JXY * pxi;
    lgrad(1,YY) = + JYX * meta - JYY * pxi;

    lgrad(2,XX) = + JXX * peta + JXY * pxi;
    lgrad(2,YY) = + JYX * peta + JYY * pxi;

    lgrad(3,XX) = - JXX * peta + JXY * mxi;
    lgrad(3,YY) = - JYX * peta + JYY * mxi;
  }
}

//////////////////////////////////////////////////////////////////////////////

void LagrangeShapeFunctionQuadP1::computeFaceJacobDetVectorAtMappedCoords(const std::vector<RealVector>& mappedCoord,
                                                                          const std::vector<Framework::Node*>& nodes,
                                                                          std::vector<RealVector>& normal)
{
  for (CFuint ip = 0; ip < mappedCoord.size(); ++ip)
  {
    RealVector& pointNormal = normal[ip];

    const CFreal xi =  mappedCoord[ip][KSI];
    const CFreal eta = mappedCoord[ip][ETA];
    const CFreal dN0dksi = -(1. - eta);
    const CFreal dN1dksi =  (1. - eta);
    const CFreal dN2dksi =  (1. + eta);
    const CFreal dN3dksi = -(1. + eta);
    const CFreal dN0deta = -(1. - xi);
    const CFreal dN1deta = -(1. + xi);
    const CFreal dN2deta =  (1. + xi);
    const CFreal dN3deta =  (1. - xi);

    // compute shape function gradient
    _3Dvec1 = (*nodes[0])*dN0dksi + (*nodes[1])*dN1dksi + (*nodes[2])*dN2dksi + (*nodes[3])*dN3dksi;
    _3Dvec2 = (*nodes[0])*dN0deta + (*nodes[1])*dN1deta + (*nodes[2])*dN2deta + (*nodes[3])*dN3deta;

    // compute face jacobian vector
    MathTools::MathFunctions::crossProd(_3Dvec1,_3Dvec2,pointNormal);
    pointNormal *= 0.0625;
  }
}

//////////////////////////////////////////////////////////////////////////////

void LagrangeShapeFunctionQuadP1::computeMappedCoordPlaneNormal(const std::vector<CFuint>& planeIdx,
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

  for (CFuint ip = 0; ip < mappedCoord.size(); ++ip)
  {
    RealVector& pointNormal = normal[ip];

    const CFreal xi =  mappedCoord[ip][KSI];
    const CFreal eta = mappedCoord[ip][ETA];

    cf_assert(planeIdx[ip] == 0 || planeIdx[ip] == 1);
    if (planeIdx[ip] == 0)
    {
      const CFreal dN0deta = -(1. - xi);
      const CFreal dN1deta = -(1. + xi);
      const CFreal dN2deta =  (1. + xi);
      const CFreal dN3deta =  (1. - xi);

      pointNormal[XX] = +0.25 * (y0*dN0deta + y1*dN1deta + y2*dN2deta + y3*dN3deta);
      pointNormal[YY] = -0.25 * (x0*dN0deta + x1*dN1deta + x2*dN2deta + x3*dN3deta);
    }
    else
    {
      const CFreal dN0dxi = -(1. - eta);
      const CFreal dN1dxi =  (1. - eta);
      const CFreal dN2dxi =  (1. + eta);
      const CFreal dN3dxi = -(1. + eta);

      pointNormal[XX] = -0.25 * (y0*dN0dxi  + y1*dN1dxi  + y2*dN2dxi  + y3*dN3dxi);
      pointNormal[YY] = +0.25 * (x0*dN0dxi  + x1*dN1dxi  + x2*dN2dxi  + x3*dN3dxi);
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void LagrangeShapeFunctionQuadP1::computeJacobian(const std::vector<Framework::Node*>& nodes,
                                                  const std::vector<RealVector>& mappedCoord,
                                                  std::vector<RealMatrix>& jacob)
{
  const CFreal x0 = (*nodes[0])[XX];
  const CFreal y0 = (*nodes[0])[YY];

  const CFreal x1 = (*nodes[1])[XX];
  const CFreal y1 = (*nodes[1])[YY];

  const CFreal x2 = (*nodes[2])[XX];
  const CFreal y2 = (*nodes[2])[YY];

  const CFreal x3 = (*nodes[3])[XX];
  const CFreal y3 = (*nodes[3])[YY];

  for (CFuint ip = 0; ip < mappedCoord.size(); ++ip) {

    RealMatrix& pointJacob = jacob[ip];

    const CFreal xi =  mappedCoord[ip][KSI];
    const CFreal eta = mappedCoord[ip][ETA];

    const CFreal a1 = 1. + xi;
    const CFreal a2 = 1. - xi;
    const CFreal b1 = 1. + eta;
    const CFreal b2 = 1. - eta;

    const CFreal dN0dxi = -b2;
    const CFreal dN1dxi =  b2;
    const CFreal dN2dxi =  b1;
    const CFreal dN3dxi = -b1;

    const CFreal dN0deta = -a2;
    const CFreal dN1deta = -a1;
    const CFreal dN2deta =  a1;
    const CFreal dN3deta =  a2;


    pointJacob(KSI,XX) = x0*dN0dxi  + x1*dN1dxi  + x2*dN2dxi  + x3*dN3dxi;
    pointJacob(ETA,XX) = x0*dN0deta + x1*dN1deta + x2*dN2deta + x3*dN3deta;

    pointJacob(KSI,YY) = y0*dN0dxi  + y1*dN1dxi  + y2*dN2dxi  + y3*dN3dxi;
    pointJacob(ETA,YY) = y0*dN0deta + y1*dN1deta + y2*dN2deta + y3*dN3deta;

    pointJacob *= 0.25;
  }
}

//////////////////////////////////////////////////////////////////////////////

void LagrangeShapeFunctionQuadP1::computeJacobianPlus1D(
         const std::vector<Framework::Node*>& nodes,
         const std::vector<RealVector>& mappedCoord,
               std::vector<RealMatrix>& jacob)
  {
    const CFreal x0x1 = (*nodes[0])[XX] - (*nodes[1])[XX];
    const CFreal x2x3 = (*nodes[2])[XX] - (*nodes[3])[XX];
    const CFreal x2x1 = (*nodes[2])[XX] - (*nodes[1])[XX];
    const CFreal x0x3 = (*nodes[0])[XX] - (*nodes[3])[XX];

    const CFreal y0y1 = (*nodes[0])[YY] - (*nodes[1])[YY];
    const CFreal y2y3 = (*nodes[2])[YY] - (*nodes[3])[YY];
    const CFreal y2y1 = (*nodes[2])[YY] - (*nodes[1])[YY];
    const CFreal y0y3 = (*nodes[0])[YY] - (*nodes[3])[YY];

    const CFreal z0z1 = (*nodes[0])[ZZ] - (*nodes[1])[ZZ];
    const CFreal z2z3 = (*nodes[2])[ZZ] - (*nodes[3])[ZZ];
    const CFreal z2z1 = (*nodes[2])[ZZ] - (*nodes[1])[ZZ];
    const CFreal z0z3 = (*nodes[0])[ZZ] - (*nodes[3])[ZZ];

    RealVector& cross = _tempVector3D;
    RealVector& jx = _3Dvec1;
    RealVector& jy = _3Dvec2;

    for (CFuint ip = 0; ip < mappedCoord.size(); ++ip) {

      const CFreal xi   = mappedCoord[ip][KSI];
      const CFreal eta  = mappedCoord[ip][ETA];

      const CFreal alfa = 0.25 * (eta - 1);
      const CFreal beta = 0.25 * (eta + 1);
      const CFreal gama = 0.25 * (xi  - 1);
      const CFreal delt = 0.25 * (xi  + 1);

      jx[XX] = alfa * x0x1 + beta * x2x3;
      jx[YY] = alfa * y0y1 + beta * y2y3;
      jx[ZZ] = alfa * z0z1 + beta * z2z3;

      jy[XX] = gama * x2x1 + delt * x0x3;
      jy[YY] = gama * y2y1 + delt * y0y3;
      jy[ZZ] = gama * z2z1 + delt * z0z3;

      jacob[ip].setRow(jx,XX);
      jacob[ip].setRow(jy,YY);

      MathTools::MathFunctions::crossProd(jx,jy,cross);

      cross.normalize();

      jacob[ip].setRow(cross,ZZ);
    }
 }

//////////////////////////////////////////////////////////////////////////////

void LagrangeShapeFunctionQuadP1::computeFaceJacobianDeterminant(
          const std::vector<RealVector>& mappedCoord,
          const std::vector<Framework::Node*>& nodes,
          const Framework::IntegratorPattern& pattern,
                std::vector<RealVector>& faceJacobian)
{
  cf_assert(pattern.nbSteps() == getNbFaces());
  cf_assert(faceJacobian.size() >= getNbFaces());

  const CFreal x0 = (*nodes[0])[XX];
  const CFreal y0 = (*nodes[0])[YY];

  const CFreal x1 = (*nodes[1])[XX];
  const CFreal y1 = (*nodes[1])[YY];

  const CFreal x2 = (*nodes[2])[XX];
  const CFreal y2 = (*nodes[2])[YY];

  const CFreal x3 = (*nodes[3])[XX];
  const CFreal y3 = (*nodes[3])[YY];

  cf_assert(pattern.nbSteps() == 4);

  const CFreal jacob0 =
    Framework::FaceJacobiansDeterminant::compute2DFaceJacobDet(0.5*(x1-x0),0.5*(y1-y0));

  const CFreal jacob1 =
    Framework::FaceJacobiansDeterminant::compute2DFaceJacobDet(0.5*(x2-x1),0.5*(y2-y1));

  const CFreal jacob2 =
    Framework::FaceJacobiansDeterminant::compute2DFaceJacobDet(0.5*(x2-x3),0.5*(y2-y3));

  const CFreal jacob3 =
    Framework::FaceJacobiansDeterminant::compute2DFaceJacobDet(0.5*(x3-x0),0.5*(y3-y0));

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

//////////////////////////////////////////////////////////////////////////////

RealVector LagrangeShapeFunctionQuadP1::computeMappedCoordinates(const RealVector& coord,
                                                                 const std::vector<Framework::Node*>& nodes)
{
  cf_assert(coord.size() == DIM_2D);
  cf_assert(nodes.size() == 4);
  cf_assert(nodes[0]->size() == DIM_2D);
  cf_assert(nodes[1]->size() == DIM_2D);
  cf_assert(nodes[2]->size() == DIM_2D);
  cf_assert(nodes[3]->size() == DIM_2D);

  RealVector& xA = (*nodes[0]);
  RealVector& xB = (*nodes[1]);
  RealVector& xC = (*nodes[2]);
  RealVector& xD = (*nodes[3]);

  ///First split the quad into 2 triangles
  ///Check in which triangle the point is
  bool isInsideABC = false;

  // Face 20
  computeFaceLineNormal(2,nodes,2,0);
  _tempVector2D = (*nodes[0]) - coord;
  CFreal inside = 0.;

  //scalar product of AP with outward pointing normal
  for(CFuint iDim=0; iDim < DIM_2D; ++iDim){
    inside += _normals[2][iDim]*_tempVector2D[iDim];
  }
  if(inside > -MathTools::MathConsts::CFrealEps()) isInsideABC = true;

  ///Interpolate in the triangle
  //if isInside triangle ABC
  if(isInsideABC){
    //look in triangle BCA
    // Node 0 = B
    // Node 1 = C
    // Node 2 = A
    _vec10 = xC-xB;
    _vec20 = xA-xB;
    _vecP0 = coord-xB;

    const CFreal overDet = 1./(_vec10[XX]*_vec20[YY] - _vec10[YY]*_vec20[XX]);

    const CFreal ksiTriag = overDet*( _vec20[YY]*_vecP0[XX] - _vec20[XX]*_vecP0[YY]);
    const CFreal etaTriag = overDet*(-_vec10[YY]*_vecP0[XX] + _vec10[XX]*_vecP0[YY]);

    _2Dvec1[KSI] = ((1.-etaTriag) * 2.) - 1.;
    _2Dvec1[ETA] = (ksiTriag * 2.) - 1.;
  }
  //else isInside triangle ACD
  else{
    //look in triangle DAC
    // Node 0 = D
    // Node 1 = A
    // Node 2 = C
    _vec10 = xA-xD;
    _vec20 = xC-xD;
    _vecP0 = coord-xD;

    const CFreal overDet = 1./(_vec10[XX]*_vec20[YY] - _vec10[YY]*_vec20[XX]);

    const CFreal ksiTriag = overDet*( _vec20[YY]*_vecP0[XX] - _vec20[XX]*_vecP0[YY]);
    const CFreal etaTriag = overDet*(-_vec10[YY]*_vecP0[XX] + _vec10[XX]*_vecP0[YY]);

    _2Dvec1[KSI] = (etaTriag * 2.) - 1.;
    _2Dvec1[ETA] = ((1.-ksiTriag) * 2.) - 1.;
  }

  computeCoordFromMappedCoord(_2Dvec1,nodes,_tempVector2D);
  _tempVector2D -= coord;
  const CFreal dist1 = _tempVector2D.norm2();

  ///Then split the quad into 2 other triangles
  ///Check in which triangle the point is
  bool isInsideABD = false;

  // Face 31
  computeFaceLineNormal(2,nodes,3,1);
  _tempVector2D = (*nodes[3]) - coord;
  inside = 0.;

  //scalar product of DP with outward pointing normal
  for(CFuint iDim=0; iDim < DIM_2D; ++iDim){
    inside += _normals[2][iDim]*_tempVector2D[iDim];
  }
  if(inside > -MathTools::MathConsts::CFrealEps()) isInsideABD = true;

  ///Interpolate in the triangle
  //if isInside triangle ABD
  if(isInsideABD){
    //look in triangle ABD
    // Node 0 = A
    // Node 1 = B
    // Node 2 = D
    _vec10 = xB-xA;
    _vec20 = xD-xA;
    _vecP0 = coord-xA;

    const CFreal overDet = 1./(_vec10[XX]*_vec20[YY] - _vec10[YY]*_vec20[XX]);

    const CFreal ksiTriag = overDet*( _vec20[YY]*_vecP0[XX] - _vec20[XX]*_vecP0[YY]);
    const CFreal etaTriag = overDet*(-_vec10[YY]*_vecP0[XX] + _vec10[XX]*_vecP0[YY]);

    _2Dvec2[KSI] = (ksiTriag * 2.) - 1.;
    _2Dvec2[ETA] = (etaTriag * 2.) - 1.;
  }
  //else isInside triangle BCD
  else{
    //look in triangle CDB
    // Node 0 = C
    // Node 1 = D
    // Node 2 = B
    _vec10 = xD-xC;
    _vec20 = xB-xC;
    _vecP0 = coord-xC;

    const CFreal overDet = 1./(_vec10[XX]*_vec20[YY] - _vec10[YY]*_vec20[XX]);

    const CFreal ksiTriag = overDet*( _vec20[YY]*_vecP0[XX] - _vec20[XX]*_vecP0[YY]);
    const CFreal etaTriag = overDet*(-_vec10[YY]*_vecP0[XX] + _vec10[XX]*_vecP0[YY]);

    _2Dvec2[KSI] = ((1.-ksiTriag) * 2.) - 1.;
    _2Dvec2[ETA] = ((1.-etaTriag) * 2.) - 1.;
  }

  computeCoordFromMappedCoord(_2Dvec2,nodes,_tempVector2D);
  _tempVector2D -= coord;
  const CFreal dist2 = _tempVector2D.norm2();

  if(dist1 < dist2){
   _mappedCoord = _2Dvec1;
//   CFout << "Approx: " << dist1 << "\n";
  }
  else{
   _mappedCoord = _2Dvec2;
//   CFout << "Approx: " << dist2 << "\n";
  }

  ///We have a first linear approximation of the mapped coordinates
  ///We build a small quad around this approximated point and
  ///we interpolate once more...to get a more accurate approximation

//   CFreal delta = 0.01;
//
// ///@todo do this ?????
// {
//   RealVector _mappedCoordA(DIM_2D);
//
//   RealVector xAn(DIM_2D);
//   RealVector xBn(DIM_2D);
//   RealVector xCn(DIM_2D);
//   RealVector xDn(DIM_2D);
//
//   _mappedCoordA[ETA] = std::max(-1.,std::min(_mappedCoord[ETA] - delta,1.));
//   _mappedCoordA[KSI] = std::max(-1.,std::min(_mappedCoord[KSI] - delta,1.));
//
//   //create a new quad around the point
//   computeCoordFromMappedCoord(_mappedCoordA,nodes,xAn);
//   RealVector _mappedCoordB = _mappedCoordA;
//   _mappedCoordB[KSI] += 2*delta;
//   _mappedCoordB[KSI] = std::max(-1.,std::min(_mappedCoordB[KSI],1.));
//   computeCoordFromMappedCoord(_mappedCoordB,nodes,xBn);
//   RealVector _mappedCoordC = _mappedCoordB;
//   _mappedCoordC[ETA] += 2*delta;
//   _mappedCoordC[ETA] = std::max(-1.,std::min(_mappedCoordC[ETA],1.));
//   computeCoordFromMappedCoord(_mappedCoordC,nodes,xCn);
//   RealVector _mappedCoordD = _mappedCoordC;
//   _mappedCoordD[KSI] -= 2*delta;
//   _mappedCoordD[KSI] = std::max(-1.,std::min(_mappedCoordD[KSI],1.));
//   computeCoordFromMappedCoord(_mappedCoordD,nodes,xDn);
//
//   ///First split the quad into 2 triangles
//   ///Check in which triangle the point is
//   bool isInside012 = false;
//
//   // Face 20
//   _normals[2][XX] = xAn[YY] - xCn[YY];
//   _normals[2][YY] = xCn[XX] - xAn[XX];
//   RealVector m_vec1 = xAn - coord;
//   CFreal inside = 0.;
//
//   //scalar product of AP with outward pointing normal
//   for(CFuint iDim=0; iDim < DIM_2D; ++iDim){
//     inside += _normals[2][iDim]*m_vec1[iDim];
//   }
//   if(inside > -MathTools::MathConsts::CFrealEps()) isInside012 = true;
//
//   ///Interpolate in the triangle
//   //if isInside triangle ABC
//   if(isInsideABC){
//     //look in triangle BCA
//     // Node 0 = B
//     // Node 1 = C
//     // Node 2 = A
//     RealVector m_v10 = xCn-xBn;
//     RealVector m_v20 = xAn-xBn;
//     RealVector m_vP0 = coord-xBn;
//
//     CFreal overDet = 1./(m_v10[XX]*m_v20[YY] - m_v10[YY]*m_v20[XX]);
//
//     CFreal ksiTriag = overDet*( m_v20[YY]*m_vP0[XX] - m_v20[XX]*m_vP0[YY]);
//     CFreal etaTriag = overDet*(-m_v10[YY]*m_vP0[XX] + m_v10[XX]*m_vP0[YY]);
//
//     CFreal deltaKsi = _mappedCoordB[KSI] - _mappedCoordA[KSI];
//     CFreal deltaEta = _mappedCoordC[ETA] - _mappedCoordB[ETA];
//     _mappedCoord[KSI] += (((1.-etaTriag) * 2.) - 1.) * deltaKsi;
//     _mappedCoord[ETA] += ((ksiTriag * 2.) - 1.) * deltaEta;
//     _mappedCoord[ETA] = std::max(-1.,std::min(_mappedCoord[ETA],1.));
//     _mappedCoord[KSI] = std::max(-1.,std::min(_mappedCoord[KSI],1.));
//
//   }
//   //else isInside triangle ACD
//   else{
//     //look in triangle DAC
//     // Node 0 = D
//     // Node 1 = A
//     // Node 2 = C
//     RealVector m_v10 = xAn-xDn;
//     RealVector m_v20 = xCn-xDn;
//     RealVector m_vP0 = coord-xDn;
//
//     CFreal overDet = 1./(m_v10[XX]*m_v20[YY] - m_v10[YY]*m_v20[XX]);
//
//     CFreal ksiTriag = overDet*( m_v20[YY]*m_vP0[XX] - m_v20[XX]*m_vP0[YY]);
//     CFreal etaTriag = overDet*(-m_v10[YY]*m_vP0[XX] + m_v10[XX]*m_vP0[YY]);
//
//     CFreal deltaKsi = _mappedCoordC[KSI] - _mappedCoordD[KSI];
//     CFreal deltaEta = _mappedCoordD[ETA] - _mappedCoordA[ETA];
//
//     _mappedCoord[KSI] += ((etaTriag * 2.) - 1.) * deltaKsi;
//     _mappedCoord[ETA] += (((1.-ksiTriag) * 2.) - 1.) * deltaEta;
//     _mappedCoord[ETA] = std::max(-1.,std::min(_mappedCoord[ETA],1.));
//     _mappedCoord[KSI] = std::max(-1.,std::min(_mappedCoord[KSI],1.));
//   }
// }
//
//   computeCoordFromMappedCoord(_mappedCoord,nodes,_tempVector2D);
// _tempVector2D -= coord;
//  CFout << "Second Approx: " << _tempVector2D.norm2()  << "\n";

  return _mappedCoord;

}
//////////////////////////////////////////////////////////////////////

    } // namespace ShapeFunctions

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
