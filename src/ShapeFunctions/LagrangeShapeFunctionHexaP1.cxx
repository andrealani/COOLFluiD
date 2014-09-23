// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "ShapeFunctions/LagrangeShapeFunctionHexaP1.hh"
#include "MathTools/MathFunctions.hh"
#include "MathTools/MatrixInverterT.hh"


//////////////////////////////////////////////////////////////////////////////

using namespace std;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

   namespace ShapeFunctions {

//////////////////////////////////////////////////////////////////////////////

CFuint LagrangeShapeFunctionHexaP1::_interpolatorID = 0;
RealVector LagrangeShapeFunctionHexaP1::_vec1 = RealVector(DIM_3D);
RealVector LagrangeShapeFunctionHexaP1::_vec2 = RealVector(DIM_3D);
RealVector LagrangeShapeFunctionHexaP1::_vec3 = RealVector(DIM_3D);
std::vector<RealVector> LagrangeShapeFunctionHexaP1::_normals = std::vector<RealVector>(6);
std::vector<RealVector> LagrangeShapeFunctionHexaP1::_gradShapFunc = std::vector<RealVector>(3,RealVector(8));

//////////////////////////////////////////////////////////////////////////////

void LagrangeShapeFunctionHexaP1::computeFaceJacobianDeterminant(
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

    const CFreal coeff = 0.250;

    const CFuint iFace0 = 0;
    for(CFuint ip = 0; ip < pattern.nbPts(iFace0); ++ip) {

      cf_assert(mappedCoord[ip].size() == DIM_3D);

      const CFreal xsi = mappedCoord[ip][KSI];
      const CFreal eta = mappedCoord[ip][ETA];
      // unused // const CFreal zta = mappedCoord[ip][ZTA];

      const CFreal jacob0 =
          Framework::FaceJacobiansDeterminant::compute3DFaceJacobDet(
            coeff * ((1.-eta)*(x1-x0)+(1.+eta)*(x2-x3)),
            coeff * ((1.-eta)*(y1-y0)+(1.+eta)*(y2-y3)),
            coeff * ((1.-eta)*(z1-z0)+(1.+eta)*(z2-z3)),
            coeff * ((1.-xsi)*(x3-x0)+(1.+xsi)*(x2-x1)),
            coeff * ((1.-xsi)*(y3-y0)+(1.+xsi)*(y2-y1)),
            coeff * ((1.-xsi)*(z3-z0)+(1.+xsi)*(z2-z1)));

      faceJacobian[iFace0][ip] = jacob0;
    }

    const CFuint iFace1 = 1;
    for(CFuint ip = 0; ip < pattern.nbPts(iFace1); ++ip) {

      cf_assert(mappedCoord[ip].size() == DIM_3D);

      const CFreal xsi = mappedCoord[ip][KSI];
      const CFreal eta = mappedCoord[ip][ETA];
      // unused // const CFreal zta = mappedCoord[ip][ZTA];

      const CFreal jacob1 =
          Framework::FaceJacobiansDeterminant::compute3DFaceJacobDet(
            coeff * ((1.-eta)*(x5-x4)+(1.+eta)*(x6-x7)),
            coeff * ((1.-eta)*(y5-y4)+(1.+eta)*(y6-y7)),
            coeff * ((1.-eta)*(z5-z4)+(1.+eta)*(z6-z7)),
            coeff * ((1.-xsi)*(x7-x4)+(1.+xsi)*(x6-x5)),
            coeff * ((1.-xsi)*(y7-y4)+(1.+xsi)*(y6-y5)),
            coeff * ((1.-xsi)*(z7-z4)+(1.+xsi)*(z6-z5)));

      faceJacobian[iFace1][ip] = jacob1;
    }

    const CFuint iFace2 = 2;
    for(CFuint ip = 0; ip < pattern.nbPts(iFace2); ++ip) {

      cf_assert(mappedCoord[ip].size() == DIM_3D);

      const CFreal xsi = mappedCoord[ip][KSI];
      // unused // const CFreal eta = mappedCoord[ip][ETA];
      const CFreal zta = mappedCoord[ip][ZTA];

      const CFreal jacob2 =
          Framework::FaceJacobiansDeterminant::compute3DFaceJacobDet(
            coeff * ((1.-zta)*(x1-x0)+(1.+zta)*(x5-x4)),
            coeff * ((1.-zta)*(y1-y0)+(1.+zta)*(y5-y4)),
            coeff * ((1.-zta)*(z1-z0)+(1.+zta)*(z5-z4)),
            coeff * ((1.-xsi)*(x4-x0)+(1.+xsi)*(x5-x1)),
            coeff * ((1.-xsi)*(y4-y0)+(1.+xsi)*(y5-y1)),
            coeff * ((1.-xsi)*(z4-z0)+(1.+xsi)*(z5-z1)));

      faceJacobian[iFace2][ip] = jacob2;
    }

    const CFuint iFace3 = 3;
    for(CFuint ip = 0; ip < pattern.nbPts(iFace3); ++ip) {

      cf_assert(mappedCoord[ip].size() == DIM_3D);

      // unused // const CFreal xsi = mappedCoord[ip][KSI];
      const CFreal eta = mappedCoord[ip][ETA];
      const CFreal zta = mappedCoord[ip][ZTA];

      const CFreal jacob3 =
          Framework::FaceJacobiansDeterminant::compute3DFaceJacobDet(
            coeff * ((1.-zta)*(x2-x1)+(1.+zta)*(x6-x5)),
            coeff * ((1.-zta)*(y2-y1)+(1.+zta)*(y6-y5)),
            coeff * ((1.-zta)*(z2-z1)+(1.+zta)*(z6-z5)),
            coeff * ((1.-eta)*(x5-x1)+(1.+eta)*(x6-x2)),
            coeff * ((1.-eta)*(y5-y1)+(1.+eta)*(y6-y2)),
            coeff * ((1.-eta)*(z5-z1)+(1.+eta)*(z6-z2)));

      faceJacobian[iFace3][ip] = jacob3;
    }

    const CFuint iFace4 = 4;
    for(CFuint ip = 0; ip < pattern.nbPts(iFace4); ++ip) {

      cf_assert(mappedCoord[ip].size() == DIM_3D);

      const CFreal xsi = mappedCoord[ip][KSI];
      // unused // const CFreal eta = mappedCoord[ip][ETA];
      const CFreal zta = mappedCoord[ip][ZTA];

      const CFreal jacob4 =
          Framework::FaceJacobiansDeterminant::compute3DFaceJacobDet(
            coeff * ((1.-zta)*(x2-x3)+(1.+zta)*(x6-x7)),
            coeff * ((1.-zta)*(y2-y3)+(1.+zta)*(y6-y7)),
            coeff * ((1.-zta)*(z2-z3)+(1.+zta)*(z6-z7)),
            coeff * ((1.-xsi)*(x7-x3)+(1.+xsi)*(x6-x2)),
            coeff * ((1.-xsi)*(y7-y3)+(1.+xsi)*(y6-y2)),
            coeff * ((1.-xsi)*(z7-z3)+(1.+xsi)*(z6-z2)));

      faceJacobian[iFace4][ip] = jacob4;
    }

    const CFuint iFace5 = 5;
    for(CFuint ip = 0; ip < pattern.nbPts(iFace5); ++ip) {

      cf_assert(mappedCoord[ip].size() == DIM_3D);

      // unused // const CFreal xsi = mappedCoord[ip][KSI];
      const CFreal eta = mappedCoord[ip][ETA];
      const CFreal zta = mappedCoord[ip][ZTA];

      const CFreal jacob5 =
          Framework::FaceJacobiansDeterminant::compute3DFaceJacobDet(
            coeff * ((1.-zta)*(x3-x0)+(1.+zta)*(x7-x4)),
            coeff * ((1.-zta)*(y3-y0)+(1.+zta)*(y7-y4)),
            coeff * ((1.-zta)*(z3-z0)+(1.+zta)*(z7-z4)),
            coeff * ((1.-eta)*(x4-x0)+(1.+eta)*(x7-x3)),
            coeff * ((1.-eta)*(y4-y0)+(1.+eta)*(y7-y3)),
            coeff * ((1.-eta)*(z4-z0)+(1.+eta)*(z7-z3)));

      faceJacobian[iFace5][ip] = jacob5;
    }
  }

//////////////////////////////////////////////////////////////////////////////

void LagrangeShapeFunctionHexaP1::computeGradientStates(
         const std::vector<RealMatrix>& jacob,
         const std::vector<RealVector>& mappedCoord,
               std::vector<RealMatrix>& grad)
{
  static RealMatrix invJ(3,3);
  MathTools::MatrixInverterT<3> inverter;

  for (CFuint ip = 0; ip < mappedCoord.size(); ++ip) {

    inverter.invert(jacob[ip], invJ);

    RealMatrix& lgrad = grad[ip];

    // coefficient is 0.25 which square will be 0.125
    const CFreal mxi   = 0.25 * (1.0 - mappedCoord[ip][KSI]);
    const CFreal meta  = 0.25 * (1.0 - mappedCoord[ip][ETA]);
    const CFreal mzeta = 0.25 * (1.0 - mappedCoord[ip][ZTA]);
    const CFreal pxi   = 0.25 * (1.0 + mappedCoord[ip][KSI]);
    const CFreal peta  = 0.25 * (1.0 + mappedCoord[ip][ETA]);
    const CFreal pzeta = 0.25 * (1.0 + mappedCoord[ip][ZTA]);

    const CFreal JXX = invJ(0,0);
    const CFreal JXY = invJ(0,1);
    const CFreal JXZ = invJ(0,2);
    const CFreal JYX = invJ(1,0);
    const CFreal JYY = invJ(1,1);
    const CFreal JYZ = invJ(1,2);
    const CFreal JZX = invJ(2,0);
    const CFreal JZY = invJ(2,1);
    const CFreal JZZ = invJ(2,2);

    lgrad(0,XX) = - JXX * meta * mzeta - JXY * mxi * mzeta - JXZ * mxi * meta;
    lgrad(0,YY) = - JYX * meta * mzeta - JYY * mxi * mzeta - JYZ * mxi * meta;
    lgrad(0,ZZ) = - JZX * meta * mzeta - JZY * mxi * mzeta - JZZ * mxi * meta;

    lgrad(1,XX) = + JXX * meta * mzeta - JXY * pxi * mzeta - JXZ * pxi * meta;
    lgrad(1,YY) = + JYX * meta * mzeta - JYY * pxi * mzeta - JYZ * pxi * meta;
    lgrad(1,ZZ) = + JZX * meta * mzeta - JZY * pxi * mzeta - JZZ * pxi * meta;

    lgrad(2,XX) = + JXX * peta * mzeta + JXY * pxi * mzeta - JXZ * pxi * peta;
    lgrad(2,YY) = + JYX * peta * mzeta + JYY * pxi * mzeta - JYZ * pxi * peta;
    lgrad(2,ZZ) = + JZX * peta * mzeta + JZY * pxi * mzeta - JZZ * pxi * peta;

    lgrad(3,XX) = - JXX * peta * mzeta + JXY * mxi * mzeta - JXZ * mxi * peta;
    lgrad(3,YY) = - JYX * peta * mzeta + JYY * mxi * mzeta - JYZ * mxi * peta;
    lgrad(3,ZZ) = - JZX * peta * mzeta + JZY * mxi * mzeta - JZZ * mxi * peta;

    lgrad(4,XX) = - JXX * meta * pzeta - JXY * mxi * pzeta + JXZ * mxi * meta;
    lgrad(4,YY) = - JYX * meta * pzeta - JYY * mxi * pzeta + JYZ * mxi * meta;
    lgrad(4,ZZ) = - JZX * meta * pzeta - JZY * mxi * pzeta + JZZ * mxi * meta;

    lgrad(5,XX) = + JXX * meta * pzeta - JXY * pxi * pzeta + JXZ * pxi * meta;
    lgrad(5,YY) = + JYX * meta * pzeta - JYY * pxi * pzeta + JYZ * pxi * meta;
    lgrad(5,ZZ) = + JZX * meta * pzeta - JZY * pxi * pzeta + JZZ * pxi * meta;

    lgrad(6,XX) = + JXX * peta * pzeta + JXY * pxi * pzeta + JXZ * pxi * peta;
    lgrad(6,YY) = + JYX * peta * pzeta + JYY * pxi * pzeta + JYZ * pxi * peta;
    lgrad(6,ZZ) = + JZX * peta * pzeta + JZY * pxi * pzeta + JZZ * pxi * peta;

    lgrad(7,XX) = - JXX * peta * pzeta + JXY * mxi * pzeta + JXZ * mxi * peta;
    lgrad(7,YY) = - JYX * peta * pzeta + JYY * mxi * pzeta + JYZ * mxi * peta;
    lgrad(7,ZZ) = - JZX * peta * pzeta + JZY * mxi * pzeta + JZZ * mxi * peta;
  }
}

//////////////////////////////////////////////////////////////////////////////

void LagrangeShapeFunctionHexaP1::computeMappedCoordPlaneNormal(const std::vector<CFuint>& planeIdx,
                                                                const std::vector<RealVector>& mappedCoord,
                                                                const std::vector<Framework::Node*>& nodes,
                                                                std::vector<RealVector>& normal)
{
  for (CFuint ip = 0; ip < mappedCoord.size(); ++ip)
  {
    RealVector& pointNormal = normal[ip];

    const CFreal xi =  mappedCoord[ip][KSI];
    const CFreal eta = mappedCoord[ip][ETA];
    const CFreal zta = mappedCoord[ip][ZTA];

    const CFreal a1 = 1. + xi;
    const CFreal a2 = 1. - xi;
    const CFreal b1 = 1. + eta;
    const CFreal b2 = 1. - eta;
    const CFreal c1 = 1. + zta;
    const CFreal c2 = 1. - zta;

    cf_assert(planeIdx[ip] == 0 || planeIdx[ip] == 1 || planeIdx[ip] == 2);
    if (planeIdx[ip] == 0)
    {
      _gradShapFunc[ETA][0] = -a2*c2;
      _gradShapFunc[ETA][1] = -a1*c2;
      _gradShapFunc[ETA][2] =  a1*c2;
      _gradShapFunc[ETA][3] =  a2*c2;
      _gradShapFunc[ETA][4] = -a2*c1;
      _gradShapFunc[ETA][5] = -a1*c1;
      _gradShapFunc[ETA][6] =  a1*c1;
      _gradShapFunc[ETA][7] =  a2*c1;

      _gradShapFunc[ZTA][0] = -a2*b2;
      _gradShapFunc[ZTA][1] = -a1*b2;
      _gradShapFunc[ZTA][2] = -a1*b1;
      _gradShapFunc[ZTA][3] = -a2*b1;
      _gradShapFunc[ZTA][4] =  b2*a2;
      _gradShapFunc[ZTA][5] =  b2*a1;
      _gradShapFunc[ZTA][6] =  b1*a1;
      _gradShapFunc[ZTA][7] =  b1*a2;

      _vec1 = _gradShapFunc[ETA][0]*(*nodes[0]);
      _vec2 = _gradShapFunc[ZTA][0]*(*nodes[0]);
      for (CFuint in = 1; in < 8; ++in)
      {
        _vec1 += _gradShapFunc[ETA][in]*(*nodes[in]);
        _vec2 += _gradShapFunc[ZTA][in]*(*nodes[in]);
      }
    }
    else if (planeIdx[ip] == 1)
    {
      _gradShapFunc[ZTA][0] = -a2*b2;
      _gradShapFunc[ZTA][1] = -a1*b2;
      _gradShapFunc[ZTA][2] = -a1*b1;
      _gradShapFunc[ZTA][3] = -a2*b1;
      _gradShapFunc[ZTA][4] =  b2*a2;
      _gradShapFunc[ZTA][5] =  b2*a1;
      _gradShapFunc[ZTA][6] =  b1*a1;
      _gradShapFunc[ZTA][7] =  b1*a2;

      _gradShapFunc[KSI][0] = -b2*c2;
      _gradShapFunc[KSI][1] =  b2*c2;
      _gradShapFunc[KSI][2] =  b1*c2;
      _gradShapFunc[KSI][3] = -b1*c2;
      _gradShapFunc[KSI][4] = -b2*c1;
      _gradShapFunc[KSI][5] =  b2*c1;
      _gradShapFunc[KSI][6] =  b1*c1;
      _gradShapFunc[KSI][7] = -b1*c1;

      _vec1 = _gradShapFunc[ZTA][0]*(*nodes[0]);
      _vec2 = _gradShapFunc[KSI][0]*(*nodes[0]);
      for (CFuint in = 1; in < 8; ++in)
      {
        _vec1 += _gradShapFunc[ZTA][in]*(*nodes[in]);
        _vec2 += _gradShapFunc[KSI][in]*(*nodes[in]);
      }
    }
    else
    {
      _gradShapFunc[KSI][0] = -b2*c2;
      _gradShapFunc[KSI][1] =  b2*c2;
      _gradShapFunc[KSI][2] =  b1*c2;
      _gradShapFunc[KSI][3] = -b1*c2;
      _gradShapFunc[KSI][4] = -b2*c1;
      _gradShapFunc[KSI][5] =  b2*c1;
      _gradShapFunc[KSI][6] =  b1*c1;
      _gradShapFunc[KSI][7] = -b1*c1;

      _gradShapFunc[ETA][0] = -a2*c2;
      _gradShapFunc[ETA][1] = -a1*c2;
      _gradShapFunc[ETA][2] =  a1*c2;
      _gradShapFunc[ETA][3] =  a2*c2;
      _gradShapFunc[ETA][4] = -a2*c1;
      _gradShapFunc[ETA][5] = -a1*c1;
      _gradShapFunc[ETA][6] =  a1*c1;
      _gradShapFunc[ETA][7] =  a2*c1;

      _vec1 = _gradShapFunc[KSI][0]*(*nodes[0]);
      _vec2 = _gradShapFunc[ETA][0]*(*nodes[0]);
      for (CFuint in = 1; in < 8; ++in)
      {
        _vec1 += _gradShapFunc[KSI][in]*(*nodes[in]);
        _vec2 += _gradShapFunc[ETA][in]*(*nodes[in]);
      }
    }

    // compute normal
    MathTools::MathFunctions::crossProd(_vec1,_vec2,pointNormal);
    pointNormal *= 0.015625;
  }
}

//////////////////////////////////////////////////////////////////////////////

void LagrangeShapeFunctionHexaP1::computeJacobian(const std::vector<Framework::Node*>& nodes,
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

    for (CFuint ip = 0; ip < mappedCoord.size(); ++ip) {

      RealMatrix& pointJacob = jacob[ip];

      const CFreal xi =  mappedCoord[ip][KSI];
      const CFreal eta = mappedCoord[ip][ETA];
      const CFreal zeta = mappedCoord[ip][ZTA];

      const CFreal a1 = 1. + xi;
      const CFreal a2 = 1. - xi;
      const CFreal b1 = 1. + eta;
      const CFreal b2 = 1. - eta;
      const CFreal c1 = 1. + zeta;
      const CFreal c2 = 1. - zeta;

      const CFreal dN0dxi = -b2*c2;
      const CFreal dN1dxi = b2*c2;
      const CFreal dN2dxi = b1*c2;
      const CFreal dN3dxi = -b1*c2;
      const CFreal dN4dxi = -b2*c1;
      const CFreal dN5dxi = b2*c1;
      const CFreal dN6dxi = b1*c1;
      const CFreal dN7dxi = -b1*c1;

      const CFreal dN0deta = -a2*c2;
      const CFreal dN1deta = -a1*c2;
      const CFreal dN2deta = a1*c2;
      const CFreal dN3deta = a2*c2;
      const CFreal dN4deta = -a2*c1;
      const CFreal dN5deta = -a1*c1;
      const CFreal dN6deta = a1*c1;
      const CFreal dN7deta = a2*c1;

      const CFreal dN0dzeta = -a2*b2;
      const CFreal dN1dzeta = -a1*b2;
      const CFreal dN2dzeta = -a1*b1;
      const CFreal dN3dzeta = -a2*b1;
      const CFreal dN4dzeta = b2*a2;
      const CFreal dN5dzeta = b2*a1;
      const CFreal dN6dzeta = b1*a1;
      const CFreal dN7dzeta = b1*a2;

      pointJacob(KSI,XX) = x0*dN0dxi   + x1*dN1dxi   + x2*dN2dxi   + x3*dN3dxi   + x4*dN4dxi   + x5*dN5dxi   + x6*dN6dxi   + x7*dN7dxi;
      pointJacob(ETA,XX) = x0*dN0deta  + x1*dN1deta  + x2*dN2deta  + x3*dN3deta  + x4*dN4deta  + x5*dN5deta  + x6*dN6deta  + x7*dN7deta;
      pointJacob(ZTA,XX) = x0*dN0dzeta + x1*dN1dzeta + x2*dN2dzeta + x3*dN3dzeta + x4*dN4dzeta + x5*dN5dzeta + x6*dN6dzeta + x7*dN7dzeta;

      pointJacob(KSI,YY) = y0*dN0dxi   + y1*dN1dxi   + y2*dN2dxi   + y3*dN3dxi   + y4*dN4dxi   + y5*dN5dxi   + y6*dN6dxi   + y7*dN7dxi;
      pointJacob(ETA,YY) = y0*dN0deta  + y1*dN1deta  + y2*dN2deta  + y3*dN3deta  + y4*dN4deta  + y5*dN5deta  + y6*dN6deta  + y7*dN7deta;
      pointJacob(ZTA,YY) = y0*dN0dzeta + y1*dN1dzeta + y2*dN2dzeta + y3*dN3dzeta + y4*dN4dzeta + y5*dN5dzeta + y6*dN6dzeta + y7*dN7dzeta;

      pointJacob(KSI,ZZ) = z0*dN0dxi   + z1*dN1dxi   + z2*dN2dxi   + z3*dN3dxi   + z4*dN4dxi   + z5*dN5dxi   + z6*dN6dxi   + z7*dN7dxi;
      pointJacob(ETA,ZZ) = z0*dN0deta  + z1*dN1deta  + z2*dN2deta  + z3*dN3deta  + z4*dN4deta  + z5*dN5deta  + z6*dN6deta  + z7*dN7deta;
      pointJacob(ZTA,ZZ) = z0*dN0dzeta + z1*dN1dzeta + z2*dN2dzeta + z3*dN3dzeta + z4*dN4dzeta + z5*dN5dzeta + z6*dN6dzeta + z7*dN7dzeta;

      pointJacob *= 0.125;
    }


/*  static RealMatrix coord(getNbNodes(),DIM_3D);
  static RealMatrix grads(DIM_3D,getNbNodes());

  // put coordinates in a matrix nodes*dim
  for (CFuint n = 0; n < getNbNodes(); ++n) {
    coord(n,XX) = (*nodes[n])[XX];
    coord(n,YY) = (*nodes[n])[YY];
    coord(n,ZZ) = (*nodes[n])[ZZ];
  }

  // compute jacobian at each quadrature point
  for (CFuint ip = 0; ip < mappedCoord.size(); ++ip) {

    // coefficient is 0.25 which square will be 0.125
    const CFreal mxi   = 0.25 * (1.0 - mappedCoord[ip][KSI]);
    const CFreal meta  = 0.25 * (1.0 - mappedCoord[ip][ETA]);
    const CFreal mzeta = 0.25 * (1.0 - mappedCoord[ip][ZTA]);
    const CFreal pxi   = 0.25 * (1.0 + mappedCoord[ip][KSI]);
    const CFreal peta  = 0.25 * (1.0 + mappedCoord[ip][ETA]);
    const CFreal pzeta = 0.25 * (1.0 + mappedCoord[ip][ZTA]);

    grads(XX,0) = - meta * mzeta;
    grads(YY,0) = - mxi  * mzeta;
    grads(ZZ,0) = - mxi  * meta ;

    grads(XX,1) = + meta * mzeta;
    grads(YY,1) = - pxi  * mzeta;
    grads(ZZ,1) = - pxi  * meta ;

    grads(XX,2) = + peta * mzeta;
    grads(YY,2) = + pxi  * mzeta;
    grads(ZZ,2) = - pxi  * peta ;

    grads(XX,3) = - peta * mzeta;
    grads(YY,3) = + mxi  * mzeta;
    grads(ZZ,3) = - mxi  * peta ;

    grads(XX,4) = - meta * pzeta;
    grads(YY,4) = - mxi  * pzeta;
    grads(ZZ,4) = + mxi  * meta ;

    grads(XX,5) = + meta * pzeta;
    grads(YY,5) = - pxi  * pzeta;
    grads(ZZ,5) = + pxi  * meta ;

    grads(XX,6) = + peta * pzeta;
    grads(YY,6) = + pxi  * pzeta;
    grads(ZZ,6) = + pxi  * peta ;

    grads(XX,7) = - peta * pzeta;
    grads(YY,7) = + mxi  * pzeta;
    grads(ZZ,7) = + mxi  * peta ;

    // results in a matrix DIM_3D*DIM_3D
    jacob[ip] = grads * coord;
  }*/
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace ShapeFunctions



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
