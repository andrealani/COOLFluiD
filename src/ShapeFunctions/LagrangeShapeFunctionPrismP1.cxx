// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "ShapeFunctions/LagrangeShapeFunctionPrismP1.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



   namespace ShapeFunctions {

//////////////////////////////////////////////////////////////////////////////

CFuint LagrangeShapeFunctionPrismP1::_interpolatorID = 0;
RealVector LagrangeShapeFunctionPrismP1::_vec1 = RealVector(DIM_3D);
RealVector LagrangeShapeFunctionPrismP1::_vec2 = RealVector(DIM_3D);
RealVector LagrangeShapeFunctionPrismP1::_vec3 = RealVector(DIM_3D);
std::vector<RealVector> LagrangeShapeFunctionPrismP1::_normals = std::vector<RealVector>(5);

//////////////////////////////////////////////////////////////////////////////

void LagrangeShapeFunctionPrismP1::computeJacobianDeterminant(
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

   const CFreal cte = 0.125 * (x5*y1*z0 - x1*y2*z0 - x4*y2*z0 + 2*x1*y3*z0 + 2*x4*y3*z0 - 2*x5*y3*z0 + x5*y4*z0 - x1*y5*z0 - x4*y5*z0 - x5*y0*z1
   + x0*y2*z1 - 2*x0*y3*z1 + x5*y3*z1 + x0*y5*z1 + x1*y0*z2 + x4*y0*z2 - x0*y1*z2 + 2*x0*y3*z2 - x1*y3*z2 - x4*y3*z2 - x0*y4*z2 -
   2*x1*y0*z3 - 2*x4*y0*z3 + 2*x5*y0*z3 + 2*x0*y1*z3 - x5*y1*z3 - 2*x0*y2*z3 + x1*y2*z3 + x4*y2*z3 + 2*x0*y4*z3 - x5*y4*z3 -
   2*x0*y5*z3 + x1*y5*z3 + x4*y5*z3 - x5*y0*z4 + x0*y2*z4 - 2*x0*y3*z4 + x5*y3*z4 + x0*y5*z4 +
   x2*(y4*z0 - y0*z1 + y1*(z0 - z3) + 2*y0*z3 - y4*z3 - y0*z4 + y3*(-2*z0 + z1 + z4)) + x1*y0*z5 + x4*y0*z5 - x0*y1*z5 +
   2*x0*y3*z5 - x1*y3*z5 - x4*y3*z5 - x0*y4*z5 + x3*
      (-2*y4*z0 + 2*y5*z0 + 2*y0*z1 - y5*z1 - 2*y0*z2 + y4*z2 + 2*y0*z4 - y5*z4 - y2*(-2*z0 + z1 + z4) - 2*y0*z5 + y4*z5 +
      y1*(-2*z0 + z2 + z5)));

   const CFreal aXi = 0.25 * (
   -(x1*y3*z0) + x2*y3*z0 + x5*y3*z0 + x1*y4*z0 - x2*y4*z0 - x5*y4*z0 + x0*y3*z1 - x2*y3*z1 - x5*y3*z1 - x0*y4*z1 + x2*y4*z1 +
   x5*y4*z1 - x0*y3*z2 + x1*y3*z2 + x0*y4*z2 - x1*y4*z2 + x1*y0*z3 - x2*y0*z3 - x5*y0*z3 - x0*y1*z3 + x2*y1*z3 + x5*y1*z3 +
   x0*y2*z3 - x1*y2*z3 - x0*y4*z3 + x1*y4*z3 + x0*y5*z3 - x1*y5*z3 - x1*y0*z4 + x2*y0*z4 + x5*y0*z4 + x0*y1*z4 - x2*y1*z4 -
   x5*y1*z4 - x0*y2*z4 + x1*y2*z4 + x0*y3*z4 - x1*y3*z4 - x0*y5*z4 + x1*y5*z4 - x0*y3*z5 + x1*y3*z5 + x0*y4*z5 - x1*y4*z5 +
   x3*(-(y2*z0) + y4*z0 - y5*z0 - y0*z1 + y2*z1 - y4*z1 + y5*z1 + y0*z2 - y0*z4 + y1*(z0 - z2 + z4 - z5) + y0*z5) +
   x4*(y2*z0 - y3*z0 + y5*z0 + y0*z1 - y2*z1 + y3*z1 - y5*z1 - y0*z2 + y0*z3 - y0*z5 + y1*(-z0 + z2 - z3 + z5)));

   const CFreal aEta =  0.25 * (
   -(x1*y3*z0) + x2*y3*z0 - x4*y3*z0 + x1*y5*z0 - x2*y5*z0 + x4*y5*z0 + x0*y3*z1 - x2*y3*z1 - x0*y5*z1 + x2*y5*z1 - x0*y3*z2 +
   x1*y3*z2 + x4*y3*z2 + x0*y5*z2 - x1*y5*z2 - x4*y5*z2 + x1*y0*z3 - x2*y0*z3 + x4*y0*z3 - x0*y1*z3 + x2*y1*z3 + x0*y2*z3 -
   x1*y2*z3 - x4*y2*z3 - x0*y4*z3 + x2*y4*z3 + x0*y5*z3 - x2*y5*z3 + x0*y3*z4 - x2*y3*z4 - x0*y5*z4 + x2*y5*z4 +
   x5*(y3*z0 - y4*z0 + y0*z1 - y0*z2 - y3*z2 + y4*z2 + y1*(-z0 + z2) - y0*z3 + y2*(z0 - z1 + z3 - z4) + y0*z4) - x1*y0*z5 +
   x2*y0*z5 - x4*y0*z5 + x0*y1*z5 - x2*y1*z5 - x0*y2*z5 + x1*y2*z5 + x4*y2*z5 - x0*y3*z5 + x2*y3*z5 + x0*y4*z5 - x2*y4*z5 +
   x3*(y4*z0 - y5*z0 - y0*z1 + y1*(z0 - z2) + y0*z2 - y4*z2 + y5*z2 - y0*z4 + y2*(-z0 + z1 + z4 - z5) + y0*z5));

   const CFreal aZeta = 0.25 * (
   x1*y2*z0 - x1*y3*z0 + x4*y3*z0 - x5*y3*z0 + x5*y4*z0 - x4*y5*z0 - x0*y2*z1 + x0*y3*z1 - x1*y0*z2 + x0*y1*z2 - x0*y3*z2 +
   x1*y3*z2 + x1*y0*z3 - x4*y0*z3 + x5*y0*z3 - x0*y1*z3 + x0*y2*z3 - x1*y2*z3 + x0*y4*z3 - x5*y4*z3 - x0*y5*z3 + x4*y5*z3 +
   x2*(y3*z0 + y0*z1 - y3*z1 - y0*z3 + y1*(-z0 + z3)) - x5*y0*z4 - x0*y3*z4 + x5*y3*z4 + x0*y5*z4 + x4*y0*z5 + x0*y3*z5 -
   x4*y3*z5 - x0*y4*z5 + x3*(y1*z0 - y2*z0 - y4*z0 + y5*z0 - y0*z1 + y2*z1 + y0*z2 - y1*z2 + y0*z4 - y5*z4 - y0*z5 + y4*z5));

   const CFreal aXiZeta = 0.25 * (
   x1*y3*z0 - x2*y3*z0 + x5*y3*z0 - x1*y4*z0 + x2*y4*z0 - x5*y4*z0 - x0*y3*z1 + x2*y3*z1 - x5*y3*z1 + x0*y4*z1 - x2*y4*z1
   + x5*y4*z1 + x0*y3*z2 - x1*y3*z2 - x0*y4*z2 + x1*y4*z2 - x1*y0*z3 + x2*y0*z3 - x5*y0*z3 + x0*y1*z3 - x2*y1*z3 + x5*y1*z3 -
   x0*y2*z3 + x1*y2*z3 - x0*y4*z3 + x1*y4*z3 + x0*y5*z3 - x1*y5*z3 + x1*y0*z4 - x2*y0*z4 + x5*y0*z4 - x0*y1*z4 + x2*y1*z4 -
   x5*y1*z4 + x0*y2*z4 - x1*y2*z4 + x0*y3*z4 - x1*y3*z4 - x0*y5*z4 + x1*y5*z4 - x0*y3*z5 + x1*y3*z5 + x0*y4*z5 - x1*y4*z5 +
   x3*(y2*z0 + y4*z0 - y5*z0 + y0*z1 - y2*z1 - y4*z1 + y5*z1 - y0*z2 - y0*z4 + y1*(-z0 + z2 + z4 - z5) + y0*z5) +
   x4*(-(y2*z0) - y3*z0 + y5*z0 - y0*z1 + y2*z1 + y3*z1 - y5*z1 + y0*z2 + y0*z3 - y0*z5 + y1*(z0 - z2 - z3 + z5)));

   const CFreal aEtaZeta = 0.25 * (
   x1*y3*z0 - x2*y3*z0 - x4*y3*z0 - x1*y5*z0 + x2*y5*z0 + x4*y5*z0 - x0*y3*z1 + x2*y3*z1 + x0*y5*z1 - x2*y5*z1 + x0*y3*z2
   - x1*y3*z2 + x4*y3*z2 - x0*y5*z2 + x1*y5*z2 - x4*y5*z2 - x1*y0*z3 + x2*y0*z3 + x4*y0*z3 + x0*y1*z3 - x2*y1*z3 - x0*y2*z3 +
   x1*y2*z3 - x4*y2*z3 - x0*y4*z3 + x2*y4*z3 + x0*y5*z3 - x2*y5*z3 + x0*y3*z4 - x2*y3*z4 - x0*y5*z4 + x2*y5*z4 +
   x5*(y3*z0 - y4*z0 - y0*z1 + y1*(z0 - z2) + y0*z2 - y3*z2 + y4*z2 - y0*z3 + y2*(-z0 + z1 + z3 - z4) + y0*z4) + x1*y0*z5 -
   x2*y0*z5 - x4*y0*z5 - x0*y1*z5 + x2*y1*z5 + x0*y2*z5 - x1*y2*z5 + x4*y2*z5 - x0*y3*z5 + x2*y3*z5 + x0*y4*z5 - x2*y4*z5 +
   x3*(y4*z0 - y5*z0 + y0*z1 - y0*z2 - y4*z2 + y5*z2 + y1*(-z0 + z2) - y0*z4 + y2*(z0 - z1 + z4 - z5) + y0*z5));

   const CFreal aZeta2 = (-(x1*y2*z0) + x4*y2*z0 + x1*y5*z0 - x4*y5*z0 + x0*y2*z1 - x3*y2*z1 - x0*y5*z1 + x3*y5*z1 + x1*y0*z2 - x4*y0*z2 - x0*y1*z2 +
   x3*y1*z2 - x1*y3*z2 + x4*y3*z2 + x0*y4*z2 - x3*y4*z2 + x1*y2*z3 - x4*y2*z3 - x1*y5*z3 + x4*y5*z3 +
   x2*(y1*(z0 - z3) + y4*(-z0 + z3) - (y0 - y3)*(z1 - z4)) + x5*(y4*(z0 - z3) + y1*(-z0 + z3) + (y0 - y3)*(z1 - z4)) -
   x0*y2*z4 + x3*y2*z4 + x0*y5*z4 - x3*y5*z4 - x1*y0*z5 + x4*y0*z5 + x0*y1*z5 - x3*y1*z5 + x1*y3*z5 - x4*y3*z5 - x0*y4*z5 +
   x3*y4*z5);

   const CFuint nbQdPts = mappedCoord.size();
   for (CFuint ip = 0; ip < nbQdPts; ++ip) {
      cf_assert(mappedCoord[ip].size() == DIM_3D);

      const CFreal xi   = mappedCoord[ip][KSI];
      const CFreal eta  = mappedCoord[ip][ETA];
      const CFreal zeta = mappedCoord[ip][ZTA];

      detJacobian[ip] = cte + aXi * xi + aEta * eta + (aXiZeta * xi + aEtaZeta * eta + aZeta  + aZeta2 * zeta) * zeta;
   }
}

//////////////////////////////////////////////////////////////////////////////

void LagrangeShapeFunctionPrismP1::computeFaceJacobianDeterminant(
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

   const CFreal coeff = 0.5;

   // Face 0
   const CFuint iFace0 = 0;
   const CFreal jacob0 =
         Framework::FaceJacobiansDeterminant::compute3DFaceJacobDet(
         (x1-x0),
         (y1-y0),
         (z1-z0),
         (x2-x0),
         (y2-y0),
         (z2-z0));

   for(CFuint ip = 0; ip < pattern.nbPts(iFace0); ++ip) {
     faceJacobian[iFace0][ip] = jacob0;
   }

   // Face 1
   const CFuint iFace1 = 1;
   const CFreal jacob1 =
         Framework::FaceJacobiansDeterminant::compute3DFaceJacobDet(
         (x4-x3),
         (y4-y3),
         (z4-z3),
         (x5-x3),
         (y5-y3),
         (z5-z3));

   for(CFuint ip = 0; ip < pattern.nbPts(iFace1); ++ip) {
     faceJacobian[iFace1][ip] = jacob1;
   }

   // Face 2
   const CFuint iFace2 = 2;
   for(CFuint ip = 0; ip < pattern.nbPts(iFace2); ++ip) {

     cf_assert(mappedCoord[ip].size() == DIM_3D);

     const CFreal xsi = mappedCoord[ip][KSI];
     // unused // const CFreal eta = mappedCoord[ip][ETA];
     const CFreal zta = mappedCoord[ip][ZTA];

     const CFreal jacob2 =
          Framework::FaceJacobiansDeterminant::compute3DFaceJacobDet(
            coeff * ((1.+zta)*(x4-x3)+(1.-zta)*(x1-x0)),
            coeff * ((1.+zta)*(y4-y3)+(1.-zta)*(y1-y0)),
            coeff * ((1.+zta)*(z4-z3)+(1.-zta)*(z1-z0)),
            coeff * (  (xsi)*(x4-x1) +(1.-xsi)*(x3-x0)),
            coeff * (  (xsi)*(y4-y1) +(1.-xsi)*(y3-y0)),
            coeff * (  (xsi)*(z4-z1) +(1.-xsi)*(z3-z0))
            );

     faceJacobian[iFace2][ip] = jacob2;
   }

   // Face 3
   const CFuint iFace3 = 3;
   const CFreal b = std::sqrt(2.0) / 2.0;
   const CFreal a = 1.0 / (2.0 * std::sqrt(2.0));
   for(CFuint ip = 0; ip < pattern.nbPts(iFace3); ++ip) {

     cf_assert(mappedCoord[ip].size() == DIM_3D);

     const CFreal s1 = 2.0 * mappedCoord[ip][KSI] * b;
     const CFreal s2 =       mappedCoord[ip][ZTA];

     const CFreal jacob3 =
          Framework::FaceJacobiansDeterminant::compute3DFaceJacobDet(
          (x1 - s2*x1 + (s2-1.0)*x2 + (1.0 + s2)*(x4 - x5))*a,
          (y1 - s2*y1 + (s2-1.0)*y2 + (1.0 + s2)*(y4 - y5))*a,
          (z1 - s2*z1 + (s2-1.0)*z2 + (1.0 + s2)*(z4 - z5))*a,
          (b*s1*(-x1 + x2 + x4 - x5) + x5 - x2) * coeff,
          (b*s1*(-y1 + y2 + y4 - y5) + y5 - y2) * coeff,
          (b*s1*(-z1 + z2 + z4 - z5) + z5 - z2) * coeff
            );

     faceJacobian[iFace3][ip] = jacob3;
   }

   // Face 4
   const CFuint iFace4 = 4;
   for(CFuint ip = 0; ip < pattern.nbPts(iFace4); ++ip) {

     cf_assert(mappedCoord[ip].size() == DIM_3D);

     // unused // const CFreal xsi = mappedCoord[ip][KSI];
     const CFreal eta = mappedCoord[ip][ETA];
     const CFreal zta = mappedCoord[ip][ZTA];

     const CFreal jacob4 =
          Framework::FaceJacobiansDeterminant::compute3DFaceJacobDet(
            coeff * ((1.+zta)*(x5-x3)+(1.-zta)*(x2-x0)),
            coeff * ((1.+zta)*(y5-y3)+(1.-zta)*(y2-y0)),
            coeff * ((1.+zta)*(z5-z3)+(1.-zta)*(z2-z0)),
            coeff * (  (eta)*(x5-x2) +(1.-eta)*(x3-x0)),
            coeff * (  (eta)*(y5-y2) +(1.-eta)*(y3-y0)),
            coeff * (  (eta)*(z5-z2) +(1.-eta)*(z3-z0))
            );

     faceJacobian[iFace4][ip] = jacob4;
   }

}

//////////////////////////////////////////////////////////////////////////////

    } // namespace ShapeFunctions



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
