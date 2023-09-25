// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "ShapeFunctions/LagrangeShapeFunctionPrismP2.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



   namespace ShapeFunctions {

//////////////////////////////////////////////////////////////////////////////

CFuint LagrangeShapeFunctionPrismP2::_interpolatorID = 0;
RealVector LagrangeShapeFunctionPrismP2::_vec1 = RealVector(DIM_3D);
RealVector LagrangeShapeFunctionPrismP2::_vec2 = RealVector(DIM_3D);
RealVector LagrangeShapeFunctionPrismP2::_vec3 = RealVector(DIM_3D);
std::vector<RealVector> LagrangeShapeFunctionPrismP2::_normals = std::vector<RealVector>(5);
std::vector<RealVector> LagrangeShapeFunctionPrismP2::_gradShapFunc = std::vector<RealVector>(3,RealVector(18));

//////////////////////////////////////////////////////////////////////////////

void LagrangeShapeFunctionPrismP2::computeShapeFunction(const RealVector& mappedCoord,
                                         RealVector& shapeFunc)
{
   const CFreal xi   = mappedCoord[KSI];
   const CFreal eta  = mappedCoord[ETA];
   const CFreal zta = mappedCoord[ZTA];


   const CFreal lbd   = 1.0 - xi - eta;

   const CFreal L1 = 0.5*zta*(zta - 1.);
   const CFreal L2 = 0.5*zta*(zta + 1.);
   const CFreal L3 = (1.0 - zta*zta);
   

   shapeFunc[0]  = lbd*(2.*lbd - 1.)*L1;
   shapeFunc[1]  = xi * (2. * xi - 1.)*L1;
   shapeFunc[2]  = eta * (2. * eta - 1.)*L1;

   shapeFunc[6]  = 4.*xi*lbd*L1;
   shapeFunc[7]  = 4.*xi*eta*L1;
   shapeFunc[8]  = 4.*eta*lbd*L1;

   shapeFunc[3]  = lbd*(2.*lbd - 1.)*L2;
   shapeFunc[4] = xi * (2. * xi - 1.)*L2;
   shapeFunc[5] = eta * (2. * eta - 1.)*L2;

   shapeFunc[15] = 4.*xi*lbd*L2;
   shapeFunc[16] = 4.*xi*eta*L2;
   shapeFunc[17] = 4.*eta*lbd*L2;

   shapeFunc[9]  = lbd*(2.*lbd - 1.)*L3;
   shapeFunc[11]  = xi * (2. * xi - 1.)*L3;
   shapeFunc[13]  = eta * (2. * eta - 1.)*L3;

   shapeFunc[10] = 4.*xi*lbd*L3;
   shapeFunc[12] = 4.*xi*eta*L3;
   shapeFunc[14] = 4.*eta*lbd*L3;

 

  
}

//////////////////////////////////////////////////////////////////////////////

void LagrangeShapeFunctionPrismP2::computeJacobian(
      const std::vector<Framework::Node*>& nodes,
      const std::vector<RealVector>& mappedCoord,
            std::vector<RealMatrix>& jacob)
 {
   for (CFuint ip = 0; ip < mappedCoord.size(); ++ip) 
   {
      RealMatrix& pointJacob = jacob[ip];
       
      const CFreal xi =  mappedCoord[ip][KSI];
      const CFreal eta = mappedCoord[ip][ETA];
      const CFreal zta = mappedCoord[ip][ZTA];

      const CFreal lbd   = 1.0 - xi - eta;

      const CFreal L1 = 0.5*zta*(zta - 1.);
      const CFreal L2 = 0.5*zta*(zta + 1.);
      const CFreal L3 = (1.0 - zta*zta);

      const CFreal dL1 = zta - 0.5;
      const CFreal dL2 = zta + 0.5;
      const CFreal dL3 = - 2.*zta;


      _gradShapFunc[KSI][0]  = (-3. + 4.*eta + 4.*xi)*L1;
      _gradShapFunc[KSI][1]  = (4.*xi - 1.)*L1;
      _gradShapFunc[KSI][2]  = 0.;

      _gradShapFunc[KSI][6]  = (4. - 8.*xi - 4.*eta)*L1;
      _gradShapFunc[KSI][7]  = (4.*eta)*L1;
      _gradShapFunc[KSI][8]  = (-4.*eta)*L1;

      _gradShapFunc[KSI][3]  = (-3. + 4.*eta + 4.*xi)*L2;
      _gradShapFunc[KSI][4] = (4.*xi - 1.)*L2;
      _gradShapFunc[KSI][5] = 0.;

      _gradShapFunc[KSI][15] = (4. - 8.*xi - 4.*eta)*L2;
      _gradShapFunc[KSI][16] = (4.*eta)*L2;
      _gradShapFunc[KSI][17] = (-4.*eta)*L2;

      _gradShapFunc[KSI][9]  = (-3. + 4.*eta + 4.*xi)*L3;
      _gradShapFunc[KSI][11]  = (4.*xi - 1.)*L3;
      _gradShapFunc[KSI][13]  = 0.;

      _gradShapFunc[KSI][10] = (4. - 8.*xi - 4.*eta)*L3;
      _gradShapFunc[KSI][12] = (4.*eta)*L3;
      _gradShapFunc[KSI][14] = (-4.*eta)*L3;




      _gradShapFunc[ETA][0]  = (-3. + 4.*eta + 4.*xi)*L1;
      _gradShapFunc[ETA][1]  = 0.;
      _gradShapFunc[ETA][2]  = (4.*eta - 1.)*L1;;

      _gradShapFunc[ETA][6]  = (-4.*xi)*L1;
      _gradShapFunc[ETA][7]  = (4.*xi)*L1;
      _gradShapFunc[ETA][8]  = (4. - 4.*xi - 8.*eta)*L1;

      _gradShapFunc[ETA][3]  = (-3. + 4.*eta + 4.*xi)*L2;
      _gradShapFunc[ETA][4] = 0.;
      _gradShapFunc[ETA][5] = (4.*eta - 1.)*L2;

      _gradShapFunc[ETA][15] = (-4.*xi)*L2;
      _gradShapFunc[ETA][16] = (4.*xi)*L2;
      _gradShapFunc[ETA][17] = (4. - 4.*xi - 8.*eta)*L2;

      _gradShapFunc[ETA][9]  = (-3. + 4.*eta + 4.*xi)*L3;
      _gradShapFunc[ETA][11]  = 0.;
      _gradShapFunc[ETA][13]  = (4.*eta - 1.)*L3;

      _gradShapFunc[ETA][10] = (-4.*xi)*L3;
      _gradShapFunc[ETA][12] = (4.*xi)*L3;
      _gradShapFunc[ETA][14] = (4. - 4.*xi - 8.*eta)*L3;



      _gradShapFunc[ZTA][0]  = lbd*(2.*lbd - 1.)*dL1;
      _gradShapFunc[ZTA][1]  = xi * (2. * xi - 1.)*dL1;
      _gradShapFunc[ZTA][2]  = eta * (2. * eta - 1.)*dL1;

      _gradShapFunc[ZTA][6]  = 4.*xi*lbd*dL1;
      _gradShapFunc[ZTA][7]  = 4.*xi*eta*dL1;
      _gradShapFunc[ZTA][8]  = 4.*eta*lbd*dL1;

      _gradShapFunc[ZTA][3]  = lbd*(2.*lbd - 1.)*dL2;
      _gradShapFunc[ZTA][4] = xi * (2. * xi - 1.)*dL2;
      _gradShapFunc[ZTA][5] = eta * (2. * eta - 1.)*dL2;

      _gradShapFunc[ZTA][15] = 4.*xi*lbd*dL2;
      _gradShapFunc[ZTA][16] = 4.*xi*eta*dL2;
      _gradShapFunc[ZTA][17] = 4.*eta*lbd*dL2;

      _gradShapFunc[ZTA][9]  = lbd*(2.*lbd - 1.)*dL3;
      _gradShapFunc[ZTA][11]  = xi * (2. * xi - 1.)*dL3;
      _gradShapFunc[ZTA][13]  = eta * (2. * eta - 1.)*dL3;

      _gradShapFunc[ZTA][10] = 4.*xi*lbd*dL3;
      _gradShapFunc[ZTA][12] = 4.*xi*eta*dL3;
      _gradShapFunc[ZTA][14] = 4.*eta*lbd*dL3;


     // compute Jacobian matrix
     const CFreal xx = (*nodes[0])[XX];
     const CFreal yy = (*nodes[0])[YY];
     const CFreal zz = (*nodes[0])[ZZ];
     const CFreal dNdksi = _gradShapFunc[KSI][0];
     const CFreal dNdeta = _gradShapFunc[ETA][0];
     const CFreal dNdzta = _gradShapFunc[ZTA][0]/2.;

     pointJacob(KSI,XX) = xx*dNdksi;
     pointJacob(ETA,XX) = xx*dNdeta;
     pointJacob(ZTA,XX) = xx*dNdzta;

     pointJacob(KSI,YY) = yy*dNdksi;
     pointJacob(ETA,YY) = yy*dNdeta;
     pointJacob(ZTA,YY) = yy*dNdzta;

     pointJacob(KSI,ZZ) = zz*dNdksi;
     pointJacob(ETA,ZZ) = zz*dNdeta;
     pointJacob(ZTA,ZZ) = zz*dNdzta;

     for (CFuint in = 1; in < 18; ++in)
     {

       const CFreal xx = (*nodes[in])[XX];
       const CFreal yy = (*nodes[in])[YY];
       const CFreal zz = (*nodes[in])[ZZ];
       const CFreal dNdksi = _gradShapFunc[KSI][in];
       const CFreal dNdeta = _gradShapFunc[ETA][in];
       const CFreal dNdzta = _gradShapFunc[ZTA][in]/2.;

       pointJacob(KSI,XX) += xx*dNdksi;
       pointJacob(ETA,XX) += xx*dNdeta;
       pointJacob(ZTA,XX) += xx*dNdzta;

       pointJacob(KSI,YY) += yy*dNdksi;
       pointJacob(ETA,YY) += yy*dNdeta;
       pointJacob(ZTA,YY) += yy*dNdzta;

       pointJacob(KSI,ZZ) += zz*dNdksi;
       pointJacob(ETA,ZZ) += zz*dNdeta;
       pointJacob(ZTA,ZZ) += zz*dNdzta;
     }
   }
 }
//////////////////////////////////////////////////////////////////////////////

void LagrangeShapeFunctionPrismP2::computeMappedCoordPlaneNormal(const std::vector<CFuint>& planeIdx,
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

      const CFreal lbd   = 1.0 - xi - eta;

      const CFreal L1 = 0.5*(zta*zta - zta);
      const CFreal L2 = 0.5*(zta*zta + zta);
      const CFreal L3 = (1.0 - zta*zta);

      const CFreal dL1 = zta - 0.5;
      const CFreal dL2 = zta + 0.5;
      const CFreal dL3 = - 2.*zta;


      _gradShapFunc[KSI][0]  = (-3. + 4.*eta + 4.*xi)*L1;
      _gradShapFunc[KSI][1]  = (4.*xi - 1.)*L1;
      _gradShapFunc[KSI][2]  = 0.;

      _gradShapFunc[KSI][6]  = (4. - 8.*xi - 4.*eta)*L1;
      _gradShapFunc[KSI][7]  = (4.*eta)*L1;
      _gradShapFunc[KSI][8]  = (-4.*eta)*L1;

      _gradShapFunc[KSI][3]  = (-3. + 4.*eta + 4.*xi)*L2;
      _gradShapFunc[KSI][4] = (4.*xi - 1.)*L2;
      _gradShapFunc[KSI][5] = 0.;

      _gradShapFunc[KSI][15] = (4. - 8.*xi - 4.*eta)*L2;
      _gradShapFunc[KSI][16] = (4.*eta)*L2;
      _gradShapFunc[KSI][17] = (-4.*eta)*L2;

      _gradShapFunc[KSI][9]  = (-3. + 4.*eta + 4.*xi)*L3;
      _gradShapFunc[KSI][11]  = (4.*xi - 1.)*L3;
      _gradShapFunc[KSI][13]  = 0.;

      _gradShapFunc[KSI][10] = (4. - 8.*xi - 4.*eta)*L3;
      _gradShapFunc[KSI][12] = (4.*eta)*L3;
      _gradShapFunc[KSI][14] = (-4.*eta)*L3;


      _gradShapFunc[ETA][0]  = (-3. + 4.*eta + 4.*xi)*L1;
      _gradShapFunc[ETA][1]  = 0.;
      _gradShapFunc[ETA][2]  = (4.*eta - 1.)*L1;;

      _gradShapFunc[ETA][6]  = (-4.*xi)*L1;
      _gradShapFunc[ETA][7]  = (4.*xi)*L1;
      _gradShapFunc[ETA][8]  = (4. - 4.*xi - 8.*eta)*L1;

      _gradShapFunc[ETA][3]  = (-3. + 4.*eta + 4.*xi)*L2;
      _gradShapFunc[ETA][4] = 0.;
      _gradShapFunc[ETA][5] = (4.*eta - 1.)*L2;

      _gradShapFunc[ETA][15] = (-4.*xi)*L2;
      _gradShapFunc[ETA][16] = (4.*xi)*L2;
      _gradShapFunc[ETA][17] = (4. - 4.*xi - 8.*eta)*L2;

      _gradShapFunc[ETA][9]  = (-3. + 4.*eta + 4.*xi)*L3;
      _gradShapFunc[ETA][11]  = 0.;
      _gradShapFunc[ETA][13]  = (4.*eta - 1.)*L3;

      _gradShapFunc[ETA][10] = (-4.*xi)*L3;
      _gradShapFunc[ETA][12] = (4.*xi)*L3;
      _gradShapFunc[ETA][14] = (4. - 4.*xi - 8.*eta)*L3;



      _gradShapFunc[ZTA][0]  = lbd*(2.*lbd - 1.)*dL1;
      _gradShapFunc[ZTA][1]  = xi * (2. * xi - 1.)*dL1;
      _gradShapFunc[ZTA][2]  = eta * (2. * eta - 1.)*dL1;

      _gradShapFunc[ZTA][6]  = 4.*xi*lbd*dL1;
      _gradShapFunc[ZTA][7]  = 4.*xi*eta*dL1;
      _gradShapFunc[ZTA][8]  = 4.*eta*lbd*dL1;

      _gradShapFunc[ZTA][3]  = lbd*(2.*lbd - 1.)*dL2;
      _gradShapFunc[ZTA][4] = xi * (2. * xi - 1.)*dL2;
      _gradShapFunc[ZTA][5] = eta * (2. * eta - 1.)*dL2;

      _gradShapFunc[ZTA][15] = 4.*xi*lbd*dL2;
      _gradShapFunc[ZTA][16] = 4.*xi*eta*dL2;
      _gradShapFunc[ZTA][17] = 4.*eta*lbd*dL2;

      _gradShapFunc[ZTA][9]  = lbd*(2.*lbd - 1.)*dL3;
      _gradShapFunc[ZTA][11]  = xi * (2. * xi - 1.)*dL3;
      _gradShapFunc[ZTA][13]  = eta * (2. * eta - 1.)*dL3;

      _gradShapFunc[ZTA][10] = 4.*xi*lbd*dL3;
      _gradShapFunc[ZTA][12] = 4.*xi*eta*dL3;
      _gradShapFunc[ZTA][14] = 4.*eta*lbd*dL3;

      cf_assert(planeIdx[ip] == 0 || planeIdx[ip] == 1 || planeIdx[ip] == 2 || planeIdx[ip] == 3);

      if (planeIdx[ip] == 0) // x-direction
    {
      _vec1 = _gradShapFunc[ETA][0]*(*nodes[0]);
      _vec2 = _gradShapFunc[ZTA][0]*(*nodes[0])/2.;
      for (CFuint in = 1; in < 18; ++in)
      {
        _vec1 += _gradShapFunc[ETA][in]*(*nodes[in]);
        _vec2 += _gradShapFunc[ZTA][in]*(*nodes[in])/2.;
      }
    }
    else if (planeIdx[ip] == 1) // y-direction
    {
      _vec1 = _gradShapFunc[ZTA][0]*(*nodes[0])/2.;
      _vec2 = _gradShapFunc[KSI][0]*(*nodes[0]);
      for (CFuint in = 1; in < 18; ++in)
      {
        _vec1 += _gradShapFunc[ZTA][in]*(*nodes[in])/2.;
        _vec2 += _gradShapFunc[KSI][in]*(*nodes[in]);
      }
    }
    else if (planeIdx[ip] == 2) // z-direction
    {
      _vec1 = _gradShapFunc[KSI][0]*(*nodes[0]);
      _vec2 = _gradShapFunc[ETA][0]*(*nodes[0]);
      for (CFuint in = 1; in < 18; ++in)
      {
        _vec1 += _gradShapFunc[KSI][in]*(*nodes[in]);
        _vec2 += _gradShapFunc[ETA][in]*(*nodes[in]);
      }
    }
    else  // oblique face noraml direction (x+y)
    {
      _vec1 = -_gradShapFunc[KSI][0]*(*nodes[0])+_gradShapFunc[ETA][0]*(*nodes[0]);
      _vec2 = -_gradShapFunc[ZTA][0]*(*nodes[0])/2.;
      for (CFuint in = 1; in < 18; ++in)
      {
        _vec1 += -_gradShapFunc[KSI][in]*(*nodes[in])+_gradShapFunc[ETA][in]*(*nodes[in]);
        _vec2 += -_gradShapFunc[ZTA][in]*(*nodes[in])/2.;
      }
    }

      // compute normal
      MathTools::MathFunctions::crossProd(_vec1,_vec2,pointNormal);
   }
}


//////////////////////////////////////////////////////////////////////////////

    } // namespace ShapeFunctions



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
