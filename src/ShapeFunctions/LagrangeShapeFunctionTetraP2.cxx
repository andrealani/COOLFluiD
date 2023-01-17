// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "ShapeFunctions/LagrangeShapeFunctionTetraP2.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



   namespace ShapeFunctions {

//////////////////////////////////////////////////////////////////////////////

CFuint LagrangeShapeFunctionTetraP2::_interpolatorID = 0;
RealMatrix LagrangeShapeFunctionTetraP2::m_invJ(3,3);
RealVector LagrangeShapeFunctionTetraP2::_vec1 = RealVector(DIM_3D);
RealVector LagrangeShapeFunctionTetraP2::_vec2 = RealVector(DIM_3D);
RealVector LagrangeShapeFunctionTetraP2::_vec3 = RealVector(DIM_3D);
RealVector LagrangeShapeFunctionTetraP2::_vec4 = RealVector(DIM_3D);
RealVector LagrangeShapeFunctionTetraP2::m_mappedCoord(DIM_3D);
RealMatrix LagrangeShapeFunctionTetraP2::m_matrix1(DIM_3D,DIM_3D);
RealMatrix LagrangeShapeFunctionTetraP2::m_matrix2(DIM_3D,DIM_3D);
std::vector<RealVector> LagrangeShapeFunctionTetraP2::_normals = std::vector<RealVector>(4,LagrangeShapeFunctionTetraP2::m_mappedCoord);
std::vector<RealVector> LagrangeShapeFunctionTetraP2::_gradShapFunc = std::vector<RealVector>(3,RealVector(10));

//////////////////////////////////////////////////////////////////////////////


void LagrangeShapeFunctionTetraP2::computeJacobian(
      const std::vector<Framework::Node*>& nodes,
      const std::vector<RealVector>& mappedCoord,
            std::vector<RealMatrix>& jacob)
 {
   for (CFuint ip = 0; ip < mappedCoord.size(); ++ip) {

     RealMatrix& pointJacob = jacob[ip];
       
       const CFreal lambda = 1.0 - mappedCoord[ip].sum();
       const CFreal xi =  mappedCoord[ip][KSI];
       const CFreal eta = mappedCoord[ip][ETA];
       const CFreal zta = mappedCoord[ip][ZTA];


       _gradShapFunc[KSI][ 0] = 1. - 4. * lambda;
       _gradShapFunc[KSI][ 1] = -1. + 4. * xi;
       _gradShapFunc[KSI][ 2] = 0.;
       _gradShapFunc[KSI][ 3] = 0.;
       _gradShapFunc[KSI][ 4] = 4.*(lambda - xi);
       _gradShapFunc[KSI][ 5] = 4. * eta;
       _gradShapFunc[KSI][ 6] = -4.* eta;
       _gradShapFunc[KSI][ 7] = 4. * zta;
       _gradShapFunc[KSI][ 8] = 0.;
       _gradShapFunc[KSI][ 9] = -4. * zta;


       _gradShapFunc[ETA][ 0] = 1. - 4. * lambda;
       _gradShapFunc[ETA][ 1] = 0.;
       _gradShapFunc[ETA][ 2] = -1. + 4. * eta;
       _gradShapFunc[ETA][ 3] = 0.;
       _gradShapFunc[ETA][ 4] = -4. * xi;
       _gradShapFunc[ETA][ 5] = 4. * xi;
       _gradShapFunc[ETA][ 6] = 4.*(lambda - eta);
       _gradShapFunc[ETA][ 7] = 0.;
       _gradShapFunc[ETA][ 8] = 4. * zta;
       _gradShapFunc[ETA][ 9] = -4. * zta;

       
       _gradShapFunc[ZTA][ 0] = 1. - 4. * lambda;
       _gradShapFunc[ZTA][ 1] = 0.;
       _gradShapFunc[ZTA][ 2] = 0.;
       _gradShapFunc[ZTA][ 3] = -1. + 4. * zta;
       _gradShapFunc[ZTA][ 4] = -4. * xi;
       _gradShapFunc[ZTA][ 5] = 0.;
       _gradShapFunc[ZTA][ 6] = -4. * eta;
       _gradShapFunc[ZTA][ 7] = 4. * xi;
       _gradShapFunc[ZTA][ 8] = 4. * eta;
       _gradShapFunc[ZTA][ 9] = 4. * (lambda - zta);


     // compute Jacobian matrix
     const CFreal xx = (*nodes[0])[XX];
     const CFreal yy = (*nodes[0])[YY];
     const CFreal zz = (*nodes[0])[ZZ];
     const CFreal dNdksi = _gradShapFunc[KSI][0];
     const CFreal dNdeta = _gradShapFunc[ETA][0];
     const CFreal dNdzta = _gradShapFunc[ZTA][0];

     pointJacob(KSI,XX) = xx*dNdksi;
     pointJacob(ETA,XX) = xx*dNdeta;
     pointJacob(ZTA,XX) = xx*dNdzta;

     pointJacob(KSI,YY) = yy*dNdksi;
     pointJacob(ETA,YY) = yy*dNdeta;
     pointJacob(ZTA,YY) = yy*dNdzta;

     pointJacob(KSI,ZZ) = zz*dNdksi;
     pointJacob(ETA,ZZ) = zz*dNdeta;
     pointJacob(ZTA,ZZ) = zz*dNdzta;

     for (CFuint in = 1; in < 10; ++in)
     {

       const CFreal xx = (*nodes[in])[XX];
       const CFreal yy = (*nodes[in])[YY];
       const CFreal zz = (*nodes[in])[ZZ];
       const CFreal dNdksi = _gradShapFunc[KSI][in];
       const CFreal dNdeta = _gradShapFunc[ETA][in];
       const CFreal dNdzta = _gradShapFunc[ZTA][in];

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

void LagrangeShapeFunctionTetraP2::computeMappedCoordPlaneNormal(const std::vector<CFuint>& planeIdx,
                                                                      const std::vector<RealVector>& mappedCoord,
                                                                      const std::vector<Framework::Node*>& nodes,
                                                                      std::vector<RealVector>& normal)
{
for (CFuint ip = 0; ip < mappedCoord.size(); ++ip)
{
  RealVector& pointNormal = normal[ip];

    const CFreal lambda = 1.0 - mappedCoord[ip].sum();
    const CFreal xi =  mappedCoord[ip][KSI];
    const CFreal eta = mappedCoord[ip][ETA];
    const CFreal zta = mappedCoord[ip][ZTA];


    _gradShapFunc[KSI][ 0] = 1. - 4. * lambda;
    _gradShapFunc[KSI][ 1] = -1. + 4. * xi;
    _gradShapFunc[KSI][ 2] = 0.;
    _gradShapFunc[KSI][ 3] = 0.;
    _gradShapFunc[KSI][ 4] = 4.*(lambda - xi);
    _gradShapFunc[KSI][ 5] = 4. * eta;
    _gradShapFunc[KSI][ 6] = -4.* eta;
    _gradShapFunc[KSI][ 7] = 4. * zta;
    _gradShapFunc[KSI][ 8] = 0.;
    _gradShapFunc[KSI][ 9] = -4. * zta;


    _gradShapFunc[ETA][ 0] = 1. - 4. * lambda;
    _gradShapFunc[ETA][ 1] = 0.;
    _gradShapFunc[ETA][ 2] = -1. + 4. * eta;
    _gradShapFunc[ETA][ 3] = 0.;
    _gradShapFunc[ETA][ 4] = -4. * xi;
    _gradShapFunc[ETA][ 5] = 4. * xi;
    _gradShapFunc[ETA][ 6] = 4.*(lambda - eta);
    _gradShapFunc[ETA][ 7] = 0.;
    _gradShapFunc[ETA][ 8] = 4. * zta;
    _gradShapFunc[ETA][ 9] = -4. * zta;

    
    _gradShapFunc[ZTA][ 0] = 1. - 4. * lambda;
    _gradShapFunc[ZTA][ 1] = 0.;
    _gradShapFunc[ZTA][ 2] = 0.;
    _gradShapFunc[ZTA][ 3] = -1. + 4. * zta;
    _gradShapFunc[ZTA][ 4] = -4. * xi;
    _gradShapFunc[ZTA][ 5] = 0.;
    _gradShapFunc[ZTA][ 6] = -4. * eta;
    _gradShapFunc[ZTA][ 7] = 4. * xi;
    _gradShapFunc[ZTA][ 8] = 4. * eta;
    _gradShapFunc[ZTA][ 9] = 4. * (lambda - zta);



  //cf_assert(planeIdx[ip] == 0 || planeIdx[ip] == 1 || planeIdx[ip] == 2);
    
  if (planeIdx[ip] == 0)
  {

      _vec1 = _gradShapFunc[ETA][0]*(*nodes[0]);
      _vec2 = _gradShapFunc[KSI][0]*(*nodes[0]);
        
      for (CFuint in = 1; in < 10; ++in)
      {
        _vec1 += _gradShapFunc[ETA][in]*(*nodes[in]);
        _vec2 += _gradShapFunc[KSI][in]*(*nodes[in]);
      }
  }
  else if (planeIdx[ip] == 1)
  {

      _vec1 = _gradShapFunc[KSI][0]*(*nodes[0]);
      _vec2 = _gradShapFunc[ZTA][0]*(*nodes[0]);
      for (CFuint in = 1; in < 10; ++in)
      {
        _vec1 += _gradShapFunc[KSI][in]*(*nodes[in]);
        _vec2 += _gradShapFunc[ZTA][in]*(*nodes[in]);
      }
  }
  else if (planeIdx[ip] == 3)
  {

      _vec1 = _gradShapFunc[ZTA][0]*(*nodes[0]);
      _vec2 = _gradShapFunc[ETA][0]*(*nodes[0]);
      for (CFuint in = 1; in < 10; ++in)
      {
        _vec1 += _gradShapFunc[ZTA][in]*(*nodes[in]);
        _vec2 += _gradShapFunc[ETA][in]*(*nodes[in]);
      }
  }
  else if (planeIdx[ip] == 2) //Face normal for equilateral triangle (oblique face)
  {
    
      _vec1 = -_gradShapFunc[ETA][0]*(*nodes[0]) + _gradShapFunc[ZTA][0]*(*nodes[0]) ;
      _vec2 = -_gradShapFunc[ETA][0]*(*nodes[0]) + _gradShapFunc[KSI][0]*(*nodes[0]);
      for (CFuint in = 1; in < 10; ++in)
      {
        _vec1 += -_gradShapFunc[ETA][in]*(*nodes[in]) + _gradShapFunc[ZTA][in]*(*nodes[in]);
        _vec2 += -_gradShapFunc[ETA][in]*(*nodes[in]) + _gradShapFunc[KSI][in]*(*nodes[in]);
      }
  }
    
  else if (planeIdx[ip] == 4) // x
  {

      _vec1 = _gradShapFunc[ETA][0]*(*nodes[0]);
      _vec2 = _gradShapFunc[ZTA][0]*(*nodes[0]);
      for (CFuint in = 1; in < 10; ++in)
      {
        _vec1 += _gradShapFunc[ETA][in]*(*nodes[in]);
        _vec2 += _gradShapFunc[ZTA][in]*(*nodes[in]);
      }
  }
  else if (planeIdx[ip] == 5) // y
  {

      _vec1 = _gradShapFunc[ZTA][0]*(*nodes[0]);
      _vec2 = _gradShapFunc[KSI][0]*(*nodes[0]);
      for (CFuint in = 1; in < 10; ++in)
      {
        _vec1 += _gradShapFunc[ZTA][in]*(*nodes[in]);
        _vec2 += _gradShapFunc[KSI][in]*(*nodes[in]);
      }
  }
  else  // z
  {

      _vec1 = _gradShapFunc[KSI][0]*(*nodes[0]);
      _vec2 = _gradShapFunc[ETA][0]*(*nodes[0]);
      for (CFuint in = 1; in < 10; ++in)
      {
        _vec1 += _gradShapFunc[KSI][in]*(*nodes[in]);
        _vec2 += _gradShapFunc[ETA][in]*(*nodes[in]);
      }
  }

  // compute normal
  MathTools::MathFunctions::crossProd(_vec1,_vec2,pointNormal);
}
}

   
    } // namespace ShapeFunctions



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
