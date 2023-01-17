// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "ShapeFunctions/LagrangeShapeFunctionTetraP1.hh"
#include "MathTools/MatrixInverterT.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



   namespace ShapeFunctions {

//////////////////////////////////////////////////////////////////////////////

CFuint LagrangeShapeFunctionTetraP1::_interpolatorID = 0;
RealVector LagrangeShapeFunctionTetraP1::_vec1 = RealVector(DIM_3D);
RealVector LagrangeShapeFunctionTetraP1::_vec2 = RealVector(DIM_3D);
RealVector LagrangeShapeFunctionTetraP1::_vec3 = RealVector(DIM_3D);
RealVector LagrangeShapeFunctionTetraP1::_vec4 = RealVector(DIM_3D);
RealVector LagrangeShapeFunctionTetraP1::m_mappedCoord(DIM_3D);
RealMatrix LagrangeShapeFunctionTetraP1::m_matrix1(DIM_3D,DIM_3D);
RealMatrix LagrangeShapeFunctionTetraP1::m_matrix2(DIM_3D,DIM_3D);
std::vector<RealVector> LagrangeShapeFunctionTetraP1::_normals = std::vector<RealVector>(4,LagrangeShapeFunctionTetraP1::m_mappedCoord);
std::vector<RealVector> LagrangeShapeFunctionTetraP1::_gradShapFunc = std::vector<RealVector>(3,RealVector(4));

//////////////////////////////////////////////////////////////////////////////

RealVector LagrangeShapeFunctionTetraP1::computeMappedCoordinates(const RealVector& coord,
                                    const std::vector<Framework::Node*>& nodes)
{
  cf_assert(coord.size() == DIM_3D);
  cf_assert(nodes.size() == 4);
  cf_assert(nodes[0]->size() ==  DIM_3D);
  cf_assert(nodes[1]->size() == DIM_3D);
  cf_assert(nodes[2]->size() == DIM_3D);
  cf_assert(nodes[3]->size() == DIM_3D);

  RealVector& xA = (*nodes[0]);
  RealVector& xB = (*nodes[1]);
  RealVector& xC = (*nodes[2]);
  RealVector& xD = (*nodes[3]);

  _vec1 = xB-xA;
  _vec2 = xC-xA;
  _vec3 = xD-xA;
  _vec4 = coord-xA;

  m_matrix1(0,0) = _vec1[XX];
  m_matrix1(0,1) = _vec2[XX];
  m_matrix1(0,2) = _vec3[XX];
  m_matrix1(1,0) = _vec1[YY];
  m_matrix1(1,1) = _vec2[YY];
  m_matrix1(1,2) = _vec3[YY];
  m_matrix1(2,0) = _vec1[ZZ];
  m_matrix1(2,1) = _vec2[ZZ];
  m_matrix1(2,2) = _vec3[ZZ];

  MathTools::MatrixInverterT<3> inverter;
  inverter.invert(m_matrix1,m_matrix2);

  m_mappedCoord = m_matrix2 * _vec4;

  return m_mappedCoord;
}
   
//////////////////////////////////////////////////////////////////////////////

void LagrangeShapeFunctionTetraP1::computeMappedCoordPlaneNormal(const std::vector<CFuint>& planeIdx,
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

     
   if (planeIdx[ip] == 0)
   {
       _gradShapFunc[KSI][ 0] = - 1.;
       _gradShapFunc[KSI][ 1] = 1. ;
       _gradShapFunc[KSI][ 2] = 0.;
       _gradShapFunc[KSI][ 3] = 0.;
       
       _gradShapFunc[ETA][ 0] = - 1.;
       _gradShapFunc[ETA][ 1] = 0.;
       _gradShapFunc[ETA][ 2] = 1.;
       _gradShapFunc[ETA][ 3] = 0.;
       
     _vec1 = _gradShapFunc[ETA][0]*(*nodes[0]);
     _vec2 = _gradShapFunc[KSI][0]*(*nodes[0]);
       
     for (CFuint in = 1; in < 4; ++in)
     {
       _vec1 += _gradShapFunc[ETA][in]*(*nodes[in]);
       _vec2 += _gradShapFunc[KSI][in]*(*nodes[in]);
     }
   }
   else if (planeIdx[ip] == 1)
   {
       
       _gradShapFunc[KSI][ 0] = - 1.;
       _gradShapFunc[KSI][ 1] = 1. ;
       _gradShapFunc[KSI][ 2] = 0.;
       _gradShapFunc[KSI][ 3] = 0.;
       
       _gradShapFunc[ZTA][ 0] = -1.;
       _gradShapFunc[ZTA][ 1] = 0.;
       _gradShapFunc[ZTA][ 2] = 0.;
       _gradShapFunc[ZTA][ 3] = 1.;
       
     _vec1 = _gradShapFunc[KSI][0]*(*nodes[0]);
     _vec2 = _gradShapFunc[ZTA][0]*(*nodes[0]);
     for (CFuint in = 1; in < 4; ++in)
     {
       _vec1 += _gradShapFunc[KSI][in]*(*nodes[in]);
       _vec2 += _gradShapFunc[ZTA][in]*(*nodes[in]);
     }
   }
   else if (planeIdx[ip] == 3)
   {
       _gradShapFunc[ETA][ 0] = - 1.;
       _gradShapFunc[ETA][ 1] = 0.;
       _gradShapFunc[ETA][ 2] = 1.;
       _gradShapFunc[ETA][ 3] = 0.;

       
       _gradShapFunc[ZTA][ 0] = -1.;
       _gradShapFunc[ZTA][ 1] = 0.;
       _gradShapFunc[ZTA][ 2] = 0.;
       _gradShapFunc[ZTA][ 3] = 1.;
              
     _vec1 = _gradShapFunc[ZTA][0]*(*nodes[0]);
     _vec2 = _gradShapFunc[ETA][0]*(*nodes[0]);
     for (CFuint in = 1; in < 4; ++in)
     {
       _vec1 += _gradShapFunc[ZTA][in]*(*nodes[in]);
       _vec2 += _gradShapFunc[ETA][in]*(*nodes[in]);
     }
   }
   else if (planeIdx[ip] == 2) //Face normal for equilateral triangle (oblique face)
   {
       _gradShapFunc[KSI][ 0] = - 1.;
       _gradShapFunc[KSI][ 1] = 1. ;
       _gradShapFunc[KSI][ 2] = 0.;
       _gradShapFunc[KSI][ 3] = 0.;
       
       _gradShapFunc[ETA][ 0] = - 1.;
       _gradShapFunc[ETA][ 1] = 0.;
       _gradShapFunc[ETA][ 2] = 1.;
       _gradShapFunc[ETA][ 3] = 0.;
       
       _gradShapFunc[ZTA][ 0] = -1.;
       _gradShapFunc[ZTA][ 1] = 0.;
       _gradShapFunc[ZTA][ 2] = 0.;
       _gradShapFunc[ZTA][ 3] = 1.;
       
     _vec1 = -_gradShapFunc[ETA][0]*(*nodes[0]) + _gradShapFunc[ZTA][0]*(*nodes[0]) ;
     _vec2 = -_gradShapFunc[ETA][0]*(*nodes[0]) + _gradShapFunc[KSI][0]*(*nodes[0]);
     for (CFuint in = 1; in < 4; ++in)
     {
       _vec1 += -_gradShapFunc[ETA][in]*(*nodes[in]) + _gradShapFunc[ZTA][in]*(*nodes[in]);
       _vec2 += -_gradShapFunc[ETA][in]*(*nodes[in]) + _gradShapFunc[KSI][in]*(*nodes[in]);
     }
   }
     
   else if (planeIdx[ip] == 4) // x
   {
       _gradShapFunc[ETA][ 0] = - 1.;
       _gradShapFunc[ETA][ 1] = 0.;
       _gradShapFunc[ETA][ 2] = 1.;
       _gradShapFunc[ETA][ 3] = 0.;

       
       _gradShapFunc[ZTA][ 0] = -1.;
       _gradShapFunc[ZTA][ 1] = 0.;
       _gradShapFunc[ZTA][ 2] = 0.;
       _gradShapFunc[ZTA][ 3] = 1.;
       
       _vec1 = _gradShapFunc[ETA][0]*(*nodes[0]);
       _vec2 = _gradShapFunc[ZTA][0]*(*nodes[0]);
       for (CFuint in = 1; in < 4; ++in)
       {
         _vec1 += _gradShapFunc[ETA][in]*(*nodes[in]);
         _vec2 += _gradShapFunc[ZTA][in]*(*nodes[in]);
       }
   }
   else if (planeIdx[ip] == 5) // y
   {
       _gradShapFunc[KSI][ 0] = - 1.;
       _gradShapFunc[KSI][ 1] = 1. ;
       _gradShapFunc[KSI][ 2] = 0.;
       _gradShapFunc[KSI][ 3] = 0.;
       
       _gradShapFunc[ZTA][ 0] = -1.;
       _gradShapFunc[ZTA][ 1] = 0.;
       _gradShapFunc[ZTA][ 2] = 0.;
       _gradShapFunc[ZTA][ 3] = 1.;
       
       _vec1 = _gradShapFunc[ZTA][0]*(*nodes[0]);
       _vec2 = _gradShapFunc[KSI][0]*(*nodes[0]);
       for (CFuint in = 1; in < 4; ++in)
       {
         _vec1 += _gradShapFunc[ZTA][in]*(*nodes[in]);
         _vec2 += _gradShapFunc[KSI][in]*(*nodes[in]);
       }
   }
   else  // z
   {
       _gradShapFunc[KSI][ 0] = - 1.;
       _gradShapFunc[KSI][ 1] = 1. ;
       _gradShapFunc[KSI][ 2] = 0.;
       _gradShapFunc[KSI][ 3] = 0.;
       
       _gradShapFunc[ETA][ 0] = - 1.;
       _gradShapFunc[ETA][ 1] = 0.;
       _gradShapFunc[ETA][ 2] = 1.;
       _gradShapFunc[ETA][ 3] = 0.;
       
       _vec1 = _gradShapFunc[KSI][0]*(*nodes[0]);
       _vec2 = _gradShapFunc[ETA][0]*(*nodes[0]);
       for (CFuint in = 1; in < 4; ++in)
       {
         _vec1 += _gradShapFunc[KSI][in]*(*nodes[in]);
         _vec2 += _gradShapFunc[ETA][in]*(*nodes[in]);
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
