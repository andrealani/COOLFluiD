// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "ShapeFunctions/LagrangeShapeFunctionLineP1.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

   namespace ShapeFunctions {

//////////////////////////////////////////////////////////////////////////////

CFuint LagrangeShapeFunctionLineP1::_interpolatorID = 0;
RealVector LagrangeShapeFunctionLineP1::_vec = RealVector(DIM_2D);
RealVector LagrangeShapeFunctionLineP1::_vec3D = RealVector(DIM_3D);
RealVector LagrangeShapeFunctionLineP1::_tmpVec1 = RealVector(DIM_2D);
RealVector LagrangeShapeFunctionLineP1::_tmpVec2 = RealVector(DIM_2D);
RealVector LagrangeShapeFunctionLineP1::_tmpVec3 = RealVector(DIM_2D);
RealVector LagrangeShapeFunctionLineP1::_tmpVec4 = RealVector(DIM_2D);
RealMatrix LagrangeShapeFunctionLineP1::_tmpMat2D = RealMatrix(DIM_2D,DIM_2D);
RealMatrix LagrangeShapeFunctionLineP1::_tmpMat3D = RealMatrix(DIM_3D,DIM_3D);
std::vector<RealVector> LagrangeShapeFunctionLineP1::_normals = std::vector<RealVector>(1);

//////////////////////////////////////////////////////////////////////////////

void LagrangeShapeFunctionLineP1::computeFaceJacobDetVectorAtMappedCoords(const std::vector<RealVector>& mappedCoord,
                                                                          const std::vector<Framework::Node*>& nodes,
                                                                          std::vector<RealVector>& normal)
{
  const CFreal x1Mx0 = 0.5 * ((*nodes[1])[XX] - (*nodes[0])[XX]);
  const CFreal y1My0 = 0.5 * ((*nodes[1])[YY] - (*nodes[0])[YY]);

  for (CFuint ip = 0; ip < mappedCoord.size(); ++ip)
  {
    RealVector& pointNormal = normal[ip];

    pointNormal[XX] = +y1My0;
    pointNormal[YY] = -x1Mx0;
  }
}

//////////////////////////////////////////////////////////////////////////////

   } // namespace ShapeFunctions

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
