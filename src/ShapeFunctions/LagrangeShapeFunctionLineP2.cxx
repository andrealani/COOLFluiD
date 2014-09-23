// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "ShapeFunctions/LagrangeShapeFunctionLineP2.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

   namespace ShapeFunctions {

//////////////////////////////////////////////////////////////////////////////

CFuint LagrangeShapeFunctionLineP2::_interpolatorID = 0;
RealVector LagrangeShapeFunctionLineP2::_vec = RealVector(DIM_2D);
std::vector<RealVector> LagrangeShapeFunctionLineP2::_normals = std::vector<RealVector>(1);

//////////////////////////////////////////////////////////////////////////////

void LagrangeShapeFunctionLineP2::computeFaceJacobDetVectorAtMappedCoords(const std::vector<RealVector>& mappedCoord,
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

    const CFreal dN0dxi = xi - 0.5;
    const CFreal dN1dxi = xi + 0.5;
    const CFreal dN2dxi = -2.*xi;

    pointNormal[XX] = + y0*dN0dxi + y1*dN1dxi + y2*dN2dxi;
    pointNormal[YY] = - x0*dN0dxi - x1*dN1dxi - x2*dN2dxi;
  }
}

//////////////////////////////////////////////////////////////////////////////

   } // namespace ShapeFunctions

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
