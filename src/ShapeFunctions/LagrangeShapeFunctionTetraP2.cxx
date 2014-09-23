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

//////////////////////////////////////////////////////////////////////////////

    } // namespace ShapeFunctions



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
