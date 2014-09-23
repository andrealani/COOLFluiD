// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "ShapeFunctions/LagrangeShapeFunctionTriagP2.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



   namespace ShapeFunctions {

//////////////////////////////////////////////////////////////////////////////

CFuint LagrangeShapeFunctionTriagP2::_interpolatorID = 0;
RealVector LagrangeShapeFunctionTriagP2::m_dPhidxi = RealVector(6);
RealVector LagrangeShapeFunctionTriagP2::m_dPhideta = RealVector(6);
RealVector LagrangeShapeFunctionTriagP2::m_vec1 = RealVector(DIM_3D);
RealVector LagrangeShapeFunctionTriagP2::m_vec2 = RealVector(DIM_3D);
RealMatrix LagrangeShapeFunctionTriagP2::m_invJ(2,2);
RealVector LagrangeShapeFunctionTriagP2::m_cross(DIM_3D);
RealVector LagrangeShapeFunctionTriagP2::m_vnormal(DIM_3D);
RealMatrix LagrangeShapeFunctionTriagP2::m_jac3d(DIM_3D,DIM_3D);
RealMatrix LagrangeShapeFunctionTriagP2::m_vmatrix(3,3);
RealVector LagrangeShapeFunctionTriagP2::m_mappedCoord(DIM_2D);
RealVector LagrangeShapeFunctionTriagP2::m_mappedCoord3D(DIM_3D);
RealVector LagrangeShapeFunctionTriagP2::m_vCA(DIM_2D);
RealVector LagrangeShapeFunctionTriagP2::m_vPA(DIM_2D);
RealVector LagrangeShapeFunctionTriagP2::m_vBA(DIM_2D);
RealVector LagrangeShapeFunctionTriagP2::m_vCA3D(DIM_3D);
RealVector LagrangeShapeFunctionTriagP2::m_vPA3D(DIM_3D);
RealVector LagrangeShapeFunctionTriagP2::m_vBA3D(DIM_3D);
RealVector LagrangeShapeFunctionTriagP2::m_temp3DVector(DIM_3D);
RealMatrix LagrangeShapeFunctionTriagP2::m_temp3DMatrix(DIM_3D,DIM_3D);
std::vector<RealVector> LagrangeShapeFunctionTriagP2::m_normals = std::vector<RealVector>(3,LagrangeShapeFunctionTriagP2::m_vBA);

//////////////////////////////////////////////////////////////////////////////

LagrangeShapeFunctionTriagP2::LagrangeShapeFunctionTriagP2()
{
}

//////////////////////////////////////////////////////////////////////////////

void LagrangeShapeFunctionTriagP2::computeFaceJacobianDeterminant(
        const std::vector<RealVector>& m_mappedCoord,
        const std::vector<Framework::Node*>& nodes,
        const Framework::IntegratorPattern& pattern,
              std::vector<RealVector>& faceJacobian)
{
    throw Common::NotImplementedException
      (FromHere(), getName()+"::computeFaceJacobianDeterminant()");
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace ShapeFunctions



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
