// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "ShapeFunctions/LagrangeShapeFunctionTriagP3.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



   namespace ShapeFunctions {

//////////////////////////////////////////////////////////////////////////////

CFuint LagrangeShapeFunctionTriagP3::_interpolatorID = 0;
RealVector LagrangeShapeFunctionTriagP3::m_dPhidxi = RealVector(10);
RealVector LagrangeShapeFunctionTriagP3::m_dPhideta = RealVector(10);
RealVector LagrangeShapeFunctionTriagP3::m_vec1 = RealVector(DIM_3D);
RealVector LagrangeShapeFunctionTriagP3::m_vec2 = RealVector(DIM_3D);
RealMatrix LagrangeShapeFunctionTriagP3::m_invJ(2,2);
RealVector LagrangeShapeFunctionTriagP3::m_cross(DIM_3D);
RealVector LagrangeShapeFunctionTriagP3::m_vnormal(DIM_3D);
RealMatrix LagrangeShapeFunctionTriagP3::m_jac3d(DIM_3D,DIM_3D);
RealMatrix LagrangeShapeFunctionTriagP3::m_vmatrix(3,3);
RealVector LagrangeShapeFunctionTriagP3::m_mappedCoord(DIM_2D);
RealVector LagrangeShapeFunctionTriagP3::m_mappedCoord3D(DIM_3D);
RealVector LagrangeShapeFunctionTriagP3::m_vCA(DIM_2D);
RealVector LagrangeShapeFunctionTriagP3::m_vPA(DIM_2D);
RealVector LagrangeShapeFunctionTriagP3::m_vBA(DIM_2D);
RealVector LagrangeShapeFunctionTriagP3::m_vCA3D(DIM_3D);
RealVector LagrangeShapeFunctionTriagP3::m_vPA3D(DIM_3D);
RealVector LagrangeShapeFunctionTriagP3::m_vBA3D(DIM_3D);
RealVector LagrangeShapeFunctionTriagP3::m_temp3DVector(DIM_3D);
RealMatrix LagrangeShapeFunctionTriagP3::m_temp3DMatrix(DIM_3D,DIM_3D);
std::vector<RealVector> LagrangeShapeFunctionTriagP3::m_normals = std::vector<RealVector>(3,LagrangeShapeFunctionTriagP3::m_vBA);

CFreal LagrangeShapeFunctionTriagP3::m_a1 = 0.797426985353;
CFreal LagrangeShapeFunctionTriagP3::m_b1 = 0.101286507323;
CFreal LagrangeShapeFunctionTriagP3::m_a2 = 0.059715871789;
CFreal LagrangeShapeFunctionTriagP3::m_b2 = 0.470142064105;
CFreal LagrangeShapeFunctionTriagP3::m_w1 = 0.225;
CFreal LagrangeShapeFunctionTriagP3::m_w2 = 0.125939180544;
CFreal LagrangeShapeFunctionTriagP3::m_w3 = 0.132394152788;
CFreal LagrangeShapeFunctionTriagP3::m_qdWeightsTriagP3[7] = { m_w1, m_w2, m_w2, m_w2, m_w3, m_w3, m_w3 };
CFreal LagrangeShapeFunctionTriagP3::m_xiQdPtsTriagP3[7] = { 1./3., m_b1, m_a1, m_b1, m_b2, m_a2, m_b2 };
CFreal LagrangeShapeFunctionTriagP3::m_etaQdPtsTriagP3[7] = { 1./3., m_b1, m_b1, m_a1, m_b2, m_b2, m_a2 };
CFuint LagrangeShapeFunctionTriagP3::m_nbQdPts = 7;

//////////////////////////////////////////////////////////////////////////////

LagrangeShapeFunctionTriagP3::LagrangeShapeFunctionTriagP3()
{
}

//////////////////////////////////////////////////////////////////////////////

void LagrangeShapeFunctionTriagP3::computeFaceJacobianDeterminant(
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
