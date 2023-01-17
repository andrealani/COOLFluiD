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
  
void LagrangeShapeFunctionTriagP2::computeFaceJacobDetVectorAtMappedCoords(const std::vector<RealVector>& mappedCoord,
  const std::vector<Framework::Node*>& nodes,
  std::vector<RealVector>& normal)                  
{
  for (CFuint ip = 0; ip < mappedCoord.size(); ++ip)
  {
    RealVector& pointNormal = normal[ip];

    const CFreal xi =  mappedCoord[ip][KSI];
    const CFreal eta = mappedCoord[ip][ETA];
    const CFreal dN0dksi = -3. + 4.*eta + 4.*xi;
    const CFreal dN1dksi = 4.*xi - 1.;
    const CFreal dN2dksi = 0.;
    const CFreal dN3dksi = 4. - 8.*xi - 4.*eta;
    const CFreal dN4dksi = 4.*eta;
    const CFreal dN5dksi = -4.*eta;

    const CFreal dN0deta = -3. + 4.*eta + 4.*xi;
    const CFreal dN1deta = 0.;
    const CFreal dN2deta = 4.*eta - 1.;
    const CFreal dN3deta = -4.*xi;
    const CFreal dN4deta = 4.*xi;
    const CFreal dN5deta = 4. - 4.*xi - 8.*eta;

    // compute shape function gradient
    m_vec1 = (*nodes[0])*dN0dksi + (*nodes[1])*dN1dksi + (*nodes[2])*dN2dksi + (*nodes[3])*dN3dksi + (*nodes[4])*dN4dksi + (*nodes[5])*dN5dksi ;
    m_vec2 = (*nodes[0])*dN0deta + (*nodes[1])*dN1deta + (*nodes[2])*dN2deta + (*nodes[3])*dN3deta + (*nodes[4])*dN4deta + (*nodes[5])*dN5deta ;

    // compute face jacobian vector
    MathTools::MathFunctions::crossProd(m_vec1,m_vec2,pointNormal);
        pointNormal*=-1.;
  }
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace ShapeFunctions



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
