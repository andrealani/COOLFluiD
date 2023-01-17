// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "ShapeFunctions/LagrangeShapeFunctionTriagP1.hh"
#include "MathTools/MatrixInverterT.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



   namespace ShapeFunctions {

//////////////////////////////////////////////////////////////////////////////

CFuint LagrangeShapeFunctionTriagP1::_interpolatorID = 0;
RealVector LagrangeShapeFunctionTriagP1::m_vec1 = RealVector(DIM_3D);
RealVector LagrangeShapeFunctionTriagP1::m_vec2 = RealVector(DIM_3D);
RealMatrix LagrangeShapeFunctionTriagP1::m_invJ(2,2);
RealVector LagrangeShapeFunctionTriagP1::m_cross(DIM_3D);
RealVector LagrangeShapeFunctionTriagP1::m_vnormal(DIM_3D);
RealMatrix LagrangeShapeFunctionTriagP1::m_jac3d(DIM_3D,DIM_3D);
RealMatrix LagrangeShapeFunctionTriagP1::m_vmatrix(3,3);
RealVector LagrangeShapeFunctionTriagP1::m_mappedCoord(DIM_2D);
RealVector LagrangeShapeFunctionTriagP1::m_mappedCoord3D(DIM_3D);
RealVector LagrangeShapeFunctionTriagP1::m_vCA(DIM_2D);
RealVector LagrangeShapeFunctionTriagP1::m_vPA(DIM_2D);
RealVector LagrangeShapeFunctionTriagP1::m_vBA(DIM_2D);
RealVector LagrangeShapeFunctionTriagP1::m_vCA3D(DIM_3D);
RealVector LagrangeShapeFunctionTriagP1::m_vPA3D(DIM_3D);
RealVector LagrangeShapeFunctionTriagP1::m_vBA3D(DIM_3D);
RealVector LagrangeShapeFunctionTriagP1::m_temp2DVector(DIM_2D);
RealVector LagrangeShapeFunctionTriagP1::m_temp3DVector(DIM_3D);
RealMatrix LagrangeShapeFunctionTriagP1::m_temp3DMatrix(DIM_3D,DIM_3D);
RealMatrix LagrangeShapeFunctionTriagP1::m_temp3DMatrix2(DIM_3D,DIM_3D);
std::vector<RealVector> LagrangeShapeFunctionTriagP1::m_normals = std::vector<RealVector>(3,LagrangeShapeFunctionTriagP1::m_vBA);

//////////////////////////////////////////////////////////////////////////////

LagrangeShapeFunctionTriagP1::LagrangeShapeFunctionTriagP1()
{
}

//////////////////////////////////////////////////////////////////////////////

void LagrangeShapeFunctionTriagP1::computeFaceJacobianDeterminant(
        const std::vector<RealVector>& m_mappedCoord,
        const std::vector<Framework::Node*>& nodes,
        const Framework::IntegratorPattern& pattern,
              std::vector<RealVector>& faceJacobian)
{
  cf_assert(pattern.nbSteps() == getNbFaces());
  cf_assert(faceJacobian.size() >= getNbFaces());

  const CFreal x0 = (*nodes[0])[XX];
  const CFreal y0 = (*nodes[0])[YY];

  const CFreal x1 = (*nodes[1])[XX];
  const CFreal y1 = (*nodes[1])[YY];

  const CFreal x2 = (*nodes[2])[XX];
  const CFreal y2 = (*nodes[2])[YY];

  CFLogDebugMax( "X0: " << x0 << "\n");
  CFLogDebugMax( "X1: " << x1 << "\n");
  CFLogDebugMax( "X2: " << x2 << "\n");
  CFLogDebugMax( "Y0: " << y0 << "\n");
  CFLogDebugMax( "Y1: " << y1 << "\n");
  CFLogDebugMax( "Y2: " << y2 << "\n");

  cf_assert(pattern.nbSteps() == 3);

  const CFreal jacob0 =
    Framework::FaceJacobiansDeterminant::compute2DFaceJacobDet(x1-x0,y1-y0);

  const CFreal jacob1 =
    Framework::FaceJacobiansDeterminant::compute2DFaceJacobDet(x1-x2,y1-y2);

  const CFreal jacob2 =
    Framework::FaceJacobiansDeterminant::compute2DFaceJacobDet(x2-x0,y2-y0);

  CFLogDebugMax( "Jacob0: " << jacob0 << "\n");
  CFLogDebugMax( "Jacob1: " << jacob1 << "\n");
  CFLogDebugMax( "Jacob2: " << jacob2 << "\n");

  const CFuint iFace0 = 0;
  for(CFuint iPoint = 0; iPoint < pattern.nbPts(iFace0); ++iPoint) {
      faceJacobian[iFace0][iPoint] = jacob0;
  }
  // the Jacobian of this face is devided by sqrt(2.0) 
  // because in the triangle used as reference 
  // The lenght of this face is sqrt(2.0) 
  const CFuint iFace1 = 1;
  for(CFuint iPoint = 0; iPoint < pattern.nbPts(iFace1); ++iPoint) {
      faceJacobian[iFace1][iPoint] = jacob1/sqrt(2.0);
  }
  const CFuint iFace2 = 2;
  for(CFuint iPoint = 0; iPoint < pattern.nbPts(iFace2); ++iPoint) {
      faceJacobian[iFace2][iPoint] = jacob2;
  }
}

//////////////////////////////////////////////////////////////////////////////

RealVector LagrangeShapeFunctionTriagP1::computeMappedCoordinates(const RealVector& coord,
                                    const std::vector<Framework::Node*>& nodes)
{
  cf_assert(coord.size() == DIM_2D);
  cf_assert(nodes.size() == 3);
  cf_assert(nodes[0]->size() == DIM_2D);
  cf_assert(nodes[1]->size() == DIM_2D);
  cf_assert(nodes[2]->size() == DIM_2D);

  RealVector& xA = (*nodes[0]);
  RealVector& xB = (*nodes[1]);
  RealVector& xC = (*nodes[2]);

  m_vBA = xB-xA;
  m_vCA = xC-xA;
  m_vPA = coord-xA;

  const CFreal overDet = 1./(m_vBA[XX]*m_vCA[YY] - m_vBA[YY]*m_vCA[XX]);


/* RealMatrix m_matrix1(2,2);
 RealMatrix m_matrix2(2,2);
///@todo improve by removing:
/// - inversion of the matrix
/// - the allocation of the matrices
/// - idem for the other functions of this type!!
  m_matrix1(0,0) = m_vBA[XX];
  m_matrix1(0,1) = m_vCA[XX];
  m_matrix1(1,0) = m_vBA[YY];
  m_matrix1(1,1) = m_vCA[YY];

  m_matrix1.invert2(m_matrix2);

  m_mappedCoord = m_matrix2 * m_vPA;
*/

  m_mappedCoord[KSI] = overDet*( m_vCA[YY]*m_vPA[XX] - m_vCA[XX]*m_vPA[YY]);
  m_mappedCoord[ETA] = overDet*(-m_vBA[YY]*m_vPA[XX] + m_vBA[XX]*m_vPA[YY]);

  return m_mappedCoord;
}

//////////////////////////////////////////////////////////////////////////////

RealVector LagrangeShapeFunctionTriagP1::computeMappedCoordinatesPlus1D(const RealVector& coord,
                                    const std::vector<Framework::Node*>& nodes)
{
  cf_assert(coord.size() == DIM_3D);
  cf_assert(nodes.size() == getNbNodes());

  RealVector& xA = (*nodes[0]);
  RealVector& xB = (*nodes[1]);
  RealVector& xC = (*nodes[2]);

  m_vBA3D = xB-xA;
  m_vCA3D = xC-xA;
  m_vPA3D = coord-xA;

  MathTools::MathFunctions::crossProd(m_vBA3D,m_vCA3D,m_temp3DVector);

  ///@todo improve by removing:
  /// - inversion of the matrix

  m_temp3DMatrix(0,0) = m_vBA3D[XX];
  m_temp3DMatrix(0,1) = m_vCA3D[XX];
  m_temp3DMatrix(0,2) = m_temp3DVector[XX];
  m_temp3DMatrix(1,0) = m_vBA3D[YY];
  m_temp3DMatrix(1,1) = m_vCA3D[YY];
  m_temp3DMatrix(1,2) = m_temp3DVector[YY];
  m_temp3DMatrix(2,0) = m_vBA3D[ZZ];
  m_temp3DMatrix(2,1) = m_vCA3D[ZZ];
  m_temp3DMatrix(2,2) = m_temp3DVector[ZZ];

  MathTools::MatrixInverterT<3> inverter;
  inverter.invert(m_temp3DMatrix,m_temp3DMatrix2);

  m_mappedCoord3D = m_temp3DMatrix2 * m_vPA3D;

  m_mappedCoord[XX] = m_mappedCoord3D[XX];
  m_mappedCoord[YY] = m_mappedCoord3D[YY];

  return m_mappedCoord;
}

//////////////////////////////////////////////////////////////////////////////
   
void LagrangeShapeFunctionTriagP1::computeFaceJacobDetVectorAtMappedCoords(const std::vector<RealVector>& mappedCoord,
   const std::vector<Framework::Node*>& nodes,
   std::vector<RealVector>& normal)                  
{
   for (CFuint ip = 0; ip < mappedCoord.size(); ++ip)
   {
     RealVector& pointNormal = normal[ip];

     const CFreal xi =  mappedCoord[ip][KSI];
     const CFreal eta = mappedCoord[ip][ETA];
     const CFreal dN0dksi = -(1.);
     const CFreal dN1dksi =  (1.);
     const CFreal dN2dksi =  (0.);

     const CFreal dN0deta = -(1.);
     const CFreal dN1deta =  (0.);
     const CFreal dN2deta =  (1.);

     // compute shape function gradient
     m_vec1 = (*nodes[0])*dN0dksi + (*nodes[1])*dN1dksi + (*nodes[2])*dN2dksi ;
     m_vec2 = (*nodes[0])*dN0deta + (*nodes[1])*dN1deta + (*nodes[2])*dN2deta ;

     // compute face jacobian vector
     MathTools::MathFunctions::crossProd(m_vec1,m_vec2,pointNormal);
     pointNormal*=-1.;
   }
}
   
//////////////////////////////////////////////////////////////////////////////

    } // namespace ShapeFunctions



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
