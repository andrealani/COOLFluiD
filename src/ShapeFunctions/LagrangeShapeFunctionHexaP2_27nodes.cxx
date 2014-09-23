// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "ShapeFunctions/LagrangeShapeFunctionHexaP2_27nodes.hh"
#include "MathTools/MatrixInverterT.hh"


//////////////////////////////////////////////////////////////////////////////

using namespace std;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



   namespace ShapeFunctions {

//////////////////////////////////////////////////////////////////////////////

CFuint LagrangeShapeFunctionHexaP2_27Nodes::_interpolatorID = 0;
RealMatrix LagrangeShapeFunctionHexaP2_27Nodes::m_invJ(3,3);
RealVector LagrangeShapeFunctionHexaP2_27Nodes::_vec1 = RealVector(DIM_3D);
RealVector LagrangeShapeFunctionHexaP2_27Nodes::_vec2 = RealVector(DIM_3D);
RealVector LagrangeShapeFunctionHexaP2_27Nodes::_vec3 = RealVector(DIM_3D);
std::vector<RealVector> LagrangeShapeFunctionHexaP2_27Nodes::_normals = std::vector<RealVector>(6);
std::vector<RealVector> LagrangeShapeFunctionHexaP2_27Nodes::_gradShapFunc = std::vector<RealVector>(3,RealVector(27));

//////////////////////////////////////////////////////////////////////////////

void LagrangeShapeFunctionHexaP2_27Nodes::computeGradientStates(
         const std::vector<RealMatrix>& jacob,
         const std::vector<RealVector>& mappedCoord,
               std::vector<RealMatrix>& grad)
{
    MathTools::MatrixInverterT<3> inverter;

    for (CFuint ip = 0; ip < mappedCoord.size(); ++ip)
    {
      inverter.invert(jacob[ip], m_invJ);

      RealMatrix& lgrad = grad[ip];

// unused //      const CFreal lambda = 1.0 - mappedCoord[ip].sum();
// unused //      const CFreal xi =  mappedCoord[ip][KSI];
// unused //      const CFreal eta = mappedCoord[ip][ETA];
// unused //      const CFreal zta = mappedCoord[ip][ZTA];
///@todo this funny task is not done yet!!
      const CFreal dN0dxi = 0.;
      const CFreal dN1dxi = 0.;
      const CFreal dN2dxi = 0.;
      const CFreal dN3dxi = 0.;
      const CFreal dN4dxi = 0.;
      const CFreal dN5dxi = 0.;
      const CFreal dN6dxi = 0.;
      const CFreal dN7dxi = 0.;
      const CFreal dN8dxi = 0.;
      const CFreal dN9dxi = 0.;
      const CFreal dN10dxi = 0.;
      const CFreal dN11dxi = 0.;
      const CFreal dN12dxi = 0.;
      const CFreal dN13dxi = 0.;
      const CFreal dN14dxi = 0.;
      const CFreal dN15dxi = 0.;
      const CFreal dN16dxi = 0.;
      const CFreal dN17dxi = 0.;
      const CFreal dN18dxi = 0.;
      const CFreal dN19dxi = 0.;
      const CFreal dN20dxi = 0.;
      const CFreal dN21dxi = 0.;
      const CFreal dN22dxi = 0.;
      const CFreal dN23dxi = 0.;
      const CFreal dN24dxi = 0.;
      const CFreal dN25dxi = 0.;
      const CFreal dN26dxi = 0.;

      const CFreal dN0deta = 0.;
      const CFreal dN1deta = 0.;
      const CFreal dN2deta = 0.;
      const CFreal dN3deta = 0.;
      const CFreal dN4deta = 0.;
      const CFreal dN5deta = 0.;
      const CFreal dN6deta = 0.;
      const CFreal dN7deta = 0.;
      const CFreal dN8deta = 0.;
      const CFreal dN9deta = 0.;
      const CFreal dN10deta = 0.;
      const CFreal dN11deta = 0.;
      const CFreal dN12deta = 0.;
      const CFreal dN13deta = 0.;
      const CFreal dN14deta = 0.;
      const CFreal dN15deta = 0.;
      const CFreal dN16deta = 0.;
      const CFreal dN17deta = 0.;
      const CFreal dN18deta = 0.;
      const CFreal dN19deta = 0.;
      const CFreal dN20deta = 0.;
      const CFreal dN21deta = 0.;
      const CFreal dN22deta = 0.;
      const CFreal dN23deta = 0.;
      const CFreal dN24deta = 0.;
      const CFreal dN25deta = 0.;
      const CFreal dN26deta = 0.;

      const CFreal dN0dzta = 0.;
      const CFreal dN1dzta = 0.;
      const CFreal dN2dzta = 0.;
      const CFreal dN3dzta = 0.;
      const CFreal dN4dzta = 0.;
      const CFreal dN5dzta = 0.;
      const CFreal dN6dzta = 0.;
      const CFreal dN7dzta = 0.;
      const CFreal dN8dzta = 0.;
      const CFreal dN9dzta = 0.;
      const CFreal dN10dzta = 0.;
      const CFreal dN11dzta = 0.;
      const CFreal dN12dzta = 0.;
      const CFreal dN13dzta = 0.;
      const CFreal dN14dzta = 0.;
      const CFreal dN15dzta = 0.;
      const CFreal dN16dzta = 0.;
      const CFreal dN17dzta = 0.;
      const CFreal dN18dzta = 0.;
      const CFreal dN19dzta = 0.;
      const CFreal dN20dzta = 0.;
      const CFreal dN21dzta = 0.;
      const CFreal dN22dzta = 0.;
      const CFreal dN23dzta = 0.;
      const CFreal dN24dzta = 0.;
      const CFreal dN25dzta = 0.;
      const CFreal dN26dzta = 0.;

      const CFreal JXX = m_invJ(0,0);
      const CFreal JXY = m_invJ(0,1);
      const CFreal JXZ = m_invJ(0,2);
      const CFreal JYX = m_invJ(1,0);
      const CFreal JYY = m_invJ(1,1);
      const CFreal JYZ = m_invJ(1,2);
      const CFreal JZX = m_invJ(2,0);
      const CFreal JZY = m_invJ(2,1);
      const CFreal JZZ = m_invJ(2,2);

      lgrad(0,XX) = dN0dxi * JXX + dN0deta * JXY + dN0dzta * JXZ;
      lgrad(0,YY) = dN0dxi * JYX + dN0deta * JYY + dN0dzta * JYZ;
      lgrad(0,ZZ) = dN0dxi * JZX + dN0deta * JZY + dN0dzta * JZZ;

      lgrad(1,XX) = dN1dxi * JXX + dN1deta * JXY + dN1dzta * JXZ;
      lgrad(1,YY) = dN1dxi * JYX + dN1deta * JYY + dN1dzta * JYZ;
      lgrad(1,ZZ) = dN1dxi * JZX + dN1deta * JZY + dN1dzta * JZZ;

      lgrad(2,XX) = dN2dxi * JXX + dN2deta * JXY + dN2dzta * JXZ;
      lgrad(2,YY) = dN2dxi * JYX + dN2deta * JYY + dN2dzta * JYZ;
      lgrad(2,ZZ) = dN2dxi * JZX + dN2deta * JZY + dN2dzta * JZZ;

      lgrad(3,XX) = dN3dxi * JXX + dN3deta * JXY + dN3dzta * JXZ;
      lgrad(3,YY) = dN3dxi * JYX + dN3deta * JYY + dN3dzta * JYZ;
      lgrad(3,ZZ) = dN3dxi * JZX + dN3deta * JZY + dN3dzta * JZZ;

      lgrad(4,XX) = dN4dxi * JXX + dN4deta * JXY + dN4dzta * JXZ;
      lgrad(4,YY) = dN4dxi * JYX + dN4deta * JYY + dN4dzta * JYZ;
      lgrad(4,ZZ) = dN4dxi * JZX + dN4deta * JZY + dN4dzta * JZZ;

      lgrad(5,XX) = dN5dxi * JXX + dN5deta * JXY + dN5dzta * JXZ;
      lgrad(5,YY) = dN5dxi * JYX + dN5deta * JYY + dN5dzta * JYZ;
      lgrad(5,ZZ) = dN5dxi * JZX + dN5deta * JZY + dN5dzta * JZZ;

      lgrad(6,XX) = dN6dxi * JXX + dN6deta * JXY + dN6dzta * JXZ;
      lgrad(6,YY) = dN6dxi * JYX + dN6deta * JYY + dN6dzta * JYZ;
      lgrad(6,ZZ) = dN6dxi * JZX + dN6deta * JZY + dN6dzta * JZZ;

      lgrad(7,XX) = dN7dxi * JXX + dN7deta * JXY + dN7dzta * JXZ;
      lgrad(7,YY) = dN7dxi * JYX + dN7deta * JYY + dN7dzta * JYZ;
      lgrad(7,ZZ) = dN7dxi * JZX + dN7deta * JZY + dN7dzta * JZZ;

      lgrad(8,XX) = dN8dxi * JXX + dN8deta * JXY + dN8dzta * JXZ;
      lgrad(8,YY) = dN8dxi * JYX + dN8deta * JYY + dN8dzta * JYZ;
      lgrad(8,ZZ) = dN8dxi * JZX + dN8deta * JZY + dN8dzta * JZZ;

      lgrad(9,XX) = dN9dxi * JXX + dN9deta * JXY + dN9dzta * JXZ;
      lgrad(9,YY) = dN9dxi * JYX + dN9deta * JYY + dN9dzta * JYZ;
      lgrad(9,ZZ) = dN9dxi * JZX + dN9deta * JZY + dN9dzta * JZZ;

      lgrad(10,XX) = dN10dxi * JXX + dN10deta * JXY + dN10dzta * JXZ;
      lgrad(10,YY) = dN10dxi * JYX + dN10deta * JYY + dN10dzta * JYZ;
      lgrad(10,ZZ) = dN10dxi * JZX + dN10deta * JZY + dN10dzta * JZZ;

      lgrad(11,XX) = dN11dxi * JXX + dN11deta * JXY + dN11dzta * JXZ;
      lgrad(11,YY) = dN11dxi * JYX + dN11deta * JYY + dN11dzta * JYZ;
      lgrad(11,ZZ) = dN11dxi * JZX + dN11deta * JZY + dN11dzta * JZZ;

      lgrad(12,XX) = dN12dxi * JXX + dN12deta * JXY + dN12dzta * JXZ;
      lgrad(12,YY) = dN12dxi * JYX + dN12deta * JYY + dN12dzta * JYZ;
      lgrad(12,ZZ) = dN12dxi * JZX + dN12deta * JZY + dN12dzta * JZZ;

      lgrad(13,XX) = dN13dxi * JXX + dN13deta * JXY + dN13dzta * JXZ;
      lgrad(13,YY) = dN13dxi * JYX + dN13deta * JYY + dN13dzta * JYZ;
      lgrad(13,ZZ) = dN13dxi * JZX + dN13deta * JZY + dN13dzta * JZZ;

      lgrad(14,XX) = dN14dxi * JXX + dN14deta * JXY + dN14dzta * JXZ;
      lgrad(14,YY) = dN14dxi * JYX + dN14deta * JYY + dN14dzta * JYZ;
      lgrad(14,ZZ) = dN14dxi * JZX + dN14deta * JZY + dN14dzta * JZZ;

      lgrad(15,XX) = dN15dxi * JXX + dN15deta * JXY + dN15dzta * JXZ;
      lgrad(15,YY) = dN15dxi * JYX + dN15deta * JYY + dN15dzta * JYZ;
      lgrad(15,ZZ) = dN15dxi * JZX + dN15deta * JZY + dN15dzta * JZZ;

      lgrad(16,XX) = dN16dxi * JXX + dN16deta * JXY + dN16dzta * JXZ;
      lgrad(16,YY) = dN16dxi * JYX + dN16deta * JYY + dN16dzta * JYZ;
      lgrad(16,ZZ) = dN16dxi * JZX + dN16deta * JZY + dN16dzta * JZZ;

      lgrad(17,XX) = dN17dxi * JXX + dN17deta * JXY + dN17dzta * JXZ;
      lgrad(17,YY) = dN17dxi * JYX + dN17deta * JYY + dN17dzta * JYZ;
      lgrad(17,ZZ) = dN17dxi * JZX + dN17deta * JZY + dN17dzta * JZZ;

      lgrad(18,XX) = dN18dxi * JXX + dN18deta * JXY + dN18dzta * JXZ;
      lgrad(18,YY) = dN18dxi * JYX + dN18deta * JYY + dN18dzta * JYZ;
      lgrad(18,ZZ) = dN18dxi * JZX + dN18deta * JZY + dN18dzta * JZZ;

      lgrad(19,XX) = dN19dxi * JXX + dN19deta * JXY + dN19dzta * JXZ;
      lgrad(19,YY) = dN19dxi * JYX + dN19deta * JYY + dN19dzta * JYZ;
      lgrad(19,ZZ) = dN19dxi * JZX + dN19deta * JZY + dN19dzta * JZZ;

      lgrad(20,XX) = dN20dxi * JXX + dN20deta * JXY + dN20dzta * JXZ;
      lgrad(20,YY) = dN20dxi * JYX + dN20deta * JYY + dN20dzta * JYZ;
      lgrad(20,ZZ) = dN20dxi * JZX + dN20deta * JZY + dN20dzta * JZZ;

      lgrad(21,XX) = dN21dxi * JXX + dN21deta * JXY + dN21dzta * JXZ;
      lgrad(21,YY) = dN21dxi * JYX + dN21deta * JYY + dN21dzta * JYZ;
      lgrad(21,ZZ) = dN21dxi * JZX + dN21deta * JZY + dN21dzta * JZZ;

      lgrad(22,XX) = dN22dxi * JXX + dN22deta * JXY + dN22dzta * JXZ;
      lgrad(22,YY) = dN22dxi * JYX + dN22deta * JYY + dN22dzta * JYZ;
      lgrad(22,ZZ) = dN22dxi * JZX + dN22deta * JZY + dN22dzta * JZZ;

      lgrad(23,XX) = dN23dxi * JXX + dN23deta * JXY + dN23dzta * JXZ;
      lgrad(23,YY) = dN23dxi * JYX + dN23deta * JYY + dN23dzta * JYZ;
      lgrad(23,ZZ) = dN23dxi * JZX + dN23deta * JZY + dN23dzta * JZZ;

      lgrad(24,XX) = dN24dxi * JXX + dN24deta * JXY + dN24dzta * JXZ;
      lgrad(24,YY) = dN24dxi * JYX + dN24deta * JYY + dN24dzta * JYZ;
      lgrad(24,ZZ) = dN24dxi * JZX + dN24deta * JZY + dN24dzta * JZZ;

      lgrad(25,XX) = dN25dxi * JXX + dN25deta * JXY + dN25dzta * JXZ;
      lgrad(25,YY) = dN25dxi * JYX + dN25deta * JYY + dN25dzta * JYZ;
      lgrad(25,ZZ) = dN25dxi * JZX + dN25deta * JZY + dN25dzta * JZZ;

      lgrad(26,XX) = dN26dxi * JXX + dN26deta * JXY + dN26dzta * JXZ;
      lgrad(26,YY) = dN26dxi * JYX + dN26deta * JYY + dN26dzta * JYZ;
      lgrad(26,ZZ) = dN26dxi * JZX + dN26deta * JZY + dN26dzta * JZZ;

    }
}

//////////////////////////////////////////////////////////////////////////////

void LagrangeShapeFunctionHexaP2_27Nodes::computeMappedCoordPlaneNormal(const std::vector<CFuint>& planeIdx,
                                                                        const std::vector<RealVector>& mappedCoord,
                                                                        const std::vector<Framework::Node*>& nodes,
                                                                        std::vector<RealVector>& normal)
{
  for (CFuint ip = 0; ip < mappedCoord.size(); ++ip)
  {
    RealVector& pointNormal = normal[ip];

    const CFreal ksi = mappedCoord[ip][KSI];
    const CFreal eta = mappedCoord[ip][ETA];
    const CFreal zta = mappedCoord[ip][ZTA];
    const CFreal ksi2 = ksi*ksi;
    const CFreal eta2 = eta*eta;
    const CFreal zta2 = zta*zta;

    const CFreal Nksi_m1 = -ksi * (1. - ksi);
    const CFreal Nksi_0  =   2. * (1. - ksi2);
    const CFreal Nksi_p1 =  ksi * (1. + ksi);

    const CFreal Neta_m1 = -eta * (1. - eta);
    const CFreal Neta_0  =   2. * (1. - eta2);
    const CFreal Neta_p1 =  eta * (1. + eta);

    const CFreal Nzeta_m1 = -zta * (1. - zta);
    const CFreal Nzeta_0  =   2. * (1. - zta2);
    const CFreal Nzeta_p1 =  zta * (1. + zta);

    cf_assert(planeIdx[ip] == 0 || planeIdx[ip] == 1 || planeIdx[ip] == 2);
    if (planeIdx[ip] == 0)
    {
      const CFreal dNeta_m1deta =  2. * eta - 1.;
      const CFreal dNeta_0deta  = -4. * eta;
      const CFreal dNeta_p1deta =  2. * eta + 1.;

      const CFreal dNzeta_m1dzeta =  2. * zta - 1.;
      const CFreal dNzeta_0dzeta  = -4. * zta;
      const CFreal dNzeta_p1dzeta =  2. * zta + 1.;

      _gradShapFunc[ETA][ 0] = Nksi_m1 * dNeta_m1deta * Nzeta_m1;
      _gradShapFunc[ETA][ 1] = Nksi_p1 * dNeta_m1deta * Nzeta_m1;
      _gradShapFunc[ETA][ 2] = Nksi_p1 * dNeta_p1deta * Nzeta_m1;
      _gradShapFunc[ETA][ 3] = Nksi_m1 * dNeta_p1deta * Nzeta_m1;
      _gradShapFunc[ETA][ 4] = Nksi_m1 * dNeta_m1deta * Nzeta_p1;
      _gradShapFunc[ETA][ 5] = Nksi_p1 * dNeta_m1deta * Nzeta_p1;
      _gradShapFunc[ETA][ 6] = Nksi_p1 * dNeta_p1deta * Nzeta_p1;
      _gradShapFunc[ETA][ 7] = Nksi_m1 * dNeta_p1deta * Nzeta_p1;
      _gradShapFunc[ETA][ 8] = Nksi_0  * dNeta_m1deta * Nzeta_m1;
      _gradShapFunc[ETA][ 9] = Nksi_p1 * dNeta_0deta  * Nzeta_m1;
      _gradShapFunc[ETA][10] = Nksi_0  * dNeta_p1deta * Nzeta_m1;
      _gradShapFunc[ETA][11] = Nksi_m1 * dNeta_0deta  * Nzeta_m1;
      _gradShapFunc[ETA][12] = Nksi_0  * dNeta_0deta  * Nzeta_m1;
      _gradShapFunc[ETA][13] = Nksi_m1 * dNeta_m1deta * Nzeta_0 ;
      _gradShapFunc[ETA][14] = Nksi_0  * dNeta_m1deta * Nzeta_0 ;
      _gradShapFunc[ETA][15] = Nksi_p1 * dNeta_m1deta * Nzeta_0 ;
      _gradShapFunc[ETA][16] = Nksi_p1 * dNeta_0deta  * Nzeta_0 ;
      _gradShapFunc[ETA][17] = Nksi_p1 * dNeta_p1deta * Nzeta_0 ;
      _gradShapFunc[ETA][18] = Nksi_0  * dNeta_p1deta * Nzeta_0 ;
      _gradShapFunc[ETA][19] = Nksi_m1 * dNeta_p1deta * Nzeta_0 ;
      _gradShapFunc[ETA][20] = Nksi_m1 * dNeta_0deta  * Nzeta_0 ;
      _gradShapFunc[ETA][21] = Nksi_0  * dNeta_0deta  * Nzeta_0 ;
      _gradShapFunc[ETA][22] = Nksi_0  * dNeta_m1deta * Nzeta_p1;
      _gradShapFunc[ETA][23] = Nksi_p1 * dNeta_0deta  * Nzeta_p1;
      _gradShapFunc[ETA][24] = Nksi_0  * dNeta_p1deta * Nzeta_p1;
      _gradShapFunc[ETA][25] = Nksi_m1 * dNeta_0deta  * Nzeta_p1;
      _gradShapFunc[ETA][26] = Nksi_0  * dNeta_0deta  * Nzeta_p1;

      _gradShapFunc[ZTA][ 0] = Nksi_m1 * Neta_m1 * dNzeta_m1dzeta;
      _gradShapFunc[ZTA][ 1] = Nksi_p1 * Neta_m1 * dNzeta_m1dzeta;
      _gradShapFunc[ZTA][ 2] = Nksi_p1 * Neta_p1 * dNzeta_m1dzeta;
      _gradShapFunc[ZTA][ 3] = Nksi_m1 * Neta_p1 * dNzeta_m1dzeta;
      _gradShapFunc[ZTA][ 4] = Nksi_m1 * Neta_m1 * dNzeta_p1dzeta;
      _gradShapFunc[ZTA][ 5] = Nksi_p1 * Neta_m1 * dNzeta_p1dzeta;
      _gradShapFunc[ZTA][ 6] = Nksi_p1 * Neta_p1 * dNzeta_p1dzeta;
      _gradShapFunc[ZTA][ 7] = Nksi_m1 * Neta_p1 * dNzeta_p1dzeta;
      _gradShapFunc[ZTA][ 8] = Nksi_0  * Neta_m1 * dNzeta_m1dzeta;
      _gradShapFunc[ZTA][ 9] = Nksi_p1 * Neta_0  * dNzeta_m1dzeta;
      _gradShapFunc[ZTA][10] = Nksi_0  * Neta_p1 * dNzeta_m1dzeta;
      _gradShapFunc[ZTA][11] = Nksi_m1 * Neta_0  * dNzeta_m1dzeta;
      _gradShapFunc[ZTA][12] = Nksi_0  * Neta_0  * dNzeta_m1dzeta;
      _gradShapFunc[ZTA][13] = Nksi_m1 * Neta_m1 * dNzeta_0dzeta ;
      _gradShapFunc[ZTA][14] = Nksi_0  * Neta_m1 * dNzeta_0dzeta ;
      _gradShapFunc[ZTA][15] = Nksi_p1 * Neta_m1 * dNzeta_0dzeta ;
      _gradShapFunc[ZTA][16] = Nksi_p1 * Neta_0  * dNzeta_0dzeta ;
      _gradShapFunc[ZTA][17] = Nksi_p1 * Neta_p1 * dNzeta_0dzeta ;
      _gradShapFunc[ZTA][18] = Nksi_0  * Neta_p1 * dNzeta_0dzeta ;
      _gradShapFunc[ZTA][19] = Nksi_m1 * Neta_p1 * dNzeta_0dzeta ;
      _gradShapFunc[ZTA][20] = Nksi_m1 * Neta_0  * dNzeta_0dzeta ;
      _gradShapFunc[ZTA][21] = Nksi_0  * Neta_0  * dNzeta_0dzeta ;
      _gradShapFunc[ZTA][22] = Nksi_0  * Neta_m1 * dNzeta_p1dzeta;
      _gradShapFunc[ZTA][23] = Nksi_p1 * Neta_0  * dNzeta_p1dzeta;
      _gradShapFunc[ZTA][24] = Nksi_0  * Neta_p1 * dNzeta_p1dzeta;
      _gradShapFunc[ZTA][25] = Nksi_m1 * Neta_0  * dNzeta_p1dzeta;
      _gradShapFunc[ZTA][26] = Nksi_0  * Neta_0  * dNzeta_p1dzeta;

      _vec1 = _gradShapFunc[ETA][0]*(*nodes[0]);
      _vec2 = _gradShapFunc[ZTA][0]*(*nodes[0]);
      for (CFuint in = 1; in < 27; ++in)
      {
        _vec1 += _gradShapFunc[ETA][in]*(*nodes[in]);
        _vec2 += _gradShapFunc[ZTA][in]*(*nodes[in]);
      }
    }
    else if (planeIdx[ip] == 1)
    {
      const CFreal dNzeta_m1dzeta =  2. * zta - 1.;
      const CFreal dNzeta_0dzeta  = -4. * zta;
      const CFreal dNzeta_p1dzeta =  2. * zta + 1.;

      const CFreal dNksi_m1dksi =  2. * ksi - 1.;
      const CFreal dNksi_0dksi  = -4. * ksi;
      const CFreal dNksi_p1dksi =  2. * ksi + 1.;

      _gradShapFunc[ZTA][ 0] = Nksi_m1 * Neta_m1 * dNzeta_m1dzeta;
      _gradShapFunc[ZTA][ 1] = Nksi_p1 * Neta_m1 * dNzeta_m1dzeta;
      _gradShapFunc[ZTA][ 2] = Nksi_p1 * Neta_p1 * dNzeta_m1dzeta;
      _gradShapFunc[ZTA][ 3] = Nksi_m1 * Neta_p1 * dNzeta_m1dzeta;
      _gradShapFunc[ZTA][ 4] = Nksi_m1 * Neta_m1 * dNzeta_p1dzeta;
      _gradShapFunc[ZTA][ 5] = Nksi_p1 * Neta_m1 * dNzeta_p1dzeta;
      _gradShapFunc[ZTA][ 6] = Nksi_p1 * Neta_p1 * dNzeta_p1dzeta;
      _gradShapFunc[ZTA][ 7] = Nksi_m1 * Neta_p1 * dNzeta_p1dzeta;
      _gradShapFunc[ZTA][ 8] = Nksi_0  * Neta_m1 * dNzeta_m1dzeta;
      _gradShapFunc[ZTA][ 9] = Nksi_p1 * Neta_0  * dNzeta_m1dzeta;
      _gradShapFunc[ZTA][10] = Nksi_0  * Neta_p1 * dNzeta_m1dzeta;
      _gradShapFunc[ZTA][11] = Nksi_m1 * Neta_0  * dNzeta_m1dzeta;
      _gradShapFunc[ZTA][12] = Nksi_0  * Neta_0  * dNzeta_m1dzeta;
      _gradShapFunc[ZTA][13] = Nksi_m1 * Neta_m1 * dNzeta_0dzeta ;
      _gradShapFunc[ZTA][14] = Nksi_0  * Neta_m1 * dNzeta_0dzeta ;
      _gradShapFunc[ZTA][15] = Nksi_p1 * Neta_m1 * dNzeta_0dzeta ;
      _gradShapFunc[ZTA][16] = Nksi_p1 * Neta_0  * dNzeta_0dzeta ;
      _gradShapFunc[ZTA][17] = Nksi_p1 * Neta_p1 * dNzeta_0dzeta ;
      _gradShapFunc[ZTA][18] = Nksi_0  * Neta_p1 * dNzeta_0dzeta ;
      _gradShapFunc[ZTA][19] = Nksi_m1 * Neta_p1 * dNzeta_0dzeta ;
      _gradShapFunc[ZTA][20] = Nksi_m1 * Neta_0  * dNzeta_0dzeta ;
      _gradShapFunc[ZTA][21] = Nksi_0  * Neta_0  * dNzeta_0dzeta ;
      _gradShapFunc[ZTA][22] = Nksi_0  * Neta_m1 * dNzeta_p1dzeta;
      _gradShapFunc[ZTA][23] = Nksi_p1 * Neta_0  * dNzeta_p1dzeta;
      _gradShapFunc[ZTA][24] = Nksi_0  * Neta_p1 * dNzeta_p1dzeta;
      _gradShapFunc[ZTA][25] = Nksi_m1 * Neta_0  * dNzeta_p1dzeta;
      _gradShapFunc[ZTA][26] = Nksi_0  * Neta_0  * dNzeta_p1dzeta;

      _gradShapFunc[KSI][ 0] = dNksi_m1dksi * Neta_m1 * Nzeta_m1;
      _gradShapFunc[KSI][ 1] = dNksi_p1dksi * Neta_m1 * Nzeta_m1;
      _gradShapFunc[KSI][ 2] = dNksi_p1dksi * Neta_p1 * Nzeta_m1;
      _gradShapFunc[KSI][ 3] = dNksi_m1dksi * Neta_p1 * Nzeta_m1;
      _gradShapFunc[KSI][ 4] = dNksi_m1dksi * Neta_m1 * Nzeta_p1;
      _gradShapFunc[KSI][ 5] = dNksi_p1dksi * Neta_m1 * Nzeta_p1;
      _gradShapFunc[KSI][ 6] = dNksi_p1dksi * Neta_p1 * Nzeta_p1;
      _gradShapFunc[KSI][ 7] = dNksi_m1dksi * Neta_p1 * Nzeta_p1;
      _gradShapFunc[KSI][ 8] = dNksi_0dksi  * Neta_m1 * Nzeta_m1;
      _gradShapFunc[KSI][ 9] = dNksi_p1dksi * Neta_0  * Nzeta_m1;
      _gradShapFunc[KSI][10] = dNksi_0dksi  * Neta_p1 * Nzeta_m1;
      _gradShapFunc[KSI][11] = dNksi_m1dksi * Neta_0  * Nzeta_m1;
      _gradShapFunc[KSI][12] = dNksi_0dksi  * Neta_0  * Nzeta_m1;
      _gradShapFunc[KSI][13] = dNksi_m1dksi * Neta_m1 * Nzeta_0 ;
      _gradShapFunc[KSI][14] = dNksi_0dksi  * Neta_m1 * Nzeta_0 ;
      _gradShapFunc[KSI][15] = dNksi_p1dksi * Neta_m1 * Nzeta_0 ;
      _gradShapFunc[KSI][16] = dNksi_p1dksi * Neta_0  * Nzeta_0 ;
      _gradShapFunc[KSI][17] = dNksi_p1dksi * Neta_p1 * Nzeta_0 ;
      _gradShapFunc[KSI][18] = dNksi_0dksi  * Neta_p1 * Nzeta_0 ;
      _gradShapFunc[KSI][19] = dNksi_m1dksi * Neta_p1 * Nzeta_0 ;
      _gradShapFunc[KSI][20] = dNksi_m1dksi * Neta_0  * Nzeta_0 ;
      _gradShapFunc[KSI][21] = dNksi_0dksi  * Neta_0  * Nzeta_0 ;
      _gradShapFunc[KSI][22] = dNksi_0dksi  * Neta_m1 * Nzeta_p1;
      _gradShapFunc[KSI][23] = dNksi_p1dksi * Neta_0  * Nzeta_p1;
      _gradShapFunc[KSI][24] = dNksi_0dksi  * Neta_p1 * Nzeta_p1;
      _gradShapFunc[KSI][25] = dNksi_m1dksi * Neta_0  * Nzeta_p1;
      _gradShapFunc[KSI][26] = dNksi_0dksi  * Neta_0  * Nzeta_p1;

      _vec1 = _gradShapFunc[ZTA][0]*(*nodes[0]);
      _vec2 = _gradShapFunc[KSI][0]*(*nodes[0]);
      for (CFuint in = 1; in < 27; ++in)
      {
        _vec1 += _gradShapFunc[ZTA][in]*(*nodes[in]);
        _vec2 += _gradShapFunc[KSI][in]*(*nodes[in]);
      }
    }
    else
    {
      const CFreal dNksi_m1dksi =  2. * ksi - 1.;
      const CFreal dNksi_0dksi  = -4. * ksi;
      const CFreal dNksi_p1dksi =  2. * ksi + 1.;

      const CFreal dNeta_m1deta =  2. * eta - 1.;
      const CFreal dNeta_0deta  = -4. * eta;
      const CFreal dNeta_p1deta =  2. * eta + 1.;

      _gradShapFunc[KSI][ 0] = dNksi_m1dksi * Neta_m1 * Nzeta_m1;
      _gradShapFunc[KSI][ 1] = dNksi_p1dksi * Neta_m1 * Nzeta_m1;
      _gradShapFunc[KSI][ 2] = dNksi_p1dksi * Neta_p1 * Nzeta_m1;
      _gradShapFunc[KSI][ 3] = dNksi_m1dksi * Neta_p1 * Nzeta_m1;
      _gradShapFunc[KSI][ 4] = dNksi_m1dksi * Neta_m1 * Nzeta_p1;
      _gradShapFunc[KSI][ 5] = dNksi_p1dksi * Neta_m1 * Nzeta_p1;
      _gradShapFunc[KSI][ 6] = dNksi_p1dksi * Neta_p1 * Nzeta_p1;
      _gradShapFunc[KSI][ 7] = dNksi_m1dksi * Neta_p1 * Nzeta_p1;
      _gradShapFunc[KSI][ 8] = dNksi_0dksi  * Neta_m1 * Nzeta_m1;
      _gradShapFunc[KSI][ 9] = dNksi_p1dksi * Neta_0  * Nzeta_m1;
      _gradShapFunc[KSI][10] = dNksi_0dksi  * Neta_p1 * Nzeta_m1;
      _gradShapFunc[KSI][11] = dNksi_m1dksi * Neta_0  * Nzeta_m1;
      _gradShapFunc[KSI][12] = dNksi_0dksi  * Neta_0  * Nzeta_m1;
      _gradShapFunc[KSI][13] = dNksi_m1dksi * Neta_m1 * Nzeta_0 ;
      _gradShapFunc[KSI][14] = dNksi_0dksi  * Neta_m1 * Nzeta_0 ;
      _gradShapFunc[KSI][15] = dNksi_p1dksi * Neta_m1 * Nzeta_0 ;
      _gradShapFunc[KSI][16] = dNksi_p1dksi * Neta_0  * Nzeta_0 ;
      _gradShapFunc[KSI][17] = dNksi_p1dksi * Neta_p1 * Nzeta_0 ;
      _gradShapFunc[KSI][18] = dNksi_0dksi  * Neta_p1 * Nzeta_0 ;
      _gradShapFunc[KSI][19] = dNksi_m1dksi * Neta_p1 * Nzeta_0 ;
      _gradShapFunc[KSI][20] = dNksi_m1dksi * Neta_0  * Nzeta_0 ;
      _gradShapFunc[KSI][21] = dNksi_0dksi  * Neta_0  * Nzeta_0 ;
      _gradShapFunc[KSI][22] = dNksi_0dksi  * Neta_m1 * Nzeta_p1;
      _gradShapFunc[KSI][23] = dNksi_p1dksi * Neta_0  * Nzeta_p1;
      _gradShapFunc[KSI][24] = dNksi_0dksi  * Neta_p1 * Nzeta_p1;
      _gradShapFunc[KSI][25] = dNksi_m1dksi * Neta_0  * Nzeta_p1;
      _gradShapFunc[KSI][26] = dNksi_0dksi  * Neta_0  * Nzeta_p1;

      _gradShapFunc[ETA][ 0] = Nksi_m1 * dNeta_m1deta * Nzeta_m1;
      _gradShapFunc[ETA][ 1] = Nksi_p1 * dNeta_m1deta * Nzeta_m1;
      _gradShapFunc[ETA][ 2] = Nksi_p1 * dNeta_p1deta * Nzeta_m1;
      _gradShapFunc[ETA][ 3] = Nksi_m1 * dNeta_p1deta * Nzeta_m1;
      _gradShapFunc[ETA][ 4] = Nksi_m1 * dNeta_m1deta * Nzeta_p1;
      _gradShapFunc[ETA][ 5] = Nksi_p1 * dNeta_m1deta * Nzeta_p1;
      _gradShapFunc[ETA][ 6] = Nksi_p1 * dNeta_p1deta * Nzeta_p1;
      _gradShapFunc[ETA][ 7] = Nksi_m1 * dNeta_p1deta * Nzeta_p1;
      _gradShapFunc[ETA][ 8] = Nksi_0  * dNeta_m1deta * Nzeta_m1;
      _gradShapFunc[ETA][ 9] = Nksi_p1 * dNeta_0deta  * Nzeta_m1;
      _gradShapFunc[ETA][10] = Nksi_0  * dNeta_p1deta * Nzeta_m1;
      _gradShapFunc[ETA][11] = Nksi_m1 * dNeta_0deta  * Nzeta_m1;
      _gradShapFunc[ETA][12] = Nksi_0  * dNeta_0deta  * Nzeta_m1;
      _gradShapFunc[ETA][13] = Nksi_m1 * dNeta_m1deta * Nzeta_0 ;
      _gradShapFunc[ETA][14] = Nksi_0  * dNeta_m1deta * Nzeta_0 ;
      _gradShapFunc[ETA][15] = Nksi_p1 * dNeta_m1deta * Nzeta_0 ;
      _gradShapFunc[ETA][16] = Nksi_p1 * dNeta_0deta  * Nzeta_0 ;
      _gradShapFunc[ETA][17] = Nksi_p1 * dNeta_p1deta * Nzeta_0 ;
      _gradShapFunc[ETA][18] = Nksi_0  * dNeta_p1deta * Nzeta_0 ;
      _gradShapFunc[ETA][19] = Nksi_m1 * dNeta_p1deta * Nzeta_0 ;
      _gradShapFunc[ETA][20] = Nksi_m1 * dNeta_0deta  * Nzeta_0 ;
      _gradShapFunc[ETA][21] = Nksi_0  * dNeta_0deta  * Nzeta_0 ;
      _gradShapFunc[ETA][22] = Nksi_0  * dNeta_m1deta * Nzeta_p1;
      _gradShapFunc[ETA][23] = Nksi_p1 * dNeta_0deta  * Nzeta_p1;
      _gradShapFunc[ETA][24] = Nksi_0  * dNeta_p1deta * Nzeta_p1;
      _gradShapFunc[ETA][25] = Nksi_m1 * dNeta_0deta  * Nzeta_p1;
      _gradShapFunc[ETA][26] = Nksi_0  * dNeta_0deta  * Nzeta_p1;

      _vec1 = _gradShapFunc[KSI][0]*(*nodes[0]);
      _vec2 = _gradShapFunc[ETA][0]*(*nodes[0]);
      for (CFuint in = 1; in < 27; ++in)
      {
        _vec1 += _gradShapFunc[KSI][in]*(*nodes[in]);
        _vec2 += _gradShapFunc[ETA][in]*(*nodes[in]);
      }
    }

    // compute normal
    MathTools::MathFunctions::crossProd(_vec1,_vec2,pointNormal);
    pointNormal *= 0.015625;
  }
}

//////////////////////////////////////////////////////////////////////////////

void LagrangeShapeFunctionHexaP2_27Nodes::computeJacobian(const std::vector<Framework::Node*>& nodes,
                                                          const std::vector<RealVector>& mappedCoord,
                                                          std::vector<RealMatrix>& jacob)
{
  for (CFuint ip = 0; ip < mappedCoord.size(); ++ip) {

    RealMatrix& pointJacob = jacob[ip];

    const CFreal ksi = mappedCoord[ip][KSI];
    const CFreal eta = mappedCoord[ip][ETA];
    const CFreal zta = mappedCoord[ip][ZTA];
    const CFreal ksi2 = ksi*ksi;
    const CFreal eta2 = eta*eta;
    const CFreal zta2 = zta*zta;

    const CFreal Nksi_m1 = -ksi * (1. - ksi);
    const CFreal Nksi_0  =   2. * (1. - ksi2);
    const CFreal Nksi_p1 =  ksi * (1. + ksi);

    const CFreal Neta_m1 = -eta * (1. - eta);
    const CFreal Neta_0  =   2. * (1. - eta2);
    const CFreal Neta_p1 =  eta * (1. + eta);

    const CFreal Nzeta_m1 = -zta * (1. - zta);
    const CFreal Nzeta_0  =   2. * (1. - zta2);
    const CFreal Nzeta_p1 =  zta * (1. + zta);

    const CFreal dNksi_m1dksi =  2. * ksi - 1.;
    const CFreal dNksi_0dksi  = -4. * ksi;
    const CFreal dNksi_p1dksi =  2. * ksi + 1.;

    const CFreal dNeta_m1deta =  2. * eta - 1.;
    const CFreal dNeta_0deta  = -4. * eta;
    const CFreal dNeta_p1deta =  2. * eta + 1.;

    const CFreal dNzeta_m1dzeta =  2. * zta - 1.;
    const CFreal dNzeta_0dzeta  = -4. * zta;
    const CFreal dNzeta_p1dzeta =  2. * zta + 1.;

    _gradShapFunc[KSI][ 0] = dNksi_m1dksi * Neta_m1 * Nzeta_m1;
    _gradShapFunc[KSI][ 1] = dNksi_p1dksi * Neta_m1 * Nzeta_m1;
    _gradShapFunc[KSI][ 2] = dNksi_p1dksi * Neta_p1 * Nzeta_m1;
    _gradShapFunc[KSI][ 3] = dNksi_m1dksi * Neta_p1 * Nzeta_m1;
    _gradShapFunc[KSI][ 4] = dNksi_m1dksi * Neta_m1 * Nzeta_p1;
    _gradShapFunc[KSI][ 5] = dNksi_p1dksi * Neta_m1 * Nzeta_p1;
    _gradShapFunc[KSI][ 6] = dNksi_p1dksi * Neta_p1 * Nzeta_p1;
    _gradShapFunc[KSI][ 7] = dNksi_m1dksi * Neta_p1 * Nzeta_p1;
    _gradShapFunc[KSI][ 8] = dNksi_0dksi  * Neta_m1 * Nzeta_m1;
    _gradShapFunc[KSI][ 9] = dNksi_p1dksi * Neta_0  * Nzeta_m1;
    _gradShapFunc[KSI][10] = dNksi_0dksi  * Neta_p1 * Nzeta_m1;
    _gradShapFunc[KSI][11] = dNksi_m1dksi * Neta_0  * Nzeta_m1;
    _gradShapFunc[KSI][12] = dNksi_0dksi  * Neta_0  * Nzeta_m1;
    _gradShapFunc[KSI][13] = dNksi_m1dksi * Neta_m1 * Nzeta_0 ;
    _gradShapFunc[KSI][14] = dNksi_0dksi  * Neta_m1 * Nzeta_0 ;
    _gradShapFunc[KSI][15] = dNksi_p1dksi * Neta_m1 * Nzeta_0 ;
    _gradShapFunc[KSI][16] = dNksi_p1dksi * Neta_0  * Nzeta_0 ;
    _gradShapFunc[KSI][17] = dNksi_p1dksi * Neta_p1 * Nzeta_0 ;
    _gradShapFunc[KSI][18] = dNksi_0dksi  * Neta_p1 * Nzeta_0 ;
    _gradShapFunc[KSI][19] = dNksi_m1dksi * Neta_p1 * Nzeta_0 ;
    _gradShapFunc[KSI][20] = dNksi_m1dksi * Neta_0  * Nzeta_0 ;
    _gradShapFunc[KSI][21] = dNksi_0dksi  * Neta_0  * Nzeta_0 ;
    _gradShapFunc[KSI][22] = dNksi_0dksi  * Neta_m1 * Nzeta_p1;
    _gradShapFunc[KSI][23] = dNksi_p1dksi * Neta_0  * Nzeta_p1;
    _gradShapFunc[KSI][24] = dNksi_0dksi  * Neta_p1 * Nzeta_p1;
    _gradShapFunc[KSI][25] = dNksi_m1dksi * Neta_0  * Nzeta_p1;
    _gradShapFunc[KSI][26] = dNksi_0dksi  * Neta_0  * Nzeta_p1;

    _gradShapFunc[ETA][ 0] = Nksi_m1 * dNeta_m1deta * Nzeta_m1;
    _gradShapFunc[ETA][ 1] = Nksi_p1 * dNeta_m1deta * Nzeta_m1;
    _gradShapFunc[ETA][ 2] = Nksi_p1 * dNeta_p1deta * Nzeta_m1;
    _gradShapFunc[ETA][ 3] = Nksi_m1 * dNeta_p1deta * Nzeta_m1;
    _gradShapFunc[ETA][ 4] = Nksi_m1 * dNeta_m1deta * Nzeta_p1;
    _gradShapFunc[ETA][ 5] = Nksi_p1 * dNeta_m1deta * Nzeta_p1;
    _gradShapFunc[ETA][ 6] = Nksi_p1 * dNeta_p1deta * Nzeta_p1;
    _gradShapFunc[ETA][ 7] = Nksi_m1 * dNeta_p1deta * Nzeta_p1;
    _gradShapFunc[ETA][ 8] = Nksi_0  * dNeta_m1deta * Nzeta_m1;
    _gradShapFunc[ETA][ 9] = Nksi_p1 * dNeta_0deta  * Nzeta_m1;
    _gradShapFunc[ETA][10] = Nksi_0  * dNeta_p1deta * Nzeta_m1;
    _gradShapFunc[ETA][11] = Nksi_m1 * dNeta_0deta  * Nzeta_m1;
    _gradShapFunc[ETA][12] = Nksi_0  * dNeta_0deta  * Nzeta_m1;
    _gradShapFunc[ETA][13] = Nksi_m1 * dNeta_m1deta * Nzeta_0 ;
    _gradShapFunc[ETA][14] = Nksi_0  * dNeta_m1deta * Nzeta_0 ;
    _gradShapFunc[ETA][15] = Nksi_p1 * dNeta_m1deta * Nzeta_0 ;
    _gradShapFunc[ETA][16] = Nksi_p1 * dNeta_0deta  * Nzeta_0 ;
    _gradShapFunc[ETA][17] = Nksi_p1 * dNeta_p1deta * Nzeta_0 ;
    _gradShapFunc[ETA][18] = Nksi_0  * dNeta_p1deta * Nzeta_0 ;
    _gradShapFunc[ETA][19] = Nksi_m1 * dNeta_p1deta * Nzeta_0 ;
    _gradShapFunc[ETA][20] = Nksi_m1 * dNeta_0deta  * Nzeta_0 ;
    _gradShapFunc[ETA][21] = Nksi_0  * dNeta_0deta  * Nzeta_0 ;
    _gradShapFunc[ETA][22] = Nksi_0  * dNeta_m1deta * Nzeta_p1;
    _gradShapFunc[ETA][23] = Nksi_p1 * dNeta_0deta  * Nzeta_p1;
    _gradShapFunc[ETA][24] = Nksi_0  * dNeta_p1deta * Nzeta_p1;
    _gradShapFunc[ETA][25] = Nksi_m1 * dNeta_0deta  * Nzeta_p1;
    _gradShapFunc[ETA][26] = Nksi_0  * dNeta_0deta  * Nzeta_p1;

    _gradShapFunc[ZTA][ 0] = Nksi_m1 * Neta_m1 * dNzeta_m1dzeta;
    _gradShapFunc[ZTA][ 1] = Nksi_p1 * Neta_m1 * dNzeta_m1dzeta;
    _gradShapFunc[ZTA][ 2] = Nksi_p1 * Neta_p1 * dNzeta_m1dzeta;
    _gradShapFunc[ZTA][ 3] = Nksi_m1 * Neta_p1 * dNzeta_m1dzeta;
    _gradShapFunc[ZTA][ 4] = Nksi_m1 * Neta_m1 * dNzeta_p1dzeta;
    _gradShapFunc[ZTA][ 5] = Nksi_p1 * Neta_m1 * dNzeta_p1dzeta;
    _gradShapFunc[ZTA][ 6] = Nksi_p1 * Neta_p1 * dNzeta_p1dzeta;
    _gradShapFunc[ZTA][ 7] = Nksi_m1 * Neta_p1 * dNzeta_p1dzeta;
    _gradShapFunc[ZTA][ 8] = Nksi_0  * Neta_m1 * dNzeta_m1dzeta;
    _gradShapFunc[ZTA][ 9] = Nksi_p1 * Neta_0  * dNzeta_m1dzeta;
    _gradShapFunc[ZTA][10] = Nksi_0  * Neta_p1 * dNzeta_m1dzeta;
    _gradShapFunc[ZTA][11] = Nksi_m1 * Neta_0  * dNzeta_m1dzeta;
    _gradShapFunc[ZTA][12] = Nksi_0  * Neta_0  * dNzeta_m1dzeta;
    _gradShapFunc[ZTA][13] = Nksi_m1 * Neta_m1 * dNzeta_0dzeta ;
    _gradShapFunc[ZTA][14] = Nksi_0  * Neta_m1 * dNzeta_0dzeta ;
    _gradShapFunc[ZTA][15] = Nksi_p1 * Neta_m1 * dNzeta_0dzeta ;
    _gradShapFunc[ZTA][16] = Nksi_p1 * Neta_0  * dNzeta_0dzeta ;
    _gradShapFunc[ZTA][17] = Nksi_p1 * Neta_p1 * dNzeta_0dzeta ;
    _gradShapFunc[ZTA][18] = Nksi_0  * Neta_p1 * dNzeta_0dzeta ;
    _gradShapFunc[ZTA][19] = Nksi_m1 * Neta_p1 * dNzeta_0dzeta ;
    _gradShapFunc[ZTA][20] = Nksi_m1 * Neta_0  * dNzeta_0dzeta ;
    _gradShapFunc[ZTA][21] = Nksi_0  * Neta_0  * dNzeta_0dzeta ;
    _gradShapFunc[ZTA][22] = Nksi_0  * Neta_m1 * dNzeta_p1dzeta;
    _gradShapFunc[ZTA][23] = Nksi_p1 * Neta_0  * dNzeta_p1dzeta;
    _gradShapFunc[ZTA][24] = Nksi_0  * Neta_p1 * dNzeta_p1dzeta;
    _gradShapFunc[ZTA][25] = Nksi_m1 * Neta_0  * dNzeta_p1dzeta;
    _gradShapFunc[ZTA][26] = Nksi_0  * Neta_0  * dNzeta_p1dzeta;

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

    for (CFuint in = 1; in < 27; ++in)
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

    pointJacob *= 0.125;
  }
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace ShapeFunctions



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
