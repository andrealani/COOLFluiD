// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "VolumeGaussLegendre8LagrangeTetra.hh"
#include "Framework/IntegratorImplProvider.hh"
#include "ShapeFunctions/LagrangeShapeFunctionTetraP0.hh"
#include "ShapeFunctions/LagrangeShapeFunctionTetraP1.hh"
#include "ShapeFunctions/LagrangeShapeFunctionTetraP2.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

   namespace ShapeFunctions {

//////////////////////////////////////////////////////////////////////////////

Framework::IntegratorImplProvider<Framework::VolumeIntegratorImpl,
                       VolumeGaussLegendre8LagrangeTetra<LagrangeShapeFunctionTetraP0> >
volGaussLegendre8LagrangeTetraP0Provider;

Framework::IntegratorImplProvider<Framework::VolumeIntegratorImpl,
                       VolumeGaussLegendre8LagrangeTetra<LagrangeShapeFunctionTetraP1> >
volGaussLegendre8LagrangeTetraP1Provider;

Framework::IntegratorImplProvider<Framework::VolumeIntegratorImpl,
                       VolumeGaussLegendre8LagrangeTetra<LagrangeShapeFunctionTetraP2> >
volGaussLegendre8LagrangeTetraP2Provider;

//////////////////////////////////////////////////////////////////////////////

template <typename INTERPOLATOR>
void VolumeGaussLegendre8LagrangeTetra<INTERPOLATOR>::setMappedCoordinates()
{
  CFLogDebugMin("VolumeGaussLegendre8LagrangeTetra<INTERPOLATOR>::setMappedCoordinates(): BEGIN\n");

  _mappedCoord.resize(getNbQuadraturePoints());
  for(CFuint i=0;i<43;i++)
  {
    _mappedCoord[i].resize(DIM_3D);
  }

  CFreal g1;
  CFreal g2;
  CFreal g3;


  _mappedCoord[0][KSI] = 0.25;
  _mappedCoord[0][ETA] = 0.25;
  _mappedCoord[0][ZTA] = 0.25;

  g1 = 0.206829931610673;
  g2 = 0.379510205167980;

  _mappedCoord[1][KSI] = g1;
  _mappedCoord[1][ETA] = g1;
  _mappedCoord[1][ZTA] = g1;

  _mappedCoord[2][KSI] = g2;
  _mappedCoord[2][ETA] = g1;
  _mappedCoord[2][ZTA] = g1;

  _mappedCoord[3][KSI] = g1;
  _mappedCoord[3][ETA] = g2;
  _mappedCoord[3][ZTA] = g1;

  _mappedCoord[4][KSI] = g1;
  _mappedCoord[4][ETA] = g1;
  _mappedCoord[4][ZTA] = g2;

  g1 = 0.082103588310547;
  g2 = 0.753689235068360;

  _mappedCoord[5][KSI] = g1;
  _mappedCoord[5][ETA] = g1;
  _mappedCoord[5][ZTA] = g1;

  _mappedCoord[6][KSI] = g2;
  _mappedCoord[6][ETA] = g1;
  _mappedCoord[6][ZTA] = g1;

  _mappedCoord[7][KSI] = g1;
  _mappedCoord[7][ETA] = g2;
  _mappedCoord[7][ZTA] = g1;

  _mappedCoord[8][KSI] = g1;
  _mappedCoord[8][ETA] = g1;
  _mappedCoord[8][ZTA] = g2;

  g1 = 0.005781950505198;
  g2 = 0.982654148484406;

  _mappedCoord[9][KSI] = g1;
  _mappedCoord[9][ETA] = g1;
  _mappedCoord[9][ZTA] = g1;

  _mappedCoord[10][KSI] = g2;
  _mappedCoord[10][ETA] = g1;
  _mappedCoord[10][ZTA] = g1;

  _mappedCoord[11][KSI] = g1;
  _mappedCoord[11][ETA] = g2;
  _mappedCoord[11][ZTA] = g1;

  _mappedCoord[12][KSI] = g1;
  _mappedCoord[12][ETA] = g1;
  _mappedCoord[12][ZTA] = g2;

  g1 = 0.050532740018894;
  g2 = 0.449467259981106;

  _mappedCoord[13][KSI] = g2;
  _mappedCoord[13][ETA] = g1;
  _mappedCoord[13][ZTA] = g1;

  _mappedCoord[14][KSI] = g1;
  _mappedCoord[14][ETA] = g2;
  _mappedCoord[14][ZTA] = g1;

  _mappedCoord[15][KSI] = g1;
  _mappedCoord[15][ETA] = g1;
  _mappedCoord[15][ZTA] = g2;

  _mappedCoord[16][KSI] = g1;
  _mappedCoord[16][ETA] = g2;
  _mappedCoord[16][ZTA] = g2;

  _mappedCoord[17][KSI] = g2;
  _mappedCoord[17][ETA] = g1;
  _mappedCoord[17][ZTA] = g2;

  _mappedCoord[18][KSI] = g2;
  _mappedCoord[18][ETA] = g2;
  _mappedCoord[18][ZTA] = g1;

  g1 = 0.229066536116811;
  g2 = 0.035639582788534;
  g3 = 0.506227344977843;

  _mappedCoord[19][KSI] = g2;
  _mappedCoord[19][ETA] = g1;
  _mappedCoord[19][ZTA] = g1;

  _mappedCoord[20][KSI] = g1;
  _mappedCoord[20][ETA] = g2;
  _mappedCoord[20][ZTA] = g1;

  _mappedCoord[21][KSI] = g1;
  _mappedCoord[21][ETA] = g1;
  _mappedCoord[21][ZTA] = g2;

  _mappedCoord[22][KSI] = g3;
  _mappedCoord[22][ETA] = g1;
  _mappedCoord[22][ZTA] = g1;

  _mappedCoord[23][KSI] = g1;
  _mappedCoord[23][ETA] = g3;
  _mappedCoord[23][ZTA] = g1;

  _mappedCoord[24][KSI] = g1;
  _mappedCoord[24][ETA] = g1;
  _mappedCoord[24][ZTA] = g3;

  _mappedCoord[25][KSI] = g2;
  _mappedCoord[25][ETA] = g3;
  _mappedCoord[25][ZTA] = g1;

  _mappedCoord[26][KSI] = g2;
  _mappedCoord[26][ETA] = g1;
  _mappedCoord[26][ZTA] = g3;

  _mappedCoord[27][KSI] = g3;
  _mappedCoord[27][ETA] = g2;
  _mappedCoord[27][ZTA] = g1;

  _mappedCoord[28][KSI] = g3;
  _mappedCoord[28][ETA] = g1;
  _mappedCoord[28][ZTA] = g2;

  _mappedCoord[29][KSI] = g1;
  _mappedCoord[29][ETA] = g2;
  _mappedCoord[29][ZTA] = g3;

  _mappedCoord[30][KSI] = g1;
  _mappedCoord[30][ETA] = g3;
  _mappedCoord[30][ZTA] = g2;

  g1 = 0.036607749553198;
  g2 = 0.190486041934633;
  g3 = 0.736298458958971;

  _mappedCoord[31][KSI] = g2;
  _mappedCoord[31][ETA] = g1;
  _mappedCoord[31][ZTA] = g1;

  _mappedCoord[32][KSI] = g1;
  _mappedCoord[32][ETA] = g2;
  _mappedCoord[32][ZTA] = g1;

  _mappedCoord[33][KSI] = g1;
  _mappedCoord[33][ETA] = g1;
  _mappedCoord[33][ZTA] = g2;

  _mappedCoord[34][KSI] = g3;
  _mappedCoord[34][ETA] = g1;
  _mappedCoord[34][ZTA] = g1;

  _mappedCoord[35][KSI] = g1;
  _mappedCoord[35][ETA] = g3;
  _mappedCoord[35][ZTA] = g1;

  _mappedCoord[36][KSI] = g1;
  _mappedCoord[36][ETA] = g1;
  _mappedCoord[36][ZTA] = g3;

  _mappedCoord[37][KSI] = g2;
  _mappedCoord[37][ETA] = g3;
  _mappedCoord[37][ZTA] = g1;

  _mappedCoord[38][KSI] = g2;
  _mappedCoord[38][ETA] = g1;
  _mappedCoord[38][ZTA] = g3;

  _mappedCoord[39][KSI] = g3;
  _mappedCoord[39][ETA] = g2;
  _mappedCoord[39][ZTA] = g1;

  _mappedCoord[40][KSI] = g3;
  _mappedCoord[40][ETA] = g1;
  _mappedCoord[40][ZTA] = g2;

  _mappedCoord[41][KSI] = g1;
  _mappedCoord[41][ETA] = g2;
  _mappedCoord[41][ZTA] = g3;

  _mappedCoord[42][KSI] = g1;
  _mappedCoord[42][ETA] = g3;
  _mappedCoord[42][ZTA] = g2;

  CFLogDebugMin("VolumeGaussLegendre8LagrangeTetra<INTERPOLATOR>::setMappedCoordinates(): END\n");
}

//////////////////////////////////////////////////////////////////////////////

template <typename INTERPOLATOR>
void VolumeGaussLegendre8LagrangeTetra<INTERPOLATOR>::setWeights()
{
  _coeff.resize(getNbQuadraturePoints());
  _coeff[0] = -0.020500188658640;

  for(CFuint i=1;i<5;i++)
    _coeff[i] = 0.014250305822867;
  for(CFuint i=5;i<9;i++)
    _coeff[i] = 0.001967033313134;
  for(CFuint i=9;i<13;i++)
    _coeff[i] = 0.000169834109093;
  for(CFuint i=13;i<19;i++)
    _coeff[i] = 0.004579683824467;
  for(CFuint i=19;i<31;i++)
    _coeff[i] = 0.005704485808682;
  for(CFuint i=31;i<43;i++)
    _coeff[i] = 0.002140519141162;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace ShapeFunctions



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
