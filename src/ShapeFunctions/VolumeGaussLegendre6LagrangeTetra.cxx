// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "VolumeGaussLegendre6LagrangeTetra.hh"
#include "Framework/IntegratorImplProvider.hh"
#include "ShapeFunctions/LagrangeShapeFunctionTetraP0.hh"
#include "ShapeFunctions/LagrangeShapeFunctionTetraP1.hh"
#include "ShapeFunctions/LagrangeShapeFunctionTetraP2.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

   namespace ShapeFunctions {

//////////////////////////////////////////////////////////////////////////////

Framework::IntegratorImplProvider<Framework::VolumeIntegratorImpl,
                       VolumeGaussLegendre6LagrangeTetra<LagrangeShapeFunctionTetraP0> >
volGaussLegendre6LagrangeTetraP0Provider;

Framework::IntegratorImplProvider<Framework::VolumeIntegratorImpl,
                       VolumeGaussLegendre6LagrangeTetra<LagrangeShapeFunctionTetraP1> >
volGaussLegendre6LagrangeTetraP1Provider;

Framework::IntegratorImplProvider<Framework::VolumeIntegratorImpl,
                       VolumeGaussLegendre6LagrangeTetra<LagrangeShapeFunctionTetraP2> >
volGaussLegendre6LagrangeTetraP2Provider;

//////////////////////////////////////////////////////////////////////////////

template <typename INTERPOLATOR>
void VolumeGaussLegendre6LagrangeTetra<INTERPOLATOR>::setMappedCoordinates()
{
  CFLogDebugMin("VolumeGaussLegendre6LagrangeTetra<INTERPOLATOR>::setMappedCoordinates(): BEGIN\n");

  _mappedCoord.resize(getNbQuadraturePoints());
  for(CFuint i=0;i<24;i++)
  {
    _mappedCoord[i].resize(DIM_3D);
  }

  CFreal g1 = 0.214602871259152;
  CFreal g2 = 0.356191386222544;
  CFreal g3;

  _mappedCoord[0][KSI] = g1;
  _mappedCoord[0][ETA] = g1;
  _mappedCoord[0][ZTA] = g1;

  _mappedCoord[1][KSI] = g2;
  _mappedCoord[1][ETA] = g1;
  _mappedCoord[1][ZTA] = g1;

  _mappedCoord[2][KSI] = g1;
  _mappedCoord[2][ETA] = g2;
  _mappedCoord[2][ZTA] = g1;

  _mappedCoord[3][KSI] = g1;
  _mappedCoord[3][ETA] = g1;
  _mappedCoord[3][ZTA] = g2;

  g1 = 0.040673958534612;
  g2 = 0.877978124396166;

  _mappedCoord[4][KSI] = g1;
  _mappedCoord[4][ETA] = g1;
  _mappedCoord[4][ZTA] = g1;

  _mappedCoord[5][KSI] = g2;
  _mappedCoord[5][ETA] = g1;
  _mappedCoord[5][ZTA] = g1;

  _mappedCoord[6][KSI] = g1;
  _mappedCoord[6][ETA] = g2;
  _mappedCoord[6][ZTA] = g1;

  _mappedCoord[7][KSI] = g1;
  _mappedCoord[7][ETA] = g1;
  _mappedCoord[7][ZTA] = g2;

  g1 = 0.322337890142275;
  g2 = 0.032986329573174;

  _mappedCoord[8][KSI] = g1;
  _mappedCoord[8][ETA] = g1;
  _mappedCoord[8][ZTA] = g1;

  _mappedCoord[9][KSI] = g2;
  _mappedCoord[9][ETA] = g1;
  _mappedCoord[9][ZTA] = g1;

  _mappedCoord[10][KSI] = g1;
  _mappedCoord[10][ETA] = g2;
  _mappedCoord[10][ZTA] = g1;

  _mappedCoord[11][KSI] = g1;
  _mappedCoord[11][ETA] = g1;
  _mappedCoord[11][ZTA] = g2;

  g1 = 0.063661001875018;
  g2 = 0.269672331458316;
  g3 = 0.603005664791649;

  _mappedCoord[12][KSI] = g2;
  _mappedCoord[12][ETA] = g1;
  _mappedCoord[12][ZTA] = g1;

  _mappedCoord[13][KSI] = g1;
  _mappedCoord[13][ETA] = g2;
  _mappedCoord[13][ZTA] = g1;

  _mappedCoord[14][KSI] = g1;
  _mappedCoord[14][ETA] = g1;
  _mappedCoord[14][ZTA] = g2;

  _mappedCoord[15][KSI] = g3;
  _mappedCoord[15][ETA] = g1;
  _mappedCoord[15][ZTA] = g1;

  _mappedCoord[16][KSI] = g1;
  _mappedCoord[16][ETA] = g3;
  _mappedCoord[16][ZTA] = g1;

  _mappedCoord[17][KSI] = g1;
  _mappedCoord[17][ETA] = g1;
  _mappedCoord[17][ZTA] = g3;

  _mappedCoord[18][KSI] = g2;
  _mappedCoord[18][ETA] = g1;
  _mappedCoord[18][ZTA] = g3;

  _mappedCoord[19][KSI] = g2;
  _mappedCoord[19][ETA] = g3;
  _mappedCoord[19][ZTA] = g1;

  _mappedCoord[20][KSI] = g3;
  _mappedCoord[20][ETA] = g1;
  _mappedCoord[20][ZTA] = g2;

  _mappedCoord[21][KSI] = g3;
  _mappedCoord[21][ETA] = g2;
  _mappedCoord[21][ZTA] = g1;

  _mappedCoord[22][KSI] = g1;
  _mappedCoord[22][ETA] = g2;
  _mappedCoord[22][ZTA] = g3;

  _mappedCoord[23][KSI] = g1;
  _mappedCoord[23][ETA] = g3;
  _mappedCoord[23][ZTA] = g2;


  CFLogDebugMin("VolumeGaussLegendre6LagrangeTetra<INTERPOLATOR>::setMappedCoordinates(): END\n");
}

//////////////////////////////////////////////////////////////////////////////

template <typename INTERPOLATOR>
void VolumeGaussLegendre6LagrangeTetra<INTERPOLATOR>::setWeights()
{
  _coeff.resize(getNbQuadraturePoints());

  for(CFuint i=0;i<4;i++)
    _coeff[i] = 0.006653791709695;
  for(CFuint i=4;i<8;i++)
    _coeff[i] = 0.001679535175887;
  for(CFuint i=8;i<12;i++)
    _coeff[i] = 0.009226196923943;
  for(CFuint i=12;i<24;i++)
    _coeff[i] = 0.008035714285714;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace ShapeFunctions



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
