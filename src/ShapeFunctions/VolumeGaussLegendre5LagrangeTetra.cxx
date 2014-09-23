// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "VolumeGaussLegendre5LagrangeTetra.hh"
#include "Framework/IntegratorImplProvider.hh"
#include "ShapeFunctions/LagrangeShapeFunctionTetraP0.hh"
#include "ShapeFunctions/LagrangeShapeFunctionTetraP1.hh"
#include "ShapeFunctions/LagrangeShapeFunctionTetraP2.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

   namespace ShapeFunctions {

//////////////////////////////////////////////////////////////////////////////

Framework::IntegratorImplProvider<Framework::VolumeIntegratorImpl,
                       VolumeGaussLegendre5LagrangeTetra<LagrangeShapeFunctionTetraP0> >
volGaussLegendre5LagrangeTetraP0Provider;

Framework::IntegratorImplProvider<Framework::VolumeIntegratorImpl,
                       VolumeGaussLegendre5LagrangeTetra<LagrangeShapeFunctionTetraP1> >
volGaussLegendre5LagrangeTetraP1Provider;

Framework::IntegratorImplProvider<Framework::VolumeIntegratorImpl,
                       VolumeGaussLegendre5LagrangeTetra<LagrangeShapeFunctionTetraP2> >
volGaussLegendre5LagrangeTetraP2Provider;

//////////////////////////////////////////////////////////////////////////////

template <typename INTERPOLATOR>
void VolumeGaussLegendre5LagrangeTetra<INTERPOLATOR>::setMappedCoordinates()
{
  CFLogDebugMin("VolumeGaussLegendre5LagrangeTetra<INTERPOLATOR>::setMappedCoordinates(): BEGIN\n");

  _mappedCoord.resize(getNbQuadraturePoints());
  for(CFuint i=0;i<14;i++)
  {
    _mappedCoord[i].resize(DIM_3D);
  }

  CFreal g1 = 0.092735250310891;
  CFreal g2 = 0.721794249067326;

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

  g1 = 0.310885919263300;
  g2 = 0.067342242210098;

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

  g1 = 0.454496295874350;
  g2 = 0.045503704125650;

  _mappedCoord[8][KSI] = g2;
  _mappedCoord[8][ETA] = g1;
  _mappedCoord[8][ZTA] = g1;


  _mappedCoord[9][KSI] = g1;
  _mappedCoord[9][ETA] = g2;
  _mappedCoord[9][ZTA] = g1;

  _mappedCoord[10][KSI] = g1;
  _mappedCoord[10][ETA] = g1;
  _mappedCoord[10][ZTA] = g2;

  _mappedCoord[11][KSI] = g1;
  _mappedCoord[11][ETA] = g2;
  _mappedCoord[11][ZTA] = g2;

  _mappedCoord[12][KSI] = g2;
  _mappedCoord[12][ETA] = g1;
  _mappedCoord[12][ZTA] = g2;

  _mappedCoord[13][KSI] = g2;
  _mappedCoord[13][ETA] = g2;
  _mappedCoord[13][ZTA] = g1;


  CFLogDebugMin("VolumeGaussLegendre5LagrangeTetra<INTERPOLATOR>::setMappedCoordinates(): END\n");
}

//////////////////////////////////////////////////////////////////////////////

template <typename INTERPOLATOR>
void VolumeGaussLegendre5LagrangeTetra<INTERPOLATOR>::setWeights()
{
  _coeff.resize(getNbQuadraturePoints());

  _coeff[0] = 0.012248840519394;
  _coeff[1] = 0.012248840519394;
  _coeff[2] = 0.012248840519394;
  _coeff[3] = 0.012248840519394;
  _coeff[4] = 0.018781320953003;
  _coeff[5] = 0.018781320953003;
  _coeff[6] = 0.018781320953003;
  _coeff[7] = 0.018781320953003;
  _coeff[8] = 0.007091003462847;
  _coeff[9] = 0.007091003462847;
  _coeff[10] = 0.007091003462847;
  _coeff[11] = 0.007091003462847;
  _coeff[12] = 0.007091003462847;
  _coeff[13] = 0.007091003462847;

}

//////////////////////////////////////////////////////////////////////////////

    } // namespace ShapeFunctions



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
