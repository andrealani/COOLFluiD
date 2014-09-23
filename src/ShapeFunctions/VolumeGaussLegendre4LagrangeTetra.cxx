// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "VolumeGaussLegendre4LagrangeTetra.hh"
#include "Framework/IntegratorImplProvider.hh"
#include "ShapeFunctions/LagrangeShapeFunctionTetraP0.hh"
#include "ShapeFunctions/LagrangeShapeFunctionTetraP1.hh"
#include "ShapeFunctions/LagrangeShapeFunctionTetraP2.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

   namespace ShapeFunctions {

//////////////////////////////////////////////////////////////////////////////

Framework::IntegratorImplProvider<Framework::VolumeIntegratorImpl,
                       VolumeGaussLegendre4LagrangeTetra<LagrangeShapeFunctionTetraP0> >
volGaussLegendre4LagrangeTetraP0Provider;

Framework::IntegratorImplProvider<Framework::VolumeIntegratorImpl,
                       VolumeGaussLegendre4LagrangeTetra<LagrangeShapeFunctionTetraP1> >
volGaussLegendre4LagrangeTetraP1Provider;

Framework::IntegratorImplProvider<Framework::VolumeIntegratorImpl,
                       VolumeGaussLegendre4LagrangeTetra<LagrangeShapeFunctionTetraP2> >
volGaussLegendre4LagrangeTetraP2Provider;

//////////////////////////////////////////////////////////////////////////////

template <typename INTERPOLATOR>
void VolumeGaussLegendre4LagrangeTetra<INTERPOLATOR>::setMappedCoordinates()
{
  CFLogDebugMin("VolumeGaussLegendre4LagrangeTetra<INTERPOLATOR>::setMappedCoordinates(): BEGIN\n");

  _mappedCoord.resize(getNbQuadraturePoints());
  _mappedCoord[0].resize(DIM_3D);
  _mappedCoord[1].resize(DIM_3D);
  _mappedCoord[2].resize(DIM_3D);
  _mappedCoord[3].resize(DIM_3D);
  _mappedCoord[4].resize(DIM_3D);
  _mappedCoord[5].resize(DIM_3D);
  _mappedCoord[6].resize(DIM_3D);
  _mappedCoord[7].resize(DIM_3D);
  _mappedCoord[8].resize(DIM_3D);
  _mappedCoord[9].resize(DIM_3D);
  _mappedCoord[10].resize(DIM_3D);

  _mappedCoord[0][KSI] = 0.25;
  _mappedCoord[0][ETA] = 0.25;
  _mappedCoord[0][ZTA] = 0.25;

  CFreal g1 = 0.071428571428572;
  CFreal g2 = 0.785714285714286;

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

  g1 = 0.399403576166799;
  g2 = 0.100596423833201;

  _mappedCoord[5][KSI] = g1;
  _mappedCoord[5][ETA] = g2;
  _mappedCoord[5][ZTA] = g2;

  _mappedCoord[6][KSI] = g2;
  _mappedCoord[6][ETA] = g1;
  _mappedCoord[6][ZTA] = g2;

  _mappedCoord[7][KSI] = g2;
  _mappedCoord[7][ETA] = g2;
  _mappedCoord[7][ZTA] = g1;

  _mappedCoord[8][KSI] = g2;
  _mappedCoord[8][ETA] = g1;
  _mappedCoord[8][ZTA] = g1;


  _mappedCoord[9][KSI] = g1;
  _mappedCoord[9][ETA] = g2;
  _mappedCoord[9][ZTA] = g1;


  _mappedCoord[10][KSI] = g1;
  _mappedCoord[10][ETA] = g1;
  _mappedCoord[10][ZTA] = g2;


  CFLogDebugMin("VolumeGaussLegendre4LagrangeTetra<INTERPOLATOR>::setMappedCoordinates(): END\n");
}

//////////////////////////////////////////////////////////////////////////////

template <typename INTERPOLATOR>
void VolumeGaussLegendre4LagrangeTetra<INTERPOLATOR>::setWeights()
{
  _coeff.resize(getNbQuadraturePoints());

  _coeff[0] = -0.013155555555556;
  _coeff[1] = 0.007622222222222;
  _coeff[2] = 0.007622222222222;
  _coeff[3] = 0.007622222222222;
  _coeff[4] = 0.007622222222222;
  _coeff[5] = 0.024888888888889;
  _coeff[6] = 0.024888888888889;
  _coeff[7] = 0.024888888888889;
  _coeff[8] = 0.024888888888889;
  _coeff[9] = 0.024888888888889;
  _coeff[10] = 0.024888888888889;

}

//////////////////////////////////////////////////////////////////////////////

    } // namespace ShapeFunctions



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
