// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "VolumeGaussLegendre3LagrangePrism.hh"
#include "Framework/IntegratorImplProvider.hh"
#include "ShapeFunctions/LagrangeShapeFunctionPrismP0.hh"
#include "ShapeFunctions/LagrangeShapeFunctionPrismP1.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



   namespace ShapeFunctions {

//////////////////////////////////////////////////////////////////////////////

Framework::IntegratorImplProvider<Framework::VolumeIntegratorImpl,
                       VolumeGaussLegendre3LagrangePrism<LagrangeShapeFunctionPrismP0> >
volGaussLegendre3LagrangePrismP0Provider;

Framework::IntegratorImplProvider<Framework::VolumeIntegratorImpl,
                       VolumeGaussLegendre3LagrangePrism<LagrangeShapeFunctionPrismP1> >
volGaussLegendre3LagrangePrismP1Provider;

//////////////////////////////////////////////////////////////////////////////

template <typename INTERPOLATOR>
void VolumeGaussLegendre3LagrangePrism<INTERPOLATOR>::setMappedCoordinates()
{
  _mappedCoord.resize(getNbQuadraturePoints());
  for (CFuint i = 0; i < getNbQuadraturePoints(); ++i) {
    _mappedCoord[i].resize(DIM_3D);
  }

  const CFreal zt = sqrt(1.0/3.0);

  _mappedCoord[0][KSI] = 0.5;
  _mappedCoord[0][ETA] = 0.5;
  _mappedCoord[0][ZTA] = zt;

  _mappedCoord[1][KSI] = 0.;
  _mappedCoord[1][ETA] = 0.5;
  _mappedCoord[0][ZTA] = zt;

  _mappedCoord[2][KSI] = 0.5;
  _mappedCoord[2][ETA] = 0.;
  _mappedCoord[0][ZTA] = zt;

  _mappedCoord[3][KSI] = 0.5;
  _mappedCoord[3][ETA] = 0.5;
  _mappedCoord[3][ZTA] = -zt;

  _mappedCoord[4][KSI] = 0.;
  _mappedCoord[4][ETA] = 0.5;
  _mappedCoord[4][ZTA] = -zt;

  _mappedCoord[5][KSI] = 0.5;
  _mappedCoord[5][ETA] = 0.;
  _mappedCoord[5][ZTA] = -zt;

}

//////////////////////////////////////////////////////////////////////////////

template <typename INTERPOLATOR>
void VolumeGaussLegendre3LagrangePrism<INTERPOLATOR>::setWeights()
{
  _coeff.resize(getNbQuadraturePoints());

  const CFreal weight = 1.0/6.0;

  _coeff[0] = weight;
  _coeff[1] = weight;
  _coeff[2] = weight;
  _coeff[3] = weight;
  _coeff[4] = weight;
  _coeff[5] = weight;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace ShapeFunctions



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
