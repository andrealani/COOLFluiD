// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "VolumeGaussLegendre1LagrangePrism.hh"
#include "Framework/IntegratorImplProvider.hh"
#include "ShapeFunctions/LagrangeShapeFunctionPrismP0.hh"
#include "ShapeFunctions/LagrangeShapeFunctionPrismP1.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



   namespace ShapeFunctions {

//////////////////////////////////////////////////////////////////////////////

Framework::IntegratorImplProvider<Framework::VolumeIntegratorImpl,
                       VolumeGaussLegendre1LagrangePrism<LagrangeShapeFunctionPrismP0> >
volGaussLegendre1LagrangePrismP0Provider;

Framework::IntegratorImplProvider<Framework::VolumeIntegratorImpl,
                       VolumeGaussLegendre1LagrangePrism<LagrangeShapeFunctionPrismP1> >
volGaussLegendre1LagrangePrismP1Provider;

//////////////////////////////////////////////////////////////////////////////

template <typename INTERPOLATOR>
void VolumeGaussLegendre1LagrangePrism<INTERPOLATOR>::setMappedCoordinates()
{
  _mappedCoord.resize(getNbQuadraturePoints());
  _mappedCoord[0].resize(DIM_3D);

  _mappedCoord[0][KSI] = 1./3.;
  _mappedCoord[0][ETA] = 1./3.;
  _mappedCoord[0][ZTA] = 0.;
}

//////////////////////////////////////////////////////////////////////////////

template <typename INTERPOLATOR>
void VolumeGaussLegendre1LagrangePrism<INTERPOLATOR>::setWeights()
{
  _coeff.resize(getNbQuadraturePoints());

  _coeff[0] = 1.0;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace ShapeFunctions



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
