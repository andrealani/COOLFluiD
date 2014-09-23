// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "VolumeGaussLegendre1LagrangeLine.hh"
#include "Framework/IntegratorImplProvider.hh"
#include "ShapeFunctions/LagrangeShapeFunctionLineP0.hh"
#include "ShapeFunctions/LagrangeShapeFunctionLineP1.hh"
#include "ShapeFunctions/LagrangeShapeFunctionLineP2.hh"
#include "ShapeFunctions/LagrangeShapeFunctionLineP3.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



   namespace ShapeFunctions {

//////////////////////////////////////////////////////////////////////////////

Framework::IntegratorImplProvider<Framework::VolumeIntegratorImpl,
                       VolumeGaussLegendre1LagrangeLine<LagrangeShapeFunctionLineP0> >
volGaussLegendre1LagrangeLineP0Provider;

Framework::IntegratorImplProvider<Framework::VolumeIntegratorImpl,
                       VolumeGaussLegendre1LagrangeLine<LagrangeShapeFunctionLineP1> >
volGaussLegendre1LagrangeLineP1Provider;

// Framework::IntegratorImplProvider<Framework::VolumeIntegratorImpl,
//                        VolumeGaussLegendre1LagrangeLine<LagrangeShapeFunctionLineP2> >
// volGaussLegendre1LagrangeLineP2Provider;

//////////////////////////////////////////////////////////////////////////////

template <typename INTERPOLATOR>
void VolumeGaussLegendre1LagrangeLine<INTERPOLATOR>::setMappedCoordinates()
{
  _mappedCoord.resize(1);
  _mappedCoord[0].resize(1);

  _mappedCoord[0][KSI] = 0.;
}

//////////////////////////////////////////////////////////////////////////////

template <typename INTERPOLATOR>
void VolumeGaussLegendre1LagrangeLine<INTERPOLATOR>::setWeights()
{
  _coeff.resize(getNbQuadraturePoints());

  _coeff[0] = 2.;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace ShapeFunctions



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
