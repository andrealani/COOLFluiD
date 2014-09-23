// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "VolumeGaussLegendre4LagrangeLine.hh"
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
                       VolumeGaussLegendre4LagrangeLine<LagrangeShapeFunctionLineP0> >
volGaussLegendre4LagrangeLineP0Provider;

Framework::IntegratorImplProvider<Framework::VolumeIntegratorImpl,
                       VolumeGaussLegendre4LagrangeLine<LagrangeShapeFunctionLineP1> >
volGaussLegendre4LagrangeLineP1Provider;

Framework::IntegratorImplProvider<Framework::VolumeIntegratorImpl,
                       VolumeGaussLegendre4LagrangeLine<LagrangeShapeFunctionLineP2> >
volGaussLegendre4LagrangeLineP2Provider;

//////////////////////////////////////////////////////////////////////////////

template <typename INTERPOLATOR>
void VolumeGaussLegendre4LagrangeLine<INTERPOLATOR>::setMappedCoordinates()
{
  _mappedCoord.resize(2);

  _mappedCoord[0].resize(1);
  _mappedCoord[1].resize(1);

  _mappedCoord[0][0] = -1./sqrt(3.);
  _mappedCoord[1][0] =  1./sqrt(3.);
}

//////////////////////////////////////////////////////////////////////////////

template <typename INTERPOLATOR>
void VolumeGaussLegendre4LagrangeLine<INTERPOLATOR>::setWeights()
{
  _coeff.resize(getNbQuadraturePoints());

  _coeff[0] = 1.;
  _coeff[1] = 1.;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace ShapeFunctions



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
