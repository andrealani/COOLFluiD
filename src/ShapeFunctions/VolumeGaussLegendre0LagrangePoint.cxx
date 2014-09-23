// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "VolumeGaussLegendre0LagrangePoint.hh"
#include "Framework/IntegratorImplProvider.hh"
#include "ShapeFunctions/LagrangeShapeFunctionPointP0.hh"
#include "ShapeFunctions/LagrangeShapeFunctionPointP1.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

   namespace ShapeFunctions {

//////////////////////////////////////////////////////////////////////////////

Framework::IntegratorImplProvider<Framework::VolumeIntegratorImpl,
                       VolumeGaussLegendre0LagrangePoint<LagrangeShapeFunctionPointP0> >
volGaussLegendre0LagrangePointP0Provider;

Framework::IntegratorImplProvider<Framework::VolumeIntegratorImpl,
                       VolumeGaussLegendre0LagrangePoint<LagrangeShapeFunctionPointP1> >
volGaussLegendre0LagrangePointP1Provider;

//////////////////////////////////////////////////////////////////////////////

template <typename INTERPOLATOR>
void VolumeGaussLegendre0LagrangePoint<INTERPOLATOR>::setMappedCoordinates()
{
  _mappedCoord.resize(1);
  _mappedCoord[0].resize(1);
  _mappedCoord[0][KSI] = 0.;
}

//////////////////////////////////////////////////////////////////////////////

template <typename INTERPOLATOR>
void VolumeGaussLegendre0LagrangePoint<INTERPOLATOR>::setWeights()
{
  _coeff.resize(getNbQuadraturePoints());
  _coeff[0] = 1.;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace ShapeFunctions



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
