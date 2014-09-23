// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "VolumeGaussLegendre1LagrangeTriag.hh"
#include "Framework/IntegratorImplProvider.hh"
#include "ShapeFunctions/LagrangeShapeFunctionTriagP0.hh"
#include "ShapeFunctions/LagrangeShapeFunctionTriagP1.hh"
#include "ShapeFunctions/LagrangeShapeFunctionTriagP2.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



   namespace ShapeFunctions {

//////////////////////////////////////////////////////////////////////////////

Framework::IntegratorImplProvider<Framework::VolumeIntegratorImpl,
                       VolumeGaussLegendre1LagrangeTriag<LagrangeShapeFunctionTriagP0> >
volGaussLegendre1LagrangeTriagP0Provider;

Framework::IntegratorImplProvider<Framework::VolumeIntegratorImpl,
                       VolumeGaussLegendre1LagrangeTriag<LagrangeShapeFunctionTriagP1> >
volGaussLegendre1LagrangeTriagP1Provider;

Framework::IntegratorImplProvider<Framework::VolumeIntegratorImpl,
                       VolumeGaussLegendre1LagrangeTriag<LagrangeShapeFunctionTriagP2> >
volGaussLegendre1LagrangeTriagP2Provider;

//////////////////////////////////////////////////////////////////////////////

template <typename INTERPOLATOR>
void VolumeGaussLegendre1LagrangeTriag<INTERPOLATOR>::setMappedCoordinates()
{
  CFLogDebugMin("VolumeGaussLegendre1LagrangeTriag<INTERPOLATOR>::setMappedCoordinates(): BEGIN\n");

  _mappedCoord.resize(1);

  _mappedCoord[0].resize(2);
  _mappedCoord[0][KSI] = 1./3.;
  _mappedCoord[0][ETA] = 1./3.;

  CFLogDebugMin("VolumeGaussLegendre1LagrangeTriag<INTERPOLATOR>::setMappedCoordinates(): BEGIN\n");
}

//////////////////////////////////////////////////////////////////////////////

template <typename INTERPOLATOR>
void VolumeGaussLegendre1LagrangeTriag<INTERPOLATOR>::setWeights()
{
  _coeff.resize(getNbQuadraturePoints());

  _coeff[0] = 0.5;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace ShapeFunctions



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
