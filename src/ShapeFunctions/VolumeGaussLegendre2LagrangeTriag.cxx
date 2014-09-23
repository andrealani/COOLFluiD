// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "VolumeGaussLegendre2LagrangeTriag.hh"
#include "Framework/IntegratorImplProvider.hh"
#include "ShapeFunctions/LagrangeShapeFunctionTriagP0.hh"
#include "ShapeFunctions/LagrangeShapeFunctionTriagP1.hh"
#include "ShapeFunctions/LagrangeShapeFunctionTriagP2.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



   namespace ShapeFunctions {

//////////////////////////////////////////////////////////////////////////////

Framework::IntegratorImplProvider<Framework::VolumeIntegratorImpl,
                       VolumeGaussLegendre2LagrangeTriag<LagrangeShapeFunctionTriagP0> >
volGaussLegendre2LagrangeTriagP0Provider;

Framework::IntegratorImplProvider<Framework::VolumeIntegratorImpl,
                       VolumeGaussLegendre2LagrangeTriag<LagrangeShapeFunctionTriagP1> >
volGaussLegendre2LagrangeTriagP1Provider;

Framework::IntegratorImplProvider<Framework::VolumeIntegratorImpl,
                       VolumeGaussLegendre2LagrangeTriag<LagrangeShapeFunctionTriagP2> >
volGaussLegendre2LagrangeTriagP2Provider;

//////////////////////////////////////////////////////////////////////////////

template <typename INTERPOLATOR>
void VolumeGaussLegendre2LagrangeTriag<INTERPOLATOR>::setMappedCoordinates()
{
  _mappedCoord.resize(3);
  for (CFuint i = 0; i < 3; ++i) {
    _mappedCoord[i].resize(2);
  }

  _mappedCoord[0][KSI] = 0.5;
  _mappedCoord[0][ETA] = 0.5;

  _mappedCoord[1][KSI] = 0.;
  _mappedCoord[1][ETA] = 0.5;

  _mappedCoord[2][KSI] = 0.5;
  _mappedCoord[2][ETA] = 0.;
}

//////////////////////////////////////////////////////////////////////////////

template <typename INTERPOLATOR>
void VolumeGaussLegendre2LagrangeTriag<INTERPOLATOR>::setWeights()
{
  _coeff.resize(getNbQuadraturePoints());

  _coeff[0] = 1./6.;
  _coeff[1] = 1./6.;
  _coeff[2] = 1./6.;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace ShapeFunctions



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
