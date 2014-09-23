// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "VolumeGaussLegendre2LagrangeQuad.hh"
#include "Framework/IntegratorImplProvider.hh"
#include "ShapeFunctions/LagrangeShapeFunctionQuadP0.hh"
#include "ShapeFunctions/LagrangeShapeFunctionQuadP1.hh"
#include "ShapeFunctions/LagrangeShapeFunctionQuadP2.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



   namespace ShapeFunctions {

//////////////////////////////////////////////////////////////////////////////

Framework::IntegratorImplProvider<Framework::VolumeIntegratorImpl,
                       VolumeGaussLegendre2LagrangeQuad<LagrangeShapeFunctionQuadP0> >
volGaussLegendre2LagrangeQuadP0Provider;

Framework::IntegratorImplProvider<Framework::VolumeIntegratorImpl,
                       VolumeGaussLegendre2LagrangeQuad<LagrangeShapeFunctionQuadP1> >
volGaussLegendre2LagrangeQuadP1Provider;

Framework::IntegratorImplProvider<Framework::VolumeIntegratorImpl,
                       VolumeGaussLegendre2LagrangeQuad<LagrangeShapeFunctionQuadP2> >
volGaussLegendre2LagrangeQuadP2Provider;

//////////////////////////////////////////////////////////////////////////////

template <typename INTERPOLATOR>
void VolumeGaussLegendre2LagrangeQuad<INTERPOLATOR>::setMappedCoordinates()
{
  CFLogDebugMin("VolumeGaussLegendre2LagrangeQuad<INTERPOLATOR>::setMappedCoordinates(): BEGIN\n");

  _mappedCoord.resize(getNbQuadraturePoints());
  _mappedCoord[0].resize(DIM_2D);
  _mappedCoord[1].resize(DIM_2D);
  _mappedCoord[2].resize(DIM_2D);
  _mappedCoord[3].resize(DIM_2D);

  const CFreal ref = 1.0 / sqrt(3.0);

  _mappedCoord[0][KSI] = - ref;
  _mappedCoord[0][ETA] = - ref;

  _mappedCoord[1][KSI] = - ref;
  _mappedCoord[1][ETA] =   ref;

  _mappedCoord[2][KSI] =   ref;
  _mappedCoord[2][ETA] = - ref;

  _mappedCoord[3][KSI] =   ref;
  _mappedCoord[3][ETA] =   ref;

  CFLogDebugMin("VolumeGaussLegendre2LagrangeQuad<INTERPOLATOR>::setMappedCoordinates(): END\n");
}

//////////////////////////////////////////////////////////////////////////////

template <typename INTERPOLATOR>
void VolumeGaussLegendre2LagrangeQuad<INTERPOLATOR>::setWeights()
{
  _coeff.resize(getNbQuadraturePoints());

  _coeff[0] = 1.;
  _coeff[1] = 1.;
  _coeff[2] = 1.;
  _coeff[3] = 1.;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace ShapeFunctions



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
