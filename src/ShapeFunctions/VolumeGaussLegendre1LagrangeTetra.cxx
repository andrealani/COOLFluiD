// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "VolumeGaussLegendre1LagrangeTetra.hh"
#include "Framework/IntegratorImplProvider.hh"
#include "ShapeFunctions/LagrangeShapeFunctionTetraP0.hh"
#include "ShapeFunctions/LagrangeShapeFunctionTetraP1.hh"
#include "ShapeFunctions/LagrangeShapeFunctionTetraP2.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



   namespace ShapeFunctions {

//////////////////////////////////////////////////////////////////////////////

Framework::IntegratorImplProvider<Framework::VolumeIntegratorImpl,
                       VolumeGaussLegendre1LagrangeTetra<LagrangeShapeFunctionTetraP0> >
volGaussLegendre1LagrangeTetraP0Provider;

Framework::IntegratorImplProvider<Framework::VolumeIntegratorImpl,
                       VolumeGaussLegendre1LagrangeTetra<LagrangeShapeFunctionTetraP1> >
volGaussLegendre1LagrangeTetraP1Provider;

Framework::IntegratorImplProvider<Framework::VolumeIntegratorImpl,
                       VolumeGaussLegendre1LagrangeTetra<LagrangeShapeFunctionTetraP2> >
volGaussLegendre1LagrangeTetraP2Provider;

//////////////////////////////////////////////////////////////////////////////

template <typename INTERPOLATOR>
void VolumeGaussLegendre1LagrangeTetra<INTERPOLATOR>::setMappedCoordinates()
{
  CFLogDebugMin("VolumeGaussLegendre1LagrangeTetra<INTERPOLATOR>::setMappedCoordinates(): BEGIN\n");

  _mappedCoord.resize(getNbQuadraturePoints());
  _mappedCoord[0].resize(DIM_3D);

  _mappedCoord[0][KSI] = 0.25;
  _mappedCoord[0][ETA] = 0.25;
  _mappedCoord[0][ZTA] = 0.25;

  CFLogDebugMin("VolumeGaussLegendre1LagrangeTetra<INTERPOLATOR>::setMappedCoordinates(): END\n");
}

//////////////////////////////////////////////////////////////////////////////

template <typename INTERPOLATOR>
void VolumeGaussLegendre1LagrangeTetra<INTERPOLATOR>::setWeights()
{
  _coeff.resize(getNbQuadraturePoints());

  _coeff[0] = 1.0/6.0;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace ShapeFunctions



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
