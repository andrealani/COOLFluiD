// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "VolumeGaussLegendre1LagrangePyram.hh"
#include "Framework/IntegratorImplProvider.hh"
#include "ShapeFunctions/LagrangeShapeFunctionPyramP0.hh"
#include "ShapeFunctions/LagrangeShapeFunctionPyramP1.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



   namespace ShapeFunctions {

//////////////////////////////////////////////////////////////////////////////

Framework::IntegratorImplProvider<Framework::VolumeIntegratorImpl,
                       VolumeGaussLegendre1LagrangePyram<LagrangeShapeFunctionPyramP0> >
volGaussLegendre1LagrangePyramP0Provider;

Framework::IntegratorImplProvider<Framework::VolumeIntegratorImpl,
                       VolumeGaussLegendre1LagrangePyram<LagrangeShapeFunctionPyramP1> >
volGaussLegendre1LagrangePyramP1Provider;

//////////////////////////////////////////////////////////////////////////////

template <typename INTERPOLATOR>
void VolumeGaussLegendre1LagrangePyram<INTERPOLATOR>::setMappedCoordinates()
{
  _mappedCoord.resize(getNbQuadraturePoints());
  _mappedCoord[0].resize(DIM_3D);

  _mappedCoord[0][KSI] = 0.;
  _mappedCoord[0][ETA] = 0.;
  _mappedCoord[0][ZTA] = 0.;
}

//////////////////////////////////////////////////////////////////////////////

template <typename INTERPOLATOR>
void VolumeGaussLegendre1LagrangePyram<INTERPOLATOR>::setWeights()
{
  _coeff.resize(getNbQuadraturePoints());

  _coeff[0] = 8.0;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace ShapeFunctions



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
