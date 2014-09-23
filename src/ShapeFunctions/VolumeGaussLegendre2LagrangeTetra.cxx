// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "VolumeGaussLegendre2LagrangeTetra.hh"
#include "Framework/IntegratorImplProvider.hh"
#include "ShapeFunctions/LagrangeShapeFunctionTetraP0.hh"
#include "ShapeFunctions/LagrangeShapeFunctionTetraP1.hh"
#include "ShapeFunctions/LagrangeShapeFunctionTetraP2.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



   namespace ShapeFunctions {

//////////////////////////////////////////////////////////////////////////////

Framework::IntegratorImplProvider<Framework::VolumeIntegratorImpl,
                       VolumeGaussLegendre2LagrangeTetra<LagrangeShapeFunctionTetraP0> >
volGaussLegendre2LagrangeTetraP0Provider;

Framework::IntegratorImplProvider<Framework::VolumeIntegratorImpl,
                       VolumeGaussLegendre2LagrangeTetra<LagrangeShapeFunctionTetraP1> >
volGaussLegendre2LagrangeTetraP1Provider;

Framework::IntegratorImplProvider<Framework::VolumeIntegratorImpl,
                       VolumeGaussLegendre2LagrangeTetra<LagrangeShapeFunctionTetraP2> >
volGaussLegendre2LagrangeTetraP2Provider;

//////////////////////////////////////////////////////////////////////////////

template <typename INTERPOLATOR>
void VolumeGaussLegendre2LagrangeTetra<INTERPOLATOR>::setMappedCoordinates()
{
  CFLogDebugMin("VolumeGaussLegendre2LagrangeTetra<INTERPOLATOR>::setMappedCoordinates(): BEGIN\n");

  _mappedCoord.resize(getNbQuadraturePoints());
  _mappedCoord[0].resize(DIM_3D);
  _mappedCoord[1].resize(DIM_3D);
  _mappedCoord[2].resize(DIM_3D);
  _mappedCoord[3].resize(DIM_3D);

  const CFreal a = (5. - sqrt(5.0))/20.;
  const CFreal b = (5. + 3.*sqrt(5.0))/20.;

  _mappedCoord[0][KSI] = a;
  _mappedCoord[0][ETA] = a;
  _mappedCoord[0][ZTA] = a;

  _mappedCoord[1][KSI] = a;
  _mappedCoord[1][ETA] = a;
  _mappedCoord[1][ZTA] = b;

  _mappedCoord[2][KSI] = a;
  _mappedCoord[2][ETA] = b;
  _mappedCoord[2][ZTA] = a;

  _mappedCoord[3][KSI] = b;
  _mappedCoord[3][ETA] = a;
  _mappedCoord[3][ZTA] = a;

  CFLogDebugMin("VolumeGaussLegendre2LagrangeTetra<INTERPOLATOR>::setMappedCoordinates(): END\n");
}

//////////////////////////////////////////////////////////////////////////////

template <typename INTERPOLATOR>
void VolumeGaussLegendre2LagrangeTetra<INTERPOLATOR>::setWeights()
{
  _coeff.resize(getNbQuadraturePoints());

  _coeff[0] = 1.0/24.0;
  _coeff[1] = 1.0/24.0;
  _coeff[2] = 1.0/24.0;
  _coeff[3] = 1.0/24.0;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace ShapeFunctions



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
