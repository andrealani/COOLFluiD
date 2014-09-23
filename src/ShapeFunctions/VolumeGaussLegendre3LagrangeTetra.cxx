// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "VolumeGaussLegendre3LagrangeTetra.hh"
#include "Framework/IntegratorImplProvider.hh"
#include "ShapeFunctions/LagrangeShapeFunctionTetraP0.hh"
#include "ShapeFunctions/LagrangeShapeFunctionTetraP1.hh"
#include "ShapeFunctions/LagrangeShapeFunctionTetraP2.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

   namespace ShapeFunctions {

//////////////////////////////////////////////////////////////////////////////

Framework::IntegratorImplProvider<Framework::VolumeIntegratorImpl,
                       VolumeGaussLegendre3LagrangeTetra<LagrangeShapeFunctionTetraP0> >
volGaussLegendre3LagrangeTetraP0Provider;

Framework::IntegratorImplProvider<Framework::VolumeIntegratorImpl,
                       VolumeGaussLegendre3LagrangeTetra<LagrangeShapeFunctionTetraP1> >
volGaussLegendre3LagrangeTetraP1Provider;

Framework::IntegratorImplProvider<Framework::VolumeIntegratorImpl,
                       VolumeGaussLegendre3LagrangeTetra<LagrangeShapeFunctionTetraP2> >
volGaussLegendre3LagrangeTetraP2Provider;

//////////////////////////////////////////////////////////////////////////////

template <typename INTERPOLATOR>
void VolumeGaussLegendre3LagrangeTetra<INTERPOLATOR>::setMappedCoordinates()
{
  CFLogDebugMin("VolumeGaussLegendre3LagrangeTetra<INTERPOLATOR>::setMappedCoordinates(): BEGIN\n");

  _mappedCoord.resize(getNbQuadraturePoints());
  _mappedCoord[0].resize(DIM_3D);
  _mappedCoord[1].resize(DIM_3D);
  _mappedCoord[2].resize(DIM_3D);
  _mappedCoord[3].resize(DIM_3D);
  _mappedCoord[4].resize(DIM_3D);

  const CFreal a = 0.25 ;
  const CFreal b = 1./6.;
  const CFreal c = 0.5;

  _mappedCoord[0][KSI] = a;
  _mappedCoord[0][ETA] = a;
  _mappedCoord[0][ZTA] = a;

  _mappedCoord[1][KSI] = b;
  _mappedCoord[1][ETA] = b;
  _mappedCoord[1][ZTA] = b;

  _mappedCoord[2][KSI] = b;
  _mappedCoord[2][ETA] = b;
  _mappedCoord[2][ZTA] = c;

  _mappedCoord[3][KSI] = b;
  _mappedCoord[3][ETA] = c;
  _mappedCoord[3][ZTA] = b;

  _mappedCoord[4][KSI] = c;
  _mappedCoord[4][ETA] = b;
  _mappedCoord[4][ZTA] = b;

  CFLogDebugMin("VolumeGaussLegendre3LagrangeTetra<INTERPOLATOR>::setMappedCoordinates(): END\n");
}

//////////////////////////////////////////////////////////////////////////////

template <typename INTERPOLATOR>
void VolumeGaussLegendre3LagrangeTetra<INTERPOLATOR>::setWeights()
{
  _coeff.resize(getNbQuadraturePoints());

  _coeff[0] = -2.0/15.0;
  _coeff[1] = 3.0/40.0;
  _coeff[2] = 3.0/40.0;
  _coeff[3] = 3.0/40.0;
  _coeff[4] = 3.0/40.0;

}

//////////////////////////////////////////////////////////////////////////////

    } // namespace ShapeFunctions



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
