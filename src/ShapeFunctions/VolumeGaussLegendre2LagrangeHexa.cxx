// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "VolumeGaussLegendre2LagrangeHexa.hh"
#include "Framework/IntegratorImplProvider.hh"
#include "ShapeFunctions/LagrangeShapeFunctionHexaP0.hh"
#include "ShapeFunctions/LagrangeShapeFunctionHexaP1.hh"
#include "ShapeFunctions/LagrangeShapeFunctionHexaP2.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



   namespace ShapeFunctions {

//////////////////////////////////////////////////////////////////////////////

Framework::IntegratorImplProvider<Framework::VolumeIntegratorImpl,
                       VolumeGaussLegendre2LagrangeHexa<LagrangeShapeFunctionHexaP0> >
volGaussLegendre2LagrangeHexaP0Provider;

Framework::IntegratorImplProvider<Framework::VolumeIntegratorImpl,
                       VolumeGaussLegendre2LagrangeHexa<LagrangeShapeFunctionHexaP1> >
volGaussLegendre2LagrangeHexaP1Provider;

Framework::IntegratorImplProvider<Framework::VolumeIntegratorImpl,
                       VolumeGaussLegendre2LagrangeHexa<LagrangeShapeFunctionHexaP2> >
volGaussLegendre2LagrangeHexaP2Provider;

//////////////////////////////////////////////////////////////////////////////

template <typename INTERPOLATOR>
void VolumeGaussLegendre2LagrangeHexa<INTERPOLATOR>::setMappedCoordinates()
{
  _mappedCoord.resize(getNbQuadraturePoints());
  for (CFuint i = 0; i <  _mappedCoord.size(); ++i) {
    _mappedCoord[i].resize(DIM_3D);
  }

  _mappedCoord[0][KSI] =   0.0000000000;
  _mappedCoord[0][ETA] =   0.8164965809;
  _mappedCoord[0][ZTA] = - 0.5773502692;

  _mappedCoord[1][KSI] =   0.0000000000;
  _mappedCoord[1][ETA] = - 0.8164965809;
  _mappedCoord[1][ZTA] = - 0.5773502692;

  _mappedCoord[2][KSI] =   0.8164965809;
  _mappedCoord[2][ETA] =   0.0000000000;
  _mappedCoord[2][ZTA] = - 0.5773502692;

  _mappedCoord[3][KSI] = - 0.8164965809;
  _mappedCoord[3][ETA] =   0.0000000000;
  _mappedCoord[3][ZTA] = - 0.5773502692;
}

//////////////////////////////////////////////////////////////////////////////

template <typename INTERPOLATOR>
void VolumeGaussLegendre2LagrangeHexa<INTERPOLATOR>::setWeights()
{
  _coeff.resize(getNbQuadraturePoints());

  _coeff[0] = 2.0;
  _coeff[1] = 2.0;
  _coeff[2] = 2.0;
  _coeff[3] = 2.0;

}

//////////////////////////////////////////////////////////////////////////////

    } // namespace ShapeFunctions



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
