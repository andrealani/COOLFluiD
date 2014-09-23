// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "VolumeGaussLegendre5LagrangeTriag.hh"
#include "Framework/IntegratorImplProvider.hh"
#include "ShapeFunctions/LagrangeShapeFunctionTriagP0.hh"
#include "ShapeFunctions/LagrangeShapeFunctionTriagP1.hh"
#include "ShapeFunctions/LagrangeShapeFunctionTriagP2.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

   namespace ShapeFunctions {

//////////////////////////////////////////////////////////////////////////////

Framework::IntegratorImplProvider<Framework::VolumeIntegratorImpl,
                       VolumeGaussLegendre5LagrangeTriag<LagrangeShapeFunctionTriagP0> >
volGaussLegendre5LagrangeTriagP0Provider;

Framework::IntegratorImplProvider<Framework::VolumeIntegratorImpl,
                       VolumeGaussLegendre5LagrangeTriag<LagrangeShapeFunctionTriagP1> >
volGaussLegendre5LagrangeTriagP1Provider;

Framework::IntegratorImplProvider<Framework::VolumeIntegratorImpl,
                       VolumeGaussLegendre5LagrangeTriag<LagrangeShapeFunctionTriagP2> >
volGaussLegendre5LagrangeTriagP2Provider;

//////////////////////////////////////////////////////////////////////////////

template <typename INTERPOLATOR>
void VolumeGaussLegendre5LagrangeTriag<INTERPOLATOR>::setMappedCoordinates()
{
  const CFuint nbPoints = getNbQuadraturePoints();
  _mappedCoord.resize(nbPoints);
  for (CFuint i = 0; i < nbPoints; ++i) {
    _mappedCoord[i].resize(DIM_2D);
  }

  const CFreal alpha  = 1./3.;
  const CFreal alpha1 = 0.0597158717;
  const CFreal alpha2 = 0.7974269853;
  const CFreal beta1  = 0.4701420641;
  const CFreal beta2  = 0.1012865073;

  _mappedCoord[0][KSI] = alpha;
  _mappedCoord[0][ETA] = alpha;

  _mappedCoord[1][KSI] = alpha1;
  _mappedCoord[1][ETA] = beta1;

  _mappedCoord[2][KSI] = beta1;
  _mappedCoord[2][ETA] = alpha1;

  _mappedCoord[3][KSI] = beta1;
  _mappedCoord[3][ETA] = beta1;

  _mappedCoord[4][KSI] = alpha2;
  _mappedCoord[4][ETA] = beta2;

  _mappedCoord[5][KSI] = beta2;
  _mappedCoord[5][ETA] = alpha2;

  _mappedCoord[6][KSI] = beta2;
  _mappedCoord[6][ETA] = beta2;
}

//////////////////////////////////////////////////////////////////////////////

template <typename INTERPOLATOR>
void VolumeGaussLegendre5LagrangeTriag<INTERPOLATOR>::setWeights()
{
  _coeff.resize(getNbQuadraturePoints());

  _coeff[0] = 0.5 * 0.225;
  _coeff[1] = 0.5 * 0.1323941527;
  _coeff[2] = 0.5 * 0.1323941527;
  _coeff[3] = 0.5 * 0.1323941527;
  _coeff[4] = 0.5 * 0.1259391805;
  _coeff[5] = 0.5 * 0.1259391805;
  _coeff[6] = 0.5 * 0.1259391805;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace ShapeFunctions



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
