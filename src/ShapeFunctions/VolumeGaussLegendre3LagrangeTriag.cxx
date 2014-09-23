// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "VolumeGaussLegendre3LagrangeTriag.hh"
#include "Framework/IntegratorImplProvider.hh"
#include "ShapeFunctions/LagrangeShapeFunctionTriagP0.hh"
#include "ShapeFunctions/LagrangeShapeFunctionTriagP1.hh"
#include "ShapeFunctions/LagrangeShapeFunctionTriagP2.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



   namespace ShapeFunctions {

//////////////////////////////////////////////////////////////////////////////

Framework::IntegratorImplProvider<Framework::VolumeIntegratorImpl,
                       VolumeGaussLegendre3LagrangeTriag<LagrangeShapeFunctionTriagP0> >
volGaussLegendre3LagrangeTriagP0Provider;

Framework::IntegratorImplProvider<Framework::VolumeIntegratorImpl,
                       VolumeGaussLegendre3LagrangeTriag<LagrangeShapeFunctionTriagP1> >
volGaussLegendre3LagrangeTriagP1Provider;

Framework::IntegratorImplProvider<Framework::VolumeIntegratorImpl,
                       VolumeGaussLegendre3LagrangeTriag<LagrangeShapeFunctionTriagP2> >
volGaussLegendre3LagrangeTriagP2Provider;

//////////////////////////////////////////////////////////////////////////////

template <typename INTERPOLATOR>
void VolumeGaussLegendre3LagrangeTriag<INTERPOLATOR>::setMappedCoordinates()
{
  const CFuint nbPoints = getNbQuadraturePoints();
  _mappedCoord.resize(nbPoints);
  for (CFuint i = 0; i < nbPoints; ++i) {
    _mappedCoord[i].resize(DIM_2D);
  }

  const CFreal t0 = 1./3.;
  const CFreal t1 = 0.2;
  const CFreal t2 = 0.6;

  _mappedCoord[0][KSI] = t0;
  _mappedCoord[0][ETA] = t0;

  _mappedCoord[1][KSI] = t1;
  _mappedCoord[1][ETA] = t1;

  _mappedCoord[2][KSI] = t1;
  _mappedCoord[2][ETA] = t2;

  _mappedCoord[3][KSI] = t2;
  _mappedCoord[3][ETA] = t1;
}

//////////////////////////////////////////////////////////////////////////////

template <typename INTERPOLATOR>
void VolumeGaussLegendre3LagrangeTriag<INTERPOLATOR>::setWeights()
{
  _coeff.resize(getNbQuadraturePoints());

  _coeff[0] = -27./96.;
  _coeff[1] =  25./96.;
  _coeff[2] =  25./96.;
  _coeff[3] =  25./96.;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace ShapeFunctions



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
