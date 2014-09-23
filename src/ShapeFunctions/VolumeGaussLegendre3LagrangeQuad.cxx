// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "VolumeGaussLegendre3LagrangeQuad.hh"
#include "Framework/IntegratorImplProvider.hh"
#include "ShapeFunctions/LagrangeShapeFunctionQuadP0.hh"
#include "ShapeFunctions/LagrangeShapeFunctionQuadP1.hh"
#include "ShapeFunctions/LagrangeShapeFunctionQuadP2.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



   namespace ShapeFunctions {

//////////////////////////////////////////////////////////////////////////////

Framework::IntegratorImplProvider<Framework::VolumeIntegratorImpl,
                       VolumeGaussLegendre3LagrangeQuad<LagrangeShapeFunctionQuadP0> >
volGaussLegendre3LagrangeQuadP0Provider;

Framework::IntegratorImplProvider<Framework::VolumeIntegratorImpl,
                       VolumeGaussLegendre3LagrangeQuad<LagrangeShapeFunctionQuadP1> >
volGaussLegendre3LagrangeQuadP1Provider;

Framework::IntegratorImplProvider<Framework::VolumeIntegratorImpl,
                       VolumeGaussLegendre3LagrangeQuad<LagrangeShapeFunctionQuadP2> >
volGaussLegendre3LagrangeQuadP2Provider;

//////////////////////////////////////////////////////////////////////////////

template <typename INTERPOLATOR>
void VolumeGaussLegendre3LagrangeQuad<INTERPOLATOR>::setMappedCoordinates()
{
  CFLogDebugMin("VolumeGaussLegendre3LagrangeQuad<INTERPOLATOR>::setMappedCoordinates(): BEGIN\n");

  _mappedCoord.resize(getNbQuadraturePoints());
  for(CFuint i = 0; i < getNbQuadraturePoints(); ++i){
    _mappedCoord[i].resize(DIM_2D);
  }

  const CFreal ref = 0.7745966692;

  _mappedCoord[0][KSI] = ref;
  _mappedCoord[0][ETA] =-ref;

  _mappedCoord[1][KSI] = ref;
  _mappedCoord[1][ETA] = 0.;

  _mappedCoord[2][KSI] = ref;
  _mappedCoord[2][ETA] = ref;

  _mappedCoord[3][KSI] = 0.;
  _mappedCoord[3][ETA] =-ref;

  _mappedCoord[4][KSI] = 0.;
  _mappedCoord[4][ETA] = 0.;

  _mappedCoord[5][KSI] = 0.;
  _mappedCoord[5][ETA] = ref;

  _mappedCoord[6][KSI] =-ref;
  _mappedCoord[6][ETA] =-ref;

  _mappedCoord[7][KSI] =-ref;
  _mappedCoord[7][ETA] = 0.;

  _mappedCoord[8][KSI] =-ref;
  _mappedCoord[8][ETA] = ref;

  CFLogDebugMin("VolumeGaussLegendre3LagrangeQuad<INTERPOLATOR>::setMappedCoordinates(): END\n");
}

//////////////////////////////////////////////////////////////////////////////

template <typename INTERPOLATOR>
void VolumeGaussLegendre3LagrangeQuad<INTERPOLATOR>::setWeights()
{
  _coeff.resize(getNbQuadraturePoints());

  const CFreal a = 0.555555555555555;
  const CFreal b = 0.888888888888888;

  _coeff[0] = a*a;
  _coeff[1] = a*b;
  _coeff[2] = a*a;
  _coeff[3] = a*b;
  _coeff[4] = b*b;
  _coeff[5] = a*b;
  _coeff[6] = a*a;
  _coeff[7] = a*b;
  _coeff[8] = a*a;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace ShapeFunctions



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
