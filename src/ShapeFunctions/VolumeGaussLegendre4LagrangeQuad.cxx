// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "VolumeGaussLegendre4LagrangeQuad.hh"
#include "Framework/IntegratorImplProvider.hh"
#include "ShapeFunctions/LagrangeShapeFunctionQuadP0.hh"
#include "ShapeFunctions/LagrangeShapeFunctionQuadP1.hh"
#include "ShapeFunctions/LagrangeShapeFunctionQuadP2.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



   namespace ShapeFunctions {

//////////////////////////////////////////////////////////////////////////////

Framework::IntegratorImplProvider<Framework::VolumeIntegratorImpl,
                       VolumeGaussLegendre4LagrangeQuad<LagrangeShapeFunctionQuadP0> >
volGaussLegendre4LagrangeQuadP0Provider;

Framework::IntegratorImplProvider<Framework::VolumeIntegratorImpl,
                       VolumeGaussLegendre4LagrangeQuad<LagrangeShapeFunctionQuadP1> >
volGaussLegendre4LagrangeQuadP1Provider;

Framework::IntegratorImplProvider<Framework::VolumeIntegratorImpl,
                       VolumeGaussLegendre4LagrangeQuad<LagrangeShapeFunctionQuadP2> >
volGaussLegendre4LagrangeQuadP2Provider;

//////////////////////////////////////////////////////////////////////////////

template <typename INTERPOLATOR>
void VolumeGaussLegendre4LagrangeQuad<INTERPOLATOR>::setMappedCoordinates()
{
  CFLogDebugMin("VolumeGaussLegendre4LagrangeQuad<INTERPOLATOR>::setMappedCoordinates(): BEGIN\n");

  _mappedCoord.resize(getNbQuadraturePoints());
  for(CFuint i = 0; i < getNbQuadraturePoints(); ++i){
    _mappedCoord[i].resize(DIM_2D);
  }

  const CFreal ref1 = 0.8611363115940526;
  const CFreal ref2 = 0.3399810435848563;

  _mappedCoord[0][KSI] = ref1;
  _mappedCoord[0][ETA] = -ref1;

  _mappedCoord[1][KSI] = ref1;
  _mappedCoord[1][ETA] = -ref2;

  _mappedCoord[2][KSI] = ref1;
  _mappedCoord[2][ETA] = ref2;

  _mappedCoord[3][KSI] = ref1;
  _mappedCoord[3][ETA] = ref1;

  _mappedCoord[4][KSI] = ref2;
  _mappedCoord[4][ETA] = -ref1;

  _mappedCoord[5][KSI] = ref2;
  _mappedCoord[5][ETA] = -ref2;

  _mappedCoord[6][KSI] = ref2;
  _mappedCoord[6][ETA] = ref2;

  _mappedCoord[7][KSI] = ref2;
  _mappedCoord[7][ETA] = ref1;

  _mappedCoord[8][KSI] = -ref2;
  _mappedCoord[8][ETA] = -ref1;

  _mappedCoord[9][KSI] = -ref2;
  _mappedCoord[9][ETA] = -ref2;

  _mappedCoord[10][KSI] = -ref2;
  _mappedCoord[10][ETA] = ref2;

  _mappedCoord[11][KSI] = -ref2;
  _mappedCoord[11][ETA] = ref1;

  _mappedCoord[12][KSI] = -ref1;
  _mappedCoord[12][ETA] = -ref1;

  _mappedCoord[13][KSI] = -ref1;
  _mappedCoord[13][ETA] = -ref2;

  _mappedCoord[14][KSI] = -ref1;
  _mappedCoord[14][ETA] = ref2;

  _mappedCoord[15][KSI] = -ref1;
  _mappedCoord[15][ETA] = ref1;

  CFLogDebugMin("VolumeGaussLegendre4LagrangeQuad<INTERPOLATOR>::setMappedCoordinates(): END\n");
}

//////////////////////////////////////////////////////////////////////////////

template <typename INTERPOLATOR>
void VolumeGaussLegendre4LagrangeQuad<INTERPOLATOR>::setWeights()
{
  _coeff.resize(getNbQuadraturePoints());

  const CFreal a = 0.3478548451374539;
  const CFreal b = 0.65211451548625461;

  _coeff[0] = a*a;
  _coeff[1] = a*b;
  _coeff[2] = a*b;
  _coeff[3] = a*a;
  _coeff[4] = a*b;
  _coeff[5] = b*b;
  _coeff[6] = b*b;
  _coeff[7] = a*b;
  _coeff[8] = a*b;
  _coeff[9] = b*b;
  _coeff[10] = b*b;
  _coeff[11] = a*b;
  _coeff[12] = a*a;
  _coeff[13] = a*b;
  _coeff[14] = a*b;
  _coeff[15] = a*a;

}

//////////////////////////////////////////////////////////////////////////////

    } // namespace ShapeFunctions



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
