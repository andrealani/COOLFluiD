// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "VolumeGaussLegendre5LagrangeQuad.hh"
#include "Framework/IntegratorImplProvider.hh"
#include "ShapeFunctions/LagrangeShapeFunctionQuadP0.hh"
#include "ShapeFunctions/LagrangeShapeFunctionQuadP1.hh"
#include "ShapeFunctions/LagrangeShapeFunctionQuadP2.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



   namespace ShapeFunctions {

//////////////////////////////////////////////////////////////////////////////

Framework::IntegratorImplProvider<Framework::VolumeIntegratorImpl,
                       VolumeGaussLegendre5LagrangeQuad<LagrangeShapeFunctionQuadP0> >
volGaussLegendre5LagrangeQuadP0Provider;

Framework::IntegratorImplProvider<Framework::VolumeIntegratorImpl,
                       VolumeGaussLegendre5LagrangeQuad<LagrangeShapeFunctionQuadP1> >
volGaussLegendre5LagrangeQuadP1Provider;

Framework::IntegratorImplProvider<Framework::VolumeIntegratorImpl,
                       VolumeGaussLegendre5LagrangeQuad<LagrangeShapeFunctionQuadP2> >
volGaussLegendre5LagrangeQuadP2Provider;

//////////////////////////////////////////////////////////////////////////////

template <typename INTERPOLATOR>
void VolumeGaussLegendre5LagrangeQuad<INTERPOLATOR>::setMappedCoordinates()
{
  CFLogDebugMin("VolumeGaussLegendre5LagrangeQuad<INTERPOLATOR>::setMappedCoordinates(): BEGIN\n");

  _mappedCoord.resize(getNbQuadraturePoints());
  for(CFuint i = 0; i < getNbQuadraturePoints(); ++i){
    _mappedCoord[i].resize(DIM_2D);
  }

  const CFreal ref1 = 0.9061798459386640;
  const CFreal ref2 = 0.5384693101056831;
  const CFreal ref3 = 0.0000000000000000;

  _mappedCoord[0][KSI] = -ref1;
  _mappedCoord[0][ETA] = -ref1;

  _mappedCoord[1][KSI] = -ref1;
  _mappedCoord[1][ETA] = -ref2;

  _mappedCoord[2][KSI] = -ref1;
  _mappedCoord[2][ETA] =  ref3;

  _mappedCoord[3][KSI] = -ref1;
  _mappedCoord[3][ETA] =  ref2;

  _mappedCoord[4][KSI] = -ref1;
  _mappedCoord[4][ETA] =  ref1;

  _mappedCoord[5][KSI] = -ref2;
  _mappedCoord[5][ETA] = -ref1;

  _mappedCoord[6][KSI] = -ref2;
  _mappedCoord[6][ETA] = -ref2;

  _mappedCoord[7][KSI] = -ref2;
  _mappedCoord[7][ETA] =  ref3;

  _mappedCoord[8][KSI] = -ref2;
  _mappedCoord[8][ETA] =  ref2;

  _mappedCoord[9][KSI] = -ref2;
  _mappedCoord[9][ETA] =  ref1;

  _mappedCoord[10][KSI] =  ref3;
  _mappedCoord[10][ETA] = -ref1;

  _mappedCoord[11][KSI] =  ref3;
  _mappedCoord[11][ETA] = -ref2;

  _mappedCoord[12][KSI] = ref3;
  _mappedCoord[12][ETA] = ref3;

  _mappedCoord[13][KSI] = ref3;
  _mappedCoord[13][ETA] = ref2;

  _mappedCoord[14][KSI] = ref3;
  _mappedCoord[14][ETA] = ref1;

  _mappedCoord[15][KSI] =  ref2;
  _mappedCoord[15][ETA] = -ref1;

  _mappedCoord[16][KSI] =  ref2;
  _mappedCoord[16][ETA] = -ref2;

  _mappedCoord[17][KSI] = ref2;
  _mappedCoord[17][ETA] = ref3;

  _mappedCoord[18][KSI] = ref2;
  _mappedCoord[18][ETA] = ref2;

  _mappedCoord[19][KSI] = ref2;
  _mappedCoord[19][ETA] = ref1;

  _mappedCoord[20][KSI] =  ref1;
  _mappedCoord[20][ETA] = -ref1;

  _mappedCoord[21][KSI] =  ref1;
  _mappedCoord[21][ETA] = -ref2;

  _mappedCoord[22][KSI] = ref1;
  _mappedCoord[22][ETA] = ref3;

  _mappedCoord[23][KSI] = ref1;
  _mappedCoord[23][ETA] = ref2;

  _mappedCoord[24][KSI] = ref1;
  _mappedCoord[24][ETA] = ref1;

  CFLogDebugMin("VolumeGaussLegendre5LagrangeQuad<INTERPOLATOR>::setMappedCoordinates(): END\n");
}

//////////////////////////////////////////////////////////////////////////////

template <typename INTERPOLATOR>
void VolumeGaussLegendre5LagrangeQuad<INTERPOLATOR>::setWeights()
{
  _coeff.resize(getNbQuadraturePoints());

  const CFreal a = 0.2369268850561891;
  const CFreal b = 0.4786286704993665;
  const CFreal c = 0.5688888888888889;

  _coeff[0] = a*a;
  _coeff[1] = a*b;
  _coeff[2] = a*c;
  _coeff[3] = a*b;
  _coeff[4] = a*a;
  _coeff[5] = b*a;
  _coeff[6] = b*b;
  _coeff[7] = b*c;
  _coeff[8] = b*b;
  _coeff[9] = b*a;
  _coeff[10] = c*a;
  _coeff[11] = c*b;
  _coeff[12] = c*c;
  _coeff[13] = c*b;
  _coeff[14] = c*a;
  _coeff[15] = b*a;
  _coeff[16] = b*b;
  _coeff[17] = b*c;
  _coeff[18] = b*b;
  _coeff[19] = b*a;
  _coeff[20] = a*a;
  _coeff[21] = a*b;
  _coeff[22] = a*c;
  _coeff[23] = a*b;
  _coeff[24] = a*a;

}

//////////////////////////////////////////////////////////////////////////////

    } // namespace ShapeFunctions



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
