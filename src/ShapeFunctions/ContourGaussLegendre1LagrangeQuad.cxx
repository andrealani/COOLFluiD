// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "ContourGaussLegendre1LagrangeQuad.hh"
#include "Framework/IntegratorImplProvider.hh"
#include "ShapeFunctions/LagrangeShapeFunctionQuadP1.hh"


//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



   namespace ShapeFunctions {

//////////////////////////////////////////////////////////////////////////////

Framework::IntegratorImplProvider<Framework::ContourIntegratorImpl,
                       ContourGaussLegendre1LagrangeQuad<LagrangeShapeFunctionQuadP1> >
conGaussLegendre1LagrangeQuadP1Provider;

//////////////////////////////////////////////////////////////////////////////

template <typename INTERPOLATOR>
void ContourGaussLegendre1LagrangeQuad<INTERPOLATOR>::setMappedCoordinates()
{
  _mappedCoord.resize(4);
  for (CFuint i = 0; i < _mappedCoord.size(); ++i) {
    _mappedCoord[i].resize(2);
  }

  _mappedCoord[0][0] = 0.;
  _mappedCoord[0][1] = -1.;

  _mappedCoord[1][0] = 1.;
  _mappedCoord[1][1] = 0.;

  _mappedCoord[2][0] = 0.;
  _mappedCoord[2][1] = 1.;

  _mappedCoord[3][0] = -1.;
  _mappedCoord[3][1] = 0.;
}

//////////////////////////////////////////////////////////////////////////////

template <typename INTERPOLATOR>
void ContourGaussLegendre1LagrangeQuad<INTERPOLATOR>::setWeights()
{
  const CFint nbFaces = 4;
  const CFint nbGaussPtsPerFace = 1;

  cf_assert(nbFaces*nbGaussPtsPerFace == getNbQuadraturePoints());

  // sum of each coefficient in face must give the
  // area of the face in the reference element, which is 2.
  const CFreal two = 2.0;

  _coeff.resize(nbFaces);

  for (CFuint iFace = 0; iFace < _coeff.size(); ++iFace) {
    _coeff[iFace].resize(nbGaussPtsPerFace);
    _coeff[iFace][0] = two;
  }
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace ShapeFunctions



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
