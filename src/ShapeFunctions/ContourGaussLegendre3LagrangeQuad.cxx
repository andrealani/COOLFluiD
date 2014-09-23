// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "ContourGaussLegendre3LagrangeQuad.hh"
#include "Framework/IntegratorImplProvider.hh"
#include "ShapeFunctions/LagrangeShapeFunctionQuadP1.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



   namespace ShapeFunctions {

//////////////////////////////////////////////////////////////////////////////

Framework::IntegratorImplProvider<Framework::ContourIntegratorImpl,
                       ContourGaussLegendre3LagrangeQuad<LagrangeShapeFunctionQuadP1> >
conGaussLegendre3LagrangeQuadP1Provider;

//////////////////////////////////////////////////////////////////////////////

template <typename INTERPOLATOR>
void ContourGaussLegendre3LagrangeQuad<INTERPOLATOR>::setMappedCoordinates()
{
  const CFuint nbGaussPts = getNbQuadraturePoints();
  _mappedCoord.resize(nbGaussPts);
  for (CFuint i = 0; i < nbGaussPts; ++i) {
    _mappedCoord[i].resize(2);
  }

  const CFreal overSq3 = 1./sqrt(3.0);

  // Face 01
  _mappedCoord[0][0] = -overSq3;
  _mappedCoord[0][1] = -1;

  _mappedCoord[1][0] = overSq3;
  _mappedCoord[1][1] = -1.;


  // Face 12
  _mappedCoord[2][0] = 1.;
  _mappedCoord[2][1] = -overSq3;

  _mappedCoord[3][0] = 1.;
  _mappedCoord[3][1] = overSq3;

  // Face 23
  _mappedCoord[4][0] = overSq3;
  _mappedCoord[4][1] = 1.;

  _mappedCoord[5][0] = -overSq3;
  _mappedCoord[5][1] = 1.;


  // Face 30
  _mappedCoord[6][0] = -1.;
  _mappedCoord[6][1] = overSq3;

  _mappedCoord[7][0] = -1.;
  _mappedCoord[7][1] = -overSq3;
}

//////////////////////////////////////////////////////////////////////////////

template <typename INTERPOLATOR>
void ContourGaussLegendre3LagrangeQuad<INTERPOLATOR>::setWeights()
{
  const CFint nbFaces = 4;
  const CFint nbGaussPtsPerFace = 2;

  cf_assert(nbFaces*nbGaussPtsPerFace == getNbQuadraturePoints());

  // sum of each coefficient in face must give the
  // area of the face in the reference element, which is 2.
  const CFreal one = 1.0;

  _coeff.resize(nbFaces);
  for (CFuint iFace = 0; iFace < _coeff.size(); ++iFace) {
    _coeff[iFace].resize(nbGaussPtsPerFace);
    _coeff[iFace][0] = one;
    _coeff[iFace][1] = one;
  }
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace ShapeFunctions



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
