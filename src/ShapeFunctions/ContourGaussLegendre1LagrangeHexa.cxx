// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "ContourGaussLegendre1LagrangeHexa.hh"
#include "Framework/IntegratorImplProvider.hh"
#include "ShapeFunctions/LagrangeShapeFunctionHexaP1.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

   namespace ShapeFunctions {

//////////////////////////////////////////////////////////////////////////////

Framework::IntegratorImplProvider<Framework::ContourIntegratorImpl,
                       ContourGaussLegendre1LagrangeHexa<LagrangeShapeFunctionHexaP1> >
conGaussLegendre1LagrangeHexaP1Provider;

//////////////////////////////////////////////////////////////////////////////

template <typename INTERPOLATOR>
void ContourGaussLegendre1LagrangeHexa<INTERPOLATOR>::setMappedCoordinates()
{
  const CFuint nbPoints = getNbQuadraturePoints();
  _mappedCoord.resize(nbPoints);
  for (CFuint i = 0; i < nbPoints; ++i) {
    _mappedCoord[i].resize(DIM_3D);
  }

  const CFreal p0 = 1.0;

  // Face 0 || Face 0321
  _mappedCoord[0][XX] = 0.;
  _mappedCoord[0][YY] = 0.;
  _mappedCoord[0][ZZ] =-p0;

  // Face 1 || Face 4567
  _mappedCoord[1][XX] = 0.;
  _mappedCoord[1][YY] = 0.;
  _mappedCoord[1][ZZ] = p0;

  // Face 2 || Face 0154
  _mappedCoord[2][XX] = 0.;
  _mappedCoord[2][YY] =-p0;
  _mappedCoord[2][ZZ] = 0.;

  // Face 3 || Face 1265
  _mappedCoord[3][XX] = p0;
  _mappedCoord[3][YY] = 0.;
  _mappedCoord[3][ZZ] = 0.;

  // Face 4 || Face 3762
  _mappedCoord[4][XX] = 0.;
  _mappedCoord[4][YY] = p0;
  _mappedCoord[4][ZZ] = 0.;

  // Face 5 || Face 0473
  _mappedCoord[5][XX] =-p0;
  _mappedCoord[5][YY] = 0.;
  _mappedCoord[5][ZZ] = 0.;
}

//////////////////////////////////////////////////////////////////////////////

template <typename INTERPOLATOR>
void ContourGaussLegendre1LagrangeHexa<INTERPOLATOR>::setWeights()
{
  const CFint nbFaces = 6;
  const CFint nbGaussPtsPerFace = 1;

  cf_assert(nbFaces*nbGaussPtsPerFace == getNbQuadraturePoints());

  // sum of each coefficient in face must give the
  // area of the face in the reference element, which is 4.
  const CFreal weight = 4.0;

  _coeff.resize(nbFaces);
  for (CFuint iFace = 0; iFace < _coeff.size(); ++iFace) {
    _coeff[iFace].resize(nbGaussPtsPerFace);
    _coeff[iFace][0] = weight;
  }
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace ShapeFunctions



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
