// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "ContourGaussLegendre3LagrangeHexa.hh"
#include "Framework/IntegratorImplProvider.hh"
#include "ShapeFunctions/LagrangeShapeFunctionHexaP1.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



   namespace ShapeFunctions {

//////////////////////////////////////////////////////////////////////////////

Framework::IntegratorImplProvider<Framework::ContourIntegratorImpl,
                       ContourGaussLegendre3LagrangeHexa<LagrangeShapeFunctionHexaP1> >
conGaussLegendre3LagrangeHexaP1Provider;

//////////////////////////////////////////////////////////////////////////////

template <typename INTERPOLATOR>
void ContourGaussLegendre3LagrangeHexa<INTERPOLATOR>::setMappedCoordinates()
{
  // resize the mapped coordinates
  _mappedCoord.resize(getNbQuadraturePoints());
  for (CFuint i = 0; i < getNbQuadraturePoints(); ++i) {
    _mappedCoord[i].resize(DIM_3D);
  }

  const CFreal delta = 0.5773502692;

  // Face 0 || Face 0321
  _mappedCoord[0][XX] =   delta;
  _mappedCoord[0][YY] =   delta;
  _mappedCoord[0][ZZ] = - 1.0;

  _mappedCoord[1][XX] =   delta;
  _mappedCoord[1][YY] = - delta;
  _mappedCoord[1][ZZ] = - 1.0;

  _mappedCoord[2][XX] = - delta;
  _mappedCoord[2][YY] =   delta;
  _mappedCoord[2][ZZ] = - 1.0;

  _mappedCoord[3][XX] = - delta;
  _mappedCoord[3][YY] = - delta;
  _mappedCoord[3][ZZ] = - 1.0;

  // Face 1 || Face 4567
  _mappedCoord[4][XX] =   delta;
  _mappedCoord[4][YY] =   delta;
  _mappedCoord[4][ZZ] =   1.0;

  _mappedCoord[5][XX] =   delta;
  _mappedCoord[5][YY] = - delta;
  _mappedCoord[5][ZZ] =   1.0;

  _mappedCoord[6][XX] = - delta;
  _mappedCoord[6][YY] =   delta;
  _mappedCoord[6][ZZ] =   1.0;

  _mappedCoord[7][XX] = - delta;
  _mappedCoord[7][YY] = - delta;
  _mappedCoord[7][ZZ] =   1.0;

  // Face 2 || Face 0154
  _mappedCoord[8][XX]  =   delta;
  _mappedCoord[8][YY]  = - 1.0;
  _mappedCoord[8][ZZ]  =   delta;

  _mappedCoord[9][XX]  = - delta;
  _mappedCoord[9][YY]  = - 1.0;
  _mappedCoord[9][ZZ]  =   delta;

  _mappedCoord[10][XX] =   delta;
  _mappedCoord[10][YY] = - 1.0;
  _mappedCoord[10][ZZ] = - delta;

  _mappedCoord[11][XX] = - delta;
  _mappedCoord[11][YY] = - 1.0;
  _mappedCoord[11][ZZ] = - delta;

  // Face 3 || Face 1265
  _mappedCoord[12][XX] =   1.0;
  _mappedCoord[12][YY] =   delta;
  _mappedCoord[12][ZZ] =   delta;

  _mappedCoord[13][XX] =   1.0;
  _mappedCoord[13][YY] =   delta;
  _mappedCoord[13][ZZ] = - delta;

  _mappedCoord[14][XX] =   1.0;
  _mappedCoord[14][YY] = - delta;
  _mappedCoord[14][ZZ] =   delta;

  _mappedCoord[15][XX] =   1.0;
  _mappedCoord[15][YY] = - delta;
  _mappedCoord[15][ZZ] = - delta;

  // Face 4 || Face 3762
  _mappedCoord[16][XX] =   delta;
  _mappedCoord[16][YY] =   1.0;
  _mappedCoord[16][ZZ] =   delta;

  _mappedCoord[17][XX] = - delta;
  _mappedCoord[17][YY] =   1.0;
  _mappedCoord[17][ZZ] =   delta;

  _mappedCoord[18][XX] =   delta;
  _mappedCoord[18][YY] =   1.0;
  _mappedCoord[18][ZZ] = - delta;

  _mappedCoord[19][XX] = - delta;
  _mappedCoord[19][YY] =   1.0;
  _mappedCoord[19][ZZ] = - delta;

  // Face 5 || Face 0473
  _mappedCoord[20][XX] = - 1.0;
  _mappedCoord[20][YY] =   delta;
  _mappedCoord[20][ZZ] =   delta;

  _mappedCoord[21][XX] = - 1.0;
  _mappedCoord[21][YY] =   delta;
  _mappedCoord[21][ZZ] = - delta;

  _mappedCoord[22][XX] = - 1.0;
  _mappedCoord[22][YY] = - delta;
  _mappedCoord[22][ZZ] =   delta;

  _mappedCoord[23][XX] = - 1.0;
  _mappedCoord[23][YY] = - delta;
  _mappedCoord[23][ZZ] = - delta;
}

//////////////////////////////////////////////////////////////////////////////

template <typename INTERPOLATOR>
void ContourGaussLegendre3LagrangeHexa<INTERPOLATOR>::setWeights()
{
  const CFint nbFaces = 6;
  const CFint nbGaussPtsPerFace = 4;

  cf_assert(nbFaces*nbGaussPtsPerFace == getNbQuadraturePoints());

  // sum of each coefficient in face must give the
  // area of the face in the reference element, which is 4.
  const CFreal one = 1.0;

  _coeff.resize(nbFaces);
  for (CFuint iFace = 0; iFace < _coeff.size(); ++iFace) {
    _coeff[iFace].resize(nbGaussPtsPerFace);
    _coeff[iFace][0] = one;
    _coeff[iFace][1] = one;
    _coeff[iFace][2] = one;
    _coeff[iFace][3] = one;
  }
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace ShapeFunctions



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
