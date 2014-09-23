// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "ContourGaussLegendre3LagrangePyram.hh"
#include "Framework/IntegratorImplProvider.hh"
#include "ShapeFunctions/LagrangeShapeFunctionPyramP1.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



   namespace ShapeFunctions {

//////////////////////////////////////////////////////////////////////////////

Framework::IntegratorImplProvider<Framework::ContourIntegratorImpl,
                       ContourGaussLegendre3LagrangePyram<LagrangeShapeFunctionPyramP1> >
conGaussLegendre3LagrangePyramP1Provider;

//////////////////////////////////////////////////////////////////////////////

template <typename INTERPOLATOR>
void ContourGaussLegendre3LagrangePyram<INTERPOLATOR>::setMappedCoordinates()
{
  CFAUTOTRACE;

  // resize the mapped coordinates
  _mappedCoord.resize(getNbQuadraturePoints());
  for (CFuint i = 0; i < getNbQuadraturePoints(); ++i) {
    _mappedCoord[i].resize(DIM_3D);
  }

  const CFreal qdelta = 0.577350269;
  const CFreal ddelta = 0.5;

  // Face 0 || Face 0321
  _mappedCoord[0][XX] =   qdelta;
  _mappedCoord[0][YY] =   qdelta;
  _mappedCoord[0][ZZ] = - 1.0;

  _mappedCoord[1][XX] =   qdelta;
  _mappedCoord[1][YY] = - qdelta;
  _mappedCoord[1][ZZ] = - 1.0;

  _mappedCoord[2][XX] = - qdelta;
  _mappedCoord[2][YY] =   qdelta;
  _mappedCoord[2][ZZ] = - 1.0;

  _mappedCoord[3][XX] = - qdelta;
  _mappedCoord[3][YY] = - qdelta;
  _mappedCoord[3][ZZ] = - 1.0;

  // Face 1 || Face 014
  _mappedCoord[4][XX] =   1.0;
  _mappedCoord[4][YY] =   0.0;
  _mappedCoord[4][ZZ] = - 1.0;

  _mappedCoord[5][XX] =   ddelta;
  _mappedCoord[5][YY] = - ddelta;
  _mappedCoord[5][ZZ] =   0.0;

  _mappedCoord[6][XX] =   ddelta;
  _mappedCoord[6][YY] =   ddelta;
  _mappedCoord[6][ZZ] =   0.0;

   // Face 2 || Face 124
  _mappedCoord[7][XX] =   0.0;
  _mappedCoord[7][YY] =   1.0;
  _mappedCoord[7][ZZ] = - 1.0;

  _mappedCoord[8][XX] = - ddelta;
  _mappedCoord[8][YY] =   ddelta;
  _mappedCoord[8][ZZ] =   0.0;

  _mappedCoord[9][XX] =   ddelta;
  _mappedCoord[9][YY] =   ddelta;
  _mappedCoord[9][ZZ] =   0.0;

  // Face 3 || Face 234
  _mappedCoord[10][XX] = - 1.0;
  _mappedCoord[10][YY] =   0.0;
  _mappedCoord[10][ZZ] = - 1.0;

  _mappedCoord[11][XX] = - ddelta;
  _mappedCoord[11][YY] =   ddelta;
  _mappedCoord[11][ZZ] =   0.0;

  _mappedCoord[12][XX] = - ddelta;
  _mappedCoord[12][YY] = - ddelta;
  _mappedCoord[12][ZZ] =   0.0;

  // Face 4 || Face 043
  _mappedCoord[13][XX] =   0.0;
  _mappedCoord[13][YY] = - 1.0;
  _mappedCoord[13][ZZ] = - 1.0;

  _mappedCoord[14][XX] = - ddelta;
  _mappedCoord[14][YY] = - ddelta;
  _mappedCoord[14][ZZ] =   0.0;

  _mappedCoord[15][XX] =   ddelta;
  _mappedCoord[15][YY] = - ddelta;
  _mappedCoord[15][ZZ] =   0.0;
}

//////////////////////////////////////////////////////////////////////////////

template <typename INTERPOLATOR>
void ContourGaussLegendre3LagrangePyram<INTERPOLATOR>::setWeights()
{
  CFAUTOTRACE;

  const CFuint nbFaces = 5;

  const CFreal wt = 1.0 / 3.0;
  const CFreal wq = 1.0 / 4.0;

  _coeff.resize(nbFaces);

  _coeff[0].resize(4);
  _coeff[0][0] = wq;
  _coeff[0][1] = wq;
  _coeff[0][2] = wq;
  _coeff[0][3] = wq;

  for (CFuint iFace = 1; iFace < nbFaces; ++iFace) {
    _coeff[iFace].resize(3);
    _coeff[iFace][0] = wt;
    _coeff[iFace][1] = wt;
    _coeff[iFace][2] = wt;
  }
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace ShapeFunctions



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
