// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "ContourGaussLegendre3LagrangeTetra.hh"
#include "Framework/IntegratorImplProvider.hh"
#include "ShapeFunctions/LagrangeShapeFunctionTetraP1.hh"
#include "ShapeFunctions/LagrangeShapeFunctionTetraP2.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



   namespace ShapeFunctions {

//////////////////////////////////////////////////////////////////////////////

Framework::IntegratorImplProvider<Framework::ContourIntegratorImpl,
                       ContourGaussLegendre3LagrangeTetra<LagrangeShapeFunctionTetraP1> >
conGaussLegendre3LagrangeTetraP1Provider;

Framework::IntegratorImplProvider<Framework::ContourIntegratorImpl,
                       ContourGaussLegendre3LagrangeTetra<LagrangeShapeFunctionTetraP2> >
conGaussLegendre3LagrangeTetraP2Provider;

//////////////////////////////////////////////////////////////////////////////

template <typename INTERPOLATOR>
void ContourGaussLegendre3LagrangeTetra<INTERPOLATOR>::setMappedCoordinates()
{

  const CFint x = 0;
  const CFint y = 1;
  const CFint z = 2;

  _mappedCoord.resize(getNbQuadraturePoints());
  for (CFuint i = 0; i < getNbQuadraturePoints(); ++i) {
    _mappedCoord[i].resize(DIM_3D);
  }

  const CFreal t0 = 1./3.;
  const CFreal t1 = 0.2;
  const CFreal t2 = 0.6;

  // first face
  _mappedCoord[0][x] = t0;
  _mappedCoord[0][y] = t0;
  _mappedCoord[0][z] = 0.;

  _mappedCoord[1][x] = t1;
  _mappedCoord[1][y] = t1;
  _mappedCoord[1][z] = 0.;

  _mappedCoord[2][x] = t1;
  _mappedCoord[2][y] = t2;
  _mappedCoord[2][z] = 0.;

  _mappedCoord[3][x] = t2;
  _mappedCoord[3][y] = t1;
  _mappedCoord[3][z] = 0.;

  // second face
  _mappedCoord[4][x] = t0;
  _mappedCoord[4][y] = 0.;
  _mappedCoord[4][z] = t0;

  _mappedCoord[5][x] = t1;
  _mappedCoord[5][y] = 0.;
  _mappedCoord[5][z] = t1;

  _mappedCoord[6][x] = t1;
  _mappedCoord[6][y] = 0.;
  _mappedCoord[6][z] = t2;

  _mappedCoord[7][x] = t2;
  _mappedCoord[7][y] = 0.;
  _mappedCoord[7][z] = t1;

  // third face
  _mappedCoord[8][x] = t0;
  _mappedCoord[8][y] = t0;
  _mappedCoord[8][z] = t0;

  _mappedCoord[9][x] = t2;
  _mappedCoord[9][y] = t1;
  _mappedCoord[9][z] = t1;

  _mappedCoord[10][x] = t1;
  _mappedCoord[10][y] = t2;
  _mappedCoord[10][z] = t1;

  _mappedCoord[11][x] = t1;
  _mappedCoord[11][y] = t1;
  _mappedCoord[11][z] = t2;

  // fourth face
  _mappedCoord[12][x] = 0.;
  _mappedCoord[12][y] = t0;
  _mappedCoord[12][z] = t0;

  _mappedCoord[13][x] = 0.;
  _mappedCoord[13][y] = t1;
  _mappedCoord[13][z] = t1;

  _mappedCoord[14][x] = 0.;
  _mappedCoord[14][y] = t1;
  _mappedCoord[14][z] = t2;

  _mappedCoord[15][x] = 0.;
  _mappedCoord[15][y] = t2;
  _mappedCoord[15][z] = t1;
}

//////////////////////////////////////////////////////////////////////////////

template <typename INTERPOLATOR>
void ContourGaussLegendre3LagrangeTetra<INTERPOLATOR>::setWeights()
{
  const CFint nbFaces = 4;
  const CFint nbGaussPtsPerFace = 4;

  cf_assert(nbFaces*nbGaussPtsPerFace == getNbQuadraturePoints());

  _coeff.resize(nbFaces);
  for (CFuint iFace = 0; iFace < _coeff.size(); ++iFace) {
    _coeff[iFace].resize(nbGaussPtsPerFace);
    _coeff[iFace][0] = - 27./48.;
    _coeff[iFace][1] =   25./48.;
    _coeff[iFace][2] =   25./48.;
    _coeff[iFace][3] =   25./48.;
  }

  // these are the areas of the reference element faces
  _coeff[0] *= 0.5;           // first face
  _coeff[1] *= 0.5;           // second face
  _coeff[2] *= sqrt(3.0)/2.0; // third face a=s^2*sqrt(3)/4
  _coeff[3] *= 0.5;           // fourth face
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace ShapeFunctions



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
