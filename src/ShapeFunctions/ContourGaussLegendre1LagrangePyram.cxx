// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "ContourGaussLegendre1LagrangePyram.hh"
#include "Framework/IntegratorImplProvider.hh"
#include "ShapeFunctions/LagrangeShapeFunctionPyramP1.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



   namespace ShapeFunctions {

//////////////////////////////////////////////////////////////////////////////

Framework::IntegratorImplProvider<Framework::ContourIntegratorImpl,
                       ContourGaussLegendre1LagrangePyram<LagrangeShapeFunctionPyramP1> >
conGaussLegendre1LagrangePyramP1Provider;

//////////////////////////////////////////////////////////////////////////////

template <typename INTERPOLATOR>
void ContourGaussLegendre1LagrangePyram<INTERPOLATOR>::setMappedCoordinates()
{
  const CFuint x = 0;
  const CFuint y = 1;
  const CFuint z = 2;

  const CFuint nbPoints = getNbQuadraturePoints();
  _mappedCoord.resize(nbPoints);
  for (CFuint i = 0; i < nbPoints; ++i) {
    _mappedCoord[i].resize(DIM_3D);
  }

  const CFreal onethrd = 1.0 / 3.0;
  const CFreal twothrd = 2.0 / 3.0;

  // FACE 014
  _mappedCoord[0][x] = twothrd;
  _mappedCoord[0][y] = 0.0;
  _mappedCoord[0][z] =-onethrd;

  // FACE 124
  _mappedCoord[1][x] = 0.0;
  _mappedCoord[1][y] = twothrd;
  _mappedCoord[1][z] =-onethrd;

  // FACE 234
  _mappedCoord[2][x] =-twothrd;
  _mappedCoord[2][y] = 0.0;
  _mappedCoord[2][z] =-onethrd;

  // FACE 043
  _mappedCoord[3][x] = 0.0;
  _mappedCoord[3][y] =-twothrd;
  _mappedCoord[3][z] =-onethrd;

  // FACE 0321
  _mappedCoord[4][x] =  0.0;
  _mappedCoord[4][y] =  0.0;
  _mappedCoord[4][z] = -1.0;
}

//////////////////////////////////////////////////////////////////////////////

template <typename INTERPOLATOR>
void ContourGaussLegendre1LagrangePyram<INTERPOLATOR>::setWeights()
{
  const CFint nbFaces = 5;
  const CFint nbGaussPtsPerFace = 1;

  cf_assert(nbFaces*nbGaussPtsPerFace == getNbQuadraturePoints());

  _coeff.resize(nbFaces);
  for (CFuint iFace = 0; iFace < _coeff.size(); ++iFace) {
    _coeff[iFace].resize(nbGaussPtsPerFace);
    _coeff[iFace][0] = 0.5;
  }
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace ShapeFunctions



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
