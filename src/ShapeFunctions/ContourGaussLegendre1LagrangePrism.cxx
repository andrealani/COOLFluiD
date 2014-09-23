// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "ContourGaussLegendre1LagrangePrism.hh"
#include "Framework/IntegratorImplProvider.hh"
#include "ShapeFunctions/LagrangeShapeFunctionPrismP1.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



   namespace ShapeFunctions {

//////////////////////////////////////////////////////////////////////////////

Framework::IntegratorImplProvider<Framework::ContourIntegratorImpl,
                       ContourGaussLegendre1LagrangePrism<LagrangeShapeFunctionPrismP1> >
conGaussLegendre1LagrangePrismP1Provider;

//////////////////////////////////////////////////////////////////////////////

template <typename INTERPOLATOR>
void ContourGaussLegendre1LagrangePrism<INTERPOLATOR>::setMappedCoordinates()
{
  const CFuint nbPoints = getNbQuadraturePoints();
  _mappedCoord.resize(nbPoints);
  for (CFuint i = 0; i < nbPoints; ++i) {
    _mappedCoord[i].resize(DIM_3D);
  }

  const CFreal p0 = 0.5;
  const CFreal p1 = 1.0 / 3.0;

  // FACE 0 || FACE 021
  _mappedCoord[0][KSI] = p1;
  _mappedCoord[0][ETA] = p1;
  _mappedCoord[0][ZTA] = -1.0;

  // FACE 1 || FACE 345
  _mappedCoord[1][KSI] = p1;
  _mappedCoord[1][ETA] = p1;
  _mappedCoord[1][ZTA] = 1.0;

  // FACE 2 || FACE 0143
  _mappedCoord[2][KSI] = p0;
  _mappedCoord[2][ETA] = 0.;
  _mappedCoord[2][ZTA] = 0.;

  // FACE 3 || FACE 1254
  _mappedCoord[3][KSI] = p0;
  _mappedCoord[3][ETA] = p0;
  _mappedCoord[3][ZTA] = 0.;

  // Face 4 || Face 0352
  _mappedCoord[4][KSI] = 0.;
  _mappedCoord[4][ETA] = p0;
  _mappedCoord[4][ZTA] = 0.;
}

//////////////////////////////////////////////////////////////////////////////

template <typename INTERPOLATOR>
void ContourGaussLegendre1LagrangePrism<INTERPOLATOR>::setWeights()
{
  const CFint nbFaces = 5;
  const CFint nbGaussPtsPerFace = 1;

  cf_assert(nbFaces*nbGaussPtsPerFace == getNbQuadraturePoints());

  _coeff.resize(nbFaces);
  for (CFuint iFace = 0; iFace < _coeff.size(); ++iFace) {
    _coeff[iFace].resize(nbGaussPtsPerFace);
  }

  // sum of each coefficient in face must give the
  // area of the face in the reference element

  _coeff[0][0] = 0.5; // face 0
  _coeff[1][0] = 0.5; // face 1
  _coeff[2][0] = 2.0; // face 2
  _coeff[3][0] = 2.0*sqrt(2.0); // face 3
  _coeff[4][0] = 2.0; // face 4
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace ShapeFunctions



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
