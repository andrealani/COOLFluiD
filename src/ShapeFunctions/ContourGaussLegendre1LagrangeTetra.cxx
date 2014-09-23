// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "ContourGaussLegendre1LagrangeTetra.hh"
#include "Framework/IntegratorImplProvider.hh"
#include "ShapeFunctions/LagrangeShapeFunctionTetraP1.hh"
#include "ShapeFunctions/LagrangeShapeFunctionTetraP2.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



   namespace ShapeFunctions {

//////////////////////////////////////////////////////////////////////////////

Framework::IntegratorImplProvider<Framework::ContourIntegratorImpl,
                       ContourGaussLegendre1LagrangeTetra<LagrangeShapeFunctionTetraP1> >
conGaussLegendre1LagrangeTetraP1Provider;

Framework::IntegratorImplProvider<Framework::ContourIntegratorImpl,
                       ContourGaussLegendre1LagrangeTetra<LagrangeShapeFunctionTetraP2> >
conGaussLegendre1LagrangeTetraP2Provider;

//////////////////////////////////////////////////////////////////////////////

template <typename INTERPOLATOR>
void ContourGaussLegendre1LagrangeTetra<INTERPOLATOR>::setMappedCoordinates()
{
  const CFuint nbPoints = getNbQuadraturePoints();
  _mappedCoord.resize(nbPoints);
  for (CFuint i = 0; i < nbPoints; ++i) {
    _mappedCoord[i].resize(DIM_3D);
  }

  const CFreal t0 = 1.0 / 3.0;

  // first face
  _mappedCoord[0][KSI] = t0;
  _mappedCoord[0][ETA] = t0;
  _mappedCoord[0][ZTA] = 0.;

  // second face
  _mappedCoord[1][KSI] = t0;
  _mappedCoord[1][ETA] = 0.;
  _mappedCoord[1][ZTA] = t0;

  // third face
  _mappedCoord[2][KSI] = t0;
  _mappedCoord[2][ETA] = t0;
  _mappedCoord[2][ZTA] = t0;

  // fourth face
  _mappedCoord[3][KSI] = 0.;
  _mappedCoord[3][ETA] = t0;
  _mappedCoord[3][ZTA] = t0;
}

//////////////////////////////////////////////////////////////////////////////

template <typename INTERPOLATOR>
void ContourGaussLegendre1LagrangeTetra<INTERPOLATOR>::setWeights()
{
  const CFuint nbFaces = INTERPOLATOR::getNbFaces();
  const CFuint nbGaussPtsPerFace = 1;

  cf_assert(nbFaces*nbGaussPtsPerFace == getNbQuadraturePoints());

  // resizing the ceff vector
  _coeff.resize(nbFaces);
  for (CFuint iFace = 0; iFace < _coeff.size(); ++iFace) {
    _coeff[iFace].resize(nbGaussPtsPerFace);
  }

  // these are the areas of the reference element faces
  _coeff[0][0] = 0.5;           // first face
  _coeff[1][0] = 0.5;           // second face
  _coeff[2][0] = sqrt(3.0)/2.0; // third face a=s^2*sqrt(3)/4
  _coeff[3][0] = 0.5;           // fourth face
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace ShapeFunctions



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
