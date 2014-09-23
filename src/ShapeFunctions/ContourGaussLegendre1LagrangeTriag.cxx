// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "ContourGaussLegendre1LagrangeTriag.hh"
#include "Framework/IntegratorImplProvider.hh"
#include "ShapeFunctions/LagrangeShapeFunctionTriagP1.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



   namespace ShapeFunctions {

//////////////////////////////////////////////////////////////////////////////

Framework::IntegratorImplProvider<Framework::ContourIntegratorImpl,
                       ContourGaussLegendre1LagrangeTriag<LagrangeShapeFunctionTriagP1> >
conGaussLegendre1LagrangeTriagP1Provider;

//////////////////////////////////////////////////////////////////////////////

template <typename INTERPOLATOR>
void ContourGaussLegendre1LagrangeTriag<INTERPOLATOR>::setMappedCoordinates()
{
  _mappedCoord.resize(getNbQuadraturePoints());
  for (CFuint i = 0; i < _mappedCoord.size(); ++i) {
    _mappedCoord[i].resize(2);
  }

  _mappedCoord[0][KSI] = 0.5;
  _mappedCoord[0][ETA] = 0.;
  _mappedCoord[1][KSI] = 0.5;
  _mappedCoord[1][ETA] = 0.5;
  _mappedCoord[2][KSI] = 0.;
  _mappedCoord[2][ETA] = 0.5;
}

//////////////////////////////////////////////////////////////////////////////

template <typename INTERPOLATOR>
void ContourGaussLegendre1LagrangeTriag<INTERPOLATOR>::setWeights()
{
  const CFint nbFaces = 3;
  const CFint nbGaussPtsPerFace = 1;

  cf_assert(nbFaces*nbGaussPtsPerFace == getNbQuadraturePoints());

  // sum of each coefficient in face must give the
  // area of the face in the reference element
  // face 1 = sqrt(2.0)
  // face 0 and face 2 = 1.0

  _coeff.resize(nbFaces);
  for (CFuint iFace = 0; iFace < _coeff.size(); ++iFace) {
    _coeff[iFace].resize(nbGaussPtsPerFace);
  }

  // face 0
  _coeff[0][0] = 1.0;

  // face 1
  _coeff[1][0] = std::sqrt(2.0);

  // face 2
  _coeff[2][0] = 1.0;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace ShapeFunctions



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
