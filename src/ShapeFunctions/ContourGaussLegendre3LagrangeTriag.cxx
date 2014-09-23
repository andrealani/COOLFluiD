// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "ContourGaussLegendre3LagrangeTriag.hh"
#include "Framework/IntegratorImplProvider.hh"
#include "ShapeFunctions/LagrangeShapeFunctionTriagP1.hh"
#include "ShapeFunctions/LagrangeShapeFunctionTriagP2.hh"

////////////////////////////////////////e//////////////////////////////

namespace COOLFluiD {



   namespace ShapeFunctions {

//////////////////////////////////////////////////////////////////////////////

Framework::IntegratorImplProvider<Framework::ContourIntegratorImpl,
    ContourGaussLegendre3LagrangeTriag<LagrangeShapeFunctionTriagP1> >
conGaussLegendre3LagrangeTriagP1Provider;

Framework::IntegratorImplProvider<Framework::ContourIntegratorImpl,
    ContourGaussLegendre3LagrangeTriag<LagrangeShapeFunctionTriagP2> >
conGaussLegendre3LagrangeTriagP2Provider;

//////////////////////////////////////////////////////////////////////////////

template <typename INTERPOLATOR>
void ContourGaussLegendre3LagrangeTriag<INTERPOLATOR>::setMappedCoordinates()
{
  const CFuint nbGaussPts = getNbQuadraturePoints();
  const CFreal sqrt3 = sqrt(3.0);
  const CFreal tdp = (3.0 + sqrt3) / 6.0;
  const CFreal tdm = (3.0 - sqrt3) / 6.0;
  
  _mappedCoord.resize(nbGaussPts);
  for (CFuint i = 0; i < nbGaussPts; ++i) {
    _mappedCoord[i].resize(DIM_2D);
  }

  // Face 01
  _mappedCoord[0][0] = tdm;
  _mappedCoord[0][1] = 0.;

  _mappedCoord[1][0] = tdp;
  _mappedCoord[1][1] = 0.;

  // Face 12
  _mappedCoord[2][0] = tdp;
  _mappedCoord[2][1] = tdm;

  _mappedCoord[3][0] = tdm;
  _mappedCoord[3][1] = tdp;

  // Face 20
  _mappedCoord[4][0] = 0.;
  _mappedCoord[4][1] = tdp;

  _mappedCoord[5][0] = 0.;
  _mappedCoord[5][1] = tdm;
}

//////////////////////////////////////////////////////////////////////////////

template <typename INTERPOLATOR>
void ContourGaussLegendre3LagrangeTriag<INTERPOLATOR>::setWeights()
{
  const CFuint nbFaces = 3;
  const CFuint nbGaussPtsPerFace = 2;
 // sum of each coefficient in face must give the
  // area of the face in the reference element
  // face 1 = sqrt(2.0)
  // face 0 and face 2 = 1.0
  _coeff.resize(nbFaces);
  for (CFuint iFace = 0; iFace < _coeff.size(); ++iFace) {
    _coeff[iFace].resize(nbGaussPtsPerFace);
    _coeff[iFace][0] = 0.5;
    _coeff[iFace][1] = 0.5;
  }

// face 0
  _coeff[0] *= 1.0;

  // face 1
  _coeff[1] *= std::sqrt(2.0);

  // face 2
  _coeff[2] *= 1.0;

}

//////////////////////////////////////////////////////////////////////////////

    } // namespace ShapeFunctions



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
