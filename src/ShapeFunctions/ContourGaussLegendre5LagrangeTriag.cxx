// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "ContourGaussLegendre5LagrangeTriag.hh"
#include "Framework/IntegratorImplProvider.hh"
#include "ShapeFunctions/LagrangeShapeFunctionTriagP1.hh"
#include "ShapeFunctions/LagrangeShapeFunctionTriagP2.hh"

////////////////////////////////////////e//////////////////////////////

namespace COOLFluiD {

   namespace ShapeFunctions {

//////////////////////////////////////////////////////////////////////////////

Framework::IntegratorImplProvider<Framework::ContourIntegratorImpl,
    ContourGaussLegendre5LagrangeTriag<LagrangeShapeFunctionTriagP1> >
conGaussLegendre5LagrangeTriagP1Provider;

Framework::IntegratorImplProvider<Framework::ContourIntegratorImpl,
    ContourGaussLegendre5LagrangeTriag<LagrangeShapeFunctionTriagP2> >
conGaussLegendre5LagrangeTriagP2Provider;

//////////////////////////////////////////////////////////////////////////////

template <typename INTERPOLATOR>
void ContourGaussLegendre5LagrangeTriag<INTERPOLATOR>::setMappedCoordinates()
{
  const CFuint nbGaussPts = getNbQuadraturePoints();
  const CFreal mid = 0.5;
  const CFreal tdp = mid + 0.5*0.77459667;
  const CFreal tdm = mid - 0.5*0.77459667;
  
  _mappedCoord.resize(nbGaussPts);
  for (CFuint i = 0; i < nbGaussPts; ++i) {
    _mappedCoord[i].resize(DIM_2D);
  }
  
  // Face 01
  _mappedCoord[0][0] = tdm;
  _mappedCoord[0][1] = 0.;
  
  _mappedCoord[1][0] = mid;
  _mappedCoord[1][1] = 0.;
  
  _mappedCoord[2][0] = tdp;
  _mappedCoord[2][1] = 0.;
  
  // Face 12
  _mappedCoord[3][0] = tdp;
  _mappedCoord[3][1] = tdm;

  _mappedCoord[4][0] = mid;
  _mappedCoord[4][1] = mid;
  
  _mappedCoord[5][0] = tdm;
  _mappedCoord[5][1] = tdp;
  
  // Face 20
  _mappedCoord[6][0] = 0.;
  _mappedCoord[6][1] = tdp;
  
  _mappedCoord[7][0] = 0.;
  _mappedCoord[7][1] = mid;
  
  _mappedCoord[8][0] = 0.;
  _mappedCoord[8][1] = tdm;
}

//////////////////////////////////////////////////////////////////////////////

template <typename INTERPOLATOR>
void ContourGaussLegendre5LagrangeTriag<INTERPOLATOR>::setWeights()
{
  const CFuint nbFaces = 3;
  const CFuint nbGaussPtsPerFace = 3;

  _coeff.resize(nbFaces);
  for (CFuint iFace = 0; iFace < _coeff.size(); ++iFace) {
    _coeff[iFace].resize(nbGaussPtsPerFace);
    _coeff[iFace][0] = 0.5*0.555555555;
    _coeff[iFace][1] = 0.5*0.888888889;
    _coeff[iFace][2] = 0.5*0.555555555;
  }
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace ShapeFunctions



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
