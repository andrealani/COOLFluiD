// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "ContourGaussLegendre3LagrangePrism.hh"
#include "Framework/IntegratorImplProvider.hh"
#include "ShapeFunctions/LagrangeShapeFunctionPrismP1.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



   namespace ShapeFunctions {

//////////////////////////////////////////////////////////////////////////////

Framework::IntegratorImplProvider<Framework::ContourIntegratorImpl,
                       ContourGaussLegendre3LagrangePrism<LagrangeShapeFunctionPrismP1> >
conGaussLegendre3LagrangePrismP1Provider;

//////////////////////////////////////////////////////////////////////////////

template <typename INTERPOLATOR>
void ContourGaussLegendre3LagrangePrism<INTERPOLATOR>::setMappedCoordinates()
{
  _mappedCoord.resize(getNbQuadraturePoints());
  for (CFuint i = 0; i < getNbQuadraturePoints(); ++i) {
    _mappedCoord[i].resize(DIM_3D);
  }

  const CFreal tdelta = 0.5;
  const CFreal qdelta = 0.5773502692;
  const CFreal q1 = 0.5 * ( 1. + qdelta);
  const CFreal q2 = 0.5 * ( 1. - qdelta);
  const CFreal scale = std::sqrt(2.0)/2.0;
  const CFreal m1 = 0.5 + ( qdelta * scale );
  const CFreal m2 = 0.5 - ( qdelta * scale );

  // Face 0 || Face 021
  _mappedCoord[0][KSI] =   tdelta;
  _mappedCoord[0][ETA] =   0.0;
  _mappedCoord[0][ZTA] = - 1.0;

  _mappedCoord[1][KSI] =   0.0;
  _mappedCoord[1][ETA] =   tdelta;
  _mappedCoord[1][ZTA] = - 1.0;

  _mappedCoord[2][KSI] =   tdelta;
  _mappedCoord[2][ETA] =   tdelta;
  _mappedCoord[2][ZTA] = - 1.0;

  // Face 1 || Face 345
  _mappedCoord[3][KSI] =   tdelta;
  _mappedCoord[3][ETA] =   0.0;
  _mappedCoord[3][ZTA] =   1.0;

  _mappedCoord[4][KSI] =   0.0;
  _mappedCoord[4][ETA] =   tdelta;
  _mappedCoord[4][ZTA] =   1.0;

  _mappedCoord[5][KSI] =   tdelta;
  _mappedCoord[5][ETA] =   tdelta;
  _mappedCoord[5][ZTA] =   1.0;

  // Face 2 || Face 0143
  _mappedCoord[6][KSI] =   q1;
  _mappedCoord[6][ETA] =   0.0;
  _mappedCoord[6][ZTA] =   qdelta;

  _mappedCoord[7][KSI] =   q2;
  _mappedCoord[7][ETA] =   0.0;
  _mappedCoord[7][ZTA] =   qdelta;

  _mappedCoord[8][KSI] =   q1;
  _mappedCoord[8][ETA] =   0.0;
  _mappedCoord[8][ZTA] = - qdelta;

  _mappedCoord[9][KSI] =   q2;
  _mappedCoord[9][ETA] =   0.0;
  _mappedCoord[9][ZTA] = - qdelta;

  // Face 3 || Face 1254
  _mappedCoord[10][KSI] =   m1;
  _mappedCoord[10][ETA] =   m2;
  _mappedCoord[10][ZTA] =   qdelta;

  _mappedCoord[11][KSI] =   m2;
  _mappedCoord[11][ETA] =   m1;
  _mappedCoord[11][ZTA] =   qdelta;

  _mappedCoord[12][KSI] =   m1;
  _mappedCoord[12][ETA] =   m2;
  _mappedCoord[12][ZTA] = - qdelta;

  _mappedCoord[13][KSI] =   m2;
  _mappedCoord[13][ETA] =   m1;
  _mappedCoord[13][ZTA] = - qdelta;

  // Face 4 || Face 0352
  _mappedCoord[14][KSI] =   0.0;
  _mappedCoord[14][ETA] =   q1;
  _mappedCoord[14][ZTA] =   qdelta;

  _mappedCoord[15][KSI] =   0.0;
  _mappedCoord[15][ETA] =   q2;
  _mappedCoord[15][ZTA] =   qdelta;

  _mappedCoord[16][KSI] =   0.0;
  _mappedCoord[16][ETA] =   q1;
  _mappedCoord[16][ZTA] = - qdelta;

  _mappedCoord[17][KSI] =   0.0;
  _mappedCoord[17][ETA] =   q2;
  _mappedCoord[17][ZTA] = - qdelta;
}

//////////////////////////////////////////////////////////////////////////////

template <typename INTERPOLATOR>
void ContourGaussLegendre3LagrangePrism<INTERPOLATOR>::setWeights()
{
  const CFuint nbFaces = 5;

  const CFreal wq  = 2.0 * 0.25;
  const CFreal wqd = 2.0*sqrt(2.0) * 0.25;
  const CFreal wt  = 0.5 * 1.0 / 3.0;

  _coeff.resize(nbFaces);

  // sum of weights in each face must equal area
  // set weights in triangular faces

  // Face 0 (Face 021) and Face 1 (Face 345)
  // Area = 0.5
  for (CFuint iFace = 0; iFace < 2; ++iFace) {
    _coeff[iFace].resize(3);
    _coeff[iFace][0] = wt;
    _coeff[iFace][1] = wt;
    _coeff[iFace][2] = wt;
  }

  // sum of weights in each face must equal area
  // set weights in quadrilateral faces

  // Face 2 || Face 0143
  // Area = 2.0
  _coeff[2].resize(4);
  _coeff[2][0] = wq;
  _coeff[2][1] = wq;
  _coeff[2][2] = wq;
  _coeff[2][3] = wq;

  // Face 3 || Face 1254 // diagonal face
  // Area = 2.0 * sqrt(2.0)
  _coeff[3].resize(4);
  _coeff[3][0] = wqd;
  _coeff[3][1] = wqd;
  _coeff[3][2] = wqd;
  _coeff[3][3] = wqd;

  // Face 2 || Face 0352
  // Area = 2.0
  _coeff[4].resize(4);
  _coeff[4][0] = wq;
  _coeff[4][1] = wq;
  _coeff[4][2] = wq;
  _coeff[4][3] = wq;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace ShapeFunctions



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
