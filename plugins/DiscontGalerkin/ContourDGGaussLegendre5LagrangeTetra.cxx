#include "ContourDGGaussLegendre5LagrangeTetra.hh"
#include "Framework/IntegratorImplProvider.hh"
#include "ShapeFunctions/LagrangeShapeFunctionTetraP1.hh"
#include "ShapeFunctions/LagrangeShapeFunctionTetraP2.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



   namespace ShapeFunctions {

//////////////////////////////////////////////////////////////////////////////

Framework::IntegratorImplProvider<Framework::ContourIntegratorImpl,
                       ContourDGGaussLegendre5LagrangeTetra<LagrangeShapeFunctionTetraP1> >
conDGGaussLegendre5LagrangeTetraP1Provider;

Framework::IntegratorImplProvider<Framework::ContourIntegratorImpl,
                       ContourDGGaussLegendre5LagrangeTetra<LagrangeShapeFunctionTetraP2> >
conDGGaussLegendre5LagrangeTetraP2Provider;

//////////////////////////////////////////////////////////////////////////////

template <typename INTERPOLATOR>
void ContourDGGaussLegendre5LagrangeTetra<INTERPOLATOR>::setMappedCoordinates()
{

  const CFint x = 0;
  const CFint y = 1;
  const CFint z = 2;

  _mappedCoord.resize(getNbQuadraturePoints());
  for (CFuint i = 0; i < getNbQuadraturePoints(); ++i) {
    _mappedCoord[i].resize(DIM_3D);
  }

  const CFreal alpha  = 1./3.;
  const CFreal alpha1 = 0.0597158717;
  const CFreal alpha2 = 0.7974269853;
  const CFreal beta1  = 0.4701420641;
  const CFreal beta2  = 0.1012865073;

  // fourth face
  _mappedCoord[0][x] = alpha;
  _mappedCoord[0][y] = alpha;
  _mappedCoord[0][z] = 0.;

  _mappedCoord[1][x] = beta1;
  _mappedCoord[1][y] = beta1;
  _mappedCoord[1][z] = 0.;

  _mappedCoord[2][x] = beta1;
  _mappedCoord[2][y] = alpha1;
  _mappedCoord[2][z] = 0.;

  _mappedCoord[3][x] = alpha1;
  _mappedCoord[3][y] = beta1;
  _mappedCoord[3][z] = 0.;

  _mappedCoord[4][x] = beta2;
  _mappedCoord[4][y] = beta2;
  _mappedCoord[4][z] = 0.;

  _mappedCoord[5][x] = beta2;
  _mappedCoord[5][y] = alpha2;
  _mappedCoord[5][z] = 0.;

  _mappedCoord[6][x] = alpha2;
  _mappedCoord[6][y] = beta2;
  _mappedCoord[6][z] = 0.;

  // third face
  _mappedCoord[21][x] = alpha;
  _mappedCoord[21][y] = 0.;
  _mappedCoord[21][z] = alpha;

  _mappedCoord[22][x] = beta1;
  _mappedCoord[22][y] = 0.;
  _mappedCoord[22][z] = beta1;

  _mappedCoord[23][x] = beta1;
  _mappedCoord[23][y] = 0.;
  _mappedCoord[23][z] = alpha1;

  _mappedCoord[24][x] = alpha1;
  _mappedCoord[24][y] = 0.;
  _mappedCoord[24][z] = beta1;

  _mappedCoord[25][x] = beta2;
  _mappedCoord[25][y] = 0.;
  _mappedCoord[25][z] = beta2;

  _mappedCoord[26][x] = beta2;
  _mappedCoord[26][y] = 0.;
  _mappedCoord[26][z] = alpha2;

  _mappedCoord[27][x] = alpha2;
  _mappedCoord[27][y] = 0.;
  _mappedCoord[27][z] = beta2;

  // first face
  _mappedCoord[7][x] = alpha;
  _mappedCoord[7][y] = alpha;
  _mappedCoord[7][z] = alpha;

  _mappedCoord[8][x] = alpha1;
  _mappedCoord[8][y] = beta1;
  _mappedCoord[8][z] = beta1;

  _mappedCoord[9][x] = beta1;
  _mappedCoord[9][y] = alpha1;
  _mappedCoord[9][z] = beta1;

  _mappedCoord[10][x] = beta1;
  _mappedCoord[10][y] = beta1;
  _mappedCoord[10][z] = alpha1;

  _mappedCoord[11][x] = alpha2;
  _mappedCoord[11][y] = beta2;
  _mappedCoord[11][z] = beta2;

  _mappedCoord[12][x] = beta2;
  _mappedCoord[12][y] = alpha2;
  _mappedCoord[12][z] = beta2;

  _mappedCoord[13][x] = beta2;
  _mappedCoord[13][y] = beta2;
  _mappedCoord[13][z] = alpha2;

  // second face
  _mappedCoord[14][x] = 0.;
  _mappedCoord[14][y] = alpha;
  _mappedCoord[14][z] = alpha;

  _mappedCoord[15][x] = 0.;
  _mappedCoord[15][y] = beta1;
  _mappedCoord[15][z] = beta1;

  _mappedCoord[16][x] = 0.;
  _mappedCoord[16][y] = beta1;
  _mappedCoord[16][z] = alpha1;

  _mappedCoord[17][x] = 0.;
  _mappedCoord[17][y] = alpha1;
  _mappedCoord[17][z] = beta1;

  _mappedCoord[18][x] = 0.;
  _mappedCoord[18][y] = beta2;
  _mappedCoord[18][z] = beta2;

  _mappedCoord[19][x] = 0.;
  _mappedCoord[19][y] = beta2;
  _mappedCoord[19][z] = alpha2;

  _mappedCoord[20][x] = 0.;
  _mappedCoord[20][y] = alpha2;
  _mappedCoord[20][z] = beta2;


}

//////////////////////////////////////////////////////////////////////////////

template <typename INTERPOLATOR>
void ContourDGGaussLegendre5LagrangeTetra<INTERPOLATOR>::setWeights()
{
  const CFint nbFaces = 4;
  const CFint nbGaussPtsPerFace = 7;

  cf_assert(nbFaces*nbGaussPtsPerFace == getNbQuadraturePoints());

  _coeff.resize(nbFaces);
  for (CFuint iFace = 0; iFace < _coeff.size(); ++iFace) {
    _coeff[iFace].resize(nbGaussPtsPerFace);
    _coeff[iFace][0] =   0.225;
    _coeff[iFace][1] =   0.1323941527;
    _coeff[iFace][2] =   0.1323941527;
    _coeff[iFace][3] =   0.1323941527;
    _coeff[iFace][4] =   0.1259391805;
    _coeff[iFace][5] =   0.1259391805;
    _coeff[iFace][6] =   0.1259391805;
  }

  // these are the areas of the reference element faces
  _coeff[0] *= 0.5;           // first face
  _coeff[1] *= 0.5;           // second face
  _coeff[2] *= 0.5;//sqrt(3.0)/2.0; // third face a=s^2*sqrt(3)/4
  _coeff[3] *= 0.5;           // fourth face
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace ShapeFunctions



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
