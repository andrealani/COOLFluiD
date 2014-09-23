#include <algorithm>
#include <numeric>

#include "Common/CFLog.hh"
#include "SpectralFV/SpectralFVElementData.hh"
#include "Common/ShouldNotBeHereException.hh"
#include "MathTools/MathConsts.hh"

//////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Framework;

namespace COOLFluiD {

  namespace SpectralFV {

//////////////////////////////////////////////////////////////////////

SpectralFVElementData::SpectralFVElementData() :
  m_dimensionality(),
  m_shape(),
  m_polyOrder(),
  m_sIntegrator(),
  m_tpIntegrator(),
  m_localNodeCoord(),
  m_localFaceNodeConn(),
  m_localCVNodeConn(),
  m_localIntFaceCVConn(),
  m_extFaceLocalNorm(),
  m_volFracCV(),
  m_invVolFracCV(),
  m_svFaceLocalExtFaceConn(),
  m_localExtFaceCVConn(),
  m_extSVFaceCVConn(),
  m_extSVFaceCVConnPerOrient(),
  m_eSVFaceLocESVFaceConnPerOrient(),
  m_svNodeCoords(),
  m_svFaceNodeCoords(),
  m_svFaceNodeConn(),
  m_svFaceNodeConnPerOrientation(),
  m_svFaceNodeConnPerOrientationNoSymm(),
  m_svFaceNodeCoordsPerOrient(),
  m_svFaceNodeCoordsPerOrientNoSymm(),
  m_extFaceNodeLocalCoords(),
  m_faceFracCV(),
  m_invVolFracFaceCV(),
  m_invVolFracFaceCVPerOrient(),
  m_polyExponents(),
  m_polyCoefSFV(),
  m_initPntsCoords(),
  m_initTransfMatrix(),
  m_intFaceQuadPntNorm(),
  m_intFaceQuadPntCoords(),
  m_intFaceQuadWheights(),
  m_intFaceQuadPntPolyVals(),
  m_intQuadPntCoords(),
  m_intQuadWheights(),
  m_intQuadPntPolyVals(),
  m_extQPntWheightCoords(),
  m_extQPntCoords(),
  m_extQWheights(),
  m_extQPntPolyVals(),
  m_extFaceQPntWheightCoords(),
  m_extFaceQuadPntCoords(),
  m_extFaceQuadWheights(),
  m_extFaceQuadPntPolyVals(),
  m_extFaceQPntCoordsPerOrient(),
  m_extFaceQWheightsPerOrient(),
  m_extFaceQPntPolyValsPerOrient(),
  m_extFaceQPntPolyValsPerOrientNoSymm(),
  m_extQPntCoordsPerOrient(),
  m_extQWheightsPerOrient(),
  m_extQPntPolyValsPerOrient(),
  m_extQPntPolyValsPerOrientNoSymm(),
  m_svFaceAvgPolyVals(),
  m_svFaceAvgPolyValsPerOrient(),
  m_fluxPolyNodeCoord(),
  m_faceFlxPntsConn(),
  m_faceFluxPolyNodeWheightCoord(),
  m_fluxPolyExponents(),
  m_fluxPolyCoefs(),
  m_solInFluxPntCoef(),
  m_volTermTensor(),
  m_solInFaceFluxPntCoef(),
  m_solInFaceFluxPntCoefPerOrient(),
  m_solInFaceFluxPntCoefPerOrientNoSymm(),
  m_cvExtFaceFluxCoef(),
  m_avgSolInSVFaceCoef(),
  m_cflConvDiffRatio(),
  m_faceOutputPntCellMappedCoords(),
  m_faceOutputPntSolPolyCoef(),
  m_faceOutputPntSolDerivCoef(),
  m_faceOutputPntConn()
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////

SpectralFVElementData::~SpectralFVElementData()
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////

void SpectralFVElementData::setPolyOrder(CFPolyOrder::Type polyOrder)
{
  cf_assert(polyOrder >= 0);
  m_polyOrder = polyOrder;
  resetSpectralFVElementData();
}

//////////////////////////////////////////////////////////////////////

void SpectralFVElementData::resetSpectralFVElementData()
{
  CFAUTOTRACE;

  computeSVNodeLocalCoords();
  createLocalNodeCoord();
  createLocalFaceNodeConn();
  createLocalCVNodeConn();
  createLocalIntFaceCVConn();
  computeVolumeFractionsOfCVs();
  computeInverseVolFracsOfCVs();
  createLocalExtFaceCVConn();
  createSVFaceLocalExtFaceConn();
  createExtFaceCVConn();
  createExtFaceCVConnPerOrient();
  computeInvVolFracsOfFaceCVs();
  computeInvVolFracsOfFaceCVsPerOrient();
  createSVFaceNodeConnectivity();
  createSVFaceNodeConnectivityPerOrient();
  createSVFaceNodeConnectivityPerOrientNoSymm();
  createSVFaceNodeCoords();
  createSVFaceNodeCoordsPerOrient();
  createSVFaceNodeCoordsPerOrientNoSymm();
  computeExtFaceLocalNormals();
  createExtFaceNodeLocalCoords();
  computeFaceFractionsOfCVs();
  createPolyExponents();
  computePolyCoefs();
  computeInitPntsCoords();
  computeInitTransfMatrix();
  computeIntFaceQuadPntsData();
  computeIntQuadPntsData();
  computeLocalFaceNormals();
  computeExtFaceQuadPntsData();
  computeExtQuadPntsData();
  computeExtFaceQuadPntsDataPerOrientation();
  computeExtFaceQuadPntsDataPerOrientationNoSymm();
  computeExtQuadPntsDataPerOrientation();
  computeExtQuadPntsDataPerOrientationNoSymm();
  computeFaceAvgPolyVals();
  computeFaceAvgPolyValsPerOrient();
  setCFLConvDiffRatio();
  createFaceOutputPntCellMappedCoords();
  createFaceOutputPntSolPolyAndDerivCoef();
  createFaceOutputPntConn();


  // data for quadrature free implementation
  if (m_dimensionality > DIM_1D) // quadrature free approach does not make sense for 1D
  {
    createFluxPolyNodeCoord();
    createFaceFluxPntsConn();
    createFaceFluxPolyNodeWheightCoord();
    createFluxPolyExponents();
    createFluxPolyCoef();
    createSolInFlxPntsCoefs();
    createVolumeTermsTensor();
    createSolInFaceFlxPntsCoefs();
    createSolInFaceFlxPntsCoefsPerOrient();
    createSolInFaceFlxPntsCoefsPerOrientNoSymm();
    createCVExtFaceFluxCoef();
    createAvgSolInSVFaceCoef();
  }
}

//////////////////////////////////////////////////////////////////////

vector< vector< CFreal > >
  SpectralFVElementData::getSVPolyValsAtNode(vector< RealVector > nodeLocalCoords)
{
  const CFuint dimUns = static_cast<CFuint>(m_dimensionality);
  cf_assert(nodeLocalCoords.size() > 0);
  cf_assert(nodeLocalCoords[0].size() == dimUns);

  // number of nodes passed
  const CFuint nbrNodes = nodeLocalCoords.size();

  // number of polynomials (== number of control volumes)
  const CFuint nbrPolys = getNbrOfCVs();

  // return variable
  vector< vector< CFreal > > polyValsAtNodes(nbrNodes);

  // compute polynomial values at given nodes
  for (CFuint iNode = 0; iNode < nbrNodes; ++iNode)
  {
    // resize
    polyValsAtNodes[iNode].resize(nbrPolys);

    // loop over polynomials
    for (CFuint iPoly = 0; iPoly < nbrPolys; ++iPoly)
    {
      // loop over terms
      for (CFuint iTerm = 0; iTerm < nbrPolys; ++iTerm)
      {
        CFreal term = m_polyCoefSFV[iPoly][iTerm];

        // loop over coordinates x-y-z
        for (CFuint iCoor = 0; iCoor < dimUns; ++iCoor)
        {
          term *= pow(nodeLocalCoords[iNode][iCoor],m_polyExponents[iTerm][iCoor]);
        }

        // add term to polynomial value
        polyValsAtNodes[iNode][iPoly] += term;
      }
    }
  }

  return polyValsAtNodes;
}

//////////////////////////////////////////////////////////////////////

void SpectralFVElementData::computeInverseVolFracsOfCVs()
{
  CFAUTOTRACE;

  // number of CVs
  const CFuint nbrCVs = getNbrOfCVs();

  // compute inverse CV volume fractions
  m_invVolFracCV.resize(nbrCVs);
  for (CFuint iCV = 0; iCV < nbrCVs; ++iCV)
  {
    m_invVolFracCV[iCV] = 1.0/m_volFracCV[iCV];
  }
}

//////////////////////////////////////////////////////////////////////

void SpectralFVElementData::computeInvVolFracsOfFaceCVs()
{
  // number of SV faces
  const CFuint nbrSVFaces = m_extSVFaceCVConn.size();

  // set inverse volume fractions of face CVs
  m_invVolFracFaceCV.resize(nbrSVFaces);
  for (CFuint iFace = 0; iFace < nbrSVFaces; ++iFace)
  {
    const CFuint nbrFaceCVs = m_extSVFaceCVConn[iFace].size();
    m_invVolFracFaceCV[iFace].resize(nbrFaceCVs);
    for (CFuint iCV = 0; iCV < nbrFaceCVs; ++iCV)
    {
      const CFuint cvID = m_extSVFaceCVConn[iFace][iCV];
      m_invVolFracFaceCV[iFace][iCV] = m_invVolFracCV[cvID];
    }
  }
}

//////////////////////////////////////////////////////////////////////

void SpectralFVElementData::computeInvVolFracsOfFaceCVsPerOrient()
{
  // number of SV faces
  const CFuint nbrOrients = m_extSVFaceCVConnPerOrient.size();

  // set inverse volume fractions of face CVs
  m_invVolFracFaceCVPerOrient.resize(nbrOrients);
  for (CFuint iOrient = 0; iOrient < nbrOrients; ++iOrient)
  {
    const CFuint nbrFaceCVs = m_extSVFaceCVConnPerOrient[iOrient].size();
    m_invVolFracFaceCVPerOrient[iOrient].resize(nbrFaceCVs);
    for (CFuint iCV = 0; iCV < nbrFaceCVs; ++iCV)
    {
      m_invVolFracFaceCVPerOrient[iOrient][iCV].resize(2);

      const CFuint cvIDL = m_extSVFaceCVConnPerOrient[iOrient][iCV][LEFT ];
      m_invVolFracFaceCVPerOrient[iOrient][iCV][LEFT ] = m_invVolFracCV[cvIDL];

      const CFuint cvIDR = m_extSVFaceCVConnPerOrient[iOrient][iCV][RIGHT];
      m_invVolFracFaceCVPerOrient[iOrient][iCV][RIGHT] = m_invVolFracCV[cvIDR];
    }
  }
}

//////////////////////////////////////////////////////////////////////

void SpectralFVElementData::createExtFaceCVConn()
{
  CFAUTOTRACE;

  // get the number of SV faces
  const CFuint nbrCellFaces = m_svFaceLocalExtFaceConn.size();

  // get the number of CVs lying at a SV face
  cf_assert(m_svFaceLocalExtFaceConn.size() > 0);
  const CFuint nbLocFacesPerCellFace = m_svFaceLocalExtFaceConn[0].size();

  // fill variable
  m_extSVFaceCVConn.resize(nbrCellFaces);
  for (CFuint iSVFace = 0; iSVFace < nbrCellFaces; ++iSVFace)
  {
// CF_DEBUG_OBJ(iSVFace);
    m_extSVFaceCVConn[iSVFace].resize(nbLocFacesPerCellFace);
    for (CFuint iExtFace = 0; iExtFace < nbLocFacesPerCellFace; ++iExtFace)
    {
      const CFuint extFaceID = m_svFaceLocalExtFaceConn[iSVFace][iExtFace];
      m_extSVFaceCVConn[iSVFace][iExtFace] = m_localExtFaceCVConn[extFaceID];
// CF_DEBUG_OBJ(m_extSVFaceCVConn[iSVFace][iExtFace]);
    }
  }
}

//////////////////////////////////////////////////////////////////////

void SpectralFVElementData::createExtFaceCVConnPerOrient()
{
  CFAUTOTRACE;

  // get the number of SV faces
  const CFuint nbrCellFaces = m_svFaceLocalExtFaceConn.size();

  // get the number of CVs lying at a SV face
  cf_assert(m_svFaceLocalExtFaceConn.size() > 0);
  const CFuint nbLocFacesPerCellFace = m_svFaceLocalExtFaceConn[0].size();

  // compute the number of CV connectivities through SV faces
  const CFuint nbrOrients = nbrCellFaces*(nbrCellFaces+1)/2;

  // resize m_extSVFaceCVConnPerOrient
  m_extSVFaceCVConnPerOrient.resize(nbrOrients);

  // resize m_eSVFaceLocESVFaceConnPerOrient
  m_eSVFaceLocESVFaceConnPerOrient.resize(nbrOrients);

  // fill m_extSVFaceCVConnPerOrient
  CFuint iOrient = 0;
  for (CFuint iSVFaceL = 0; iSVFaceL < nbrCellFaces; ++iSVFaceL)
  {
    for (CFuint iSVFaceR = iSVFaceL; iSVFaceR < nbrCellFaces; ++iSVFaceR, ++iOrient)
    {
      // resize m_eSVFaceLocESVFaceConnPerOrient[iOrient]
      m_eSVFaceLocESVFaceConnPerOrient[iOrient].resize(2);

      // set local SV face idxs
      m_eSVFaceLocESVFaceConnPerOrient[iOrient][LEFT ] = iSVFaceL ;
      m_eSVFaceLocESVFaceConnPerOrient[iOrient][RIGHT] = iSVFaceR;

      // resize m_extSVFaceCVConnPerOrient[iOrient]
      m_extSVFaceCVConnPerOrient[iOrient].resize(nbLocFacesPerCellFace);

      // create external face CV connectivity
      for (CFuint iExtFace = 0; iExtFace < nbLocFacesPerCellFace; ++iExtFace)
      {
        // resize the connectivity
        m_extSVFaceCVConnPerOrient[iOrient][iExtFace].resize(2);

        // get local external face IDs
        const CFuint localExtFace1ID = m_svFaceLocalExtFaceConn[iSVFaceL][iExtFace];
        const CFuint localExtFace2ID =
                               m_svFaceLocalExtFaceConn[iSVFaceR][nbLocFacesPerCellFace-1-iExtFace];

        // put neighbouring cv IDs in connectivity
        m_extSVFaceCVConnPerOrient[iOrient][iExtFace][LEFT ] = m_localExtFaceCVConn[localExtFace1ID];
        m_extSVFaceCVConnPerOrient[iOrient][iExtFace][RIGHT] = m_localExtFaceCVConn[localExtFace2ID];
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////

void SpectralFVElementData::computeInitPntsCoords()
{
  CFAUTOTRACE;

  // set nodal set
  setInterpolationNodeSet(m_polyOrder,m_initPntsCoords);
}

//////////////////////////////////////////////////////////////////////

void SpectralFVElementData::computeInitTransfMatrix()
{
  CFAUTOTRACE;

  // number of CVs
  const CFuint nbrCVs = getNbrOfCVs();

  // evaluate SV polynomials at initialization points
  vector< vector< CFreal > > lhs;
  lhs.resize(nbrCVs);
  for (CFuint iCV = 0; iCV < nbrCVs; ++iCV)
  {
    lhs[iCV].resize(nbrCVs);
  }
  lhs = getSVPolyValsAtNode(m_initPntsCoords);


  // fill in left hand side matrix
  RealMatrix lhsRM(nbrCVs,nbrCVs);
  for (CFuint iNode = 0; iNode < nbrCVs; ++iNode)
  {
    for (CFuint iCV = 0; iCV < nbrCVs; ++iCV)
    {
      lhsRM(iNode,iCV) = lhs[iNode][iCV];
    }
  }

  // compute values in transformation matrix
  m_initTransfMatrix.resize(nbrCVs,nbrCVs);
  InvertMatrix(lhsRM,m_initTransfMatrix);
}

//////////////////////////////////////////////////////////////////////

void SpectralFVElementData::createSVFaceNodeCoords()
{
  CFAUTOTRACE;

  // number of cell faces
  const CFuint nbrCellFaces = m_svFaceNodeConn.size();

  // create SV face-node coordinates
  m_svFaceNodeCoords.resize(nbrCellFaces);
  for (CFuint iFace = 0; iFace < nbrCellFaces; ++iFace)
  {
    // number of face nodes
    const CFuint nbrFaceNodes = m_svFaceNodeConn[iFace].size();

    m_svFaceNodeCoords[iFace].resize(nbrFaceNodes);
    for (CFuint iNode = 0; iNode < nbrFaceNodes; ++iNode)
    {
      m_svFaceNodeCoords[iFace][iNode].resize(m_dimensionality);

      // node ID
      const CFuint nodeID = m_svFaceNodeConn[iFace][iNode];

      // set coordinate
      m_svFaceNodeCoords[iFace][iNode] = m_svNodeCoords[nodeID];
    }
  }
}

//////////////////////////////////////////////////////////////////////

void SpectralFVElementData::createSVFaceNodeCoordsPerOrient()
{
  CFAUTOTRACE;

  // number of orientations
  const CFuint nbrOrients = m_svFaceNodeConnPerOrientation.size();

  // number of cell faces
  const CFuint nbrCellFaces = m_svFaceNodeCoords.size();

  // number of face nodes
  cf_assert(nbrCellFaces > 0);
  const CFuint nbrFaceNodes = m_svFaceNodeCoords[0].size();

  // create SV face-node coordinates per face connectivity orientation
  m_svFaceNodeCoordsPerOrient.resize(nbrOrients);
  for (CFuint iOrient = 0; iOrient < nbrOrients; ++iOrient)
  {
    // resize
    m_svFaceNodeCoordsPerOrient[iOrient].resize(2);
    for (CFuint iSide = 0; iSide < 2; ++iSide)
    {
      m_svFaceNodeCoordsPerOrient[iOrient][iSide].resize(nbrFaceNodes);
      for (CFuint iNode = 0; iNode < nbrFaceNodes; ++iNode)
      {
        m_svFaceNodeCoordsPerOrient[iOrient][iSide][iNode].resize(m_dimensionality);
      }
    }

    // fill variable
    for (CFuint iNode = 0; iNode < nbrFaceNodes; ++iNode)
    {
      // node IDs
      const CFuint nodeIDL = m_svFaceNodeConnPerOrientation[iOrient][LEFT ][iNode];
      const CFuint nodeIDR = m_svFaceNodeConnPerOrientation[iOrient][RIGHT][iNode];

      // set coordinates
      m_svFaceNodeCoordsPerOrient[iOrient][LEFT ][iNode] = m_svNodeCoords[nodeIDL];
      m_svFaceNodeCoordsPerOrient[iOrient][RIGHT][iNode] = m_svNodeCoords[nodeIDR];
    }
  }
}

//////////////////////////////////////////////////////////////////////

void SpectralFVElementData::createSVFaceNodeCoordsPerOrientNoSymm()
{
  CFAUTOTRACE;

  // number of orientations
  const CFuint nbrOrients = m_svFaceNodeConnPerOrientationNoSymm.size();

  // number of cell faces
  const CFuint nbrCellFaces = m_svFaceNodeCoords.size();

  // number of face nodes
  cf_assert(nbrCellFaces > 0);
  const CFuint nbrFaceNodes = m_svFaceNodeCoords[0].size();

  // create SV face-node coordinates per face connectivity orientation
  m_svFaceNodeCoordsPerOrientNoSymm.resize(nbrOrients);
  for (CFuint iOrient = 0; iOrient < nbrOrients; ++iOrient)
  {
    // resize
    m_svFaceNodeCoordsPerOrientNoSymm[iOrient].resize(2);
    for (CFuint iSide = 0; iSide < 2; ++iSide)
    {
      m_svFaceNodeCoordsPerOrientNoSymm[iOrient][iSide].resize(nbrFaceNodes);
      for (CFuint iNode = 0; iNode < nbrFaceNodes; ++iNode)
      {
        m_svFaceNodeCoordsPerOrientNoSymm[iOrient][iSide][iNode].resize(m_dimensionality);
      }
    }

    // fill variable
    for (CFuint iNode = 0; iNode < nbrFaceNodes; ++iNode)
    {
      // node IDs
      const CFuint nodeIDL = m_svFaceNodeConnPerOrientationNoSymm[iOrient][LEFT ][iNode];
      const CFuint nodeIDR = m_svFaceNodeConnPerOrientationNoSymm[iOrient][RIGHT][iNode];

      // set coordinates
      m_svFaceNodeCoordsPerOrientNoSymm[iOrient][LEFT ][iNode] = m_svNodeCoords[nodeIDL];
      m_svFaceNodeCoordsPerOrientNoSymm[iOrient][RIGHT][iNode] = m_svNodeCoords[nodeIDR];
    }
  }
}

//////////////////////////////////////////////////////////////////////

void SpectralFVElementData::computeIntFaceQuadPntsData()
{
  const CFuint dimUns = static_cast<CFuint>(m_dimensionality);

  CFAUTOTRACE;

  // set the dimensionality and order of the simplex integrator
  m_sIntegrator.setDimensionality(static_cast<CFDim>(m_dimensionality-1));
  m_sIntegrator.setIntegratorOrder(m_polyOrder);

  // number of internal and external faces
  const CFuint nbrIntFaces = m_localIntFaceCVConn.size();
  const CFuint nbrExtFaces = m_localFaceNodeConn.size() - nbrIntFaces;

  // create quadrature point data
  m_intFaceQuadWheights   .resize(nbrIntFaces);
  m_intFaceQuadPntCoords  .resize(nbrIntFaces);
  m_intFaceQuadPntPolyVals.resize(nbrIntFaces);
  for (CFuint iIntFace = 0; iIntFace < nbrIntFaces; ++iIntFace)
  {
    // face index
    const CFuint faceIdx = iIntFace+nbrExtFaces;

    // get face node coordinates
    const CFuint nbrFaceNodes = m_localFaceNodeConn[faceIdx].size();
    vector< RealVector > nodeCoord(nbrFaceNodes);
    for (CFuint iNode = 0; iNode < nbrFaceNodes; ++iNode)
    {
      const CFuint nodeID = m_localFaceNodeConn[faceIdx][iNode];
      nodeCoord[iNode].resize(m_dimensionality);
      nodeCoord[iNode] = m_localNodeCoord[nodeID];
    }

    // number of simplices in CV face
    const CFuint nbrSimplices = nbrFaceNodes - m_dimensionality + 1;

    // get face quadrature data
    for (CFuint iSimplex = 0; iSimplex < nbrSimplices; ++iSimplex)
    {
      // get node coordinates of simplex
      vector< RealVector > simplexNodeCoord(m_dimensionality);
      simplexNodeCoord[0].resize(m_dimensionality);
      simplexNodeCoord[0] = nodeCoord[0];
      for (CFuint iNode = 1; iNode < dimUns; ++iNode)
      {
        simplexNodeCoord[iNode].resize(m_dimensionality);
        simplexNodeCoord[iNode] = nodeCoord[iSimplex+iNode];
      }

      // get simplex quadrature nodes and wheights
      vector< RealVector > qNodeCoordSimplex = m_sIntegrator.getQuadPntsCoordsPlus1D  (simplexNodeCoord);
      vector< CFreal >     qWheightsSimplex  = m_sIntegrator.getQuadPntsWheightsPlus1D(simplexNodeCoord);

      // add simplex quadrature nodes and wheights to global list
      m_intFaceQuadPntCoords[iIntFace]
            .insert(m_intFaceQuadPntCoords[iIntFace].end(),qNodeCoordSimplex.begin(),qNodeCoordSimplex.end());
      m_intFaceQuadWheights [iIntFace]
            .insert(m_intFaceQuadWheights [iIntFace].end(),qWheightsSimplex .begin(),qWheightsSimplex .end());
    }

//     // The sum of these quadrature wheights should be one,
//     // because the face surface is included in the normal vector!!!
//     const CFreal invQWheightsSum = 1.0/accumulate(m_intFaceQuadWheights[iIntFace].begin(),
//                                                   m_intFaceQuadWheights[iIntFace].end()  ,0.0);
//     const CFuint nbrQPntsFace = m_intFaceQuadWheights[iIntFace].size();
//     for (CFuint iQPnt = 0; iQPnt < nbrQPntsFace; ++iQPnt)
//     {
//       m_intFaceQuadWheights[iIntFace][iQPnt] *= invQWheightsSum;
//     }

    const CFuint nbrQPntsFace = m_intFaceQuadWheights[iIntFace].size();
    // compute quadrature point polynomial values
    for (CFuint iQPnt = 0; iQPnt < nbrQPntsFace; ++iQPnt)
    {
      m_intFaceQuadPntPolyVals[iIntFace] = getSVPolyValsAtNode(m_intFaceQuadPntCoords[iIntFace]);
    }
  }
}

//////////////////////////////////////////////////////////////////////

void SpectralFVElementData::computeIntQuadPntsData()
{
  CFAUTOTRACE;

  // number of local internal faces
  const CFuint localIntFaces = getNbrOfLocalIntFaces();

  // number of local internal quadrature points
  CFuint nbrIntQPnts = 0;
  for (CFuint iFace = 0; iFace < localIntFaces; ++iFace)
  {
    nbrIntQPnts += m_intFaceQuadWheights[iFace].size();
  }

  // number of local CVs
  const CFuint nbrCVs = getNbrOfCVs();

  // resize variables
  m_intQuadWheights.resize(nbrIntQPnts);
  m_intQuadPntCoords.resize(nbrIntQPnts);
  m_intQuadPntPolyVals.resize(nbrIntQPnts);
  for (CFuint iIntQPnt = 0; iIntQPnt < nbrIntQPnts; ++iIntQPnt)
  {
    m_intQuadPntCoords  [iIntQPnt].resize(m_dimensionality);
    m_intQuadPntPolyVals[iIntQPnt].resize(nbrCVs);
  }

  // copy quadrature data to new variable
  CFuint iIntQPnt = 0;
  for (CFuint iFace = 0; iFace < localIntFaces; ++iFace)
  {
    // number of quadrature points in this face
    const CFuint nbrQPnts = m_intFaceQuadWheights[iFace].size();

    for (CFuint iQPnt = 0; iQPnt < nbrQPnts; ++iQPnt, ++iIntQPnt)
    {
      m_intQuadWheights   [iIntQPnt] = m_intFaceQuadWheights   [iFace][iQPnt];
      m_intQuadPntCoords  [iIntQPnt] = m_intFaceQuadPntCoords  [iFace][iQPnt];
      m_intQuadPntPolyVals[iIntQPnt] = m_intFaceQuadPntPolyVals[iFace][iQPnt];
    }
  }
}

//////////////////////////////////////////////////////////////////////

void SpectralFVElementData::computeExtFaceQuadPntsData()
{
  CFAUTOTRACE;

  const CFuint dimUns = static_cast<CFuint>(m_dimensionality);

  // get number of cell faces
  const CFuint nbrCellFaces = m_svFaceNodeCoords.size();

  // get the number of CVs lying at a SV face
  cf_assert(m_svFaceLocalExtFaceConn.size() > 0);
  const CFuint nbLocFacesPerCellFace = m_svFaceLocalExtFaceConn[0].size();

  // get the number of polynomials
  const CFuint nbrPolys = getNbrOfCVs();

  // dimensionality -1
  CFDim dimM1 = static_cast<CFDim>(m_dimensionality-1);

  // simplex integrator
  m_sIntegrator.setDimensionality(dimM1);
  m_sIntegrator.setIntegratorOrder(m_polyOrder);

  // compute wheight coordinates of quadrature points on a SV face (dimensionality of the face + 1 coords)
  // use the coordinates local to the SV face for this, (should be created before)
  // these wheight coordinates must be compatible with the order of the nodes in svFaceNodeCoords!!!
  m_extFaceQPntWheightCoords.resize(nbLocFacesPerCellFace);
  vector< vector< CFreal > >      qNodeQuadWheights (nbLocFacesPerCellFace);
  for (CFuint iExtFace = 0; iExtFace < nbLocFacesPerCellFace; ++iExtFace)
  {
    if (m_dimensionality > DIM_1D) /// @note is this really necessary? verify
    {
      // number of nodes to this face
      const CFuint nbrFaceNodes =  m_extFaceNodeLocalCoords[iExtFace].size();

      // number of simplices in this face
      const CFuint nbrSimplices = nbrFaceNodes - dimM1;

      // loop over simplices
      vector< RealVector >  qNodeCoord;
      vector< CFreal >      qWheights;
      for (CFuint iSimplex = 0; iSimplex < nbrSimplices; ++iSimplex)
      {
        // create vector of simplex node coordinates
        vector< RealVector > simpNodeCoord(m_dimensionality);
        simpNodeCoord[0].resize(m_dimensionality-1);
        simpNodeCoord[0] = m_extFaceNodeLocalCoords[iExtFace][0];
        for (CFuint iNode = 1; iNode < dimUns; ++iNode)
        {
          simpNodeCoord[iNode].resize(m_dimensionality-1);
          simpNodeCoord[iNode] = m_extFaceNodeLocalCoords[iExtFace][iNode+iSimplex];
        }

        // get simplex quadrature nodes
        vector< RealVector > qNodeCoordSimp = m_sIntegrator.getQuadPntsCoords  (simpNodeCoord);
        vector< CFreal >     qWheightsSimp  = m_sIntegrator.getQuadPntsWheights(simpNodeCoord);

        // add simplex quadrature nodes to global list
        qNodeCoord.insert(qNodeCoord.end(),qNodeCoordSimp.begin(),qNodeCoordSimp.end());
        qWheights .insert(qWheights .end(),qWheightsSimp .begin(),qWheightsSimp .end());
      }

      // total number of quadrature nodes in this face
      const CFuint nbrQNodesFace = qNodeCoord.size();

      // resize qNodeWheightCoords[iExtFace] and qNodeQuadWheights[iExtFace]
      m_extFaceQPntWheightCoords[iExtFace].resize(nbrQNodesFace);
      qNodeQuadWheights         [iExtFace].resize(nbrQNodesFace);

      // compute quadrature node wheight coordinates and store quadrature wheights
      for (CFuint iQNode = 0; iQNode < nbrQNodesFace; ++iQNode)
      {
        // resize
        m_extFaceQPntWheightCoords[iExtFace][iQNode].resize(m_dimensionality);

        // compute wheight coordinates
        m_extFaceQPntWheightCoords[iExtFace][iQNode][0] = 1.0 - qNodeCoord[iQNode].sum();
        // remaining wheight coordinates
        for (CFuint iCoor = 1; iCoor < dimUns; ++iCoor)
        {
          m_extFaceQPntWheightCoords[iExtFace][iQNode][iCoor] = qNodeCoord[iQNode][iCoor-1];
        }

        // store quadrature wheights
        qNodeQuadWheights[iExtFace][iQNode] = qWheights[iQNode];
      }
    }
    else
    {
      // resize
      m_extFaceQPntWheightCoords[iExtFace].resize(1);
      m_extFaceQPntWheightCoords[iExtFace][0].resize(1);

      // quadrature node wheight coordinate
      m_extFaceQPntWheightCoords[iExtFace][0][0] = 1.;

      // resize
      qNodeQuadWheights[iExtFace].resize(1);

      // quadrature node quadrature wheight
      qNodeQuadWheights[iExtFace][0] = 1.;
    }
  }

  // resize variables
  m_extFaceQuadWheights   .resize(nbrCellFaces);
  m_extFaceQuadPntCoords  .resize(nbrCellFaces);
  m_extFaceQuadPntPolyVals.resize(nbrCellFaces);
  for (CFuint iCellFace = 0; iCellFace < nbrCellFaces; ++iCellFace)
  {
    m_extFaceQuadWheights   [iCellFace].resize(nbLocFacesPerCellFace);
    m_extFaceQuadPntCoords  [iCellFace].resize(nbLocFacesPerCellFace);
    m_extFaceQuadPntPolyVals[iCellFace].resize(nbLocFacesPerCellFace);
    for (CFuint iExtFace = 0; iExtFace < nbLocFacesPerCellFace; ++iExtFace)
    {
      // number of quadrature points on this face
      const CFuint nbrFaceQPnts = m_extFaceQPntWheightCoords[iExtFace].size();

      m_extFaceQuadWheights   [iCellFace][iExtFace].resize(nbrFaceQPnts);
      m_extFaceQuadPntCoords  [iCellFace][iExtFace].resize(nbrFaceQPnts);
      m_extFaceQuadPntPolyVals[iCellFace][iExtFace].resize(nbrFaceQPnts);

      for (CFuint iQNode = 0; iQNode < nbrFaceQPnts; ++iQNode)
      {
        m_extFaceQuadPntCoords  [iCellFace][iExtFace][iQNode].resize(m_dimensionality);
        m_extFaceQuadPntPolyVals[iCellFace][iExtFace][iQNode].resize(nbrPolys);
      }
    }
  }

  // create quadrature data for SV faces
  for (CFuint iCellFace = 0; iCellFace < nbrCellFaces; ++iCellFace)
  {
    // vector for face node coordinates
    vector< RealVector > faceNodeCoor(m_dimensionality);
    for (CFuint iNode = 0; iNode < dimUns; ++iNode)
    {
      faceNodeCoor[iNode].resize(m_dimensionality);
      faceNodeCoor[iNode] = m_svFaceNodeCoords[iCellFace][iNode];
    }

    // compute quadrature data for this face
    for (CFuint iExtFace = 0; iExtFace < nbLocFacesPerCellFace; ++iExtFace)
    {
      // number of quadrature nodes for this subface
      const CFuint nbrFaceQNodes = qNodeQuadWheights[iExtFace].size();

      // set quadrature wheights
      m_extFaceQuadWheights[iCellFace][iExtFace] = qNodeQuadWheights[iExtFace];
      const CFreal scaleFactor =
        m_faceFracCV[iCellFace][iExtFace]/accumulate(m_extFaceQuadWheights[iCellFace][iExtFace].begin(),
                                          m_extFaceQuadWheights[iCellFace][iExtFace].end(),0.0);

      // scale the quadrature wheights
      for (CFuint iQNode = 0; iQNode < nbrFaceQNodes; ++iQNode)
      {
        m_extFaceQuadWheights[iCellFace][iExtFace][iQNode] *= scaleFactor;
      }

      // quadrature node coordinates
      for (CFuint iQNode = 0; iQNode < nbrFaceQNodes; ++iQNode)
      {
        m_extFaceQuadPntCoords[iCellFace][iExtFace][iQNode] = 0.0;
        for (CFuint iFaceNode = 0; iFaceNode < dimUns; ++iFaceNode)
        {
          m_extFaceQuadPntCoords[iCellFace][iExtFace][iQNode] +=
            m_extFaceQPntWheightCoords[iExtFace][iQNode][iFaceNode]*faceNodeCoor[iFaceNode];
        }
      }

      // quadrature node polynomial values
      m_extFaceQuadPntPolyVals[iCellFace][iExtFace] =
        getSVPolyValsAtNode(m_extFaceQuadPntCoords[iCellFace][iExtFace]);
    }
  }

/*  for (CFuint iCellFace = 0; iCellFace < nbrCellFaces; ++iCellFace)
  {
    CF_DEBUG_OBJ(iCellFace);
    for (CFuint iExtFace = 0; iExtFace < nbLocFacesPerCellFace; ++iExtFace)
    {
      // number of quadrature nodes for this subface
      const CFuint nbrFaceQNodes = m_extFaceQuadWheights[iCellFace][iExtFace].size();

      CF_DEBUG_OBJ(iExtFace);
      for (CFuint iQNode = 0; iQNode < nbrFaceQNodes; ++iQNode)
      {
        CF_DEBUG_OBJ(iQNode);
        CF_DEBUG_OBJ(m_extFaceQuadWheights[iCellFace][iExtFace][iQNode]);
        CF_DEBUG_OBJ(m_extFaceQuadPntCoords[iCellFace][iExtFace][iQNode]);
        for (CFuint iPoly = 0; iPoly < nbrPolys; ++iPoly)
        {
          CF_DEBUG_OBJ(m_extFaceQuadPntPolyVals[iCellFace][iExtFace][iQNode][iPoly]);
        }
      }
    }
  }*/
}

//////////////////////////////////////////////////////////////////////

void SpectralFVElementData::computeExtQuadPntsData()
{
  CFAUTOTRACE;

  // get number of cell faces
  const CFuint nbrCellFaces = m_svFaceNodeCoords.size();

  // get the number of CVs lying at a SV face
  cf_assert(m_svFaceLocalExtFaceConn.size() > 0);
  const CFuint nbLocFacesPerCellFace = m_svFaceLocalExtFaceConn[0].size();

  // get the number of polynomials
  const CFuint nbrCVs = getNbrOfCVs();

  // number of external quadrature points
  CFuint nbrExtQPnts = 0;
  for (CFuint iFace = 0; iFace < nbLocFacesPerCellFace; ++iFace)
  {
    nbrExtQPnts += m_extFaceQuadWheights[0][iFace].size();
  }

  // resize variables
  m_extQPntWheightCoords.resize(nbrExtQPnts);
  for (CFuint iQPnt = 0; iQPnt < nbrExtQPnts; ++iQPnt)
  {
    m_extQPntWheightCoords[iQPnt].resize(m_dimensionality);
  }

  m_extQWheights.resize(nbrCellFaces);
  m_extQPntCoords.resize(nbrCellFaces);
  m_extQPntPolyVals.resize(nbrCellFaces);
  for (CFuint iFace = 0; iFace < nbrCellFaces; ++iFace)
  {
    m_extQWheights   [iFace].resize(nbrExtQPnts);
    m_extQPntCoords  [iFace].resize(nbrExtQPnts);
    m_extQPntPolyVals[iFace].resize(nbrExtQPnts);

    for (CFuint iQPnt = 0; iQPnt < nbrExtQPnts; ++iQPnt)
    {
      m_extQPntCoords  [iFace][iQPnt].resize(m_dimensionality);
      m_extQPntPolyVals[iFace][iQPnt].resize(nbrCVs);
    }
  }

  // copy quadrature data
  CFuint iQPnt = 0;
  for (CFuint iExtFace = 0; iExtFace < nbLocFacesPerCellFace; ++iExtFace)
  {
    // number of quadrature nodes for this subface
    const CFuint nbrFaceQNodes = m_extFaceQPntWheightCoords[iExtFace].size();
    for (CFuint iQNode = 0; iQNode < nbrFaceQNodes; ++iQNode, ++iQPnt)
    {
      m_extQPntWheightCoords[iQPnt] = m_extFaceQPntWheightCoords[iExtFace][iQNode];
    }
  }

  for (CFuint iFace = 0; iFace < nbrCellFaces; ++iFace)
  {
    CFuint iQPnt = 0;
    for (CFuint iExtFace = 0; iExtFace < nbLocFacesPerCellFace; ++iExtFace)
    {
      // number of quadrature nodes for this subface
      const CFuint nbrFaceQNodes = m_extFaceQuadWheights[iFace][iExtFace].size();
      for (CFuint iQNode = 0; iQNode < nbrFaceQNodes; ++iQNode, ++iQPnt)
      {
        m_extQWheights   [iFace][iQPnt] = m_extFaceQuadWheights   [iFace][iExtFace][iQNode];
        m_extQPntCoords  [iFace][iQPnt] = m_extFaceQuadPntCoords  [iFace][iExtFace][iQNode];
        m_extQPntPolyVals[iFace][iQPnt] = m_extFaceQuadPntPolyVals[iFace][iExtFace][iQNode];
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////

void SpectralFVElementData::computeExtFaceQuadPntsDataPerOrientation()
{
  CFAUTOTRACE;

  const CFuint dimUns = static_cast<CFuint>(m_dimensionality);

  // get number of orientations
  const CFuint nbrOrients = m_svFaceNodeCoordsPerOrient.size();

  // get the number of CVs lying at a SV face
  cf_assert(m_svFaceLocalExtFaceConn.size() > 0);
  const CFuint nbLocFacesPerCellFace = m_svFaceLocalExtFaceConn[0].size();

  // get the number of polynomials
  const CFuint nbrPolys = getNbrOfCVs();

  // dimensionality -1
  CFDim dimM1 = static_cast<CFDim>(m_dimensionality-1);

  // simplex integrator
  m_sIntegrator.setDimensionality(dimM1);
  m_sIntegrator.setIntegratorOrder(m_polyOrder);

  // resize variables
  m_extFaceQWheightsPerOrient.resize(nbrOrients);
  m_extFaceQPntCoordsPerOrient.resize(nbrOrients);
  m_extFaceQPntPolyValsPerOrient.resize(nbrOrients);
  for (CFuint iOrient = 0; iOrient < nbrOrients; ++iOrient)
  {
    m_extFaceQWheightsPerOrient[iOrient].resize(nbLocFacesPerCellFace);
    for (CFuint iExtFace = 0; iExtFace < nbLocFacesPerCellFace; ++iExtFace)
    {
      // number of quadrature points on this face
      cf_assert(m_extFaceQuadWheights.size() > 0);
      const CFuint nbrFaceQPnts = m_extFaceQuadWheights[0][iExtFace].size();

      m_extFaceQWheightsPerOrient[iOrient][iExtFace].resize(nbrFaceQPnts);
    }

    m_extFaceQPntCoordsPerOrient[iOrient].resize(2);
    m_extFaceQPntPolyValsPerOrient[iOrient].resize(2);
    for (CFuint iSide = 0; iSide < 2; ++iSide)
    {
      m_extFaceQPntCoordsPerOrient[iOrient][iSide].resize(nbLocFacesPerCellFace);
      m_extFaceQPntPolyValsPerOrient[iOrient][iSide].resize(nbLocFacesPerCellFace);

      for (CFuint iExtFace = 0; iExtFace < nbLocFacesPerCellFace; ++iExtFace)
      {
        // number of quadrature points on this face
        cf_assert(m_extFaceQuadWheights.size() > 0);
        const CFuint nbrFaceQPnts = m_extFaceQuadWheights[0][iExtFace].size();

        m_extFaceQPntCoordsPerOrient[iOrient][iSide][iExtFace].resize(nbrFaceQPnts);
        m_extFaceQPntPolyValsPerOrient[iOrient][iSide][iExtFace].resize(nbrFaceQPnts);
        for (CFuint iQNode = 0; iQNode < nbrFaceQPnts; ++iQNode)
        {
          m_extFaceQPntCoordsPerOrient[iOrient][iSide][iExtFace][iQNode].resize(m_dimensionality);
          m_extFaceQPntPolyValsPerOrient[iOrient][iSide][iExtFace][iQNode].resize(nbrPolys);
        }
      }
    }
  }

  // create quadrature data
  for (CFuint iOrient = 0; iOrient < nbrOrients; ++iOrient)
  {
    for (CFuint iSide = 0; iSide < 2; ++iSide)
    {
      // vector for face node coordinates
      vector< RealVector > faceNodeCoor(m_dimensionality);
      for (CFuint iNode = 0; iNode < dimUns; ++iNode)
      {
        faceNodeCoor[iNode].resize(m_dimensionality);
        faceNodeCoor[iNode] = m_svFaceNodeCoordsPerOrient[iOrient][iSide][iNode];
      }

      // compute quadrature data for this face
      for (CFuint iExtFace = 0; iExtFace < nbLocFacesPerCellFace; ++iExtFace)
      {
        // number of quadrature nodes for this subface
        const CFuint nbrFaceQNodes = m_extFaceQuadWheights[0][iExtFace].size();

        // set quadrature wheights
        if (iSide == 0)
        {
          m_extFaceQWheightsPerOrient[iOrient][iExtFace] = m_extFaceQuadWheights[0][iExtFace];
        }

        // quadrature node coordinates
        for (CFuint iQNode = 0; iQNode < nbrFaceQNodes; ++iQNode)
        {
          m_extFaceQPntCoordsPerOrient[iOrient][iSide][iExtFace][iQNode] = 0.0;
          for (CFuint iFaceNode = 0; iFaceNode < dimUns; ++iFaceNode)
          {
            m_extFaceQPntCoordsPerOrient[iOrient][iSide][iExtFace][iQNode] +=
              m_extFaceQPntWheightCoords[iExtFace][iQNode][iFaceNode]*faceNodeCoor[iFaceNode];
          }
        }

        // quadrature node polynomial values
        m_extFaceQPntPolyValsPerOrient[iOrient][iSide][iExtFace] =
          getSVPolyValsAtNode(m_extFaceQPntCoordsPerOrient[iOrient][iSide][iExtFace]);
      }
    }
  }

/*  for (CFuint iOrient = 0; iOrient < nbrOrients; ++iOrient)
  {
    CF_DEBUG_OBJ(iOrient);
    for (CFuint iExtFace = 0; iExtFace < nbLocFacesPerCellFace; ++iExtFace)
    {
      // number of quadrature nodes for this subface
      const CFuint nbrFaceQNodes = m_extFaceQWheightsPerOrient[iOrient][iExtFace].size();

      CF_DEBUG_OBJ(iExtFace);
      for (CFuint iQNode = 0; iQNode < nbrFaceQNodes; ++iQNode)
      {
        CF_DEBUG_OBJ(iQNode);
        CF_DEBUG_OBJ(m_extFaceQWheightsPerOrient[iOrient][iExtFace][iQNode]);
        CF_DEBUG_OBJ(m_extFaceQPntCoordsPerOrient[iOrient][LEFT ][iExtFace][iQNode]);
        CF_DEBUG_OBJ(m_extFaceQPntCoordsPerOrient[iOrient][RIGHT][iExtFace][iQNode]);
        for (CFuint iPoly = 0; iPoly < nbrPolys; ++iPoly)
        {
          CF_DEBUG_OBJ(m_extFaceQPntPolyValsPerOrient[iOrient][LEFT ][iExtFace][iQNode][iPoly]);
        }
        for (CFuint iPoly = 0; iPoly < nbrPolys; ++iPoly)
        {
          CF_DEBUG_OBJ(m_extFaceQPntPolyValsPerOrient[iOrient][RIGHT][iExtFace][iQNode][iPoly]);
        }
      }
    }
  }*/
}

//////////////////////////////////////////////////////////////////////

void SpectralFVElementData::computeExtFaceQuadPntsDataPerOrientationNoSymm()
{
  CFAUTOTRACE;

  const CFuint dimUns = static_cast<CFuint>(m_dimensionality);

  // get number of orientations
  const CFuint nbrOrients = m_svFaceNodeCoordsPerOrientNoSymm.size();

  // get the number of CVs lying at a SV face
  cf_assert(m_svFaceLocalExtFaceConn.size() > 0);
  const CFuint nbLocFacesPerCellFace = m_svFaceLocalExtFaceConn[0].size();

  // get the number of polynomials
  const CFuint nbrPolys = getNbrOfCVs();

  // dimensionality -1
  CFDim dimM1 = static_cast<CFDim>(m_dimensionality-1);

  // simplex integrator
  m_sIntegrator.setDimensionality(dimM1);
  m_sIntegrator.setIntegratorOrder(m_polyOrder);

  // resize variables
  m_extFaceQPntPolyValsPerOrientNoSymm.resize(nbrOrients);
  for (CFuint iOrient = 0; iOrient < nbrOrients; ++iOrient)
  {
    m_extFaceQPntPolyValsPerOrientNoSymm[iOrient].resize(2);
    for (CFuint iSide = 0; iSide < 2; ++iSide)
    {
      m_extFaceQPntPolyValsPerOrientNoSymm[iOrient][iSide].resize(nbLocFacesPerCellFace);
      for (CFuint iExtFace = 0; iExtFace < nbLocFacesPerCellFace; ++iExtFace)
      {
        // number of quadrature points on this face
        cf_assert(m_extFaceQuadWheights.size() > 0);
        const CFuint nbrFaceQPnts = m_extFaceQuadWheights[0][iExtFace].size();
        m_extFaceQPntPolyValsPerOrientNoSymm[iOrient][iSide][iExtFace].resize(nbrFaceQPnts);
        for (CFuint iQNode = 0; iQNode < nbrFaceQPnts; ++iQNode)
        {
          m_extFaceQPntPolyValsPerOrientNoSymm[iOrient][iSide][iExtFace][iQNode].resize(nbrPolys);
        }
      }
    }
  }

  // create quadrature data
  for (CFuint iOrient = 0; iOrient < nbrOrients; ++iOrient)
  {
    for (CFuint iSide = 0; iSide < 2; ++iSide)
    {
      // vector for face node coordinates
      vector< RealVector > faceNodeCoor(m_dimensionality);
      for (CFuint iNode = 0; iNode < dimUns; ++iNode)
      {
        faceNodeCoor[iNode].resize(m_dimensionality);
        faceNodeCoor[iNode] = m_svFaceNodeCoordsPerOrientNoSymm[iOrient][iSide][iNode];
      }

      // compute quadrature data for this face
      for (CFuint iExtFace = 0; iExtFace < nbLocFacesPerCellFace; ++iExtFace)
      {
        // number of quadrature nodes for this subface
        const CFuint nbrFaceQNodes = m_extFaceQuadWheights[0][iExtFace].size();

        // quadrature node coordinates
        vector< RealVector > qNodeCoords(nbrFaceQNodes);
        for (CFuint iQNode = 0; iQNode < nbrFaceQNodes; ++iQNode)
        {
          qNodeCoords[iQNode].resize(m_dimensionality);
          qNodeCoords[iQNode] = 0.0;
          for (CFuint iFaceNode = 0; iFaceNode < dimUns; ++iFaceNode)
          {
            qNodeCoords[iQNode] +=
                m_extFaceQPntWheightCoords[iExtFace][iQNode][iFaceNode]*faceNodeCoor[iFaceNode];
          }
        }

        // quadrature node polynomial values
        m_extFaceQPntPolyValsPerOrientNoSymm[iOrient][iSide][iExtFace] =
            getSVPolyValsAtNode(qNodeCoords);
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////

void SpectralFVElementData::computeExtQuadPntsDataPerOrientation()
{
  CFAUTOTRACE;

  // get number of orientations
  const CFuint nbrOrients = m_svFaceNodeCoordsPerOrient.size();

  // get the number of CVs lying at a SV face
  cf_assert(m_svFaceLocalExtFaceConn.size() > 0);
  const CFuint nbLocFacesPerCellFace = m_svFaceLocalExtFaceConn[0].size();

  // get the number of polynomials
  const CFuint nbrCVs = getNbrOfCVs();

  // number of external quadrature points
  CFuint nbrExtQPnts = 0;
  for (CFuint iFace = 0; iFace < nbLocFacesPerCellFace; ++iFace)
  {
    nbrExtQPnts += m_extFaceQWheightsPerOrient[0][iFace].size();
  }

  // resize variables
  m_extQWheightsPerOrient.resize(nbrOrients);
  m_extQPntCoordsPerOrient.resize(nbrOrients);
  m_extQPntPolyValsPerOrient.resize(nbrOrients);
  for (CFuint iOrient = 0; iOrient < nbrOrients; ++iOrient)
  {
    m_extQWheightsPerOrient[iOrient].resize(nbrExtQPnts);

    m_extQPntCoordsPerOrient  [iOrient].resize(2);
    m_extQPntPolyValsPerOrient[iOrient].resize(2);
    for (CFuint iSide = 0; iSide < 2; ++iSide)
    {
      m_extQPntCoordsPerOrient  [iOrient][iSide].resize(nbrExtQPnts);
      m_extQPntPolyValsPerOrient[iOrient][iSide].resize(nbrExtQPnts);

      for (CFuint iQPnt = 0; iQPnt < nbrExtQPnts; ++iQPnt)
      {
        m_extQPntCoordsPerOrient  [iOrient][iSide][iQPnt].resize(m_dimensionality);
        m_extQPntPolyValsPerOrient[iOrient][iSide][iQPnt].resize(nbrCVs);
      }
    }
  }

  // copy quadrature data
  for (CFuint iOrient = 0; iOrient < nbrOrients; ++iOrient)
  {
    CFuint iQPnt = 0;
    for (CFuint iExtFace = 0; iExtFace < nbLocFacesPerCellFace; ++iExtFace)
    {
      // number of quadrature nodes for this subface
      const CFuint nbrFaceQNodes = m_extFaceQWheightsPerOrient[iOrient][iExtFace].size();
      for (CFuint iQNode = 0; iQNode < nbrFaceQNodes; ++iQNode, ++iQPnt)
      {
        m_extQWheightsPerOrient[iOrient][iQPnt] = m_extFaceQWheightsPerOrient[iOrient][iExtFace][iQNode];
        for (CFuint iSide = 0; iSide < 2; ++iSide)
        {
          m_extQPntCoordsPerOrient  [iOrient][iSide][iQPnt] =
                                            m_extFaceQPntCoordsPerOrient[iOrient][iSide][iExtFace][iQNode];
          m_extQPntPolyValsPerOrient[iOrient][iSide][iQPnt] =
                                            m_extFaceQPntPolyValsPerOrient[iOrient][iSide][iExtFace][iQNode];
        }
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////

void SpectralFVElementData::computeExtQuadPntsDataPerOrientationNoSymm()
{
  CFAUTOTRACE;

  // get number of orientations
  const CFuint nbrOrients = m_svFaceNodeCoordsPerOrientNoSymm.size();

  // get the number of CVs lying at a SV face
  cf_assert(m_svFaceLocalExtFaceConn.size() > 0);
  const CFuint nbLocFacesPerCellFace = m_svFaceLocalExtFaceConn[0].size();

  // get the number of polynomials
  const CFuint nbrCVs = getNbrOfCVs();

  // number of external quadrature points
  CFuint nbrExtQPnts = 0;
  for (CFuint iFace = 0; iFace < nbLocFacesPerCellFace; ++iFace)
  {
    nbrExtQPnts += m_extFaceQWheightsPerOrient[0][iFace].size();
  }

  // resize variables
  m_extQPntPolyValsPerOrientNoSymm.resize(nbrOrients);
  for (CFuint iOrient = 0; iOrient < nbrOrients; ++iOrient)
  {
    m_extQPntPolyValsPerOrientNoSymm[iOrient].resize(2);
    for (CFuint iSide = 0; iSide < 2; ++iSide)
    {
      m_extQPntPolyValsPerOrientNoSymm[iOrient][iSide].resize(nbrExtQPnts);
      for (CFuint iQPnt = 0; iQPnt < nbrExtQPnts; ++iQPnt)
      {
        m_extQPntPolyValsPerOrientNoSymm[iOrient][iSide][iQPnt].resize(nbrCVs);
      }
    }
  }

  // copy quadrature data
  for (CFuint iOrient = 0; iOrient < nbrOrients; ++iOrient)
  {
    CFuint iQPnt = 0;
    for (CFuint iExtFace = 0; iExtFace < nbLocFacesPerCellFace; ++iExtFace)
    {
      // number of quadrature nodes for this subface
      const CFuint nbrFaceQNodes =
          m_extFaceQPntPolyValsPerOrientNoSymm[iOrient][0][iExtFace].size();
      cf_assert(nbrFaceQNodes == m_extFaceQPntPolyValsPerOrientNoSymm[iOrient][1][iExtFace].size());
      for (CFuint iQNode = 0; iQNode < nbrFaceQNodes; ++iQNode, ++iQPnt)
      {
        for (CFuint iSide = 0; iSide < 2; ++iSide)
        {
          m_extQPntPolyValsPerOrientNoSymm[iOrient][iSide][iQPnt] =
              m_extFaceQPntPolyValsPerOrientNoSymm[iOrient][iSide][iExtFace][iQNode];
        }
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////

void SpectralFVElementData::computeFaceAvgPolyVals()
{
  CFAUTOTRACE;

  const CFuint dimUns = static_cast<CFuint>(m_dimensionality);

  // get number of cell faces
  const CFuint nbrCellFaces = m_svFaceNodeCoords.size();

  // get the number of polynomials
  const CFuint nbrPolys = getNbrOfCVs();

  // resize m_svFaceAvgPolyVals
  m_svFaceAvgPolyVals.resize(nbrCellFaces);
  for (CFuint iFace = 0; iFace < nbrCellFaces; ++iFace)
  {
    m_svFaceAvgPolyVals[iFace].resize(nbrPolys);
  }

  // set the dimensionality and order of the simplex integrator
  const CFDim dimM1 = static_cast<CFDim>(m_dimensionality-1);
  m_sIntegrator.setDimensionality(dimM1);
  m_sIntegrator.setIntegratorOrder(m_polyOrder);

  // loop over the faces
  for (CFuint iFace = 0; iFace < nbrCellFaces; ++iFace)
  {
    // quadrature node coordinates and wheights
    vector< RealVector > qNodeCoord = m_sIntegrator.getQuadPntsCoordsPlus1D  (m_svFaceNodeCoords[iFace]);
    vector< CFreal >     qWheights  = m_sIntegrator.getQuadPntsWheightsPlus1D(m_svFaceNodeCoords[iFace]);

    // number of quadrature nodes
    const CFuint nbrQNodes = qWheights.size();

    // surface of face (in local coordinates)
    const CFreal faceSurf = accumulate(qWheights.begin(),qWheights.end(),0.0);

    // loop over polynomials
    for (CFuint iPoly = 0; iPoly < nbrPolys; ++iPoly)
    {
      m_svFaceAvgPolyVals[iFace][iPoly] = 0.0;
      for (CFuint iQNode = 0; iQNode < nbrQNodes; ++iQNode)
      {
        for (CFuint iTerm = 0; iTerm < nbrPolys; ++iTerm)
        {
          CFreal term = m_polyCoefSFV[iPoly][iTerm];
          for (CFuint iCoor = 0; iCoor < dimUns; ++iCoor)
          {
            term *= pow(qNodeCoord[iQNode][iCoor],m_polyExponents[iTerm][iCoor]);
          }
          m_svFaceAvgPolyVals[iFace][iPoly] += qWheights[iQNode]*term;
        }
      }

      // divide by face surface
      m_svFaceAvgPolyVals[iFace][iPoly] /= faceSurf;
    }
  }

/*  for (CFuint iFace = 0; iFace < nbrCellFaces; ++iFace)
  {
    CF_DEBUG_OBJ(iFace);
    for (CFuint iPoly = 0; iPoly < nbrPolys; ++iPoly)
    {
      CF_DEBUG_OBJ(m_svFaceAvgPolyVals[iFace][iPoly]);
    }
  }*/
}

//////////////////////////////////////////////////////////////////////

void SpectralFVElementData::computeFaceAvgPolyValsPerOrient()
{
  const CFuint dimUns = static_cast<CFuint>(m_dimensionality);

  // get number of orientations
  const CFuint nbrOrients = m_svFaceNodeCoordsPerOrient.size();

  // get the number of polynomials
  const CFuint nbrPolys = getNbrOfCVs();

  // resize m_svFaceAvgPolyValsPerOrient
  m_svFaceAvgPolyValsPerOrient.resize(nbrOrients);
  for (CFuint iOrient = 0; iOrient < nbrOrients; ++iOrient)
  {
    m_svFaceAvgPolyValsPerOrient[iOrient].resize(2);
    for (CFuint iSide = 0; iSide < 2; ++iSide)
    {
      m_svFaceAvgPolyValsPerOrient[iOrient][iSide].resize(nbrPolys);
    }
  }

  // set the dimensionality and order of the simplex integrator
  const CFDim dimM1 = static_cast<CFDim>(m_dimensionality-1);
  m_sIntegrator.setDimensionality(dimM1);
  m_sIntegrator.setIntegratorOrder(m_polyOrder);

  // loop over the orientations
  for (CFuint iOrient = 0; iOrient < nbrOrients; ++iOrient)
  {
    // number of face nodes
    const CFuint nbrFaceNodes = m_svFaceNodeCoordsPerOrient[iOrient][LEFT].size();
    cf_assert(nbrFaceNodes == m_svFaceNodeCoordsPerOrient[iOrient][RIGHT].size());

    for (CFuint iSide = 0; iSide < 2; ++iSide)
    {
      // quadrature node coordinates and wheights
      vector< RealVector > qNodeCoord =
                          m_sIntegrator.getQuadPntsCoordsPlus1D  (m_svFaceNodeCoordsPerOrient[iOrient][iSide]);
      vector< CFreal >     qWheights  =
                          m_sIntegrator.getQuadPntsWheightsPlus1D(m_svFaceNodeCoordsPerOrient[iOrient][iSide]);

      // number of quadrature nodes
      const CFuint nbrQNodes = qWheights.size();

      // surface of face (in local coordinates)
      const CFreal faceSurf = accumulate(qWheights.begin(),qWheights.end(),0.0);

      // loop over polynomials
      for (CFuint iPoly = 0; iPoly < nbrPolys; ++iPoly)
      {
        m_svFaceAvgPolyValsPerOrient[iOrient][iSide][iPoly] = 0.0;
        for (CFuint iQNode = 0; iQNode < nbrQNodes; ++iQNode)
        {
          for (CFuint iTerm = 0; iTerm < nbrPolys; ++iTerm)
          {
            CFreal term = m_polyCoefSFV[iPoly][iTerm];
            for (CFuint iCoor = 0; iCoor < dimUns; ++iCoor)
            {
              term *= pow(qNodeCoord[iQNode][iCoor],m_polyExponents[iTerm][iCoor]);
            }
            m_svFaceAvgPolyValsPerOrient[iOrient][iSide][iPoly] += qWheights[iQNode]*term;
          }
        }

        // divide by face surface
        m_svFaceAvgPolyValsPerOrient[iOrient][iSide][iPoly] /= faceSurf;
      }
    }
  }

/*  for (CFuint iOrient = 0; iOrient < nbrOrients; ++iOrient)
  {
    CF_DEBUG_OBJ(iOrient);
    for (CFuint iSide = 0; iSide < 2; ++iSide)
    {
      CF_DEBUG_OBJ(iSide);
      for (CFuint iPoly = 0; iPoly < nbrPolys; ++iPoly)
      {
        CF_DEBUG_OBJ(m_svFaceAvgPolyValsPerOrient[iOrient][iSide][iPoly]);
      }
    }
  }*/
}

//////////////////////////////////////////////////////////////////////

void SpectralFVElementData::createFluxPolyNodeCoord()
{
  CFAUTOTRACE;

  // set the node coordinates
  setInterpolationNodeSet(static_cast<CFPolyOrder::Type>(m_polyOrder+1),m_fluxPolyNodeCoord);
}

//////////////////////////////////////////////////////////////////////

void SpectralFVElementData::createFluxPolyCoef()
{
  CFAUTOTRACE;

  const CFuint dimUns = static_cast<CFuint>(m_dimensionality);

  // number of flux polynomials
  const CFuint nbrFlxPolys = getNbrOfFlxPnts();

  // resize m_fluxPolyCoefs
  m_fluxPolyCoefs.resize(nbrFlxPolys);
  for (CFuint iPoly = 0; iPoly < nbrFlxPolys; ++iPoly)
  {
    m_fluxPolyCoefs[iPoly].resize(nbrFlxPolys);
  }

  // variable for LHS of linear system
  RealMatrix lhs(nbrFlxPolys,nbrFlxPolys);
  RealMatrix lhsInv(nbrFlxPolys,nbrFlxPolys);

  // create LHS matrix
  for (CFuint iPoly = 0; iPoly < nbrFlxPolys; ++iPoly)
  {
    for (CFuint iTerm = 0; iTerm < nbrFlxPolys; ++iTerm)
    {
      lhs(iPoly,iTerm) = 1.0;
      for (CFuint iCoor = 0; iCoor < dimUns; ++ iCoor)
      {
        lhs(iPoly,iTerm) *= pow(m_fluxPolyNodeCoord[iPoly][iCoor],m_fluxPolyExponents[iTerm][iCoor]);
      }
    }
  }

  // invert the LHS matrix
  InvertMatrix(lhs,lhsInv);

  // store flux polynomial coefficients
  for (CFuint iPoly = 0; iPoly < nbrFlxPolys; ++iPoly)
  {
    for (CFuint iTerm = 0; iTerm < nbrFlxPolys; ++iTerm)
    {
      m_fluxPolyCoefs[iPoly][iTerm] = lhsInv(iTerm,iPoly);
    }
  }
}

//////////////////////////////////////////////////////////////////////

void SpectralFVElementData::createSolInFlxPntsCoefs()
{
  CFAUTOTRACE;

  // compute solution reconstruction coefficients
  m_solInFluxPntCoef = getSVPolyValsAtNode(m_fluxPolyNodeCoord);

/*  for (CFuint iFlx = 0; iFlx < m_solInFluxPntCoef.size(); ++iFlx)
  {
    CF_DEBUG_OBJ(iFlx);
    CFreal sum = 0.0;
    for (CFuint iCV = 0; iCV < m_solInFluxPntCoef[iFlx].size(); ++iCV)
    {
      // CF_DEBUG_OBJ(m_solInFluxPntCoef[iFlx][iCV]);
      sum += m_solInFluxPntCoef[iFlx][iCV];
    }
    CF_DEBUG_OBJ(sum);
  }*/
}

//////////////////////////////////////////////////////////////////////

void SpectralFVElementData::createVolumeTermsTensor()
// this should also work in 3D, if the local face node connectivity is ordered in the correct way
// and if the faces are subdivided into triangles
{
  CFAUTOTRACE;

  const CFuint dimUns = static_cast<CFuint>(m_dimensionality);

  // number of local faces
  const CFuint nbrLocalFaces = getNbrOfLocalFaces();

  // number of local internal faces
  const CFuint nbrLocalIntFaces = getNbrOfLocalIntFaces();

  // number of local external faces
  const CFuint nbrLocalExtFaces = nbrLocalFaces - nbrLocalIntFaces;

  // number of CVs
  const CFuint nbrCVs = getNbrOfCVs();

  // number of flux points
  const CFuint nbrFlxPnts = getNbrOfFlxPnts();

  // resize and initialize  m_volTermTensor
  m_volTermTensor.resize(nbrCVs);
  for (CFuint iCV = 0; iCV < nbrCVs; ++iCV)
  {
    m_volTermTensor[iCV].resize(nbrFlxPnts);
    for (CFuint iFlx = 0; iFlx < nbrFlxPnts; ++iFlx)
    {
      m_volTermTensor[iCV][iFlx].resize(m_dimensionality);
      for (CFuint iCoor = 0; iCoor < dimUns; ++iCoor)
      {
        m_volTermTensor[iCV][iFlx][iCoor] = 0.0;
      }
    }
  }

  // set the dimensionality and order of the simplex integrator
  const CFDim dimM1 = static_cast<CFDim>(m_dimensionality-1);
  m_sIntegrator.setDimensionality(dimM1);
  const CFPolyOrder::Type polyOrderP1 = static_cast<CFPolyOrder::Type>(m_polyOrder+1);
  m_sIntegrator.setIntegratorOrder(polyOrderP1);

  // loop over internal faces
  for (CFuint iFace = 0; iFace < nbrLocalIntFaces; ++iFace)
  {
    // local face ID
    const CFuint faceID = iFace + nbrLocalExtFaces;

    // face unit normal
    /// @warning the normal in the first quadrature point is taken here!!!
    RealVector unitNormal = m_intFaceQuadPntNorm[iFace][0]/m_intFaceQuadPntNorm[iFace][0].norm2();

    // number of control volume nodes
    const CFuint nbrFaceNodes = m_localFaceNodeConn[faceID].size();

    // number of simplices in face
    const CFuint nbrSimplices = nbrFaceNodes - dimM1;

    // get node coordinates
    vector< RealVector > faceNodeCoord(nbrFaceNodes);
    for (CFuint iNode = 0; iNode < nbrFaceNodes; ++iNode)
    {
      faceNodeCoord[iNode].resize(m_dimensionality);
      const CFuint nodeID = m_localFaceNodeConn[faceID][iNode];
      faceNodeCoord[iNode] = m_localNodeCoord[nodeID];
    }

    // get quadrature nodes and wheights
    vector< RealVector > qNodeCoord;
    vector< CFreal > qWheights;
    for (CFuint iSimplex = 0; iSimplex < nbrSimplices; ++iSimplex)
    {
      // get node coordinates of simplex
      vector< RealVector > simpNodeCoord(m_dimensionality);
      simpNodeCoord[0].resize(m_dimensionality);
      simpNodeCoord[0] = faceNodeCoord[0];
      for (CFuint iNode = 1; iNode < dimUns; ++iNode)
      {
        simpNodeCoord[iNode].resize(m_dimensionality);
        simpNodeCoord[iNode] = faceNodeCoord[iSimplex+iNode];
      }

      // get simplex quadrature nodes and wheights
      vector< RealVector > qNodeCoordSimp = m_sIntegrator.getQuadPntsCoordsPlus1D  (simpNodeCoord);
      vector< CFreal >     qWheightsSimp  = m_sIntegrator.getQuadPntsWheightsPlus1D(simpNodeCoord);

      // add simplex quadrature nodes and wheights to global list
      qNodeCoord.insert(qNodeCoord.end(),qNodeCoordSimp.begin(),qNodeCoordSimp.end());
      qWheights.insert(qWheights.end(),qWheightsSimp.begin(),qWheightsSimp.end());
    }

    // compute integral of flux polynomial over face
    const CFuint nbrQNodes = qWheights.size();
    for (CFuint iFlx = 0; iFlx < nbrFlxPnts; ++iFlx)
    {
      CFreal integral = 0.0;
      for (CFuint iQNode = 0; iQNode < nbrQNodes; ++iQNode)
      {
        for (CFuint iTerm = 0; iTerm < nbrFlxPnts; ++iTerm)
        {
          CFreal term = m_fluxPolyCoefs[iFlx][iTerm];
          for (CFuint iCoor = 0; iCoor < dimUns; ++iCoor)
          {
            term *= pow(qNodeCoord[iQNode][iCoor],m_fluxPolyExponents[iTerm][iCoor]);
          }
          integral += qWheights[iQNode]*term;
        }
      }

      // add contribution to volume term tensor
      for (CFuint iCoor = 0; iCoor < dimUns; ++iCoor)
      {
        m_volTermTensor[m_localIntFaceCVConn[iFace][LEFT ]][iFlx][iCoor] -= integral*unitNormal[iCoor];
        m_volTermTensor[m_localIntFaceCVConn[iFace][RIGHT]][iFlx][iCoor] += integral*unitNormal[iCoor];
      }
    }
  }

/*  for (CFuint iCV = 0; iCV < nbrCVs; ++iCV)
  {
    CF_DEBUG_OBJ(iCV);
    for (CFuint iFlx = 0; iFlx < nbrFlxPnts; ++iFlx)
    {
      for (CFuint iCoor = 0; iCoor < m_dimensionality; ++iCoor)
      {
      CF_DEBUG_OBJ(m_volTermTensor[iCV][iFlx][iCoor]);
      }
    }
  }*/
}

//////////////////////////////////////////////////////////////////////

void SpectralFVElementData::createSolInFaceFlxPntsCoefs()
{
  CFAUTOTRACE;

  // get the number of cell faces
  const CFuint nbrCellFaces = m_faceFlxPntsConn.size();

  // get the number of CVs
  const CFuint nbrCVs = getNbrOfCVs();

  // get the number of face flux points
  const CFuint nbrFaceFlxPnts = getNbrOfFaceFlxPnts();

  // resize m_solInFaceFluxPntCoef
  m_solInFaceFluxPntCoef.resize(nbrCellFaces);

  // set coefficients
  for (CFuint iFace = 0; iFace < nbrCellFaces; ++iFace)
  {
    // resize
    m_solInFaceFluxPntCoef[iFace].resize(nbrFaceFlxPnts);
    for (CFuint iFlx = 0; iFlx < nbrFaceFlxPnts; ++iFlx)
    {
      // resize
      m_solInFaceFluxPntCoef[iFace][iFlx].resize(nbrCVs);

      // get the flux point ID
      const CFuint flxPntID = m_faceFlxPntsConn[iFace][iFlx];

      // set the coefficients
      m_solInFaceFluxPntCoef[iFace][iFlx] = m_solInFluxPntCoef[flxPntID];
    }
  }

/*  for (CFuint iFace = 0; iFace < m_solInFaceFluxPntCoef.size(); ++iFace)
  {
    CF_DEBUG_OBJ(iFace);
    for (CFuint iFlx = 0; iFlx < m_solInFaceFluxPntCoef[iFace].size(); ++iFlx)
    {
      CF_DEBUG_OBJ(iFlx);
      for (CFuint iCV = 0; iCV < m_solInFaceFluxPntCoef[iFace][iFlx].size(); ++iCV)
      {
        CF_DEBUG_OBJ(m_solInFaceFluxPntCoef[iFace][iFlx][iCV]);
      }
    }
  }*/
}

//////////////////////////////////////////////////////////////////////

void SpectralFVElementData::createSolInFaceFlxPntsCoefsPerOrient()
{
  CFAUTOTRACE;

  const CFuint dimUns = static_cast<CFuint>(m_dimensionality);

  // get number of orientations
  const CFuint nbrOrients = m_svFaceNodeCoordsPerOrient.size();

  // get the number of CVs
  const CFuint nbrCVs = getNbrOfCVs();

  // get the number of face flux points
  const CFuint nbrFaceFlxPnts = getNbrOfFaceFlxPnts();

  // resize m_solInFaceFluxPntCoefPerOrient
  m_solInFaceFluxPntCoefPerOrient.resize(nbrOrients);
  for (CFuint iOrient = 0; iOrient < nbrOrients; ++iOrient)
  {
    m_solInFaceFluxPntCoefPerOrient[iOrient].resize(2);
    for (CFuint iSide = 0; iSide < 2; ++iSide)
    {
      m_solInFaceFluxPntCoefPerOrient[iOrient][iSide].resize(nbrFaceFlxPnts);
      for (CFuint iFlx = 0; iFlx < nbrFaceFlxPnts; ++iFlx)
      {
        m_solInFaceFluxPntCoefPerOrient[iOrient][iSide][iFlx].resize(nbrCVs);
      }
    }
  }

  // compute face flux point solution reconstruction coefficients per face connectivity orientation
  for (CFuint iOrient = 0; iOrient < nbrOrients; ++iOrient)
  {
    for (CFuint iSide = 0; iSide < 2; ++iSide)
    {
      // vector for face node coordinates
      vector< RealVector > faceNodeCoor(m_dimensionality);
      for (CFuint iNode = 0; iNode < dimUns; ++iNode)
      {
        faceNodeCoor[iNode].resize(m_dimensionality);
        faceNodeCoor[iNode] = m_svFaceNodeCoordsPerOrient[iOrient][iSide][iNode];
      }

      // face flux point local coordinates
      vector< RealVector > flxPntsLocCoords(nbrFaceFlxPnts);
      for (CFuint iFlx = 0; iFlx < nbrFaceFlxPnts; ++iFlx)
      {
        flxPntsLocCoords[iFlx].resize(m_dimensionality);

        flxPntsLocCoords[iFlx] = 0.0;
        for (CFuint iNode = 0; iNode < dimUns; ++iNode)
        {
          flxPntsLocCoords[iFlx] += m_faceFluxPolyNodeWheightCoord[iFlx][iNode]*faceNodeCoor[iNode];
        }
      }

      // compute solution reconstruction coefficients
      m_solInFaceFluxPntCoefPerOrient[iOrient][iSide] = getSVPolyValsAtNode(flxPntsLocCoords);
    }
  }

/*  for (CFuint iOrient = 0; iOrient < m_solInFaceFluxPntCoefPerOrient.size(); ++iOrient)
  {
    CF_DEBUG_OBJ(iOrient);
    for (CFuint iFlx = 0; iFlx < m_solInFaceFluxPntCoefPerOrient[iOrient][LEFT ].size(); ++iFlx)
    {
      CF_DEBUG_OBJ(iFlx);
      for (CFuint iCV = 0; iCV < m_solInFaceFluxPntCoefPerOrient[iOrient][LEFT ][iFlx].size(); ++iCV)
      {
        CF_DEBUG_OBJ(m_solInFaceFluxPntCoefPerOrient[iOrient][LEFT ][iFlx][iCV]);
      }
      for (CFuint iCV = 0; iCV < m_solInFaceFluxPntCoefPerOrient[iOrient][LEFT ][iFlx].size(); ++iCV)
      {
        CF_DEBUG_OBJ(m_solInFaceFluxPntCoefPerOrient[iOrient][RIGHT][iFlx][iCV]);
      }
    }
  }*/
}

//////////////////////////////////////////////////////////////////////

void SpectralFVElementData::createSolInFaceFlxPntsCoefsPerOrientNoSymm()
{
  CFAUTOTRACE;

  const CFuint dimUns = static_cast<CFuint>(m_dimensionality);

  // get number of orientations
  const CFuint nbrOrients = m_svFaceNodeCoordsPerOrientNoSymm.size();

  // get the number of CVs
  const CFuint nbrCVs = getNbrOfCVs();

  // get the number of face flux points
  const CFuint nbrFaceFlxPnts = getNbrOfFaceFlxPnts();

  // resize m_solInFaceFluxPntCoefPerOrient
  m_solInFaceFluxPntCoefPerOrientNoSymm.resize(nbrOrients);
  for (CFuint iOrient = 0; iOrient < nbrOrients; ++iOrient)
  {
    m_solInFaceFluxPntCoefPerOrientNoSymm[iOrient].resize(2);
    for (CFuint iSide = 0; iSide < 2; ++iSide)
    {
      m_solInFaceFluxPntCoefPerOrientNoSymm[iOrient][iSide].resize(nbrFaceFlxPnts);
      for (CFuint iFlx = 0; iFlx < nbrFaceFlxPnts; ++iFlx)
      {
        m_solInFaceFluxPntCoefPerOrientNoSymm[iOrient][iSide][iFlx].resize(nbrCVs);
      }
    }
  }

  // compute face flux point solution reconstruction coefficients per face connectivity orientation
  for (CFuint iOrient = 0; iOrient < nbrOrients; ++iOrient)
  {
    for (CFuint iSide = 0; iSide < 2; ++iSide)
    {
      // vector for face node coordinates
      vector< RealVector > faceNodeCoor(m_dimensionality);
      for (CFuint iNode = 0; iNode < dimUns; ++iNode)
      {
        faceNodeCoor[iNode].resize(m_dimensionality);
        faceNodeCoor[iNode] = m_svFaceNodeCoordsPerOrientNoSymm[iOrient][iSide][iNode];
      }

      // face flux point local coordinates
      vector< RealVector > flxPntsLocCoords(nbrFaceFlxPnts);
      for (CFuint iFlx = 0; iFlx < nbrFaceFlxPnts; ++iFlx)
      {
        flxPntsLocCoords[iFlx].resize(m_dimensionality);

        flxPntsLocCoords[iFlx] = 0.0;
        for (CFuint iNode = 0; iNode < dimUns; ++iNode)
        {
          flxPntsLocCoords[iFlx] += m_faceFluxPolyNodeWheightCoord[iFlx][iNode]*faceNodeCoor[iNode];
        }
      }

      // compute solution reconstruction coefficients
      m_solInFaceFluxPntCoefPerOrientNoSymm[iOrient][iSide] = getSVPolyValsAtNode(flxPntsLocCoords);
    }
  }
}

//////////////////////////////////////////////////////////////////////

void SpectralFVElementData::createCVExtFaceFluxCoef()
{
  CFAUTOTRACE;

  const CFuint dimUns = static_cast<CFuint>(m_dimensionality);

  // number of CV faces at SV face
  /// @todo change the quadrature free implementation to take into account elements with different kinds of faces
  const CFuint nbrCVFacesAtSVFace = m_extSVFaceCVConn[0].size();

  // number of flux points
  const CFuint nbrFlxPnts = getNbrOfFlxPnts();

  // get the number of face flux points
  const CFuint nbrFaceFlxPnts = getNbrOfFaceFlxPnts();

  // set the dimensionality and order of the simplex integrator
  const CFDim dimM1 = static_cast<CFDim>(m_dimensionality-1);
  m_sIntegrator.setDimensionality(dimM1);
  const CFPolyOrder::Type polyOrderP1 = static_cast<CFPolyOrder::Type>(m_polyOrder+1);
  m_sIntegrator.setIntegratorOrder(polyOrderP1);

  // compute coefficients, using the first SV face (for which the last local coordinate should be zero)
  // resize m_cvExtFaceFluxCoef
  m_cvExtFaceFluxCoef.resize(nbrCVFacesAtSVFace);
  for (CFuint iCVFace = 0; iCVFace < nbrCVFacesAtSVFace; ++iCVFace)
  {
    // resize m_cvExtFaceFluxCoef[iCVFace]
    m_cvExtFaceFluxCoef[iCVFace].resize(nbrFaceFlxPnts);

    // number of CV face nodes
    const CFuint nbrFaceNodes = m_localFaceNodeConn[iCVFace].size();

    // number of simplices in face
    const CFuint nbrSimplices = nbrFaceNodes - dimM1;

    // get node coordinates
    vector< RealVector > faceNodeCoord(nbrFaceNodes);
    for (CFuint iNode = 0; iNode < nbrFaceNodes; ++iNode)
    {
      faceNodeCoord[iNode].resize(m_dimensionality);
      const CFuint nodeID = m_localFaceNodeConn[iCVFace][iNode];
      faceNodeCoord[iNode] = m_localNodeCoord[nodeID];
    }

    // get quadrature nodes and wheights
    vector< RealVector > qNodeCoord;
    vector< CFreal > qWheights;
    for (CFuint iSimplex = 0; iSimplex < nbrSimplices; ++iSimplex)
    {
      // get node coordinates of simplex
      vector< RealVector > simpNodeCoord(m_dimensionality);
      simpNodeCoord[0].resize(m_dimensionality);
      simpNodeCoord[0] = faceNodeCoord[0];
      for (CFuint iNode = 1; iNode < dimUns; ++iNode)
      {
        simpNodeCoord[iNode].resize(m_dimensionality);
        simpNodeCoord[iNode] = faceNodeCoord[iSimplex+iNode];
      }

      // get simplex quadrature nodes and wheights
      vector< RealVector > qNodeCoordSimp = m_sIntegrator.getQuadPntsCoordsPlus1D(simpNodeCoord);
      vector< CFreal >     qWheightsSimp  = m_sIntegrator.getQuadPntsWheightsPlus1D(simpNodeCoord);

      // add simplex quadrature nodes and wheights to global list
      qNodeCoord.insert(qNodeCoord.end(),qNodeCoordSimp.begin(),qNodeCoordSimp.end());
      qWheights .insert(qWheights .end(),qWheightsSimp .begin(),qWheightsSimp .end());
    }

    // compute integral of flux polynomial over face
    const CFuint nbrQNodes = qWheights.size();
    for (CFuint iFlx = 0; iFlx < nbrFaceFlxPnts; ++iFlx)
    {
      // flux point ID
      cf_assert(m_faceFlxPntsConn.size() > 0);
      const CFuint flxID = m_faceFlxPntsConn[0][iFlx];

      m_cvExtFaceFluxCoef[iCVFace][iFlx] = 0.0;
      for (CFuint iQNode = 0; iQNode < nbrQNodes; ++iQNode)
      {
        for (CFuint iTerm = 0; iTerm < nbrFlxPnts; ++iTerm)
        {
          CFreal term = m_fluxPolyCoefs[flxID][iTerm];
          for (CFuint iCoor = 0; iCoor < dimUns; ++iCoor)
          {
            term *= pow(qNodeCoord[iQNode][iCoor],m_fluxPolyExponents[iTerm][iCoor]);
          }
          m_cvExtFaceFluxCoef[iCVFace][iFlx] += qWheights[iQNode]*term;
        }
      }

      // sum of coefficients should be equal to one (face surface is taken into account seperately)
      for (CFuint iDim = 2; iDim < dimUns; ++iDim)
      {
        m_cvExtFaceFluxCoef[iCVFace][iFlx] *= iDim;
      }
    }
  }

/*  for (CFuint iCVFace = 0; iCVFace < nbrCVFacesAtSVFace; ++iCVFace)
  {
    CF_DEBUG_OBJ(iCVFace);
    for (CFuint iFlx = 0; iFlx < nbrFaceFlxPnts; ++iFlx)
    {
      CF_DEBUG_OBJ(m_cvExtFaceFluxCoef[iCVFace][iFlx]);
    }
  }*/
}

//////////////////////////////////////////////////////////////////////

void SpectralFVElementData::createAvgSolInSVFaceCoef()
{
  CFAUTOTRACE;

  const CFuint dimUns = static_cast<CFuint>(m_dimensionality);

  // number of CV faces at SV face
  /// @todo change the quadrature free implementation to take into account elements with different kinds of faces
  const CFuint nbrCVFacesAtSVFace = m_extSVFaceCVConn[0].size();

  // get the number of face flux points
  const CFuint nbrFaceFlxPnts = getNbrOfFaceFlxPnts();

  // resize m_avgSolInSVFaceCoef
  m_avgSolInSVFaceCoef.resize(nbrFaceFlxPnts);
  for (CFuint iFlx = 0; iFlx < nbrFaceFlxPnts; ++iFlx)
  {
    m_avgSolInSVFaceCoef[iFlx] = 0.0;
    for (CFuint iCVFace = 0; iCVFace < nbrCVFacesAtSVFace; ++iCVFace)
    {
      m_avgSolInSVFaceCoef[iFlx] += m_cvExtFaceFluxCoef[iCVFace][iFlx];
    }

    for (CFuint iDim = 1; iDim < dimUns; ++iDim)
    {
      m_avgSolInSVFaceCoef[iFlx] *= iDim;
    }
// CF_DEBUG_OBJ(m_avgSolInSVFaceCoef[iFlx]);
  }
}

//////////////////////////////////////////////////////////////////////

void SpectralFVElementData::createFaceOutputPntSolPolyAndDerivCoef()
{
  CFAUTOTRACE;

  // get number of cell faces
  const CFuint nbrCellFaces = getNbrSVFaces();
  cf_assert(m_faceOutputPntCellMappedCoords.size() == nbrCellFaces);

  // compute polynomial values in cell face output points
  m_faceOutputPntSolPolyCoef .resize(nbrCellFaces);
  m_faceOutputPntSolDerivCoef.resize(nbrCellFaces);
  for (CFuint iFace = 0; iFace < nbrCellFaces; ++iFace)
  {
    m_faceOutputPntSolPolyCoef [iFace] = getSVPolyValsAtNode   (m_faceOutputPntCellMappedCoords[iFace]);
//    m_faceOutputPntSolDerivCoef[iFace] = getSolPolyDerivsAtNode(m_faceOutputPntCellMappedCoords[iFace]);
  }
}

//////////////////////////////////////////////////////////////////////

void SpectralFVElementData::InvertMatrix(RealMatrix A, RealMatrix& AI)
{
  cf_assert(A.nbRows() == AI.nbRows());
  cf_assert(A.nbCols() == AI.nbCols());

  const CFuint n = A.nbRows();

  for (CFuint i = 0; i < n; ++i)
  {
    for (CFuint j = 0; j < n; ++j)
    {
      AI(i,j) = 0.;
    }
    AI(i,i) = 1.;
  }

  for (CFuint i = 0; i < n; ++i)
  {
    CFreal fac = fabs(A(i,i));
    CFuint k = i;
    for (CFuint j = i+1; j < n; ++j)
    {
      if (fabs(A(j,i)) > fac)
      {
        fac = fabs(A(j,i));
        k = j;
      }
    }

    if (fac < MathTools::MathConsts::CFrealEps()) throw ZeroDeterminantException (FromHere(),"Matrix is singular to working precision!!!");

    SwapRows(A ,i,k);
    SwapRows(AI,i,k);

    fac = 1./A(i,i);
    A(i,i) = 1.;
    for (CFuint j = i+1; j < n; ++j)
    {
      A(i,j) = A(i,j)*fac;
    }
    for (CFuint j = 0; j < n; ++j)
    {
      AI(i,j) = AI(i,j)*fac;
    }
    for (CFuint k = i+1; k < n; ++k)
    {
      fac = A(k,i);
      for (CFuint j = i+1; j < n; ++j)
      {
        A(k,j) = A(k,j) - A(i,j)*fac;
      }
      for (CFuint j = 0; j < n; ++j)
      {
        AI(k,j) = AI(k,j) - AI(i,j)*fac;
      }
    }
  }

  for (CFuint i = 0; i < n-1; ++i)
  {
    const CFuint ii = n-2-i;
    for (CFuint j = 0; j < n; ++j)
    {
      for (CFuint k = ii+1; k < n; ++k)
      {
        AI(ii,j) = AI(ii,j) - AI(k,j)*A(ii,k);
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////

void SpectralFVElementData::SwapRows(RealMatrix& A, CFuint row1, CFuint row2)
{
  RealVector swapRow = A.getRow<RealVector>(row1);
  A.setRow(A.getRow<RealVector>(row2),row1);
  A.setRow(swapRow,row2);
}

//////////////////////////////////////////////////////////////////////

  } // namespace SpectralFV

} // namespace COOLFluiD
