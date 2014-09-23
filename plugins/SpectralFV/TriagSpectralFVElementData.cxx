#include <fstream>


#include "Common/CFLog.hh"
#include "Environment/DirPaths.hh"
#include "Common/NotImplementedException.hh"
#include "MathTools/RealMatrix.hh"
#include "SpectralFV/TriagSpectralFVElementData.hh"
#include "MathTools/MathConsts.hh"

//////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Common;

namespace COOLFluiD {

  namespace SpectralFV {

//////////////////////////////////////////////////////////////////////

TriagSpectralFVElementData::TriagSpectralFVElementData() :
  SpectralFVElementData()
{
  m_shape = CFGeoShape::TRIAG;
  m_dimensionality = DIM_2D;
}

//////////////////////////////////////////////////////////////////////

TriagSpectralFVElementData::TriagSpectralFVElementData(CFPolyOrder::Type polyOrder)
{
  m_shape = CFGeoShape::TRIAG;
  m_dimensionality = DIM_2D;
  m_polyOrder = polyOrder;

  resetSpectralFVElementData();
}

//////////////////////////////////////////////////////////////////////

TriagSpectralFVElementData::~TriagSpectralFVElementData()
{
}

//////////////////////////////////////////////////////////////////////

void TriagSpectralFVElementData::computeSVNodeLocalCoords()
{
  m_svNodeCoords.resize(3);

  // first node
  m_svNodeCoords[0].resize(2);
  m_svNodeCoords[0][KSI] = 0.0;
  m_svNodeCoords[0][ETA] = 0.0;

  // second node
  m_svNodeCoords[1].resize(2);
  m_svNodeCoords[1][KSI] = 1.0;
  m_svNodeCoords[1][ETA] = 0.0;

  // third node
  m_svNodeCoords[2].resize(2);
  m_svNodeCoords[2][KSI] = 0.0;
  m_svNodeCoords[2][ETA] = 1.0;
}

//////////////////////////////////////////////////////////////////////

void TriagSpectralFVElementData::createLocalNodeCoord()
{
  CFAUTOTRACE;

  // Path where to read files from
  std::string filename;
  const std::string append = Environment::DirPaths::getInstance().getWorkingDir().string() + "/";

  // number of local nodes
  const CFuint nbrLocalNodes = (m_polyOrder+1)*(m_polyOrder+1) + m_polyOrder+2;

  // resize the vectors containing the node coordinates
  m_localNodeCoord.resize(nbrLocalNodes);
  for (CFuint iNode = 0; iNode < nbrLocalNodes; ++iNode)
  {
    m_localNodeCoord[iNode].resize(2);
  }

  switch (m_polyOrder)
  {
    case CFPolyOrder::ORDER0:
    {
      // fill the matrix containing the node coordinates
      m_localNodeCoord[0][KSI] = 0.0;
      m_localNodeCoord[0][ETA] = 0.0;

      m_localNodeCoord[1][KSI] = 1.0;
      m_localNodeCoord[1][ETA] = 0.0;

      m_localNodeCoord[2][KSI] = 0.0;
      m_localNodeCoord[2][ETA] = 1.0;

    } break;
    case CFPolyOrder::ORDER1:
    {
      // fill the matrix containing the node coordinates
      m_localNodeCoord[0][KSI] = 0.0;
      m_localNodeCoord[0][ETA] = 0.0;

      m_localNodeCoord[1][KSI] = 0.5;
      m_localNodeCoord[1][ETA] = 0.0;

      m_localNodeCoord[2][KSI] = 1.0/3.0;
      m_localNodeCoord[2][ETA] = 1.0/3.0;

      m_localNodeCoord[3][KSI] = 0.0;
      m_localNodeCoord[3][ETA] = 0.5;

      m_localNodeCoord[4][KSI] = 1.0;
      m_localNodeCoord[4][ETA] = 0.0;

      m_localNodeCoord[5][KSI] = 0.5;
      m_localNodeCoord[5][ETA] = 0.5;

      m_localNodeCoord[6][KSI] = 0.0;
      m_localNodeCoord[6][ETA] = 1.0;

    } break;
    case CFPolyOrder::ORDER2:
    {
      // partition parameters
      CFreal alpha, beta;

      // read defining parameters from file
      filename = append + "SV3TRIAGDEF.DAT";
      ifstream inputFile;
      inputFile.open(filename.c_str(), ios::in);

      if (inputFile.is_open())
      {
        inputFile >> alpha;
        inputFile >> beta;
      }
      else
      {
        alpha = 0.091;
        beta  = 0.18;
/*        alpha = 1.0/4.0;
        beta  = 1.999/3.0;*/
/*        alpha = 1.0/4.0;
        beta  = 1.0/4.0;*/
/*        alpha = 1.0/4.0;
        beta  = 1.0/3.0;*/
      }
      inputFile.close();

      CF_DEBUG_OBJ(alpha);
      CF_DEBUG_OBJ(beta);

      // fill the matrix containing the node coordinates
      m_localNodeCoord[0][KSI] = 0.0;
      m_localNodeCoord[0][ETA] = 0.0;

      m_localNodeCoord[1][KSI] = alpha;
      m_localNodeCoord[1][ETA] = 0.0;

      m_localNodeCoord[2][KSI] = 0.5*beta;
      m_localNodeCoord[2][ETA] = 0.5*beta;

      m_localNodeCoord[3][KSI] = 0.0;
      m_localNodeCoord[3][ETA] = alpha;

      m_localNodeCoord[4][KSI] = 1.0 - alpha;
      m_localNodeCoord[4][ETA] = 0.0;

      m_localNodeCoord[5][KSI] = 1.0 - beta;
      m_localNodeCoord[5][ETA] = 0.5*beta;

      m_localNodeCoord[6][KSI] = 1.0/3.0;
      m_localNodeCoord[6][ETA] = 1.0/3.0;

      m_localNodeCoord[7][KSI] = 0.5*beta;
      m_localNodeCoord[7][ETA] = 1.0 - beta;

      m_localNodeCoord[8][KSI] = 0.0;
      m_localNodeCoord[8][ETA] = 1.0 - alpha;

      m_localNodeCoord[9][KSI] = 1.0;
      m_localNodeCoord[9][ETA] = 0.0;

      m_localNodeCoord[10][KSI] = 1.0 - alpha;
      m_localNodeCoord[10][ETA] = alpha;

      m_localNodeCoord[11][KSI] = alpha;
      m_localNodeCoord[11][ETA] = 1.0 - alpha;

      m_localNodeCoord[12][KSI] = 0.0;
      m_localNodeCoord[12][ETA] = 1.0;
    } break;
    case CFPolyOrder::ORDER3:
    {
      // partition parameters
      CFreal alpha, beta, gamma, delta;

      // read defining parameters from file
      filename = append + "SV4TRIAGDEF.DAT";
      ifstream inputFile;
      inputFile.open(filename.c_str(), ios::in);

      if (inputFile.is_open())
      {
        inputFile >> alpha;
        inputFile >> beta;
        inputFile >> gamma;
        inputFile >> delta;
      }
      else
      {
        alpha = 0.078;
        beta  = 0.104;
        gamma = 0.052;
        delta = 0.351;
/*        alpha = 1.0/15.0;
        beta  = 2.00000/15.0;
        gamma = 1.0/15.0;
        delta = 2.00001/15.0;*/
      }
      inputFile.close();

      // fill the matrix containing the node coordinates
      m_localNodeCoord[0][KSI] = 0.0;
      m_localNodeCoord[0][ETA] = 0.0;

      m_localNodeCoord[1][KSI] = alpha;
      m_localNodeCoord[1][ETA] = 0.0;

      m_localNodeCoord[2][KSI] = 0.5*beta;
      m_localNodeCoord[2][ETA] = 0.5*beta;

      m_localNodeCoord[3][KSI] = 0.0;
      m_localNodeCoord[3][ETA] = alpha;

      m_localNodeCoord[4][KSI] = 0.5;
      m_localNodeCoord[4][ETA] = 0.0;

      m_localNodeCoord[5][KSI] = 0.5 - 0.5*gamma;
      m_localNodeCoord[5][ETA] = gamma;

      m_localNodeCoord[6][KSI] = 0.5*delta;
      m_localNodeCoord[6][ETA] = 0.5*delta;

      m_localNodeCoord[7][KSI] = gamma;
      m_localNodeCoord[7][ETA] = 0.5 - 0.5*gamma;

      m_localNodeCoord[8][KSI] = 0.0;
      m_localNodeCoord[8][ETA] = 0.5;

      m_localNodeCoord[9][KSI] = 1.0 - alpha;
      m_localNodeCoord[9][ETA] = 0.0;

      m_localNodeCoord[10][KSI] = 1.0 - beta;
      m_localNodeCoord[10][ETA] = 0.5*beta;

      m_localNodeCoord[11][KSI] = 1.0 - delta;
      m_localNodeCoord[11][ETA] = 0.5*delta;

      m_localNodeCoord[12][KSI] = 0.5 - 0.5*gamma;
      m_localNodeCoord[12][ETA] = 0.5 - 0.5*gamma;

      m_localNodeCoord[13][KSI] = 0.5*delta;
      m_localNodeCoord[13][ETA] = 1.0 - delta;

      m_localNodeCoord[14][KSI] = 0.5*beta;
      m_localNodeCoord[14][ETA] = 1.0 - beta;

      m_localNodeCoord[15][KSI] = 0.0;
      m_localNodeCoord[15][ETA] = 1.0 - alpha;

      m_localNodeCoord[16][KSI] = 1.0;
      m_localNodeCoord[16][ETA] = 0.0;

      m_localNodeCoord[17][KSI] = 1.0 - alpha;
      m_localNodeCoord[17][ETA] = alpha;

      m_localNodeCoord[18][KSI] = 0.5;
      m_localNodeCoord[18][ETA] = 0.5;

      m_localNodeCoord[19][KSI] = alpha;
      m_localNodeCoord[19][ETA] = 1.0 - alpha;

      m_localNodeCoord[20][KSI] = 0.0;
      m_localNodeCoord[20][ETA] = 1.0;
    } break;
    default:
    {
      throw Common::NotImplementedException (FromHere(),"Only solution orders up to 3 have been implemented for spectral FV!");
    }
  }
}

//////////////////////////////////////////////////////////////////////

void TriagSpectralFVElementData::createLocalFaceNodeConn()
{
  CFAUTOTRACE;

  /*    .
        .
        .    CV5   .
        | \      .
        |   \  .
        |     \
        |  CV2  \   CV4    .
        |\       /\      .
        |  \   /    \  .
        |    \   CV1  \   CV3
        | CV0  \        \
        |________\________\.  .  .
  */

  // compute total number of faces
  const CFuint nbrLocalFaces = 3*(m_polyOrder+2)*(m_polyOrder+1)/2;

  // helper variable
  const CFuint solOrderP1 = m_polyOrder+1;

  // resize the local face-node connectivity
  m_localFaceNodeConn.resize(nbrLocalFaces);
  for (CFuint iFace = 0; iFace < nbrLocalFaces; ++iFace)
  {
    m_localFaceNodeConn[iFace].resize(2);
  }

  // local face ID variable
  CFuint faceID = 0;

  // create local connectivity

  // external faces
  // first SV face
  for (CFuint iFace = 0; iFace < solOrderP1; ++iFace, ++faceID)
  {
    CFuint node0 = (iFace+1)*(iFace+1);
    CFuint node1 = iFace*iFace;
    m_localFaceNodeConn[faceID][0] = node0;
    m_localFaceNodeConn[faceID][1] = node1;
  }
  // second SV face
  CFuint nodeID = solOrderP1*solOrderP1;
  for (CFuint iFace = 0; iFace < solOrderP1; ++iFace, ++faceID)
  {
    m_localFaceNodeConn[faceID][0] = nodeID+iFace+1;
    m_localFaceNodeConn[faceID][1] = nodeID+iFace;
  }
  // third SV face
  m_localFaceNodeConn[faceID][0] = nodeID - 1;
  m_localFaceNodeConn[faceID][1] = nodeID + solOrderP1;
  ++faceID;
  for (CFuint iFace = m_polyOrder; iFace > 0; --iFace, ++faceID)
  {
    CFuint node0 = iFace*iFace - 1;
    CFuint node1 = (iFace+1)*(iFace+1) - 1;
    m_localFaceNodeConn[faceID][0] = node0;
    m_localFaceNodeConn[faceID][1] = node1;
  }

  // internal faces
  // first for (positively) diagonal faces in sketch
  const CFuint polyOrderUns = static_cast<CFuint>(m_polyOrder);
  for (CFuint iRow = 0; iRow < polyOrderUns; ++iRow)
  {
    CFuint upperNodeID = iRow*iRow + 1;
    CFuint lowerNodeID = (iRow+1)*(iRow+1) + 2;
    for (CFuint iFace = 0; iFace < iRow; ++iFace, ++faceID)
    {
      m_localFaceNodeConn[faceID][0] = upperNodeID;
      m_localFaceNodeConn[faceID][1] = lowerNodeID;
      upperNodeID += 2;
      lowerNodeID += 2;
    }
  }
  CFuint upperNodeID = m_polyOrder*m_polyOrder + 1;
  CFuint lowerNodeID = solOrderP1*solOrderP1 + 1;
  for (CFuint iFace = 0; iFace < polyOrderUns; ++iFace, ++faceID)
  {
    m_localFaceNodeConn[faceID][0] = upperNodeID;
    m_localFaceNodeConn[faceID][1] = lowerNodeID;
    upperNodeID += 2;
    lowerNodeID += 1;
  }

  // then for (negatively) diagonal faces in sketch
  nodeID = 1;
  CFuint nbrFacesInRow = 2;
  for (CFuint iRow = 0; iRow < polyOrderUns; ++iRow, nbrFacesInRow += 2, ++nodeID)
  {
    for (CFuint iFace = 0; iFace < nbrFacesInRow; ++iFace, ++faceID, ++nodeID)
    {
      m_localFaceNodeConn[faceID][0] = nodeID;
      m_localFaceNodeConn[faceID][1] = nodeID+1;
    }
  }
}

//////////////////////////////////////////////////////////////////////

void TriagSpectralFVElementData::createLocalCVNodeConn()
{
  CFAUTOTRACE;

  // number of CVs
  const CFuint nbrCVs = (m_polyOrder+1)*(m_polyOrder+2)/2;

  // resize local CV-node connectivity
  m_localCVNodeConn.resize(nbrCVs);

  // create connectivity
  if (nbrCVs == 1)
  {
    // resize local CV-node connectivity
    m_localCVNodeConn[0].resize(3);

    // fill connectivity
    m_localCVNodeConn[0][0] = 0;
    m_localCVNodeConn[0][1] = 1;
    m_localCVNodeConn[0][2] = 2;
  }
  else
  {
    // control volume ID
    CFuint cvID = 0;

    // first nodes
    CFuint upperNode = 0;
    CFuint lowerNode = 1;

    // first row (CV 0)
    // resize local CV-node connectivity
    m_localCVNodeConn[cvID].resize(4);
    // fill connectivity
    m_localCVNodeConn[cvID][0] = lowerNode;
    m_localCVNodeConn[cvID][1] = lowerNode+1;
    m_localCVNodeConn[cvID][2] = lowerNode+2;
    m_localCVNodeConn[cvID][3] = upperNode;
    ++cvID;
    ++upperNode;
    lowerNode += 3;

    // intermediate rows
    const CFuint polyOrderUns = static_cast<CFuint>(m_polyOrder);
    for (CFuint iRow = 1; iRow < polyOrderUns; ++iRow)
    {
      // first element of row
      // resize local CV-node connectivity
      m_localCVNodeConn[cvID].resize(5);
      // fill connectivity
      m_localCVNodeConn[cvID][0] = lowerNode;
      m_localCVNodeConn[cvID][1] = lowerNode+1;
      m_localCVNodeConn[cvID][2] = lowerNode+2;
      m_localCVNodeConn[cvID][3] = upperNode+1;
      m_localCVNodeConn[cvID][4] = upperNode;
      ++cvID;
      lowerNode += 2;
      ++upperNode;

      // intermediate elements in row
      for (CFuint iElem = 1; iElem < iRow; ++iElem, ++cvID, lowerNode += 2, upperNode += 2)
      {
        // resize local CV-node connectivity
        m_localCVNodeConn[cvID].resize(6);
        // fill connectivity
        m_localCVNodeConn[cvID][0] = lowerNode;
        m_localCVNodeConn[cvID][1] = lowerNode+1;
        m_localCVNodeConn[cvID][2] = lowerNode+2;
        m_localCVNodeConn[cvID][3] = upperNode+2;
        m_localCVNodeConn[cvID][4] = upperNode+1;
        m_localCVNodeConn[cvID][5] = upperNode;
      }

      // last element of row
      // resize local CV-node connectivity
      m_localCVNodeConn[cvID].resize(5);
      // fill connectivity
      m_localCVNodeConn[cvID][0] = lowerNode;
      m_localCVNodeConn[cvID][1] = lowerNode+1;
      m_localCVNodeConn[cvID][2] = lowerNode+2;
      m_localCVNodeConn[cvID][3] = upperNode+1;
      m_localCVNodeConn[cvID][4] = upperNode;
      ++cvID;
      lowerNode += 3;
      upperNode += 2;
    }

    // last row
    // first element of row
    // resize local CV-node connectivity
    m_localCVNodeConn[cvID].resize(4);
    // fill connectivity
    m_localCVNodeConn[cvID][0] = lowerNode;
    m_localCVNodeConn[cvID][1] = lowerNode+1;
    m_localCVNodeConn[cvID][2] = upperNode+1;
    m_localCVNodeConn[cvID][3] = upperNode;
    ++cvID;
    ++lowerNode;
    ++upperNode;

    // intermediate elements in row
    for (CFuint iElem = 1; iElem < polyOrderUns; ++iElem, ++cvID, ++lowerNode, upperNode += 2)
    {
      // resize local CV-node connectivity
      m_localCVNodeConn[cvID].resize(5);
      // fill connectivity
      m_localCVNodeConn[cvID][0] = lowerNode;
      m_localCVNodeConn[cvID][1] = lowerNode+1;
      m_localCVNodeConn[cvID][2] = upperNode+2;
      m_localCVNodeConn[cvID][3] = upperNode+1;
      m_localCVNodeConn[cvID][4] = upperNode;
    }

    // last element of row
    // resize local CV-node connectivity
    m_localCVNodeConn[cvID].resize(4);
    // fill connectivity
    m_localCVNodeConn[cvID][0] = lowerNode;
    m_localCVNodeConn[cvID][1] = lowerNode+1;
    m_localCVNodeConn[cvID][2] = upperNode+1;
    m_localCVNodeConn[cvID][3] = upperNode;
  }
}

//////////////////////////////////////////////////////////////////////

void TriagSpectralFVElementData::computeLocalFaceNormals()
{
  CFAUTOTRACE;

  // get number of local internal faces
  const CFuint nbrLocalFaces    = m_localFaceNodeConn.size();
  const CFuint nbrLocalExtFaces = 3*(m_polyOrder+1);
  const CFuint nbrLocalIntFaces = nbrLocalFaces - nbrLocalExtFaces;


  // compute local internal face normals
  m_intFaceQuadPntNorm.resize(nbrLocalIntFaces);
  for (CFuint iIntFace = 0; iIntFace < nbrLocalIntFaces; ++iIntFace)
  {
    // number of quadrature points on this face
    const CFuint nbrQPnts = m_intFaceQuadWheights[iIntFace].size();
    m_intFaceQuadPntNorm[iIntFace].resize(nbrQPnts,RealVector(2));

    const CFuint node0 = m_localFaceNodeConn[iIntFace+nbrLocalExtFaces][0];
    const CFuint node1 = m_localFaceNodeConn[iIntFace+nbrLocalExtFaces][1];

    RealVector normal(2);
    normal[KSI] =  (m_localNodeCoord[node1][ETA] - m_localNodeCoord[node0][ETA]);
    normal[ETA] = -(m_localNodeCoord[node1][KSI] - m_localNodeCoord[node0][KSI]);
    normal /= normal.norm2(); // divide by size to create a unit normal

    for (CFuint iQPnt = 0; iQPnt < nbrQPnts; ++iQPnt)
    {
      m_intFaceQuadPntNorm[iIntFace][iQPnt] = normal;
    }
//     CF_DEBUG_OBJ(m_intFaceQuadPntNorm[iIntFace]);
  }
}

//////////////////////////////////////////////////////////////////////

void TriagSpectralFVElementData::computeExtFaceLocalNormals()
{
  CFAUTOTRACE;

  // number of SV faces
  const CFuint nbrSVFaces = 3;

  // compute the normals
  m_extFaceLocalNorm.resize(nbrSVFaces);
  for (CFuint iSVFace = 0; iSVFace < nbrSVFaces; ++iSVFace)
  {
    // resize the variable
    m_extFaceLocalNorm[iSVFace].resize(m_dimensionality);

    // compute normal
    m_extFaceLocalNorm[iSVFace][KSI] =  (m_svFaceNodeCoords[iSVFace][1][ETA] - m_svFaceNodeCoords[iSVFace][0][ETA]);
    m_extFaceLocalNorm[iSVFace][ETA] = -(m_svFaceNodeCoords[iSVFace][1][KSI] - m_svFaceNodeCoords[iSVFace][0][KSI]);

//     CF_DEBUG_OBJ(m_extFaceLocalNorm[iSVFace]);
  }
}

//////////////////////////////////////////////////////////////////////

void TriagSpectralFVElementData::createLocalIntFaceCVConn()
{
  CFAUTOTRACE;

  /*    .
        .
        .    CV5   .
        | \      .
        |   \  .
        |     \
        |  CV2  \   CV4    .
        |\       /\      .
        |  \   /    \  .
        |    \   CV1  \   CV3
        | CV0  \        \
        |________\________\.  .  .
  */

  const CFuint solOrderP1 = m_polyOrder+1;

  // compute number of internal faces
  const CFuint nbrIntFaces = 3*m_polyOrder*(m_polyOrder+1)/2;

  // resize the local internal face-CV connectivity
  m_localIntFaceCVConn.resize(nbrIntFaces);
  for (CFuint iFace = 0; iFace < nbrIntFaces; ++iFace)
  {
    m_localIntFaceCVConn[iFace].resize(2);
  }

  // create local connectivity
  CFuint faceID = 0;
  // first for (positively) diagonal faces in sketch
  CFint cvID = 0;
  for (CFuint iRow = 0; iRow < solOrderP1; ++iRow, ++cvID)
    for (CFuint iFace = 0; iFace < iRow; ++iFace, ++faceID, ++cvID)
    {
      m_localIntFaceCVConn[faceID][0] = cvID+1;
      m_localIntFaceCVConn[faceID][1] = cvID;
    }

  // then for (negatively) diagonal faces in sketch above
  CFuint upperCVID = 0;
  CFuint lowerCVID = 1;
  CFuint nbrFacesInRow = 2;
  const CFuint polyOrderUns = static_cast<CFuint>(m_polyOrder);
  for (CFuint iRow = 0; iRow < polyOrderUns; ++iRow, nbrFacesInRow += 2, ++lowerCVID)
  {
    for (CFuint iFace = 0; iFace < nbrFacesInRow; ++iFace, ++faceID)
    {
      m_localIntFaceCVConn[faceID][LEFT ] = upperCVID;
      m_localIntFaceCVConn[faceID][RIGHT] = lowerCVID;

      //increment CV IDs
      lowerCVID += ((iFace+1)%2);
      upperCVID += (iFace    %2);
    }
  }

/*  for (CFuint iFace = 0; iFace < nbrIntFaces; ++iFace)
  {
    CF_DEBUG_OBJ(m_localIntFaceCVConn[iFace][LEFT ]);
    CF_DEBUG_OBJ(m_localIntFaceCVConn[iFace][RIGHT]);
  }*/
}

//////////////////////////////////////////////////////////////////////

void TriagSpectralFVElementData::computeVolumeFractionsOfCVs()
{
  CFAUTOTRACE;

  // number of CVs in SV
  const CFuint nbrCVs = m_localCVNodeConn.size();

  // resize the vector containing the CV face fractions
  m_volFracCV.resize(nbrCVs);

  // compute CV volume fractions
  for (CFuint iCV = 0; iCV < nbrCVs; ++iCV)
  {
    const CFuint nbrNodesInCV = m_localCVNodeConn[iCV].size();

    // Create vector for node coordinates of CV
    vector< RealVector > cvNodesCoord(nbrNodesInCV);
    for (CFuint iNode = 0; iNode < nbrNodesInCV; ++iNode)
    {
      cvNodesCoord[iNode].resize(2);
      const CFuint localNodeID = m_localCVNodeConn[iCV][iNode];
      cvNodesCoord[iNode] = m_localNodeCoord[localNodeID];
    }

    // Compute the volume fraction of the CV
    m_volFracCV[iCV] = computePolygonSurface(cvNodesCoord)/0.5;
  }
}

//////////////////////////////////////////////////////////////////////

void TriagSpectralFVElementData::createLocalExtFaceCVConn()
{
  CFAUTOTRACE;

  // number of external faces at SV face
  const CFuint nbrCVsAtSVFace = m_polyOrder + 1;

  // total number of external faces
  const CFuint nbrExtFaces = 3*nbrCVsAtSVFace;

  // resize m_localExtFaceCVConn
  m_localExtFaceCVConn.resize(nbrExtFaces);

  // create local external face - CV connectivity
  CFuint extFaceID = 0;
  // first SV face
  for (CFuint iFace = 0; iFace < nbrCVsAtSVFace; ++iFace, ++extFaceID)
  {
    m_localExtFaceCVConn[extFaceID] = iFace*(iFace+1)/2;
  }

  // second SV face
  CFuint cvID = m_polyOrder*(m_polyOrder+1)/2;
  for (CFuint iFace = 0; iFace < nbrCVsAtSVFace; ++iFace, ++cvID, ++extFaceID)
  {
    m_localExtFaceCVConn[extFaceID] = cvID;
  }

  // third SV face
  for (CFuint iFace = nbrCVsAtSVFace; iFace > 0; --iFace, ++extFaceID)
  {
    m_localExtFaceCVConn[extFaceID] = iFace*(iFace+1)/2 - 1;
  }
}

//////////////////////////////////////////////////////////////////////

void TriagSpectralFVElementData::createSVFaceLocalExtFaceConn()
{
  CFAUTOTRACE;

  // number of external faces at SV face
  const CFuint nbrCVsAtSVFace = m_polyOrder + 1;

  // resize m_svFaceLocalExtFaceConn
  m_svFaceLocalExtFaceConn.resize(3);
  for (CFuint iSVFace = 0; iSVFace < 3; ++iSVFace)
  {
    m_svFaceLocalExtFaceConn[iSVFace].resize(nbrCVsAtSVFace);
  }

  // create SV face - local external face connectivity
  CFuint extFaceID = 0;
  for (CFuint iSVFace = 0; iSVFace < 3; ++iSVFace)
  {
    for (CFuint iExtFace = 0; iExtFace < nbrCVsAtSVFace; ++iExtFace, ++extFaceID)
    {
      m_svFaceLocalExtFaceConn[iSVFace][iExtFace] = extFaceID;
    }
  }
}

//////////////////////////////////////////////////////////////////////

void TriagSpectralFVElementData::createExtFaceNodeLocalCoords()
{
  CFAUTOTRACE;

  // number of CVs at SV face
  const CFuint nbrCVsAtSVFace = m_polyOrder + 1;

  // create node coordinates local to face
  m_extFaceNodeLocalCoords.resize(nbrCVsAtSVFace);
  for (CFuint iExtFace = 0; iExtFace < nbrCVsAtSVFace; ++iExtFace)
  {
    // number of nodes to this face
    const CFuint nbrFaceNodes = m_localFaceNodeConn[iExtFace].size();

    m_extFaceNodeLocalCoords[iExtFace].resize(nbrFaceNodes);
    for (CFuint iNode = 0; iNode < nbrFaceNodes; ++iNode)
    {
      const CFuint nodeID = m_localFaceNodeConn[iExtFace][iNode];
      m_extFaceNodeLocalCoords[iExtFace][iNode].resize(m_dimensionality-1);
      m_extFaceNodeLocalCoords[iExtFace][iNode][KSI] = m_localNodeCoord[nodeID][KSI];
    }
  }
}

//////////////////////////////////////////////////////////////////////

void TriagSpectralFVElementData::createSVFaceNodeConnectivity()
{
  CFAUTOTRACE;

  // number of SV faces
  const CFuint nbrSVFaces = 3;

  // resize the variable
  m_svFaceNodeConn.resize(nbrSVFaces);
  for (CFuint iSVFace = 0; iSVFace < nbrSVFaces; ++iSVFace)
  {
    m_svFaceNodeConn[iSVFace].resize(2);
  }

  // fill the variable
  // first SV face
  m_svFaceNodeConn[0][0] = 0;
  m_svFaceNodeConn[0][1] = 1;

  // second SV face
  m_svFaceNodeConn[1][0] = 1;
  m_svFaceNodeConn[1][1] = 2;

  // third SV face
  m_svFaceNodeConn[2][0] = 2;
  m_svFaceNodeConn[2][1] = 0;
}

//////////////////////////////////////////////////////////////////////

void TriagSpectralFVElementData::createSVFaceNodeConnectivityPerOrient()
{
  CFAUTOTRACE;

  // number of SV faces
  const CFuint nbrSVFaces = m_svFaceNodeConn.size();

  // number of possible orientations
  const CFuint nbrOrient = 6;

  // resize the variable
  m_svFaceNodeConnPerOrientation.resize(nbrOrient);
  for (CFuint iOrient = 0; iOrient < nbrOrient; ++iOrient)
  {
    m_svFaceNodeConnPerOrientation[iOrient].resize(2);
    for (CFuint iCell = 0; iCell < 2; ++iCell)
      m_svFaceNodeConnPerOrientation[iOrient][iCell].resize(2);
  }

  // fill the variable
  CFuint iOrient = 0;
  for (CFuint iSVFaceL = 0; iSVFaceL < nbrSVFaces; ++iSVFaceL)
  {
    for (CFuint iSVFaceR = iSVFaceL; iSVFaceR < nbrSVFaces; ++iSVFaceR, ++iOrient)
    {
      for (CFuint iNode = 0; iNode < 2; ++iNode)
      {
        m_svFaceNodeConnPerOrientation[iOrient][LEFT ][iNode] = m_svFaceNodeConn[iSVFaceL][iNode  ];
        m_svFaceNodeConnPerOrientation[iOrient][RIGHT][iNode] = m_svFaceNodeConn[iSVFaceR][1-iNode];
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////

void TriagSpectralFVElementData::createSVFaceNodeConnectivityPerOrientNoSymm()
{
  CFAUTOTRACE;

  // number of SV faces
  const CFuint nbrSVFaces = m_svFaceNodeConn.size();

  // number of possible orientations
  const CFuint nbrOrient = 9;

  // resize the variable
  m_svFaceNodeConnPerOrientationNoSymm.resize(nbrOrient);
  for (CFuint iOrient = 0; iOrient < nbrOrient; ++iOrient)
  {
    m_svFaceNodeConnPerOrientationNoSymm[iOrient].resize(2);
    for (CFuint iCell = 0; iCell < 2; ++iCell)
      m_svFaceNodeConnPerOrientationNoSymm[iOrient][iCell].resize(2);
  }

  // fill the variable
  CFuint iOrient = 0;
  for (CFuint iSVFaceL = 0; iSVFaceL < nbrSVFaces; ++iSVFaceL)
  {
    for (CFuint iSVFaceR = 0; iSVFaceR < nbrSVFaces; ++iSVFaceR, ++iOrient)
    {
      for (CFuint iNode = 0; iNode < 2; ++iNode)
      {
        m_svFaceNodeConnPerOrientationNoSymm[iOrient][LEFT ][iNode] = m_svFaceNodeConn[iSVFaceL][iNode  ];
        m_svFaceNodeConnPerOrientationNoSymm[iOrient][RIGHT][iNode] = m_svFaceNodeConn[iSVFaceR][1-iNode];
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////

void TriagSpectralFVElementData::computeFaceFractionsOfCVs()
{
  CFAUTOTRACE;

  // number of SV faces
  const CFuint nbrSVFaces = 3;

  // number of CVs at a SV face
  const CFuint nbrCVsAtSVFace = m_polyOrder + 1;

  // resize the vector containing the CV face fractions
  m_faceFracCV.resize(nbrSVFaces);
  for (CFuint iFace = 0; iFace < nbrSVFaces; ++iFace)
  {
    m_faceFracCV[iFace].resize(nbrCVsAtSVFace);
  }

  // fill in CV fractions first SV face (first nbrCVsAtSVFace faces in m_localFaceNodeConn make up one SV face)
  for (CFuint iCV = 0; iCV < nbrCVsAtSVFace; ++iCV)
  {
    const CFuint node1 = m_localFaceNodeConn[iCV][0];
    const CFuint node2 = m_localFaceNodeConn[iCV][1];
    m_faceFracCV[0][iCV] = m_localNodeCoord[node1][KSI] - m_localNodeCoord[node2][KSI];
  }

  // fill in CV fractions of other SV faces
  for (CFuint iFace = 1; iFace < nbrSVFaces; ++iFace)
  {
    for (CFuint iCV = 0; iCV < nbrCVsAtSVFace; ++iCV)
    {
      m_faceFracCV[iFace][iCV] = m_faceFracCV[0][iCV];
    }
  }
}

//////////////////////////////////////////////////////////////////////

void TriagSpectralFVElementData::createPolyExponents()
{
  CFAUTOTRACE;

  // helper variable
  const CFuint polyOrderP1 = m_polyOrder + 1;

  // number of polynomial terms
  const CFuint nbrPolyTerms = (m_polyOrder + 1)*(m_polyOrder + 2)/2;

  // resize the variable
  m_polyExponents.resize(nbrPolyTerms);
  for (CFuint iTerm = 0; iTerm < nbrPolyTerms; ++iTerm)
  {
    m_polyExponents[iTerm].resize(2);
  }

  // define exponents
  CFuint iTerm = 0;
  for (CFuint iP = 0; iP < polyOrderP1; ++iP)
  {
    for (CFuint iY = 0; iY < iP+1; ++iY, ++iTerm)
    {
      m_polyExponents[iTerm][KSI] = iP-iY;
      m_polyExponents[iTerm][ETA] = iY;
    }
  }
}

//////////////////////////////////////////////////////////////////////

void TriagSpectralFVElementData::computePolyCoefs()
{
  CFAUTOTRACE;

  // set the dimensionality and order of the simplex integrator
  m_sIntegrator.setDimensionality(m_dimensionality);
  m_sIntegrator.setIntegratorOrder(m_polyOrder);

  // number of control volumes in a SV
  const CFuint nbrOfCVs = m_localCVNodeConn.size();

  // matrix for linear system
  RealMatrix polyCoefSystem(nbrOfCVs,nbrOfCVs);

  // fill the linear system matrix
  for (CFuint iCV = 0; iCV < nbrOfCVs; ++iCV)
  {
    // number of control volume nodes
    const CFuint nbrCVNodes = m_localCVNodeConn[iCV].size();

    // number of triangles in polygon
    const CFuint nbrTriangles = nbrCVNodes - 2;

    // get node coordinates
    vector< RealVector > cvNodeCoord(nbrCVNodes);
    for (CFuint iNode = 0; iNode < nbrCVNodes; ++iNode)
    {
      cvNodeCoord[iNode].resize(m_dimensionality);
      const CFuint nodeID = m_localCVNodeConn[iCV][iNode];
      cvNodeCoord[iNode] = m_localNodeCoord[nodeID];
    }

    // get quadrature nodes and wheights
    vector< RealVector > qNodeCoord;
    vector< CFreal > qWheights;
    for (CFuint iTriangle = 0; iTriangle < nbrTriangles; ++iTriangle)
    {
      // get node coordinates of triangle
      vector< RealVector > triagNodeCoord(3);
      for (CFuint iNode = 0; iNode < 3; ++iNode)
      {
        triagNodeCoord[iNode].resize(m_dimensionality);
      }
      triagNodeCoord[0] = cvNodeCoord[0];
      triagNodeCoord[1] = cvNodeCoord[iTriangle+1];
      triagNodeCoord[2] = cvNodeCoord[iTriangle+2];

      // get triangle quadrature nodes and wheights
      vector< RealVector > qNodeCoordTriangle = m_sIntegrator.getQuadPntsCoords  (triagNodeCoord);
      vector< CFreal >     qWheightsTriangle  = m_sIntegrator.getQuadPntsWheights(triagNodeCoord);

      // add triangle quadrature nodes and wheights to global list
      qNodeCoord.insert(qNodeCoord.end(),qNodeCoordTriangle.begin(),qNodeCoordTriangle.end());
      qWheights .insert(qWheights .end(),qWheightsTriangle .begin(),qWheightsTriangle .end());
    }

    // compute left hand side of linear system
    const CFuint nbrQNodes = qWheights.size();
    for (CFuint iTerm = 0; iTerm < nbrOfCVs; ++iTerm)
    {
      polyCoefSystem(iCV,iTerm) = 0;

      for (CFuint iQNode = 0; iQNode < nbrQNodes; ++iQNode)
      {
        polyCoefSystem(iCV,iTerm) += qWheights[iQNode]
                                        *pow(qNodeCoord[iQNode][KSI],m_polyExponents[iTerm][KSI])
                                        *pow(qNodeCoord[iQNode][ETA],m_polyExponents[iTerm][ETA]);
      }
    }
  }

  // invert the linear system matrix
  RealMatrix invLinSysMatrix(nbrOfCVs,nbrOfCVs);
  InvertMatrix(polyCoefSystem,invLinSysMatrix);

  // right hand side of linear system
  RealMatrix rhsLinSys(nbrOfCVs,nbrOfCVs,0.0);
  for (CFuint iCV = 0; iCV < nbrOfCVs; ++iCV)
  {
    rhsLinSys(iCV,iCV) = 0.5*m_volFracCV[iCV];
  }

  // multiply imverted linear system matrix with rhs
  RealMatrix polyCoefMatrix(nbrOfCVs,nbrOfCVs);
  polyCoefMatrix = invLinSysMatrix*rhsLinSys;

  // store polynomial coefficients in m_polyCoefSFV
  m_polyCoefSFV.resize(nbrOfCVs);
  for (CFuint iPoly = 0; iPoly < nbrOfCVs; ++iPoly)
  {
//    CF_DEBUG_OBJ(iPoly);
    m_polyCoefSFV[iPoly].resize(nbrOfCVs);
    for (CFuint iTerm = 0; iTerm < nbrOfCVs; ++iTerm)
    {
      m_polyCoefSFV[iPoly][iTerm] = polyCoefMatrix(iTerm,iPoly);
//      cout << m_polyCoefSFV[iPoly][iTerm] << "*ksi" << m_polyExponents[iTerm][KSI] << "eta" << m_polyExponents[iTerm][ETA] << " + ";
    }
//    cout << endl;
  }
}

//////////////////////////////////////////////////////////////////////

void TriagSpectralFVElementData::createFaceFluxPntsConn()
{
 CFAUTOTRACE;

  // number of flux points at a face
  const CFuint nbrFaceFlxPnts = m_polyOrder +2;

  // resize m_faceFlxPntsConn (nbr of faces)
  m_faceFlxPntsConn.resize(3);

  // create connectivity
  // face 0
  m_faceFlxPntsConn[0].resize(nbrFaceFlxPnts);
  CFuint flxPntID = 0;
  for (CFuint iFlx = 0; iFlx < nbrFaceFlxPnts; ++iFlx)
  {
    m_faceFlxPntsConn[0][iFlx] = flxPntID;
    flxPntID += iFlx+1;
  }

  // face 1
  m_faceFlxPntsConn[1].resize(nbrFaceFlxPnts);
  flxPntID = m_faceFlxPntsConn[0][nbrFaceFlxPnts-1];
  for (CFuint iFlx = 0; iFlx < nbrFaceFlxPnts; ++iFlx, ++flxPntID)
  {
    m_faceFlxPntsConn[1][iFlx] = flxPntID;
  }

  // face 2
  m_faceFlxPntsConn[2].resize(nbrFaceFlxPnts);
  flxPntID = m_faceFlxPntsConn[1][nbrFaceFlxPnts-1];
  for (CFuint iFlx = 0; iFlx < nbrFaceFlxPnts; ++iFlx)
  {
    m_faceFlxPntsConn[2][iFlx] = flxPntID;
    flxPntID -= nbrFaceFlxPnts-iFlx;
  }

/*  for (CFuint iFace = 0; iFace < 3; ++iFace)
  {
    CF_DEBUG_OBJ(iFace);
    for (CFuint iFlx = 0; iFlx < nbrFaceFlxPnts; ++iFlx)
    {
      CF_DEBUG_OBJ(m_faceFlxPntsConn[iFace][iFlx]);
    }
  }*/
}

//////////////////////////////////////////////////////////////////////

void TriagSpectralFVElementData::createFaceFluxPolyNodeWheightCoord()
{
  CFAUTOTRACE;

  // get the number of face flux points
  const CFuint nbrFaceFlxPnts = getNbrOfFaceFlxPnts();

  // resize m_faceFluxPolyNodeWheightCoord
  m_faceFluxPolyNodeWheightCoord.resize(nbrFaceFlxPnts);

  // dimensionality -1
  const CFuint dimM1 = m_dimensionality - 1;

  // use the first SV face (for which the last local coordinate should be zero)
  for (CFuint iFlx = 0; iFlx < nbrFaceFlxPnts; ++iFlx)
  {
    // get node ID
    const CFuint nodeID = m_faceFlxPntsConn[0][iFlx];

    // resize
    m_faceFluxPolyNodeWheightCoord[iFlx].resize(m_dimensionality);

    // compute wheight coordinates
    m_faceFluxPolyNodeWheightCoord[iFlx][KSI] = 1.0 - m_fluxPolyNodeCoord[nodeID].sum();
    for (CFuint iCoor = 0; iCoor < dimM1; ++iCoor)
    {
      m_faceFluxPolyNodeWheightCoord[iFlx][iCoor+1] = m_fluxPolyNodeCoord[nodeID][iCoor];
    }
  }
}

//////////////////////////////////////////////////////////////////////

void TriagSpectralFVElementData::createFluxPolyExponents()
{
  CFAUTOTRACE;

  // helper variable
  const CFuint polyOrderP2 = m_polyOrder + 2;

  // number of polynomial terms
  const CFuint nbrPolyTerms = (m_polyOrder + 2)*(m_polyOrder + 3)/2;

  // resize the variable
  m_fluxPolyExponents.resize(nbrPolyTerms);
  for (CFuint iTerm = 0; iTerm < nbrPolyTerms; ++iTerm)
  {
    m_fluxPolyExponents[iTerm].resize(2);
  }

  // define exponents
  CFuint iTerm = 0;
  for (CFuint iP = 0; iP < polyOrderP2; ++iP)
  {
    for (CFuint iY = 0; iY < iP+1; ++iY, ++iTerm)
    {
      m_fluxPolyExponents[iTerm][KSI] = iP-iY;
      m_fluxPolyExponents[iTerm][ETA] = iY;
    }
  }
}

//////////////////////////////////////////////////////////////////////

void TriagSpectralFVElementData::setInterpolationNodeSet(const CFPolyOrder::Type order,
                                                         vector< RealVector >& nodalSet)
{
  // number of nodes
  const CFuint nbrNodes = (order+1)*(order+2)/2;

  // return variable
  nodalSet.resize(nbrNodes);
  for (CFuint iNode = 0; iNode < nbrNodes; ++iNode)
  {
    nodalSet[iNode].resize(2);
  }

  // fill the vector containing the flux polynomial node coordinates
  // Legendre-Gauss-Lobatto nodes from Hesthaven, FROM ELECTROSTATICS TO ALMOST OPTIMAL NODAL SETS
  // FOR POLYNOMIAL INTERPOLATION IN A SIMPLEX. SIAM J. Numer. Anal. (1998); 35(2):655-676.
  switch (order)
  {
    case CFPolyOrder::ORDER0:
    {
      nodalSet[0][KSI] = 1.0/3.0;
      nodalSet[0][ETA] = 1.0/3.0;
    } break;
    case CFPolyOrder::ORDER1:
    {
      nodalSet[0][KSI] = 0.0;
      nodalSet[0][ETA] = 0.0;

      nodalSet[1][KSI] = 1.0;
      nodalSet[1][ETA] = 0.0;

      nodalSet[2][KSI] = 0.0;
      nodalSet[2][ETA] = 1.0;
    } break;
    case CFPolyOrder::ORDER2:
    {
      nodalSet[0][KSI] = 0.0;
      nodalSet[0][ETA] = 0.0;

      nodalSet[1][KSI] = 0.5;
      nodalSet[1][ETA] = 0.0;

      nodalSet[2][KSI] = 0.0;
      nodalSet[2][ETA] = 0.5;

      nodalSet[3][KSI] = 1.0;
      nodalSet[3][ETA] = 0.0;

      nodalSet[4][KSI] = 0.5;
      nodalSet[4][ETA] = 0.5;

      nodalSet[5][KSI] = 0.0;
      nodalSet[5][ETA] = 1.0;
    } break;
    case CFPolyOrder::ORDER3:
    {
      nodalSet[0][KSI] = 0.0;
      nodalSet[0][ETA] = 0.0;

      nodalSet[1][KSI] = 0.25;
      nodalSet[1][ETA] = 0.0;

      nodalSet[2][KSI] = 0.0;
      nodalSet[2][ETA] = 0.25;

      nodalSet[3][KSI] = 0.75;
      nodalSet[3][ETA] = 0.0;

      nodalSet[4][KSI] = 1.0/3.0;
      nodalSet[4][ETA] = 1.0/3.0;

      nodalSet[5][KSI] = 0.0;
      nodalSet[5][ETA] = 0.75;

      nodalSet[6][KSI] = 1.0;
      nodalSet[6][ETA] = 0.0;

      nodalSet[7][KSI] = 0.75;
      nodalSet[7][ETA] = 0.25;

      nodalSet[8][KSI] = 0.25;
      nodalSet[8][ETA] = 0.75;

      nodalSet[9][KSI] = 0.0;
      nodalSet[9][ETA] = 1.0;
    } break;
    case CFPolyOrder::ORDER4:
    {
      nodalSet[0][KSI] = 0.0;
      nodalSet[0][ETA] = 0.0;

      nodalSet[1][KSI] = 0.5*(1-cos(MathTools::MathConsts::CFrealPi()/4.0));
      nodalSet[1][ETA] = 0.0;

      nodalSet[2][KSI] = 0.0;
      nodalSet[2][ETA] = 0.5*(1-cos(MathTools::MathConsts::CFrealPi()/4.0));

      nodalSet[3][KSI] = 0.5;
      nodalSet[3][ETA] = 0.0;

      nodalSet[4][KSI] = 0.2371200168;
      nodalSet[4][ETA] = 0.2371200168;

      nodalSet[5][KSI] = 0.0;
      nodalSet[5][ETA] = 0.5;

      nodalSet[6][KSI] = 0.5*(1+cos(MathTools::MathConsts::CFrealPi()/4.0));
      nodalSet[6][ETA] = 0.0;

      nodalSet[7][KSI] = 0.5257599664;
      nodalSet[7][ETA] = 0.2371200168;

      nodalSet[8][KSI] = 0.2371200168;
      nodalSet[8][ETA] = 0.5257599664;

      nodalSet[9][KSI] = 0.0;
      nodalSet[9][ETA] = 0.5*(1+cos(MathTools::MathConsts::CFrealPi()/4.0));

      nodalSet[10][KSI] = 1.0;
      nodalSet[10][ETA] = 0.0;

      nodalSet[11][KSI] = 0.5*(1+cos(MathTools::MathConsts::CFrealPi()/4.0));
      nodalSet[11][ETA] = 0.5*(1-cos(MathTools::MathConsts::CFrealPi()/4.0));

      nodalSet[12][KSI] = 0.5;
      nodalSet[12][ETA] = 0.5;

      nodalSet[13][KSI] = 0.5*(1-cos(MathTools::MathConsts::CFrealPi()/4.0));
      nodalSet[13][ETA] = 0.5*(1+cos(MathTools::MathConsts::CFrealPi()/4.0));

      nodalSet[14][KSI] = 0.0;
      nodalSet[14][ETA] = 1.0;
    } break;
    default:
    {
      throw Common::NotImplementedException (FromHere(),"Higher-order nodal set not defined!");
    }
  }
}

//////////////////////////////////////////////////////////////////////

void TriagSpectralFVElementData::setCFLConvDiffRatio()
{
  CFAUTOTRACE;

  switch(m_polyOrder)
  {
    case CFPolyOrder::ORDER0:
    {
      m_cflConvDiffRatio = 4.0; // check this!
    } break;
    case CFPolyOrder::ORDER1:
    {
      m_cflConvDiffRatio = 6.5; // check this!
    } break;
    case CFPolyOrder::ORDER2:
    {
      m_cflConvDiffRatio = 17.0; // check this!
    } break;
    case CFPolyOrder::ORDER3:
    {
      m_cflConvDiffRatio = 50.0; // check this!
    } break;
    default:
    {
      throw Common::NotImplementedException (FromHere(),"Higher-order triangular SV element not defined!");
    }
  }
}

//////////////////////////////////////////////////////////////////////

void TriagSpectralFVElementData::createFaceOutputPntCellMappedCoords()
{
  // number of points on a face
  CFuint nbrFacePnts;
  CFreal dKsi;
  if (m_polyOrder == CFPolyOrder::ORDER0)
  {
    nbrFacePnts = 2;
    dKsi = 2.0;
  }
  else
  {
    nbrFacePnts = m_polyOrder + 1;
    dKsi = 2.0/m_polyOrder;
  }

  // face mapped coordinates of uniform distribution of points
  CFreal ksi = -1.0;
  m_faceOutputPntFaceMappedCoords.resize(0);
  for (CFuint iPnt = 0; iPnt < nbrFacePnts; ++iPnt, ksi += dKsi)
  {
    RealVector mapCoord(1);
    mapCoord[KSI] = ksi;
    m_faceOutputPntFaceMappedCoords.push_back(mapCoord);
  }

  // compute cell mapped coordinates for distribution on each face
  const CFuint nbrCellFaces = getNbrSVFaces();
  m_faceOutputPntCellMappedCoords.resize(nbrCellFaces);
  for (CFuint iFace = 0; iFace < nbrCellFaces; ++iFace)
  {
    // current face node coordinates
    const vector<RealVector>& faceNodeCoords = m_svFaceNodeCoords[iFace];
    m_faceOutputPntCellMappedCoords[iFace].resize(0);
    for (CFuint iPnt = 0; iPnt < nbrFacePnts; ++iPnt)
    {
      const CFreal fun0 = 0.5*(1.0-m_faceOutputPntFaceMappedCoords[iPnt][KSI]);
      const CFreal fun1 = 0.5*(1.0+m_faceOutputPntFaceMappedCoords[iPnt][KSI]);
      m_faceOutputPntCellMappedCoords[iFace].push_back(fun0*faceNodeCoords[0]+fun1*faceNodeCoords[1]);
    }
  }
}

//////////////////////////////////////////////////////////////////////

void TriagSpectralFVElementData::createFaceOutputPntConn()
{
  CFAUTOTRACE;

  m_faceOutputPntConn.resize(0);
  for (CFuint iCell = 0; iCell < (CFuint) m_polyOrder; ++iCell)
  {
    vector<CFuint> cellNode(2);
    cellNode[0] = iCell;
    cellNode[1] = iCell+1;
    m_faceOutputPntConn.push_back(cellNode);
  }
}

//////////////////////////////////////////////////////////////////////

  } // namespace SpectralFV

} // namespace COOLFluiD

