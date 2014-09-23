#include <fstream>


#include "Common/CFLog.hh"
#include "Environment/DirPaths.hh"
#include "Framework/BadFormatException.hh"
#include "Common/NotImplementedException.hh"
#include "MathTools/RealMatrix.hh"
#include "SpectralFV/TetraSpectralFVElementData.hh"

//////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

namespace COOLFluiD {

  namespace SpectralFV {

//////////////////////////////////////////////////////////////////////

TetraSpectralFVElementData::TetraSpectralFVElementData() :
  SpectralFVElementData()
{
  m_shape = CFGeoShape::TETRA;
  m_dimensionality = DIM_3D;
}

//////////////////////////////////////////////////////////////////////

TetraSpectralFVElementData::TetraSpectralFVElementData(CFPolyOrder::Type polyOrder)
{
  m_shape = CFGeoShape::TETRA;
  m_dimensionality = DIM_3D;
  m_polyOrder = polyOrder;

  resetSpectralFVElementData();
}

//////////////////////////////////////////////////////////////////////

TetraSpectralFVElementData::~TetraSpectralFVElementData()
{
}

//////////////////////////////////////////////////////////////////////

void TetraSpectralFVElementData::computeSVNodeLocalCoords()
{
  m_svNodeCoords.resize(4);

  // first node
  m_svNodeCoords[0].resize(3);
  m_svNodeCoords[0][KSI] = 0.0;
  m_svNodeCoords[0][ETA] = 0.0;
  m_svNodeCoords[0][ZTA] = 0.0;

  // second node
  m_svNodeCoords[1].resize(3);
  m_svNodeCoords[1][KSI] = 1.0;
  m_svNodeCoords[1][ETA] = 0.0;
  m_svNodeCoords[1][ZTA] = 0.0;

  // third node
  m_svNodeCoords[2].resize(3);
  m_svNodeCoords[2][KSI] = 0.0;
  m_svNodeCoords[2][ETA] = 1.0;
  m_svNodeCoords[2][ZTA] = 0.0;

  // fourth node
  m_svNodeCoords[3].resize(3);
  m_svNodeCoords[3][KSI] = 0.0;
  m_svNodeCoords[3][ETA] = 0.0;
  m_svNodeCoords[3][ZTA] = 1.0;
}

//////////////////////////////////////////////////////////////////////

void TetraSpectralFVElementData::createLocalNodeCoord()
{
  CFAUTOTRACE;

  // Path where to read files from
  std::string filename;
  const std::string append = Environment::DirPaths::getInstance().getWorkingDir().string() + "/";

//   // number of local nodes
//   const CFuint nbrLocalNodes = (m_polyOrder+1)*(m_polyOrder+1)*(m_polyOrder+1) +
//                                (m_polyOrder+1)*(m_polyOrder+1) + m_polyOrder+2;
//
//   // resize the vectors containing the node coordinates
//   m_localNodeCoord.resize(nbrLocalNodes);
//   for (CFuint iNode = 0; iNode < nbrLocalNodes; ++iNode)
//   {
//     m_localNodeCoord[iNode].resize(3);
//   }

  switch (m_polyOrder)
  {
    case CFPolyOrder::ORDER0:
    {
      // number of local nodes
      const CFuint nbrLocalNodes = 4;

      m_localNodeCoord.resize(nbrLocalNodes);
      for (CFuint iNode = 0; iNode < nbrLocalNodes; ++iNode)
      {
        m_localNodeCoord[iNode].resize(3);
      }

      // fill the matrix containing the node coordinates
      m_localNodeCoord[0][KSI] = 0.0;
      m_localNodeCoord[0][ETA] = 0.0;
      m_localNodeCoord[0][ZTA] = 0.0;

      m_localNodeCoord[1][KSI] = 1.0;
      m_localNodeCoord[1][ETA] = 0.0;
      m_localNodeCoord[1][ZTA] = 0.0;

      m_localNodeCoord[2][KSI] = 0.0;
      m_localNodeCoord[2][ETA] = 1.0;
      m_localNodeCoord[2][ZTA] = 0.0;

      m_localNodeCoord[3][KSI] = 0.0;
      m_localNodeCoord[3][ETA] = 0.0;
      m_localNodeCoord[3][ZTA] = 1.0;
    } break;
    case CFPolyOrder::ORDER1:
    {
      // number of local nodes
      const CFuint nbrLocalNodes = 15;

      m_localNodeCoord.resize(nbrLocalNodes);
      for (CFuint iNode = 0; iNode < nbrLocalNodes; ++iNode)
      {
        m_localNodeCoord[iNode].resize(3);
      }

      // fill the matrix containing the node coordinates
      m_localNodeCoord[ 0][KSI] = 0.0;
      m_localNodeCoord[ 0][ETA] = 0.0;
      m_localNodeCoord[ 0][ZTA] = 0.0;


      m_localNodeCoord[ 1][KSI] = 0.5;
      m_localNodeCoord[ 1][ETA] = 0.0;
      m_localNodeCoord[ 1][ZTA] = 0.0;

      m_localNodeCoord[ 2][KSI] = 1.0/3.0;
      m_localNodeCoord[ 2][ETA] = 1.0/3.0;
      m_localNodeCoord[ 2][ZTA] = 0.0;

      m_localNodeCoord[ 3][KSI] = 0.0;
      m_localNodeCoord[ 3][ETA] = 0.5;
      m_localNodeCoord[ 3][ZTA] = 0.0;

      m_localNodeCoord[ 4][KSI] = 1.0/3.0;
      m_localNodeCoord[ 4][ETA] = 0.0;
      m_localNodeCoord[ 4][ZTA] = 1.0/3.0;

      m_localNodeCoord[ 5][KSI] = 1.0/4.0;
      m_localNodeCoord[ 5][ETA] = 1.0/4.0;
      m_localNodeCoord[ 5][ZTA] = 1.0/4.0;

      m_localNodeCoord[ 6][KSI] = 0.0;
      m_localNodeCoord[ 6][ETA] = 1.0/3.0;
      m_localNodeCoord[ 6][ZTA] = 1.0/3.0;

      m_localNodeCoord[ 7][KSI] = 0.0;
      m_localNodeCoord[ 7][ETA] = 0.0;
      m_localNodeCoord[ 7][ZTA] = 0.5;


      m_localNodeCoord[ 8][KSI] = 1.0;
      m_localNodeCoord[ 8][ETA] = 0.0;
      m_localNodeCoord[ 8][ZTA] = 0.0;

      m_localNodeCoord[ 9][KSI] = 0.5;
      m_localNodeCoord[ 9][ETA] = 0.5;
      m_localNodeCoord[ 9][ZTA] = 0.0;

      m_localNodeCoord[10][KSI] = 0.0;
      m_localNodeCoord[10][ETA] = 1.0;
      m_localNodeCoord[10][ZTA] = 0.0;

      m_localNodeCoord[11][KSI] = 0.5;
      m_localNodeCoord[11][ETA] = 0.0;
      m_localNodeCoord[11][ZTA] = 0.5;

      m_localNodeCoord[12][KSI] = 1.0/3.0;
      m_localNodeCoord[12][ETA] = 1.0/3.0;
      m_localNodeCoord[12][ZTA] = 1.0/3.0;

      m_localNodeCoord[13][KSI] = 0.0;
      m_localNodeCoord[13][ETA] = 0.5;
      m_localNodeCoord[13][ZTA] = 0.5;

      m_localNodeCoord[14][KSI] = 0.0;
      m_localNodeCoord[14][ETA] = 0.0;
      m_localNodeCoord[14][ZTA] = 1.0;
    } break;
    case CFPolyOrder::ORDER2:
    {
      // partition parameters
      CFreal alpha, beta, gamma;

      // read defining parameters from file
      filename = append + "SV3TETRADEF.DAT";
      ifstream inputFile;
      inputFile.open(filename.c_str(), ios::in);

      if (inputFile.is_open())
      {
        inputFile >> alpha;
        inputFile >> beta;
        gamma = -1.0;
        if (!inputFile.eof())
        {
          inputFile >> gamma;
        }
      }
      else
      {
        alpha = 0.1025;
        beta  = 0.1875;
        gamma = -1.0;
/*        alpha = 0.1093621117;
        beta  = 0.1730022492;*/
/*        alpha = 1.0/4.0;
        beta  = 1.999/3.0;*/
/*        alpha = 1.0/4.0;
        beta  = 1.0/4.0;*/
/*        alpha = 1.0/4.0;
        beta  = 1.0/3.0;*/
      }
      inputFile.close();

      if (gamma == -1.0)
      {
        // determine gamma parameter as a function of the other two,
        // to ensure that the four nodes of the corner CV faces are coplanar
        gamma = 3.0*alpha*beta/(4.0*alpha - beta);
        if (gamma < 0.0)
        {
          CF_DEBUG_OBJ(gamma);
          gamma = 0.001;
        }
        if (gamma > 3.0/4.0)
        {
          CF_DEBUG_OBJ(gamma);
          gamma = 3.0/4.0;
        }
      }

      CF_DEBUG_OBJ(alpha);
      CF_DEBUG_OBJ(beta);
      CF_DEBUG_OBJ(gamma);

      // number of local nodes
      const CFuint nbrLocalNodes = 37;

      m_localNodeCoord.resize(nbrLocalNodes);
      for (CFuint iNode = 0; iNode < nbrLocalNodes; ++iNode)
      {
        m_localNodeCoord[iNode].resize(3);
      }

      // fill the matrix containing the node coordinates
      m_localNodeCoord[ 0][KSI] = 0.0;
      m_localNodeCoord[ 0][ETA] = 0.0;
      m_localNodeCoord[ 0][ZTA] = 0.0;


      m_localNodeCoord[ 1][KSI] = alpha;
      m_localNodeCoord[ 1][ETA] = 0.0;
      m_localNodeCoord[ 1][ZTA] = 0.0;

      m_localNodeCoord[ 2][KSI] = 0.5*beta;
      m_localNodeCoord[ 2][ETA] = 0.5*beta;
      m_localNodeCoord[ 2][ZTA] = 0.0;

      m_localNodeCoord[ 3][KSI] = 0.0;
      m_localNodeCoord[ 3][ETA] = alpha;
      m_localNodeCoord[ 3][ZTA] = 0.0;

      m_localNodeCoord[ 4][KSI] = 0.5*beta;
      m_localNodeCoord[ 4][ETA] = 0.0;
      m_localNodeCoord[ 4][ZTA] = 0.5*beta;

      m_localNodeCoord[ 5][KSI] = 1.0/3.0*gamma;
      m_localNodeCoord[ 5][ETA] = 1.0/3.0*gamma;
      m_localNodeCoord[ 5][ZTA] = 1.0/3.0*gamma;

      m_localNodeCoord[ 6][KSI] = 0.0;
      m_localNodeCoord[ 6][ETA] = 0.5*beta;
      m_localNodeCoord[ 6][ZTA] = 0.5*beta;

      m_localNodeCoord[ 7][KSI] = 0.0;
      m_localNodeCoord[ 7][ETA] = 0.0;
      m_localNodeCoord[ 7][ZTA] = alpha;


      m_localNodeCoord[ 8][KSI] = 1.0 - alpha;
      m_localNodeCoord[ 8][ETA] = 0.0;
      m_localNodeCoord[ 8][ZTA] = 0.0;

      m_localNodeCoord[ 9][KSI] = 1.0 - beta;
      m_localNodeCoord[ 9][ETA] = 0.5*beta;
      m_localNodeCoord[ 9][ZTA] = 0.0;

      m_localNodeCoord[10][KSI] = 1.0/3.0;
      m_localNodeCoord[10][ETA] = 1.0/3.0;
      m_localNodeCoord[10][ZTA] = 0.0;

      m_localNodeCoord[11][KSI] = 0.5*beta;
      m_localNodeCoord[11][ETA] = 1.0 - beta;
      m_localNodeCoord[11][ZTA] = 0.0;

      m_localNodeCoord[12][KSI] = 0.0;
      m_localNodeCoord[12][ETA] = 1.0 - alpha;
      m_localNodeCoord[12][ZTA] = 0.0;

      m_localNodeCoord[13][KSI] = 1.0 - beta;
      m_localNodeCoord[13][ETA] = 0.0;
      m_localNodeCoord[13][ZTA] = 0.5*beta;

      m_localNodeCoord[14][KSI] = 1.0 - gamma;
      m_localNodeCoord[14][ETA] = 1.0/3.0*gamma;
      m_localNodeCoord[14][ZTA] = 1.0/3.0*gamma;

      m_localNodeCoord[15][KSI] = 1.0/4.0;
      m_localNodeCoord[15][ETA] = 1.0/4.0;
      m_localNodeCoord[15][ZTA] = 1.0/4.0;

      m_localNodeCoord[16][KSI] = 1.0/3.0*gamma;
      m_localNodeCoord[16][ETA] = 1.0 - gamma;
      m_localNodeCoord[16][ZTA] = 1.0/3.0*gamma;

      m_localNodeCoord[17][KSI] = 0.0;
      m_localNodeCoord[17][ETA] = 1.0 - beta;
      m_localNodeCoord[17][ZTA] = 0.5*beta;

      m_localNodeCoord[18][KSI] = 1.0/3.0;
      m_localNodeCoord[18][ETA] = 0.0;
      m_localNodeCoord[18][ZTA] = 1.0/3.0;

      m_localNodeCoord[19][KSI] = 0.0;
      m_localNodeCoord[19][ETA] = 1.0/3.0;
      m_localNodeCoord[19][ZTA] = 1.0/3.0;

      m_localNodeCoord[20][KSI] = 0.5*beta;
      m_localNodeCoord[20][ETA] = 0.0;
      m_localNodeCoord[20][ZTA] = 1.0 - beta;

      m_localNodeCoord[21][KSI] = 1.0/3.0*gamma;
      m_localNodeCoord[21][ETA] = 1.0/3.0*gamma;
      m_localNodeCoord[21][ZTA] = 1.0 - gamma;

      m_localNodeCoord[22][KSI] = 0.0;
      m_localNodeCoord[22][ETA] = 0.5*beta;
      m_localNodeCoord[22][ZTA] = 1.0 - beta;

      m_localNodeCoord[23][KSI] = 0.0;
      m_localNodeCoord[23][ETA] = 0.0;
      m_localNodeCoord[23][ZTA] = 1.0 - alpha;


      m_localNodeCoord[24][KSI] = 1.0;
      m_localNodeCoord[24][ETA] = 0.0;
      m_localNodeCoord[24][ZTA] = 0.0;

      m_localNodeCoord[25][KSI] = 1.0 - alpha;
      m_localNodeCoord[25][ETA] = alpha;
      m_localNodeCoord[25][ZTA] = 0.0;

      m_localNodeCoord[26][KSI] = alpha;
      m_localNodeCoord[26][ETA] = 1.0 - alpha;
      m_localNodeCoord[26][ZTA] = 0.0;

      m_localNodeCoord[27][KSI] = 0.0;
      m_localNodeCoord[27][ETA] = 1.0;
      m_localNodeCoord[27][ZTA] = 0.0;

      m_localNodeCoord[28][KSI] = 1.0 - alpha;
      m_localNodeCoord[28][ETA] = 0.0;
      m_localNodeCoord[28][ZTA] = alpha;

      m_localNodeCoord[29][KSI] = 1.0 - beta;
      m_localNodeCoord[29][ETA] = 0.5*beta;
      m_localNodeCoord[29][ZTA] = 0.5*beta;

      m_localNodeCoord[30][KSI] = 1.0/3.0;
      m_localNodeCoord[30][ETA] = 1.0/3.0;
      m_localNodeCoord[30][ZTA] = 1.0/3.0;

      m_localNodeCoord[31][KSI] = 0.5*beta;
      m_localNodeCoord[31][ETA] = 1.0 - beta;
      m_localNodeCoord[31][ZTA] = 0.5*beta;

      m_localNodeCoord[32][KSI] = 0.0;
      m_localNodeCoord[32][ETA] = 1.0 - alpha;
      m_localNodeCoord[32][ZTA] = alpha;

      m_localNodeCoord[33][KSI] = alpha;
      m_localNodeCoord[33][ETA] = 0.0;
      m_localNodeCoord[33][ZTA] = 1.0 - alpha;

      m_localNodeCoord[34][KSI] = 0.5*beta;
      m_localNodeCoord[34][ETA] = 0.5*beta;
      m_localNodeCoord[34][ZTA] = 1.0 - beta;

      m_localNodeCoord[35][KSI] = 0.0;
      m_localNodeCoord[35][ETA] = alpha;
      m_localNodeCoord[35][ZTA] = 1.0 - alpha;

      m_localNodeCoord[36][KSI] = 0.0;
      m_localNodeCoord[36][ETA] = 0.0;
      m_localNodeCoord[36][ZTA] = 1.0;
    } break;
    case CFPolyOrder::ORDER3:
    {
/*      // partition parameters
      CFreal alpha, beta, gamma, delta;

      // read defining parameters from file
      filename = append + "SV4TETRADEF.DAT";
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
      }
      inputFile.close();*/
      throw Common::NotImplementedException (FromHere(),"4th order accurate tetrahedral SV cells have not yet been implemented!");
    } break;
    default:
    {
      throw Common::NotImplementedException (FromHere(),"Only solution orders up to 3 have been implemented for spectral FV!");
    }
  }
}

//////////////////////////////////////////////////////////////////////

void TetraSpectralFVElementData::createLocalFaceNodeConn()
{
  CFAUTOTRACE;


  switch (m_polyOrder)
  {
    case CFPolyOrder::ORDER0:
    {
      // number of local faces
      const CFuint nbrLocalFaces = 4;

      // create the local face-node connectivity
      m_localFaceNodeConn.resize(nbrLocalFaces);

      // index for faces
      CFuint faceIdx = 0;

      // index for nodes
      CFuint nodeIdx = 0;

      // EXTERNAL FACES
      // first SV face
      m_localFaceNodeConn[faceIdx].resize(3);
      m_localFaceNodeConn[faceIdx][nodeIdx] = 0; ++nodeIdx;
      m_localFaceNodeConn[faceIdx][nodeIdx] = 1; ++nodeIdx;
      m_localFaceNodeConn[faceIdx][nodeIdx] = 2;
      ++faceIdx; nodeIdx = 0;

      // second SV face
      m_localFaceNodeConn[faceIdx].resize(3);
      m_localFaceNodeConn[faceIdx][nodeIdx] = 0; ++nodeIdx;
      m_localFaceNodeConn[faceIdx][nodeIdx] = 3; ++nodeIdx;
      m_localFaceNodeConn[faceIdx][nodeIdx] = 1;
      ++faceIdx; nodeIdx = 0;

      // third SV face
      m_localFaceNodeConn[faceIdx].resize(3);
      m_localFaceNodeConn[faceIdx][nodeIdx] = 1; ++nodeIdx;
      m_localFaceNodeConn[faceIdx][nodeIdx] = 3; ++nodeIdx;
      m_localFaceNodeConn[faceIdx][nodeIdx] = 2;
      ++faceIdx; nodeIdx = 0;

      // fourth SV face
      m_localFaceNodeConn[faceIdx].resize(3);
      m_localFaceNodeConn[faceIdx][nodeIdx] = 0; ++nodeIdx;
      m_localFaceNodeConn[faceIdx][nodeIdx] = 2; ++nodeIdx;
      m_localFaceNodeConn[faceIdx][nodeIdx] = 3;

      // INTERNAL FACES
    } break;
    case CFPolyOrder::ORDER1:
    {
      // number of local faces
      const CFuint nbrLocalFaces = 18;

      // create the local face-node connectivity
      m_localFaceNodeConn.resize(nbrLocalFaces);

      // index for faces
      CFuint faceIdx = 0;

      // index for nodes
      CFuint nodeIdx = 0;

      // EXTERNAL FACES
      // first SV face
      m_localFaceNodeConn[faceIdx].resize(4);
      m_localFaceNodeConn[faceIdx][nodeIdx] = 2; ++nodeIdx;
      m_localFaceNodeConn[faceIdx][nodeIdx] = 3; ++nodeIdx;
      m_localFaceNodeConn[faceIdx][nodeIdx] = 0; ++nodeIdx;
      m_localFaceNodeConn[faceIdx][nodeIdx] = 1;
      ++faceIdx; nodeIdx = 0;

      m_localFaceNodeConn[faceIdx].resize(4);
      m_localFaceNodeConn[faceIdx][nodeIdx] = 2; ++nodeIdx;
      m_localFaceNodeConn[faceIdx][nodeIdx] = 1; ++nodeIdx;
      m_localFaceNodeConn[faceIdx][nodeIdx] = 8; ++nodeIdx;
      m_localFaceNodeConn[faceIdx][nodeIdx] = 9;
      ++faceIdx; nodeIdx = 0;

      m_localFaceNodeConn[faceIdx].resize(4);
      m_localFaceNodeConn[faceIdx][nodeIdx] = 2; ++nodeIdx;
      m_localFaceNodeConn[faceIdx][nodeIdx] = 9; ++nodeIdx;
      m_localFaceNodeConn[faceIdx][nodeIdx] = 10; ++nodeIdx;
      m_localFaceNodeConn[faceIdx][nodeIdx] = 3;
      ++faceIdx; nodeIdx = 0;

      // second SV face
      m_localFaceNodeConn[faceIdx].resize(4);
      m_localFaceNodeConn[faceIdx][nodeIdx] = 4; ++nodeIdx;
      m_localFaceNodeConn[faceIdx][nodeIdx] = 1; ++nodeIdx;
      m_localFaceNodeConn[faceIdx][nodeIdx] = 0; ++nodeIdx;
      m_localFaceNodeConn[faceIdx][nodeIdx] = 7;
      ++faceIdx; nodeIdx = 0;

      m_localFaceNodeConn[faceIdx].resize(4);
      m_localFaceNodeConn[faceIdx][nodeIdx] = 4; ++nodeIdx;
      m_localFaceNodeConn[faceIdx][nodeIdx] = 7; ++nodeIdx;
      m_localFaceNodeConn[faceIdx][nodeIdx] = 14; ++nodeIdx;
      m_localFaceNodeConn[faceIdx][nodeIdx] = 11;
      ++faceIdx; nodeIdx = 0;

      m_localFaceNodeConn[faceIdx].resize(4);
      m_localFaceNodeConn[faceIdx][nodeIdx] = 4; ++nodeIdx;
      m_localFaceNodeConn[faceIdx][nodeIdx] = 11; ++nodeIdx;
      m_localFaceNodeConn[faceIdx][nodeIdx] = 8; ++nodeIdx;
      m_localFaceNodeConn[faceIdx][nodeIdx] = 1;
      ++faceIdx; nodeIdx = 0;

      // third SV face
      m_localFaceNodeConn[faceIdx].resize(4);
      m_localFaceNodeConn[faceIdx][nodeIdx] = 12; ++nodeIdx;
      m_localFaceNodeConn[faceIdx][nodeIdx] = 9; ++nodeIdx;
      m_localFaceNodeConn[faceIdx][nodeIdx] = 8; ++nodeIdx;
      m_localFaceNodeConn[faceIdx][nodeIdx] = 11;
      ++faceIdx; nodeIdx = 0;

      m_localFaceNodeConn[faceIdx].resize(4);
      m_localFaceNodeConn[faceIdx][nodeIdx] = 12; ++nodeIdx;
      m_localFaceNodeConn[faceIdx][nodeIdx] = 11; ++nodeIdx;
      m_localFaceNodeConn[faceIdx][nodeIdx] = 14; ++nodeIdx;
      m_localFaceNodeConn[faceIdx][nodeIdx] = 13;
      ++faceIdx; nodeIdx = 0;

      m_localFaceNodeConn[faceIdx].resize(4);
      m_localFaceNodeConn[faceIdx][nodeIdx] = 12; ++nodeIdx;
      m_localFaceNodeConn[faceIdx][nodeIdx] = 13; ++nodeIdx;
      m_localFaceNodeConn[faceIdx][nodeIdx] = 10; ++nodeIdx;
      m_localFaceNodeConn[faceIdx][nodeIdx] = 9;
      ++faceIdx; nodeIdx = 0;

      // fourth SV face
      m_localFaceNodeConn[faceIdx].resize(4);
      m_localFaceNodeConn[faceIdx][nodeIdx] = 6; ++nodeIdx;
      m_localFaceNodeConn[faceIdx][nodeIdx] = 7; ++nodeIdx;
      m_localFaceNodeConn[faceIdx][nodeIdx] = 0; ++nodeIdx;
      m_localFaceNodeConn[faceIdx][nodeIdx] = 3;
      ++faceIdx; nodeIdx = 0;

      m_localFaceNodeConn[faceIdx].resize(4);
      m_localFaceNodeConn[faceIdx][nodeIdx] = 6; ++nodeIdx;
      m_localFaceNodeConn[faceIdx][nodeIdx] = 3; ++nodeIdx;
      m_localFaceNodeConn[faceIdx][nodeIdx] = 10; ++nodeIdx;
      m_localFaceNodeConn[faceIdx][nodeIdx] = 13;
      ++faceIdx; nodeIdx = 0;

      m_localFaceNodeConn[faceIdx].resize(4);
      m_localFaceNodeConn[faceIdx][nodeIdx] = 6; ++nodeIdx;
      m_localFaceNodeConn[faceIdx][nodeIdx] = 13; ++nodeIdx;
      m_localFaceNodeConn[faceIdx][nodeIdx] = 14; ++nodeIdx;
      m_localFaceNodeConn[faceIdx][nodeIdx] = 7;
      ++faceIdx; nodeIdx = 0;

      // INTERNAL FACES
      m_localFaceNodeConn[faceIdx].resize(4);
      m_localFaceNodeConn[faceIdx][nodeIdx] = 5; ++nodeIdx;
      m_localFaceNodeConn[faceIdx][nodeIdx] = 4; ++nodeIdx;
      m_localFaceNodeConn[faceIdx][nodeIdx] = 1; ++nodeIdx;
      m_localFaceNodeConn[faceIdx][nodeIdx] = 2;
      ++faceIdx; nodeIdx = 0;

      m_localFaceNodeConn[faceIdx].resize(4);
      m_localFaceNodeConn[faceIdx][nodeIdx] = 5; ++nodeIdx;
      m_localFaceNodeConn[faceIdx][nodeIdx] = 2; ++nodeIdx;
      m_localFaceNodeConn[faceIdx][nodeIdx] = 9; ++nodeIdx;
      m_localFaceNodeConn[faceIdx][nodeIdx] = 12;
      ++faceIdx; nodeIdx = 0;

      m_localFaceNodeConn[faceIdx].resize(4);
      m_localFaceNodeConn[faceIdx][nodeIdx] = 5; ++nodeIdx;
      m_localFaceNodeConn[faceIdx][nodeIdx] = 2; ++nodeIdx;
      m_localFaceNodeConn[faceIdx][nodeIdx] = 3; ++nodeIdx;
      m_localFaceNodeConn[faceIdx][nodeIdx] = 6;
      ++faceIdx; nodeIdx = 0;

      m_localFaceNodeConn[faceIdx].resize(4);
      m_localFaceNodeConn[faceIdx][nodeIdx] = 5; ++nodeIdx;
      m_localFaceNodeConn[faceIdx][nodeIdx] = 4; ++nodeIdx;
      m_localFaceNodeConn[faceIdx][nodeIdx] = 11; ++nodeIdx;
      m_localFaceNodeConn[faceIdx][nodeIdx] = 12;
      ++faceIdx; nodeIdx = 0;

      m_localFaceNodeConn[faceIdx].resize(4);
      m_localFaceNodeConn[faceIdx][nodeIdx] = 5; ++nodeIdx;
      m_localFaceNodeConn[faceIdx][nodeIdx] = 12; ++nodeIdx;
      m_localFaceNodeConn[faceIdx][nodeIdx] = 13; ++nodeIdx;
      m_localFaceNodeConn[faceIdx][nodeIdx] = 6;
      ++faceIdx; nodeIdx = 0;

      m_localFaceNodeConn[faceIdx].resize(4);
      m_localFaceNodeConn[faceIdx][nodeIdx] = 5; ++nodeIdx;
      m_localFaceNodeConn[faceIdx][nodeIdx] = 6; ++nodeIdx;
      m_localFaceNodeConn[faceIdx][nodeIdx] = 7; ++nodeIdx;
      m_localFaceNodeConn[faceIdx][nodeIdx] = 4;
    } break;
    case CFPolyOrder::ORDER2:
    {
      // number of local faces
//       const CFuint nbrLocalFaces = 60;// for first and second partition family with triangular faces at corner CVs
      const CFuint nbrLocalFaces = 48;// for partition family with quadrilateral faces at corner CVs

      // create the local face-node connectivity
      m_localFaceNodeConn.resize(nbrLocalFaces);

      // index for faces
      CFuint faceIdx = 0;

      // index for nodes
      CFuint nodeIdx = 0;

      // EXTERNAL FACES
      // first SV face
      m_localFaceNodeConn[faceIdx].resize(4);
      m_localFaceNodeConn[faceIdx][nodeIdx] = 0; ++nodeIdx;
      m_localFaceNodeConn[faceIdx][nodeIdx] = 1; ++nodeIdx;
      m_localFaceNodeConn[faceIdx][nodeIdx] = 2; ++nodeIdx;
      m_localFaceNodeConn[faceIdx][nodeIdx] = 3;
      ++faceIdx; nodeIdx = 0;

      m_localFaceNodeConn[faceIdx].resize(5);
      m_localFaceNodeConn[faceIdx][nodeIdx] = 10; ++nodeIdx;
      m_localFaceNodeConn[faceIdx][nodeIdx] = 2; ++nodeIdx;
      m_localFaceNodeConn[faceIdx][nodeIdx] = 1; ++nodeIdx;
      m_localFaceNodeConn[faceIdx][nodeIdx] = 8; ++nodeIdx;
      m_localFaceNodeConn[faceIdx][nodeIdx] = 9;
      ++faceIdx; nodeIdx = 0;

      m_localFaceNodeConn[faceIdx].resize(5);
      m_localFaceNodeConn[faceIdx][nodeIdx] = 10; ++nodeIdx;
      m_localFaceNodeConn[faceIdx][nodeIdx] = 11; ++nodeIdx;
      m_localFaceNodeConn[faceIdx][nodeIdx] = 12; ++nodeIdx;
      m_localFaceNodeConn[faceIdx][nodeIdx] = 3; ++nodeIdx;
      m_localFaceNodeConn[faceIdx][nodeIdx] = 2;
      ++faceIdx; nodeIdx = 0;

      m_localFaceNodeConn[faceIdx].resize(4);
      m_localFaceNodeConn[faceIdx][nodeIdx] = 24; ++nodeIdx;
      m_localFaceNodeConn[faceIdx][nodeIdx] = 25; ++nodeIdx;
      m_localFaceNodeConn[faceIdx][nodeIdx] = 9; ++nodeIdx;
      m_localFaceNodeConn[faceIdx][nodeIdx] = 8;
      ++faceIdx; nodeIdx = 0;

      m_localFaceNodeConn[faceIdx].resize(5);
      m_localFaceNodeConn[faceIdx][nodeIdx] = 10; ++nodeIdx;
      m_localFaceNodeConn[faceIdx][nodeIdx] = 9; ++nodeIdx;
      m_localFaceNodeConn[faceIdx][nodeIdx] = 25; ++nodeIdx;
      m_localFaceNodeConn[faceIdx][nodeIdx] = 26; ++nodeIdx;
      m_localFaceNodeConn[faceIdx][nodeIdx] = 11;
      ++faceIdx; nodeIdx = 0;

      m_localFaceNodeConn[faceIdx].resize(4);
      m_localFaceNodeConn[faceIdx][nodeIdx] = 27; ++nodeIdx;
      m_localFaceNodeConn[faceIdx][nodeIdx] = 12; ++nodeIdx;
      m_localFaceNodeConn[faceIdx][nodeIdx] = 11; ++nodeIdx;
      m_localFaceNodeConn[faceIdx][nodeIdx] = 26;
      ++faceIdx; nodeIdx = 0;

      // second SV face
      m_localFaceNodeConn[faceIdx].resize(4);
      m_localFaceNodeConn[faceIdx][nodeIdx] = 0; ++nodeIdx;
      m_localFaceNodeConn[faceIdx][nodeIdx] = 7; ++nodeIdx;
      m_localFaceNodeConn[faceIdx][nodeIdx] = 4; ++nodeIdx;
      m_localFaceNodeConn[faceIdx][nodeIdx] = 1;
      ++faceIdx; nodeIdx = 0;

      m_localFaceNodeConn[faceIdx].resize(5);
      m_localFaceNodeConn[faceIdx][nodeIdx] = 18; ++nodeIdx;
      m_localFaceNodeConn[faceIdx][nodeIdx] = 4; ++nodeIdx;
      m_localFaceNodeConn[faceIdx][nodeIdx] = 7; ++nodeIdx;
      m_localFaceNodeConn[faceIdx][nodeIdx] = 23; ++nodeIdx;
      m_localFaceNodeConn[faceIdx][nodeIdx] = 20;
      ++faceIdx; nodeIdx = 0;

      m_localFaceNodeConn[faceIdx].resize(5);
      m_localFaceNodeConn[faceIdx][nodeIdx] = 18; ++nodeIdx;
      m_localFaceNodeConn[faceIdx][nodeIdx] = 13; ++nodeIdx;
      m_localFaceNodeConn[faceIdx][nodeIdx] = 8; ++nodeIdx;
      m_localFaceNodeConn[faceIdx][nodeIdx] = 1; ++nodeIdx;
      m_localFaceNodeConn[faceIdx][nodeIdx] = 4;
      ++faceIdx; nodeIdx = 0;

      m_localFaceNodeConn[faceIdx].resize(4);
      m_localFaceNodeConn[faceIdx][nodeIdx] = 36; ++nodeIdx;
      m_localFaceNodeConn[faceIdx][nodeIdx] = 33; ++nodeIdx;
      m_localFaceNodeConn[faceIdx][nodeIdx] = 20; ++nodeIdx;
      m_localFaceNodeConn[faceIdx][nodeIdx] = 23;
      ++faceIdx; nodeIdx = 0;

      m_localFaceNodeConn[faceIdx].resize(5);
      m_localFaceNodeConn[faceIdx][nodeIdx] = 18; ++nodeIdx;
      m_localFaceNodeConn[faceIdx][nodeIdx] = 20; ++nodeIdx;
      m_localFaceNodeConn[faceIdx][nodeIdx] = 33; ++nodeIdx;
      m_localFaceNodeConn[faceIdx][nodeIdx] = 28; ++nodeIdx;
      m_localFaceNodeConn[faceIdx][nodeIdx] = 13;
      ++faceIdx; nodeIdx = 0;

      m_localFaceNodeConn[faceIdx].resize(4);
      m_localFaceNodeConn[faceIdx][nodeIdx] = 24; ++nodeIdx;
      m_localFaceNodeConn[faceIdx][nodeIdx] = 8; ++nodeIdx;
      m_localFaceNodeConn[faceIdx][nodeIdx] = 13; ++nodeIdx;
      m_localFaceNodeConn[faceIdx][nodeIdx] = 28;
      ++faceIdx; nodeIdx = 0;

      // third SV face
      m_localFaceNodeConn[faceIdx].resize(4);
      m_localFaceNodeConn[faceIdx][nodeIdx] = 24; ++nodeIdx;
      m_localFaceNodeConn[faceIdx][nodeIdx] = 28; ++nodeIdx;
      m_localFaceNodeConn[faceIdx][nodeIdx] = 29; ++nodeIdx;
      m_localFaceNodeConn[faceIdx][nodeIdx] = 25;
      ++faceIdx; nodeIdx = 0;

      m_localFaceNodeConn[faceIdx].resize(5);
      m_localFaceNodeConn[faceIdx][nodeIdx] = 30; ++nodeIdx;
      m_localFaceNodeConn[faceIdx][nodeIdx] = 29; ++nodeIdx;
      m_localFaceNodeConn[faceIdx][nodeIdx] = 28; ++nodeIdx;
      m_localFaceNodeConn[faceIdx][nodeIdx] = 33; ++nodeIdx;
      m_localFaceNodeConn[faceIdx][nodeIdx] = 34;
      ++faceIdx; nodeIdx = 0;

      m_localFaceNodeConn[faceIdx].resize(5);
      m_localFaceNodeConn[faceIdx][nodeIdx] = 30; ++nodeIdx;
      m_localFaceNodeConn[faceIdx][nodeIdx] = 31; ++nodeIdx;
      m_localFaceNodeConn[faceIdx][nodeIdx] = 26; ++nodeIdx;
      m_localFaceNodeConn[faceIdx][nodeIdx] = 25; ++nodeIdx;
      m_localFaceNodeConn[faceIdx][nodeIdx] = 29;
      ++faceIdx; nodeIdx = 0;

      m_localFaceNodeConn[faceIdx].resize(4);
      m_localFaceNodeConn[faceIdx][nodeIdx] = 36; ++nodeIdx;
      m_localFaceNodeConn[faceIdx][nodeIdx] = 35; ++nodeIdx;
      m_localFaceNodeConn[faceIdx][nodeIdx] = 34; ++nodeIdx;
      m_localFaceNodeConn[faceIdx][nodeIdx] = 33;
      ++faceIdx; nodeIdx = 0;

      m_localFaceNodeConn[faceIdx].resize(5);
      m_localFaceNodeConn[faceIdx][nodeIdx] = 30; ++nodeIdx;
      m_localFaceNodeConn[faceIdx][nodeIdx] = 34; ++nodeIdx;
      m_localFaceNodeConn[faceIdx][nodeIdx] = 35; ++nodeIdx;
      m_localFaceNodeConn[faceIdx][nodeIdx] = 32; ++nodeIdx;
      m_localFaceNodeConn[faceIdx][nodeIdx] = 31;
      ++faceIdx; nodeIdx = 0;

      m_localFaceNodeConn[faceIdx].resize(4);
      m_localFaceNodeConn[faceIdx][nodeIdx] = 27; ++nodeIdx;
      m_localFaceNodeConn[faceIdx][nodeIdx] = 26; ++nodeIdx;
      m_localFaceNodeConn[faceIdx][nodeIdx] = 31; ++nodeIdx;
      m_localFaceNodeConn[faceIdx][nodeIdx] = 32;
      ++faceIdx; nodeIdx = 0;

      // fourth SV face
      m_localFaceNodeConn[faceIdx].resize(4);
      m_localFaceNodeConn[faceIdx][nodeIdx] = 0; ++nodeIdx;
      m_localFaceNodeConn[faceIdx][nodeIdx] = 3; ++nodeIdx;
      m_localFaceNodeConn[faceIdx][nodeIdx] = 6; ++nodeIdx;
      m_localFaceNodeConn[faceIdx][nodeIdx] = 7;
      ++faceIdx; nodeIdx = 0;

      m_localFaceNodeConn[faceIdx].resize(5);
      m_localFaceNodeConn[faceIdx][nodeIdx] = 19; ++nodeIdx;
      m_localFaceNodeConn[faceIdx][nodeIdx] = 6; ++nodeIdx;
      m_localFaceNodeConn[faceIdx][nodeIdx] = 3; ++nodeIdx;
      m_localFaceNodeConn[faceIdx][nodeIdx] = 12; ++nodeIdx;
      m_localFaceNodeConn[faceIdx][nodeIdx] = 17;
      ++faceIdx; nodeIdx = 0;

      m_localFaceNodeConn[faceIdx].resize(5);
      m_localFaceNodeConn[faceIdx][nodeIdx] = 19; ++nodeIdx;
      m_localFaceNodeConn[faceIdx][nodeIdx] = 22; ++nodeIdx;
      m_localFaceNodeConn[faceIdx][nodeIdx] = 23; ++nodeIdx;
      m_localFaceNodeConn[faceIdx][nodeIdx] = 7; ++nodeIdx;
      m_localFaceNodeConn[faceIdx][nodeIdx] = 6;
      ++faceIdx; nodeIdx = 0;

      m_localFaceNodeConn[faceIdx].resize(4);
      m_localFaceNodeConn[faceIdx][nodeIdx] = 27; ++nodeIdx;
      m_localFaceNodeConn[faceIdx][nodeIdx] = 32; ++nodeIdx;
      m_localFaceNodeConn[faceIdx][nodeIdx] = 17; ++nodeIdx;
      m_localFaceNodeConn[faceIdx][nodeIdx] = 12;
      ++faceIdx; nodeIdx = 0;

      m_localFaceNodeConn[faceIdx].resize(5);
      m_localFaceNodeConn[faceIdx][nodeIdx] = 19; ++nodeIdx;
      m_localFaceNodeConn[faceIdx][nodeIdx] = 17; ++nodeIdx;
      m_localFaceNodeConn[faceIdx][nodeIdx] = 32; ++nodeIdx;
      m_localFaceNodeConn[faceIdx][nodeIdx] = 35; ++nodeIdx;
      m_localFaceNodeConn[faceIdx][nodeIdx] = 22;
      ++faceIdx; nodeIdx = 0;

      m_localFaceNodeConn[faceIdx].resize(4);
      m_localFaceNodeConn[faceIdx][nodeIdx] = 36; ++nodeIdx;
      m_localFaceNodeConn[faceIdx][nodeIdx] = 23; ++nodeIdx;
      m_localFaceNodeConn[faceIdx][nodeIdx] = 22; ++nodeIdx;
      m_localFaceNodeConn[faceIdx][nodeIdx] = 35;
      ++faceIdx; nodeIdx = 0;

      // INTERNAL FACES
      m_localFaceNodeConn[faceIdx].resize(4);
      m_localFaceNodeConn[faceIdx][nodeIdx] = 15; ++nodeIdx;
      m_localFaceNodeConn[faceIdx][nodeIdx] = 5; ++nodeIdx;
      m_localFaceNodeConn[faceIdx][nodeIdx] = 2; ++nodeIdx;
      m_localFaceNodeConn[faceIdx][nodeIdx] = 10;
      ++faceIdx; nodeIdx = 0;

      m_localFaceNodeConn[faceIdx].resize(4);
      m_localFaceNodeConn[faceIdx][nodeIdx] = 15; ++nodeIdx;
      m_localFaceNodeConn[faceIdx][nodeIdx] = 14; ++nodeIdx;
      m_localFaceNodeConn[faceIdx][nodeIdx] = 9; ++nodeIdx;
      m_localFaceNodeConn[faceIdx][nodeIdx] = 10;
      ++faceIdx; nodeIdx = 0;

      m_localFaceNodeConn[faceIdx].resize(4);
      m_localFaceNodeConn[faceIdx][nodeIdx] = 15; ++nodeIdx;
      m_localFaceNodeConn[faceIdx][nodeIdx] = 16; ++nodeIdx;
      m_localFaceNodeConn[faceIdx][nodeIdx] = 11; ++nodeIdx;
      m_localFaceNodeConn[faceIdx][nodeIdx] = 10;
      ++faceIdx; nodeIdx = 0;


      m_localFaceNodeConn[faceIdx].resize(4);
      m_localFaceNodeConn[faceIdx][nodeIdx] = 15; ++nodeIdx;
      m_localFaceNodeConn[faceIdx][nodeIdx] = 5; ++nodeIdx;
      m_localFaceNodeConn[faceIdx][nodeIdx] = 4; ++nodeIdx;
      m_localFaceNodeConn[faceIdx][nodeIdx] = 18;
      ++faceIdx; nodeIdx = 0;

      m_localFaceNodeConn[faceIdx].resize(4);
      m_localFaceNodeConn[faceIdx][nodeIdx] = 15; ++nodeIdx;
      m_localFaceNodeConn[faceIdx][nodeIdx] = 21; ++nodeIdx;
      m_localFaceNodeConn[faceIdx][nodeIdx] = 20; ++nodeIdx;
      m_localFaceNodeConn[faceIdx][nodeIdx] = 18;
      ++faceIdx; nodeIdx = 0;

      m_localFaceNodeConn[faceIdx].resize(4);
      m_localFaceNodeConn[faceIdx][nodeIdx] = 15; ++nodeIdx;
      m_localFaceNodeConn[faceIdx][nodeIdx] = 14; ++nodeIdx;
      m_localFaceNodeConn[faceIdx][nodeIdx] = 13; ++nodeIdx;
      m_localFaceNodeConn[faceIdx][nodeIdx] = 18;
      ++faceIdx; nodeIdx = 0;


      m_localFaceNodeConn[faceIdx].resize(4);
      m_localFaceNodeConn[faceIdx][nodeIdx] = 15; ++nodeIdx;
      m_localFaceNodeConn[faceIdx][nodeIdx] = 14; ++nodeIdx;
      m_localFaceNodeConn[faceIdx][nodeIdx] = 29; ++nodeIdx;
      m_localFaceNodeConn[faceIdx][nodeIdx] = 30;
      ++faceIdx; nodeIdx = 0;

      m_localFaceNodeConn[faceIdx].resize(4);
      m_localFaceNodeConn[faceIdx][nodeIdx] = 15; ++nodeIdx;
      m_localFaceNodeConn[faceIdx][nodeIdx] = 21; ++nodeIdx;
      m_localFaceNodeConn[faceIdx][nodeIdx] = 34; ++nodeIdx;
      m_localFaceNodeConn[faceIdx][nodeIdx] = 30;
      ++faceIdx; nodeIdx = 0;

      m_localFaceNodeConn[faceIdx].resize(4);
      m_localFaceNodeConn[faceIdx][nodeIdx] = 15; ++nodeIdx;
      m_localFaceNodeConn[faceIdx][nodeIdx] = 16; ++nodeIdx;
      m_localFaceNodeConn[faceIdx][nodeIdx] = 31; ++nodeIdx;
      m_localFaceNodeConn[faceIdx][nodeIdx] = 30;
      ++faceIdx; nodeIdx = 0;


      m_localFaceNodeConn[faceIdx].resize(4);
      m_localFaceNodeConn[faceIdx][nodeIdx] = 15; ++nodeIdx;
      m_localFaceNodeConn[faceIdx][nodeIdx] = 5; ++nodeIdx;
      m_localFaceNodeConn[faceIdx][nodeIdx] = 6; ++nodeIdx;
      m_localFaceNodeConn[faceIdx][nodeIdx] = 19;
      ++faceIdx; nodeIdx = 0;

      m_localFaceNodeConn[faceIdx].resize(4);
      m_localFaceNodeConn[faceIdx][nodeIdx] = 15; ++nodeIdx;
      m_localFaceNodeConn[faceIdx][nodeIdx] = 16; ++nodeIdx;
      m_localFaceNodeConn[faceIdx][nodeIdx] = 17; ++nodeIdx;
      m_localFaceNodeConn[faceIdx][nodeIdx] = 19;
      ++faceIdx; nodeIdx = 0;

      m_localFaceNodeConn[faceIdx].resize(4);
      m_localFaceNodeConn[faceIdx][nodeIdx] = 15; ++nodeIdx;
      m_localFaceNodeConn[faceIdx][nodeIdx] = 21; ++nodeIdx;
      m_localFaceNodeConn[faceIdx][nodeIdx] = 22; ++nodeIdx;
      m_localFaceNodeConn[faceIdx][nodeIdx] = 19;
      ++faceIdx; nodeIdx = 0;


//       // following faces are for SV3_3D schemes of first family with triangular corner CV faces
//       m_localFaceNodeConn[faceIdx].resize(3);
//       m_localFaceNodeConn[faceIdx][nodeIdx] = 7; ++nodeIdx;
//       m_localFaceNodeConn[faceIdx][nodeIdx] = 6; ++nodeIdx;
//       m_localFaceNodeConn[faceIdx][nodeIdx] = 4;
//       ++faceIdx; nodeIdx = 0;
//
//       m_localFaceNodeConn[faceIdx].resize(3);
//       m_localFaceNodeConn[faceIdx][nodeIdx] = 5; ++nodeIdx;
//       m_localFaceNodeConn[faceIdx][nodeIdx] = 4; ++nodeIdx;
//       m_localFaceNodeConn[faceIdx][nodeIdx] = 6;
//       ++faceIdx; nodeIdx = 0;
//
//       m_localFaceNodeConn[faceIdx].resize(3);
//       m_localFaceNodeConn[faceIdx][nodeIdx] = 3; ++nodeIdx;
//       m_localFaceNodeConn[faceIdx][nodeIdx] = 2; ++nodeIdx;
//       m_localFaceNodeConn[faceIdx][nodeIdx] = 6;
//       ++faceIdx; nodeIdx = 0;
//
//       m_localFaceNodeConn[faceIdx].resize(3);
//       m_localFaceNodeConn[faceIdx][nodeIdx] = 5; ++nodeIdx;
//       m_localFaceNodeConn[faceIdx][nodeIdx] = 6; ++nodeIdx;
//       m_localFaceNodeConn[faceIdx][nodeIdx] = 2;
//       ++faceIdx; nodeIdx = 0;
//
//       m_localFaceNodeConn[faceIdx].resize(3);
//       m_localFaceNodeConn[faceIdx][nodeIdx] = 1; ++nodeIdx;
//       m_localFaceNodeConn[faceIdx][nodeIdx] = 4; ++nodeIdx;
//       m_localFaceNodeConn[faceIdx][nodeIdx] = 2;
//       ++faceIdx; nodeIdx = 0;
//
//       m_localFaceNodeConn[faceIdx].resize(3);
//       m_localFaceNodeConn[faceIdx][nodeIdx] = 5; ++nodeIdx;
//       m_localFaceNodeConn[faceIdx][nodeIdx] = 2; ++nodeIdx;
//       m_localFaceNodeConn[faceIdx][nodeIdx] = 4;
//       ++faceIdx; nodeIdx = 0;
//
//
//       m_localFaceNodeConn[faceIdx].resize(3);
//       m_localFaceNodeConn[faceIdx][nodeIdx] = 8; ++nodeIdx;
//       m_localFaceNodeConn[faceIdx][nodeIdx] = 9; ++nodeIdx;
//       m_localFaceNodeConn[faceIdx][nodeIdx] = 13;
//       ++faceIdx; nodeIdx = 0;
//
//       m_localFaceNodeConn[faceIdx].resize(3);
//       m_localFaceNodeConn[faceIdx][nodeIdx] = 14; ++nodeIdx;
//       m_localFaceNodeConn[faceIdx][nodeIdx] = 13; ++nodeIdx;
//       m_localFaceNodeConn[faceIdx][nodeIdx] = 9;
//       ++faceIdx; nodeIdx = 0;
//
//       m_localFaceNodeConn[faceIdx].resize(3);
//       m_localFaceNodeConn[faceIdx][nodeIdx] = 25; ++nodeIdx;
//       m_localFaceNodeConn[faceIdx][nodeIdx] = 29; ++nodeIdx;
//       m_localFaceNodeConn[faceIdx][nodeIdx] = 9;
//       ++faceIdx; nodeIdx = 0;
//
//       m_localFaceNodeConn[faceIdx].resize(3);
//       m_localFaceNodeConn[faceIdx][nodeIdx] = 14; ++nodeIdx;
//       m_localFaceNodeConn[faceIdx][nodeIdx] = 9; ++nodeIdx;
//       m_localFaceNodeConn[faceIdx][nodeIdx] = 29;
//       ++faceIdx; nodeIdx = 0;
//
//       m_localFaceNodeConn[faceIdx].resize(3);
//       m_localFaceNodeConn[faceIdx][nodeIdx] = 28; ++nodeIdx;
//       m_localFaceNodeConn[faceIdx][nodeIdx] = 13; ++nodeIdx;
//       m_localFaceNodeConn[faceIdx][nodeIdx] = 29;
//       ++faceIdx; nodeIdx = 0;
//
//       m_localFaceNodeConn[faceIdx].resize(3);
//       m_localFaceNodeConn[faceIdx][nodeIdx] = 14; ++nodeIdx;
//       m_localFaceNodeConn[faceIdx][nodeIdx] = 29; ++nodeIdx;
//       m_localFaceNodeConn[faceIdx][nodeIdx] = 13;
//       ++faceIdx; nodeIdx = 0;
//
//
//       m_localFaceNodeConn[faceIdx].resize(3);
//       m_localFaceNodeConn[faceIdx][nodeIdx] = 26; ++nodeIdx;
//       m_localFaceNodeConn[faceIdx][nodeIdx] = 11; ++nodeIdx;
//       m_localFaceNodeConn[faceIdx][nodeIdx] = 31;
//       ++faceIdx; nodeIdx = 0;
//
//       m_localFaceNodeConn[faceIdx].resize(3);
//       m_localFaceNodeConn[faceIdx][nodeIdx] = 16; ++nodeIdx;
//       m_localFaceNodeConn[faceIdx][nodeIdx] = 31; ++nodeIdx;
//       m_localFaceNodeConn[faceIdx][nodeIdx] = 11;
//       ++faceIdx; nodeIdx = 0;
//
//       m_localFaceNodeConn[faceIdx].resize(3);
//       m_localFaceNodeConn[faceIdx][nodeIdx] = 12; ++nodeIdx;
//       m_localFaceNodeConn[faceIdx][nodeIdx] = 17; ++nodeIdx;
//       m_localFaceNodeConn[faceIdx][nodeIdx] = 11;
//       ++faceIdx; nodeIdx = 0;
//
//       m_localFaceNodeConn[faceIdx].resize(3);
//       m_localFaceNodeConn[faceIdx][nodeIdx] = 16; ++nodeIdx;
//       m_localFaceNodeConn[faceIdx][nodeIdx] = 11; ++nodeIdx;
//       m_localFaceNodeConn[faceIdx][nodeIdx] = 17;
//       ++faceIdx; nodeIdx = 0;
//
//       m_localFaceNodeConn[faceIdx].resize(3);
//       m_localFaceNodeConn[faceIdx][nodeIdx] = 32; ++nodeIdx;
//       m_localFaceNodeConn[faceIdx][nodeIdx] = 31; ++nodeIdx;
//       m_localFaceNodeConn[faceIdx][nodeIdx] = 17;
//       ++faceIdx; nodeIdx = 0;
//
//       m_localFaceNodeConn[faceIdx].resize(3);
//       m_localFaceNodeConn[faceIdx][nodeIdx] = 16; ++nodeIdx;
//       m_localFaceNodeConn[faceIdx][nodeIdx] = 17; ++nodeIdx;
//       m_localFaceNodeConn[faceIdx][nodeIdx] = 31;
//       ++faceIdx; nodeIdx = 0;
//
//
//       m_localFaceNodeConn[faceIdx].resize(3);
//       m_localFaceNodeConn[faceIdx][nodeIdx] = 23; ++nodeIdx;
//       m_localFaceNodeConn[faceIdx][nodeIdx] = 20; ++nodeIdx;
//       m_localFaceNodeConn[faceIdx][nodeIdx] = 22;
//       ++faceIdx; nodeIdx = 0;
//
//       m_localFaceNodeConn[faceIdx].resize(3);
//       m_localFaceNodeConn[faceIdx][nodeIdx] = 21; ++nodeIdx;
//       m_localFaceNodeConn[faceIdx][nodeIdx] = 22; ++nodeIdx;
//       m_localFaceNodeConn[faceIdx][nodeIdx] = 20;
//       ++faceIdx; nodeIdx = 0;
//
//       m_localFaceNodeConn[faceIdx].resize(3);
//       m_localFaceNodeConn[faceIdx][nodeIdx] = 35; ++nodeIdx;
//       m_localFaceNodeConn[faceIdx][nodeIdx] = 22; ++nodeIdx;
//       m_localFaceNodeConn[faceIdx][nodeIdx] = 34;
//       ++faceIdx; nodeIdx = 0;
//
//       m_localFaceNodeConn[faceIdx].resize(3);
//       m_localFaceNodeConn[faceIdx][nodeIdx] = 21; ++nodeIdx;
//       m_localFaceNodeConn[faceIdx][nodeIdx] = 34; ++nodeIdx;
//       m_localFaceNodeConn[faceIdx][nodeIdx] = 22;
//       ++faceIdx; nodeIdx = 0;
//
//       m_localFaceNodeConn[faceIdx].resize(3);
//       m_localFaceNodeConn[faceIdx][nodeIdx] = 33; ++nodeIdx;
//       m_localFaceNodeConn[faceIdx][nodeIdx] = 34; ++nodeIdx;
//       m_localFaceNodeConn[faceIdx][nodeIdx] = 20;
//       ++faceIdx; nodeIdx = 0;
//
//       m_localFaceNodeConn[faceIdx].resize(3);
//       m_localFaceNodeConn[faceIdx][nodeIdx] = 21; ++nodeIdx;
//       m_localFaceNodeConn[faceIdx][nodeIdx] = 20; ++nodeIdx;
//       m_localFaceNodeConn[faceIdx][nodeIdx] = 34;
//       ++faceIdx; nodeIdx = 0;


//       // following faces are for SV3_3D schemes of second family with triangular corner CV faces
//       m_localFaceNodeConn[faceIdx].resize(3);
//       m_localFaceNodeConn[faceIdx][nodeIdx] = 7; ++nodeIdx;
//       m_localFaceNodeConn[faceIdx][nodeIdx] = 5; ++nodeIdx;
//       m_localFaceNodeConn[faceIdx][nodeIdx] = 4;
//       ++faceIdx; nodeIdx = 0;
//
//       m_localFaceNodeConn[faceIdx].resize(3);
//       m_localFaceNodeConn[faceIdx][nodeIdx] = 7; ++nodeIdx;
//       m_localFaceNodeConn[faceIdx][nodeIdx] = 6; ++nodeIdx;
//       m_localFaceNodeConn[faceIdx][nodeIdx] = 5;
//       ++faceIdx; nodeIdx = 0;
//
//       m_localFaceNodeConn[faceIdx].resize(3);
//       m_localFaceNodeConn[faceIdx][nodeIdx] = 3; ++nodeIdx;
//       m_localFaceNodeConn[faceIdx][nodeIdx] = 5; ++nodeIdx;
//       m_localFaceNodeConn[faceIdx][nodeIdx] = 6;
//       ++faceIdx; nodeIdx = 0;
//
//       m_localFaceNodeConn[faceIdx].resize(3);
//       m_localFaceNodeConn[faceIdx][nodeIdx] = 3; ++nodeIdx;
//       m_localFaceNodeConn[faceIdx][nodeIdx] = 2; ++nodeIdx;
//       m_localFaceNodeConn[faceIdx][nodeIdx] = 5;
//       ++faceIdx; nodeIdx = 0;
//
//       m_localFaceNodeConn[faceIdx].resize(3);
//       m_localFaceNodeConn[faceIdx][nodeIdx] = 1; ++nodeIdx;
//       m_localFaceNodeConn[faceIdx][nodeIdx] = 5; ++nodeIdx;
//       m_localFaceNodeConn[faceIdx][nodeIdx] = 2;
//       ++faceIdx; nodeIdx = 0;
//
//       m_localFaceNodeConn[faceIdx].resize(3);
//       m_localFaceNodeConn[faceIdx][nodeIdx] = 1; ++nodeIdx;
//       m_localFaceNodeConn[faceIdx][nodeIdx] = 4; ++nodeIdx;
//       m_localFaceNodeConn[faceIdx][nodeIdx] = 5;
//       ++faceIdx; nodeIdx = 0;
//
//
//       m_localFaceNodeConn[faceIdx].resize(3);
//       m_localFaceNodeConn[faceIdx][nodeIdx] = 8; ++nodeIdx;
//       m_localFaceNodeConn[faceIdx][nodeIdx] = 9; ++nodeIdx;
//       m_localFaceNodeConn[faceIdx][nodeIdx] = 14;
//       ++faceIdx; nodeIdx = 0;
//
//       m_localFaceNodeConn[faceIdx].resize(3);
//       m_localFaceNodeConn[faceIdx][nodeIdx] = 8; ++nodeIdx;
//       m_localFaceNodeConn[faceIdx][nodeIdx] = 14; ++nodeIdx;
//       m_localFaceNodeConn[faceIdx][nodeIdx] = 13;
//       ++faceIdx; nodeIdx = 0;
//
//       m_localFaceNodeConn[faceIdx].resize(3);
//       m_localFaceNodeConn[faceIdx][nodeIdx] = 25; ++nodeIdx;
//       m_localFaceNodeConn[faceIdx][nodeIdx] = 29; ++nodeIdx;
//       m_localFaceNodeConn[faceIdx][nodeIdx] = 14;
//       ++faceIdx; nodeIdx = 0;
//
//       m_localFaceNodeConn[faceIdx].resize(3);
//       m_localFaceNodeConn[faceIdx][nodeIdx] = 25; ++nodeIdx;
//       m_localFaceNodeConn[faceIdx][nodeIdx] = 14; ++nodeIdx;
//       m_localFaceNodeConn[faceIdx][nodeIdx] = 9;
//       ++faceIdx; nodeIdx = 0;
//
//       m_localFaceNodeConn[faceIdx].resize(3);
//       m_localFaceNodeConn[faceIdx][nodeIdx] = 28; ++nodeIdx;
//       m_localFaceNodeConn[faceIdx][nodeIdx] = 13; ++nodeIdx;
//       m_localFaceNodeConn[faceIdx][nodeIdx] = 14;
//       ++faceIdx; nodeIdx = 0;
//
//       m_localFaceNodeConn[faceIdx].resize(3);
//       m_localFaceNodeConn[faceIdx][nodeIdx] = 28; ++nodeIdx;
//       m_localFaceNodeConn[faceIdx][nodeIdx] = 14; ++nodeIdx;
//       m_localFaceNodeConn[faceIdx][nodeIdx] = 29;
//       ++faceIdx; nodeIdx = 0;
//
//
//       m_localFaceNodeConn[faceIdx].resize(3);
//       m_localFaceNodeConn[faceIdx][nodeIdx] = 26; ++nodeIdx;
//       m_localFaceNodeConn[faceIdx][nodeIdx] = 11; ++nodeIdx;
//       m_localFaceNodeConn[faceIdx][nodeIdx] = 16;
//       ++faceIdx; nodeIdx = 0;
//
//       m_localFaceNodeConn[faceIdx].resize(3);
//       m_localFaceNodeConn[faceIdx][nodeIdx] = 26; ++nodeIdx;
//       m_localFaceNodeConn[faceIdx][nodeIdx] = 16; ++nodeIdx;
//       m_localFaceNodeConn[faceIdx][nodeIdx] = 31;
//       ++faceIdx; nodeIdx = 0;
//
//       m_localFaceNodeConn[faceIdx].resize(3);
//       m_localFaceNodeConn[faceIdx][nodeIdx] = 12; ++nodeIdx;
//       m_localFaceNodeConn[faceIdx][nodeIdx] = 17; ++nodeIdx;
//       m_localFaceNodeConn[faceIdx][nodeIdx] = 16;
//       ++faceIdx; nodeIdx = 0;
//
//       m_localFaceNodeConn[faceIdx].resize(3);
//       m_localFaceNodeConn[faceIdx][nodeIdx] = 12; ++nodeIdx;
//       m_localFaceNodeConn[faceIdx][nodeIdx] = 16; ++nodeIdx;
//       m_localFaceNodeConn[faceIdx][nodeIdx] = 11;
//       ++faceIdx; nodeIdx = 0;
//
//       m_localFaceNodeConn[faceIdx].resize(3);
//       m_localFaceNodeConn[faceIdx][nodeIdx] = 32; ++nodeIdx;
//       m_localFaceNodeConn[faceIdx][nodeIdx] = 31; ++nodeIdx;
//       m_localFaceNodeConn[faceIdx][nodeIdx] = 16;
//       ++faceIdx; nodeIdx = 0;
//
//       m_localFaceNodeConn[faceIdx].resize(3);
//       m_localFaceNodeConn[faceIdx][nodeIdx] = 32; ++nodeIdx;
//       m_localFaceNodeConn[faceIdx][nodeIdx] = 16; ++nodeIdx;
//       m_localFaceNodeConn[faceIdx][nodeIdx] = 17;
//       ++faceIdx; nodeIdx = 0;
//
//
//       m_localFaceNodeConn[faceIdx].resize(3);
//       m_localFaceNodeConn[faceIdx][nodeIdx] = 23; ++nodeIdx;
//       m_localFaceNodeConn[faceIdx][nodeIdx] = 20; ++nodeIdx;
//       m_localFaceNodeConn[faceIdx][nodeIdx] = 21;
//       ++faceIdx; nodeIdx = 0;
//
//       m_localFaceNodeConn[faceIdx].resize(3);
//       m_localFaceNodeConn[faceIdx][nodeIdx] = 23; ++nodeIdx;
//       m_localFaceNodeConn[faceIdx][nodeIdx] = 21; ++nodeIdx;
//       m_localFaceNodeConn[faceIdx][nodeIdx] = 22;
//       ++faceIdx; nodeIdx = 0;
//
//       m_localFaceNodeConn[faceIdx].resize(3);
//       m_localFaceNodeConn[faceIdx][nodeIdx] = 35; ++nodeIdx;
//       m_localFaceNodeConn[faceIdx][nodeIdx] = 22; ++nodeIdx;
//       m_localFaceNodeConn[faceIdx][nodeIdx] = 21;
//       ++faceIdx; nodeIdx = 0;
//
//       m_localFaceNodeConn[faceIdx].resize(3);
//       m_localFaceNodeConn[faceIdx][nodeIdx] = 35; ++nodeIdx;
//       m_localFaceNodeConn[faceIdx][nodeIdx] = 21; ++nodeIdx;
//       m_localFaceNodeConn[faceIdx][nodeIdx] = 34;
//       ++faceIdx; nodeIdx = 0;
//
//       m_localFaceNodeConn[faceIdx].resize(3);
//       m_localFaceNodeConn[faceIdx][nodeIdx] = 33; ++nodeIdx;
//       m_localFaceNodeConn[faceIdx][nodeIdx] = 34; ++nodeIdx;
//       m_localFaceNodeConn[faceIdx][nodeIdx] = 21;
//       ++faceIdx; nodeIdx = 0;
//
//       m_localFaceNodeConn[faceIdx].resize(3);
//       m_localFaceNodeConn[faceIdx][nodeIdx] = 33; ++nodeIdx;
//       m_localFaceNodeConn[faceIdx][nodeIdx] = 21; ++nodeIdx;
//       m_localFaceNodeConn[faceIdx][nodeIdx] = 20;
//       ++faceIdx; nodeIdx = 0;

      // following faces are for SV3_3D schemes with quadrilateral corner CV faces
      m_localFaceNodeConn[faceIdx].resize(4);
      m_localFaceNodeConn[faceIdx][nodeIdx] = 7; ++nodeIdx;
      m_localFaceNodeConn[faceIdx][nodeIdx] = 6; ++nodeIdx;
      m_localFaceNodeConn[faceIdx][nodeIdx] = 5; ++nodeIdx;
      m_localFaceNodeConn[faceIdx][nodeIdx] = 4;
      ++faceIdx; nodeIdx = 0;

      m_localFaceNodeConn[faceIdx].resize(4);
      m_localFaceNodeConn[faceIdx][nodeIdx] = 3; ++nodeIdx;
      m_localFaceNodeConn[faceIdx][nodeIdx] = 2; ++nodeIdx;
      m_localFaceNodeConn[faceIdx][nodeIdx] = 5; ++nodeIdx;
      m_localFaceNodeConn[faceIdx][nodeIdx] = 6;
      ++faceIdx; nodeIdx = 0;

      m_localFaceNodeConn[faceIdx].resize(4);
      m_localFaceNodeConn[faceIdx][nodeIdx] = 1; ++nodeIdx;
      m_localFaceNodeConn[faceIdx][nodeIdx] = 4; ++nodeIdx;
      m_localFaceNodeConn[faceIdx][nodeIdx] = 5; ++nodeIdx;
      m_localFaceNodeConn[faceIdx][nodeIdx] = 2;
      ++faceIdx; nodeIdx = 0;


      m_localFaceNodeConn[faceIdx].resize(4);
      m_localFaceNodeConn[faceIdx][nodeIdx] = 28; ++nodeIdx;
      m_localFaceNodeConn[faceIdx][nodeIdx] = 13; ++nodeIdx;
      m_localFaceNodeConn[faceIdx][nodeIdx] = 14; ++nodeIdx;
      m_localFaceNodeConn[faceIdx][nodeIdx] = 29;
      ++faceIdx; nodeIdx = 0;

      m_localFaceNodeConn[faceIdx].resize(4);
      m_localFaceNodeConn[faceIdx][nodeIdx] = 25; ++nodeIdx;
      m_localFaceNodeConn[faceIdx][nodeIdx] = 29; ++nodeIdx;
      m_localFaceNodeConn[faceIdx][nodeIdx] = 14; ++nodeIdx;
      m_localFaceNodeConn[faceIdx][nodeIdx] = 9;
      ++faceIdx; nodeIdx = 0;

      m_localFaceNodeConn[faceIdx].resize(4);
      m_localFaceNodeConn[faceIdx][nodeIdx] = 8; ++nodeIdx;
      m_localFaceNodeConn[faceIdx][nodeIdx] = 9; ++nodeIdx;
      m_localFaceNodeConn[faceIdx][nodeIdx] = 14; ++nodeIdx;
      m_localFaceNodeConn[faceIdx][nodeIdx] = 13;
      ++faceIdx; nodeIdx = 0;


      m_localFaceNodeConn[faceIdx].resize(4);
      m_localFaceNodeConn[faceIdx][nodeIdx] = 32; ++nodeIdx;
      m_localFaceNodeConn[faceIdx][nodeIdx] = 31; ++nodeIdx;
      m_localFaceNodeConn[faceIdx][nodeIdx] = 16; ++nodeIdx;
      m_localFaceNodeConn[faceIdx][nodeIdx] = 17;
      ++faceIdx; nodeIdx = 0;

      m_localFaceNodeConn[faceIdx].resize(4);
      m_localFaceNodeConn[faceIdx][nodeIdx] = 26; ++nodeIdx;
      m_localFaceNodeConn[faceIdx][nodeIdx] = 11; ++nodeIdx;
      m_localFaceNodeConn[faceIdx][nodeIdx] = 16; ++nodeIdx;
      m_localFaceNodeConn[faceIdx][nodeIdx] = 31;
      ++faceIdx; nodeIdx = 0;

      m_localFaceNodeConn[faceIdx].resize(4);
      m_localFaceNodeConn[faceIdx][nodeIdx] = 12; ++nodeIdx;
      m_localFaceNodeConn[faceIdx][nodeIdx] = 17; ++nodeIdx;
      m_localFaceNodeConn[faceIdx][nodeIdx] = 16; ++nodeIdx;
      m_localFaceNodeConn[faceIdx][nodeIdx] = 11;
      ++faceIdx; nodeIdx = 0;


      m_localFaceNodeConn[faceIdx].resize(4);
      m_localFaceNodeConn[faceIdx][nodeIdx] = 23; ++nodeIdx;
      m_localFaceNodeConn[faceIdx][nodeIdx] = 20; ++nodeIdx;
      m_localFaceNodeConn[faceIdx][nodeIdx] = 21; ++nodeIdx;
      m_localFaceNodeConn[faceIdx][nodeIdx] = 22;
      ++faceIdx; nodeIdx = 0;

      m_localFaceNodeConn[faceIdx].resize(4);
      m_localFaceNodeConn[faceIdx][nodeIdx] = 33; ++nodeIdx;
      m_localFaceNodeConn[faceIdx][nodeIdx] = 34; ++nodeIdx;
      m_localFaceNodeConn[faceIdx][nodeIdx] = 21; ++nodeIdx;
      m_localFaceNodeConn[faceIdx][nodeIdx] = 20;
      ++faceIdx; nodeIdx = 0;

      m_localFaceNodeConn[faceIdx].resize(4);
      m_localFaceNodeConn[faceIdx][nodeIdx] = 35; ++nodeIdx;
      m_localFaceNodeConn[faceIdx][nodeIdx] = 22; ++nodeIdx;
      m_localFaceNodeConn[faceIdx][nodeIdx] = 21; ++nodeIdx;
      m_localFaceNodeConn[faceIdx][nodeIdx] = 34;
      ++faceIdx; nodeIdx = 0;

      cf_assert(faceIdx == m_localFaceNodeConn.size());

    } break;
    case CFPolyOrder::ORDER3:
    {
      throw Common::NotImplementedException (FromHere(),"4th order accurate tetrahedral SV cells have not yet been implemented!");
    } break;
    default:
    {
      throw Common::NotImplementedException (FromHere(),"Only solution orders up to 3 have been implemented for spectral FV!");
    }
  }

/*  for (CFuint iFace = 0; iFace < m_localFaceNodeConn.size(); ++iFace)
  {
CF_DEBUG_OBJ(iFace);
    for (CFuint iNode = 0; iNode < m_localFaceNodeConn[iFace].size(); ++iNode)
    {
CF_DEBUG_OBJ(m_localFaceNodeConn[iFace][iNode]);
    }
  }*/
}

//////////////////////////////////////////////////////////////////////

void TetraSpectralFVElementData::createLocalCVNodeConn()
{
  CFAUTOTRACE;

  // number of CVs
  const CFuint nbrCVs = (m_polyOrder+1)*(m_polyOrder+2)*(m_polyOrder+3)/6;

  // resize local CV-node connectivity
  m_localCVNodeConn.resize(nbrCVs);

  // create connectivity
  if (nbrCVs == 1)
  {
    // resize local CV-node connectivity
    m_localCVNodeConn[0].resize(4);

    // fill connectivity
    m_localCVNodeConn[0][0] = 0;
    m_localCVNodeConn[0][1] = 1;
    m_localCVNodeConn[0][2] = 2;
    m_localCVNodeConn[0][3] = 3;
  }
  else
  {
    switch (m_polyOrder)
    {
      case CFPolyOrder::ORDER1:
      {
        // index for CVs
        CFuint cvIdx = 0;

        // index for nodes
        CFuint nodeIdx = 0;

        // first CV
        m_localCVNodeConn[cvIdx].resize(8);
        m_localCVNodeConn[cvIdx][nodeIdx] = 0; ++nodeIdx;
        m_localCVNodeConn[cvIdx][nodeIdx] = 1; ++nodeIdx;
        m_localCVNodeConn[cvIdx][nodeIdx] = 2; ++nodeIdx;
        m_localCVNodeConn[cvIdx][nodeIdx] = 3; ++nodeIdx;
        m_localCVNodeConn[cvIdx][nodeIdx] = 7; ++nodeIdx;
        m_localCVNodeConn[cvIdx][nodeIdx] = 4; ++nodeIdx;
        m_localCVNodeConn[cvIdx][nodeIdx] = 5; ++nodeIdx;
        m_localCVNodeConn[cvIdx][nodeIdx] = 6;
        ++cvIdx; nodeIdx = 0;

        // second CV
        m_localCVNodeConn[cvIdx].resize(8);
        m_localCVNodeConn[cvIdx][nodeIdx] = 8; ++nodeIdx;
        m_localCVNodeConn[cvIdx][nodeIdx] = 9; ++nodeIdx;
        m_localCVNodeConn[cvIdx][nodeIdx] = 2; ++nodeIdx;
        m_localCVNodeConn[cvIdx][nodeIdx] = 1; ++nodeIdx;
        m_localCVNodeConn[cvIdx][nodeIdx] = 11; ++nodeIdx;
        m_localCVNodeConn[cvIdx][nodeIdx] = 12; ++nodeIdx;
        m_localCVNodeConn[cvIdx][nodeIdx] = 5; ++nodeIdx;
        m_localCVNodeConn[cvIdx][nodeIdx] = 4;
        ++cvIdx; nodeIdx = 0;

        // third CV
        m_localCVNodeConn[cvIdx].resize(8);
        m_localCVNodeConn[cvIdx][nodeIdx] = 10; ++nodeIdx;
        m_localCVNodeConn[cvIdx][nodeIdx] = 3; ++nodeIdx;
        m_localCVNodeConn[cvIdx][nodeIdx] = 2; ++nodeIdx;
        m_localCVNodeConn[cvIdx][nodeIdx] = 9; ++nodeIdx;
        m_localCVNodeConn[cvIdx][nodeIdx] = 13; ++nodeIdx;
        m_localCVNodeConn[cvIdx][nodeIdx] = 6; ++nodeIdx;
        m_localCVNodeConn[cvIdx][nodeIdx] = 5; ++nodeIdx;
        m_localCVNodeConn[cvIdx][nodeIdx] = 12;
        ++cvIdx; nodeIdx = 0;

        // fourth CV
        m_localCVNodeConn[cvIdx].resize(8);
        m_localCVNodeConn[cvIdx][nodeIdx] = 14; ++nodeIdx;
        m_localCVNodeConn[cvIdx][nodeIdx] = 13; ++nodeIdx;
        m_localCVNodeConn[cvIdx][nodeIdx] = 12; ++nodeIdx;
        m_localCVNodeConn[cvIdx][nodeIdx] = 11; ++nodeIdx;
        m_localCVNodeConn[cvIdx][nodeIdx] = 7; ++nodeIdx;
        m_localCVNodeConn[cvIdx][nodeIdx] = 6; ++nodeIdx;
        m_localCVNodeConn[cvIdx][nodeIdx] = 5; ++nodeIdx;
        m_localCVNodeConn[cvIdx][nodeIdx] = 4;
      } break;
      case CFPolyOrder::ORDER2:
      {
        // index for CVs
        CFuint cvIdx = 0;

        // index for nodes
        CFuint nodeIdx = 0;

        // first CV
        m_localCVNodeConn[cvIdx].resize(8);
        m_localCVNodeConn[cvIdx][nodeIdx] = 0; ++nodeIdx;
        m_localCVNodeConn[cvIdx][nodeIdx] = 1; ++nodeIdx;
        m_localCVNodeConn[cvIdx][nodeIdx] = 2; ++nodeIdx;
        m_localCVNodeConn[cvIdx][nodeIdx] = 3; ++nodeIdx;
        m_localCVNodeConn[cvIdx][nodeIdx] = 7; ++nodeIdx;
        m_localCVNodeConn[cvIdx][nodeIdx] = 4; ++nodeIdx;
        m_localCVNodeConn[cvIdx][nodeIdx] = 5; ++nodeIdx;
        m_localCVNodeConn[cvIdx][nodeIdx] = 6;
        ++cvIdx; nodeIdx = 0;

        // second CV
        m_localCVNodeConn[cvIdx].resize(11);
        m_localCVNodeConn[cvIdx][nodeIdx] = 2; ++nodeIdx;
        m_localCVNodeConn[cvIdx][nodeIdx] = 5; ++nodeIdx;
        m_localCVNodeConn[cvIdx][nodeIdx] = 4; ++nodeIdx;
        m_localCVNodeConn[cvIdx][nodeIdx] = 1; ++nodeIdx;
        m_localCVNodeConn[cvIdx][nodeIdx] = 10; ++nodeIdx;
        m_localCVNodeConn[cvIdx][nodeIdx] = 15; ++nodeIdx;
        m_localCVNodeConn[cvIdx][nodeIdx] = 18; ++nodeIdx;
        m_localCVNodeConn[cvIdx][nodeIdx] = 8; ++nodeIdx;
        m_localCVNodeConn[cvIdx][nodeIdx] = 9; ++nodeIdx;
        m_localCVNodeConn[cvIdx][nodeIdx] = 14; ++nodeIdx;
        m_localCVNodeConn[cvIdx][nodeIdx] = 13;
        ++cvIdx; nodeIdx = 0;

        // third CV
        m_localCVNodeConn[cvIdx].resize(11);
        m_localCVNodeConn[cvIdx][nodeIdx] = 6; ++nodeIdx;
        m_localCVNodeConn[cvIdx][nodeIdx] = 5; ++nodeIdx;
        m_localCVNodeConn[cvIdx][nodeIdx] = 2; ++nodeIdx;
        m_localCVNodeConn[cvIdx][nodeIdx] = 3; ++nodeIdx;
        m_localCVNodeConn[cvIdx][nodeIdx] = 19; ++nodeIdx;
        m_localCVNodeConn[cvIdx][nodeIdx] = 15; ++nodeIdx;
        m_localCVNodeConn[cvIdx][nodeIdx] = 10; ++nodeIdx;
        m_localCVNodeConn[cvIdx][nodeIdx] = 12; ++nodeIdx;
        m_localCVNodeConn[cvIdx][nodeIdx] = 17; ++nodeIdx;
        m_localCVNodeConn[cvIdx][nodeIdx] = 16; ++nodeIdx;
        m_localCVNodeConn[cvIdx][nodeIdx] = 11;
        ++cvIdx; nodeIdx = 0;

        // fourth CV
        m_localCVNodeConn[cvIdx].resize(11);
        m_localCVNodeConn[cvIdx][nodeIdx] = 4; ++nodeIdx;
        m_localCVNodeConn[cvIdx][nodeIdx] = 5; ++nodeIdx;
        m_localCVNodeConn[cvIdx][nodeIdx] = 6; ++nodeIdx;
        m_localCVNodeConn[cvIdx][nodeIdx] = 7; ++nodeIdx;
        m_localCVNodeConn[cvIdx][nodeIdx] = 18; ++nodeIdx;
        m_localCVNodeConn[cvIdx][nodeIdx] = 15; ++nodeIdx;
        m_localCVNodeConn[cvIdx][nodeIdx] = 19; ++nodeIdx;
        m_localCVNodeConn[cvIdx][nodeIdx] = 23; ++nodeIdx;
        m_localCVNodeConn[cvIdx][nodeIdx] = 20; ++nodeIdx;
        m_localCVNodeConn[cvIdx][nodeIdx] = 21; ++nodeIdx;
        m_localCVNodeConn[cvIdx][nodeIdx] = 22;
        ++cvIdx; nodeIdx = 0;

        // fifth CV
        m_localCVNodeConn[cvIdx].resize(8);
        m_localCVNodeConn[cvIdx][nodeIdx] = 24; ++nodeIdx;
        m_localCVNodeConn[cvIdx][nodeIdx] = 25; ++nodeIdx;
        m_localCVNodeConn[cvIdx][nodeIdx] = 9; ++nodeIdx;
        m_localCVNodeConn[cvIdx][nodeIdx] = 8; ++nodeIdx;
        m_localCVNodeConn[cvIdx][nodeIdx] = 28; ++nodeIdx;
        m_localCVNodeConn[cvIdx][nodeIdx] = 29; ++nodeIdx;
        m_localCVNodeConn[cvIdx][nodeIdx] = 14; ++nodeIdx;
        m_localCVNodeConn[cvIdx][nodeIdx] = 13;
        ++cvIdx; nodeIdx = 0;

        // sixth CV
        m_localCVNodeConn[cvIdx].resize(11);
        m_localCVNodeConn[cvIdx][nodeIdx] = 9; ++nodeIdx;
        m_localCVNodeConn[cvIdx][nodeIdx] = 14; ++nodeIdx;
        m_localCVNodeConn[cvIdx][nodeIdx] = 29; ++nodeIdx;
        m_localCVNodeConn[cvIdx][nodeIdx] = 25; ++nodeIdx;
        m_localCVNodeConn[cvIdx][nodeIdx] = 10; ++nodeIdx;
        m_localCVNodeConn[cvIdx][nodeIdx] = 15; ++nodeIdx;
        m_localCVNodeConn[cvIdx][nodeIdx] = 30; ++nodeIdx;
        m_localCVNodeConn[cvIdx][nodeIdx] = 26; ++nodeIdx;
        m_localCVNodeConn[cvIdx][nodeIdx] = 11; ++nodeIdx;
        m_localCVNodeConn[cvIdx][nodeIdx] = 16; ++nodeIdx;
        m_localCVNodeConn[cvIdx][nodeIdx] = 31;
        ++cvIdx; nodeIdx = 0;

        // seventh CV
        m_localCVNodeConn[cvIdx].resize(8);
        m_localCVNodeConn[cvIdx][nodeIdx] = 27; ++nodeIdx;
        m_localCVNodeConn[cvIdx][nodeIdx] = 12; ++nodeIdx;
        m_localCVNodeConn[cvIdx][nodeIdx] = 11; ++nodeIdx;
        m_localCVNodeConn[cvIdx][nodeIdx] = 26; ++nodeIdx;
        m_localCVNodeConn[cvIdx][nodeIdx] = 32; ++nodeIdx;
        m_localCVNodeConn[cvIdx][nodeIdx] = 17; ++nodeIdx;
        m_localCVNodeConn[cvIdx][nodeIdx] = 16; ++nodeIdx;
        m_localCVNodeConn[cvIdx][nodeIdx] = 31;
        ++cvIdx; nodeIdx = 0;

        // eighth CV
        m_localCVNodeConn[cvIdx].resize(11);
        m_localCVNodeConn[cvIdx][nodeIdx] = 29; ++nodeIdx;
        m_localCVNodeConn[cvIdx][nodeIdx] = 14; ++nodeIdx;
        m_localCVNodeConn[cvIdx][nodeIdx] = 13; ++nodeIdx;
        m_localCVNodeConn[cvIdx][nodeIdx] = 28; ++nodeIdx;
        m_localCVNodeConn[cvIdx][nodeIdx] = 30; ++nodeIdx;
        m_localCVNodeConn[cvIdx][nodeIdx] = 15; ++nodeIdx;
        m_localCVNodeConn[cvIdx][nodeIdx] = 18; ++nodeIdx;
        m_localCVNodeConn[cvIdx][nodeIdx] = 33; ++nodeIdx;
        m_localCVNodeConn[cvIdx][nodeIdx] = 34; ++nodeIdx;
        m_localCVNodeConn[cvIdx][nodeIdx] = 21; ++nodeIdx;
        m_localCVNodeConn[cvIdx][nodeIdx] = 20;
        ++cvIdx; nodeIdx = 0;

        // ninth CV
        m_localCVNodeConn[cvIdx].resize(11);
        m_localCVNodeConn[cvIdx][nodeIdx] = 17; ++nodeIdx;
        m_localCVNodeConn[cvIdx][nodeIdx] = 16; ++nodeIdx;
        m_localCVNodeConn[cvIdx][nodeIdx] = 31; ++nodeIdx;
        m_localCVNodeConn[cvIdx][nodeIdx] = 32; ++nodeIdx;
        m_localCVNodeConn[cvIdx][nodeIdx] = 19; ++nodeIdx;
        m_localCVNodeConn[cvIdx][nodeIdx] = 15; ++nodeIdx;
        m_localCVNodeConn[cvIdx][nodeIdx] = 30; ++nodeIdx;
        m_localCVNodeConn[cvIdx][nodeIdx] = 35; ++nodeIdx;
        m_localCVNodeConn[cvIdx][nodeIdx] = 22; ++nodeIdx;
        m_localCVNodeConn[cvIdx][nodeIdx] = 21; ++nodeIdx;
        m_localCVNodeConn[cvIdx][nodeIdx] = 34;
        ++cvIdx; nodeIdx = 0;

        // tenth CV
        m_localCVNodeConn[cvIdx].resize(8);
        m_localCVNodeConn[cvIdx][nodeIdx] = 36; ++nodeIdx;
        m_localCVNodeConn[cvIdx][nodeIdx] = 35; ++nodeIdx;
        m_localCVNodeConn[cvIdx][nodeIdx] = 34; ++nodeIdx;
        m_localCVNodeConn[cvIdx][nodeIdx] = 33; ++nodeIdx;
        m_localCVNodeConn[cvIdx][nodeIdx] = 23; ++nodeIdx;
        m_localCVNodeConn[cvIdx][nodeIdx] = 22; ++nodeIdx;
        m_localCVNodeConn[cvIdx][nodeIdx] = 21; ++nodeIdx;
        m_localCVNodeConn[cvIdx][nodeIdx] = 20;
        ++cvIdx; nodeIdx = 0;
      } break;
      case CFPolyOrder::ORDER3:
      {
        throw Common::NotImplementedException (FromHere(),"4th order accurate tetrahedral SV cells have not yet been implemented!");
      } break;
      default:
      {
        throw Common::NotImplementedException (FromHere(),"Only solution orders up to 3 have been implemented for spectral FV!");
      }
    }
  }

  // write file with node coordinates and edge node connectivity for each CV
  ofstream outFile("svPlotData.dat");

  // number of nodes
  outFile << m_localNodeCoord.size() << "\n";

  // write node coordinates
  for (CFuint iNode = 0; iNode < m_localNodeCoord.size(); ++iNode)
  {
    for (CFuint iCoor = 0; iCoor < m_localNodeCoord[iNode].size(); ++iCoor)
    {
      outFile << m_localNodeCoord[iNode][iCoor] << " ";
    }
    outFile << "\n";
  }

  // number of CVs
  outFile << m_localCVNodeConn.size() << "\n";

  // write connectivity for each CV
  for (CFuint iCV = 0; iCV < m_localCVNodeConn.size(); ++iCV)
  {
    // get edge node connectivity for this CV
    vector< vector< CFuint > > edgeToNode;
    setPolyHedronEdgeNodeConn(m_localCVNodeConn[iCV].size(), edgeToNode);

    // number of edges in this CV
    outFile << edgeToNode.size() << "\n";

    // write edge node connectivity
    for (CFuint iEdge = 0; iEdge < edgeToNode.size(); ++iEdge)
    {
      const CFuint node0 = m_localCVNodeConn[iCV][edgeToNode[iEdge][0]];
      const CFuint node1 = m_localCVNodeConn[iCV][edgeToNode[iEdge][1]];
      outFile << node0 << " " << node1 << "\n";
    }
  }

  // number of faces
  outFile << m_localFaceNodeConn.size() << "\n";

  // write connectivity for each face
  for (CFuint iFace = 0; iFace < m_localFaceNodeConn.size(); ++iFace)
  {
    // number of nodes in this face
    outFile << m_localFaceNodeConn[iFace].size() << "\n";

    // write face node connectivity
    for (CFuint iNode = 0; iNode < m_localFaceNodeConn[iFace].size(); ++iNode)
    {
      outFile << m_localFaceNodeConn[iFace][iNode] << " ";
    }
    outFile << "\n";
  }

  //close the file
  outFile.close();
}

//////////////////////////////////////////////////////////////////////

void TetraSpectralFVElementData::computeLocalFaceNormals()
{
  CFAUTOTRACE;

  // get number of local internal faces
  const CFuint nbrLocalFaces    = m_localFaceNodeConn.size();
  const CFuint nbrLocalExtFaces = 2*(m_polyOrder+1)*(m_polyOrder+2);
  const CFuint nbrLocalIntFaces = nbrLocalFaces - nbrLocalExtFaces;

//   // compute local internal face normals
//   m_intFaceQuadPntNorm.resize(nbrLocalIntFaces);
//   for (CFuint iIntFace = 0; iIntFace < nbrLocalIntFaces; ++iIntFace)
//   {
// //     CF_DEBUG_OBJ(iIntFace);
//     // number of quadrature points on this face
//     const CFuint nbrQPnts = m_intFaceQuadWheights[iIntFace].size();
//     m_intFaceQuadPntNorm[iIntFace].resize(nbrQPnts,RealVector(3));
//
//     // number of face nodes
//     const CFuint nbrFaceNodes = m_localFaceNodeConn[iIntFace+nbrLocalExtFaces].size();
//
//     // vector containing face node coordinates
//     vector< RealVector > faceNodes;
//     faceNodes.resize(nbrFaceNodes);
//     for (CFuint iNode = 0; iNode < nbrFaceNodes; ++iNode)
//     {
//       const CFuint nodeID = m_localFaceNodeConn[iIntFace+nbrLocalExtFaces][iNode];
//       faceNodes[iNode].resize(m_dimensionality);
//       faceNodes[iNode] = m_localNodeCoord[nodeID];
//     }
//
//     // compute normal
//     RealVector normal(3);
//     normal = computePolygonNormal(faceNodes);
//     normal /= normal.norm2(); // divide by size to create a unit normal
//
//     for (CFuint iQPnt = 0; iQPnt < nbrQPnts; ++iQPnt)
//     {
//       m_intFaceQuadPntNorm[iIntFace][iQPnt] = normal;
// //       CF_DEBUG_OBJ(m_intFaceQuadPntNorm[iIntFace][iQPnt]);
//     }
//   }

  // set the dimensionality and order of the tensor product integrator
  m_tpIntegrator.setDimensionality(static_cast<CFDim>(m_dimensionality-1));
  m_tpIntegrator.setIntegratorOrder(m_polyOrder);

  // compute local internal face normals
  m_intFaceQuadPntNorm.resize(nbrLocalIntFaces);
  for (CFuint iIntFace = 0; iIntFace < nbrLocalIntFaces; ++iIntFace)
  {
//     CF_DEBUG_OBJ(iIntFace);
    // number of face nodes
    const CFuint nbrFaceNodes = m_localFaceNodeConn[iIntFace+nbrLocalExtFaces].size();

    // vector containing face node coordinates
    vector< RealVector > faceNodes(nbrFaceNodes,RealVector(m_dimensionality));
    faceNodes.resize(nbrFaceNodes);
    for (CFuint iNode = 0; iNode < nbrFaceNodes; ++iNode)
    {
      const CFuint nodeID = m_localFaceNodeConn[iIntFace+nbrLocalExtFaces][iNode];
      faceNodes[iNode] = m_localNodeCoord[nodeID];
    }

    // number of quadrilaterals in CV face
    vector< vector< RealVector > > quadDecomp;
    setPolyGonQuadrilateralDecomposition(faceNodes,quadDecomp);
    const CFuint nbrQuads = quadDecomp.size();

    // get face local normals at quadrature points
    for (CFuint iQuad = 0; iQuad < nbrQuads; ++iQuad)
    {
      // get quadrilateral quadrature nodes and wheights
      vector< RealVector > qNodeNormalQuad = m_tpIntegrator.getQuadPntsUnitNormalsPlus1D(quadDecomp[iQuad]);

      // add simplex quadrature nodes and wheights to global list
      m_intFaceQuadPntNorm[iIntFace]
          .insert(m_intFaceQuadPntNorm[iIntFace].end(),qNodeNormalQuad.begin(),qNodeNormalQuad.end());
    }

//     const CFuint nbrQPnts = m_intFaceQuadPntNorm[iIntFace].size();
//     for (CFuint iQPnt = 0; iQPnt < nbrQPnts; ++iQPnt)
//     {
//       CF_DEBUG_OBJ(m_intFaceQuadPntNorm[iIntFace][iQPnt]);
//     }
  }
//   CF_DEBUG_EXIT;
}

//////////////////////////////////////////////////////////////////////

void TetraSpectralFVElementData::computeExtFaceLocalNormals()
{
  CFAUTOTRACE;

  // number of SV faces
  const CFuint nbrSVFaces = 4;

  // compute the normals
  m_extFaceLocalNorm.resize(nbrSVFaces);
  for (CFuint iSVFace = 0; iSVFace < nbrSVFaces; ++iSVFace)
  {
    // resize the variable
    m_extFaceLocalNorm[iSVFace].resize(m_dimensionality);

    // compute normal
    m_extFaceLocalNorm[iSVFace] = computeTriangleNormal(m_svFaceNodeCoords[iSVFace][0],
                                                        m_svFaceNodeCoords[iSVFace][1],
                                                        m_svFaceNodeCoords[iSVFace][2]);

//     CF_DEBUG_OBJ(m_extFaceLocalNorm[iSVFace]);
  }
}

//////////////////////////////////////////////////////////////////////

void TetraSpectralFVElementData::createLocalIntFaceCVConn()
{
  CFAUTOTRACE;

  // get number of CVs
  const CFuint nbrCVs = getNbrOfCVs();

  // get number of local internal faces
  const CFuint nbrLocalFaces    = m_localFaceNodeConn.size();
  const CFuint nbrLocalExtFaces = 2*(m_polyOrder+1)*(m_polyOrder+2);
  const CFuint nbrLocalIntFaces = nbrLocalFaces - nbrLocalExtFaces;

  // create the local internal face-CV connectivity
  m_localIntFaceCVConn.resize(nbrLocalIntFaces);
  for (CFuint iFace = 0; iFace < nbrLocalIntFaces; ++iFace)
  {
    // face ID
    const CFuint faceID = iFace+nbrLocalExtFaces;

    // number of face nodes
    const CFuint nbrFaceNodes = m_localFaceNodeConn[faceID].size();

    // resize
    m_localIntFaceCVConn[iFace].resize(2);

    // loop over CVs, to find face-neighbours
    CFuint iNeighbour = 0;
    for (CFuint iCV = 0; iCV < nbrCVs && iNeighbour < 2; ++iCV)
    {
      // number of CV nodes
      const CFuint nbrCVNodes = m_localCVNodeConn[iCV].size();

      // check if CV is face-neigbour
      bool isNeighbour = true;
      for (CFuint iFaceNode = 0; iFaceNode < nbrFaceNodes && isNeighbour; ++iFaceNode)
      {
        bool foundNode = false;
        for (CFuint iCVNode = 0; iCVNode < nbrCVNodes && !foundNode; ++iCVNode)
        {
          foundNode = (m_localFaceNodeConn[faceID][iFaceNode] == m_localCVNodeConn[iCV][iCVNode]);
        }
        isNeighbour = foundNode;
      }

      if (isNeighbour)
      {
        m_localIntFaceCVConn[iFace][iNeighbour] = iCV;
        ++iNeighbour;
      }
    }

    cf_assert(iNeighbour == 2);

    // correct connectivity according to face orientation
    // (normal points into the RIGHT CV)
    // compute face ~center and average face normal
    RealVector faceCenter(0.0,m_dimensionality);
    vector< RealVector > faceNodes(nbrFaceNodes,RealVector(3));
    for (CFuint iFaceNode = 0; iFaceNode < nbrFaceNodes; ++iFaceNode)
    {
      const CFuint nodeID = m_localFaceNodeConn[faceID][iFaceNode];
      faceCenter += m_localNodeCoord[nodeID];
      faceNodes[iFaceNode] = m_localNodeCoord[nodeID];
    }
    faceCenter /= nbrFaceNodes;
    RealVector avgNormal = computePolygonNormal(faceNodes);

    // computer cell center
    const CFuint cvID = m_localIntFaceCVConn[iFace][LEFT];
    const CFuint nbrCVNodes = m_localCVNodeConn[cvID].size();
    RealVector cvCenter(0.0,m_dimensionality);
    for (CFuint iCVNode = 0; iCVNode < nbrCVNodes; ++iCVNode)
    {
      const CFuint nodeID = m_localCVNodeConn[cvID][iCVNode];
      cvCenter += m_localNodeCoord[nodeID];
    }
    cvCenter /= nbrCVNodes;

    // check orientation
    RealVector face2CV = cvCenter-faceCenter;
    if (face2CV[XX]*avgNormal[XX] + face2CV[YY]*avgNormal[YY] + face2CV[ZZ]*avgNormal[ZZ] > 0.0)
    {
      CFuint swap = m_localIntFaceCVConn[iFace][LEFT ];
      m_localIntFaceCVConn[iFace][LEFT ] = m_localIntFaceCVConn[iFace][RIGHT];
      m_localIntFaceCVConn[iFace][RIGHT] = swap;
    }
  }

/*  for (CFuint iFace = 0; iFace < m_localIntFaceCVConn.size(); ++iFace)
  {
    CFreal normSize = 0.0;
    for (CFuint iCoor = 0; iCoor < m_dimensionality; ++iCoor)
    {
      normSize += m_intFaceQuadPntNorm[iFace][iCoor]*m_intFaceQuadPntNorm[iFace][iCoor];
    }
    normSize = sqrt(normSize);
CF_DEBUG_OBJ(iFace);
CF_DEBUG_OBJ(m_intFaceQuadPntNorm[iFace]);
CF_DEBUG_OBJ(normSize);
CF_DEBUG_OBJ(m_localIntFaceCVConn[iFace][LEFT ]);
CF_DEBUG_OBJ(m_localIntFaceCVConn[iFace][RIGHT]);
  }*/
//   CF_DEBUG_EXIT;
}

//////////////////////////////////////////////////////////////////////

void TetraSpectralFVElementData::computeVolumeFractionsOfCVs()
{
  CFAUTOTRACE;

  // number of CVs in SV
  const CFuint nbrCVs = getNbrOfCVs();

  // resize the vector containing the CV face fractions
  m_volFracCV.resize(nbrCVs);

  // setup integrator
  m_tpIntegrator.setDimensionality(m_dimensionality);
  m_tpIntegrator.setIntegratorOrder(CFPolyOrder::ORDER0);

  // compute CV volume fractions
//   CFreal sum = 0.0;
  for (CFuint iCV = 0; iCV < nbrCVs; ++iCV)
  {
    const CFuint nbrNodesInCV = m_localCVNodeConn[iCV].size();

    // Create vector for node coordinates of CV
    vector< RealVector > cvNodesCoord(nbrNodesInCV);
    for (CFuint iNode = 0; iNode < nbrNodesInCV; ++iNode)
    {
      cvNodesCoord[iNode].resize(m_dimensionality);
      const CFuint localNodeID = m_localCVNodeConn[iCV][iNode];
      cvNodesCoord[iNode] = m_localNodeCoord[localNodeID];
    }

//     // Compute the volume fraction of the CV
//     m_volFracCV[iCV] = computePolyHedronVolume(cvNodesCoord)/(1.0/6.0);
//     CF_DEBUG_OBJ(m_volFracCV[iCV]);

    // Compute the volume fraction of the CV
    // (based on the tensorproductguassintegrator, to be able to handle curved CVs)
    m_volFracCV[iCV] = 0.0;
    vector< vector< RealVector > > boxDecomp;
    setPolyHedronBoxDecomposition(cvNodesCoord,boxDecomp);
    const CFuint nbrBoxes = boxDecomp.size();
    for (CFuint iBox = 0; iBox < nbrBoxes; ++iBox)
    {
      vector< CFreal > quadWheights = m_tpIntegrator.getQuadPntsWheights(boxDecomp[iBox]);
      const CFuint nbrQPnts = quadWheights.size();
      for (CFuint iQPnt = 0; iQPnt < nbrQPnts; ++iQPnt)
      {
        m_volFracCV[iCV] += quadWheights[iQPnt];
      }
    }
    m_volFracCV[iCV] *= 6.0;
//     CF_DEBUG_OBJ(m_volFracCV[iCV]);
//     sum += m_volFracCV[iCV];
    if (m_volFracCV[iCV] < 0.0)
    {
      throw BadFormatException (FromHere(),"Invalid SV partition, results in CVs with negative volume...");
    }
  }
//   CF_DEBUG_OBJ(sum);
}

//////////////////////////////////////////////////////////////////////

void TetraSpectralFVElementData::createLocalExtFaceCVConn()
{
  CFAUTOTRACE;

  // number of external faces at SV face
  const CFuint nbrCVsAtSVFace = (m_polyOrder + 1)*(m_polyOrder + 2)/2;

  // total number of external faces
  const CFuint nbrExtFaces = 4*nbrCVsAtSVFace;

  // resize m_localExtFaceCVConn
  m_localExtFaceCVConn.resize(nbrExtFaces);

  // create local external face - CV connectivity
  CFuint extFaceID = 0;

  // first SV face
  CFuint iFace = 0;
  for (CFuint iCVLayer = 0; iFace < nbrCVsAtSVFace; ++iCVLayer)
  {
    // number of CVs at face in this layer
    const CFuint nbrFaceCVsInLayer = iCVLayer+1;

    // set neighbouring CVs
    CFuint cvIdx = iCVLayer*(iCVLayer+1)*(iCVLayer+2)/6;
    for (CFuint iFaceCV = 0; iFaceCV < nbrFaceCVsInLayer; ++iFaceCV, ++iFace, ++extFaceID)
    {
      m_localExtFaceCVConn[extFaceID] = cvIdx; ++cvIdx;
    }
  }

  // second SV face
  iFace = 0;
  for (CFuint iCVLayer = 0; iFace < nbrCVsAtSVFace; ++iCVLayer)
  {
    // number of CVs at face in this layer
    const CFuint nbrFaceCVsInLayer = iCVLayer+1;

    // set neighbouring CVs
    CFuint cvIdx = (iCVLayer+1)*(iCVLayer+2)*(iCVLayer+3)/6 - 1;
    for (CFuint iFaceCV = 0; iFaceCV < nbrFaceCVsInLayer; ++iFaceCV, ++iFace, ++extFaceID)
    {
      m_localExtFaceCVConn[extFaceID] = cvIdx; cvIdx -= iFaceCV + 2;
    }
  }

  // third SV face
  // CV rows in other direction
  vector< CFuint > cvRow2(m_polyOrder+1);
  CFuint cvID = m_polyOrder*(m_polyOrder+1)*(m_polyOrder+2)/6;
  const CFuint polyOrderP1 = m_polyOrder+1;
  for (CFuint iRow2 = 0; iRow2 < polyOrderP1; ++iRow2)
  {
    cvRow2[iRow2] = cvID; cvID += polyOrderP1-iRow2;
  }

  for (CFuint iCVRow = 0; iCVRow < polyOrderP1; ++iCVRow)
  {
    for (CFuint iFaceCV = 0; iFaceCV < iCVRow+1; ++iFaceCV, ++extFaceID)
    {
      m_localExtFaceCVConn[extFaceID] = cvRow2[iCVRow-iFaceCV]; ++cvRow2[iCVRow-iFaceCV];
    }
  }

  // fourth SV face
  iFace = 0;
  for (CFuint iCVLayer = 0; iFace < nbrCVsAtSVFace; ++iCVLayer)
  {
    // number of CVs at face in this layer
    const CFuint nbrFaceCVsInLayer = iCVLayer+1;

    // set neighbouring CVs
    CFuint cvIdx = iCVLayer*(iCVLayer+1)*(iCVLayer+2)/6 + nbrFaceCVsInLayer - 1;
    for (CFuint iFaceCV = 0; iFaceCV < nbrFaceCVsInLayer; ++iFaceCV, ++iFace, ++extFaceID)
    {
      m_localExtFaceCVConn[extFaceID] = cvIdx; cvIdx += nbrFaceCVsInLayer - 1 - iFaceCV;
    }
  }

/*  for (CFuint iFace = 0; iFace < m_localExtFaceCVConn.size(); ++iFace)
  {
CF_DEBUG_OBJ(m_localExtFaceCVConn[iFace]);
  }
  CF_DEBUG_EXIT;*/
}

//////////////////////////////////////////////////////////////////////

void TetraSpectralFVElementData::createSVFaceLocalExtFaceConn()
{
  CFAUTOTRACE;

  // number of external faces at SV face
  const CFuint nbrCVsAtSVFace = (m_polyOrder + 1)*(m_polyOrder + 2)/2;

  // create SV face - local external face connectivity
  CFuint extFaceID = 0;
  m_svFaceLocalExtFaceConn.resize(4);
  for (CFuint iSVFace = 0; iSVFace < 4; ++iSVFace)
  {
// CF_DEBUG_OBJ(iSVFace);
    m_svFaceLocalExtFaceConn[iSVFace].resize(nbrCVsAtSVFace);
    for (CFuint iExtFace = 0; iExtFace < nbrCVsAtSVFace; ++iExtFace, ++extFaceID)
    {
      m_svFaceLocalExtFaceConn[iSVFace][iExtFace] = extFaceID;
// CF_DEBUG_OBJ(m_svFaceLocalExtFaceConn[iSVFace][iExtFace]);
    }
  }
//   CF_DEBUG_EXIT;
}

//////////////////////////////////////////////////////////////////////

void TetraSpectralFVElementData::createExtFaceNodeLocalCoords()
{
  CFAUTOTRACE;

  // number of CVs at SV face
  const CFuint nbrCVsAtSVFace = (m_polyOrder+1)*(m_polyOrder+2)/2;

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
      // !!! KSI and ETA axis switched between cell coordinate system and face coordinate system,
      // to ensure outward pointing normal
      m_extFaceNodeLocalCoords[iExtFace][iNode][KSI] = m_localNodeCoord[nodeID][ETA];
      m_extFaceNodeLocalCoords[iExtFace][iNode][ETA] = m_localNodeCoord[nodeID][KSI];
    }
  }
}

//////////////////////////////////////////////////////////////////////

void TetraSpectralFVElementData::createExtFaceCVConnPerOrient()
{
  CFAUTOTRACE;

  // get the number of SV faces
  const CFuint nbrCellFaces = m_svFaceLocalExtFaceConn.size();

  // get the number of CVs lying at a SV face
  cf_assert(m_svFaceLocalExtFaceConn.size() > 0);
  const CFuint nbLocFacesPerCellFace = m_svFaceLocalExtFaceConn[0].size();

  // number CV rows
  const CFuint nbrCVRows = m_polyOrder + 1;

  // face indexes for right face, inverted order
  vector< CFuint > invFaceIdxes(nbLocFacesPerCellFace);
  for (CFuint iCVRow = 0; iCVRow < nbrCVRows; ++iCVRow)
  {
    // index of first face in this row
    const CFuint firstFaceIdx = iCVRow*(iCVRow+1)/2;

    // first put the faces
    for (CFuint iFace = 0; iFace < iCVRow+1; ++iFace)
    {
      invFaceIdxes[firstFaceIdx+iFace] = firstFaceIdx+iFace;
    }

    // then switch the outlying faces
    for (CFuint iFace = 0; iFace < (iCVRow+1)/2; ++iFace)
    {
      CFuint swap = invFaceIdxes[firstFaceIdx+iFace];
      invFaceIdxes[firstFaceIdx+iFace] = invFaceIdxes[firstFaceIdx+iCVRow-iFace];
      invFaceIdxes[firstFaceIdx+iCVRow-iFace] = swap;
    }
  }

//   for (CFuint iFace = 0; iFace < invFaceIdxes.size(); ++iFace)
//   {
// CF_DEBUG_OBJ(invFaceIdxes[iFace]);
//   }

  // arrange faces in rows with three different orientations
  vector< vector< vector< CFuint > > > faceRowArrang(3);

  // first row arrangement
  {
    faceRowArrang[0].resize(nbrCVRows);
    CFuint faceIdx = 0;
    for (CFuint iRow = 0; iRow < nbrCVRows; ++iRow)
    {
      faceRowArrang[0][iRow].resize(iRow+1);
      for (CFuint iFace = 0; iFace < iRow+1; ++iFace, ++faceIdx)
      {
        faceRowArrang[0][iRow][iFace] = faceIdx;
      }
    }
  }

  // second row arrangement
  {
    faceRowArrang[1].resize(nbrCVRows);
    for (CFuint iRow = 0; iRow < nbrCVRows; ++iRow)
    {
      faceRowArrang[1][iRow].resize(iRow+1);
      for (CFuint iFace = 0; iFace < iRow+1; ++iFace)
      {
        faceRowArrang[1][iRow][iFace] = faceRowArrang[0][nbrCVRows-1-iRow+iFace][nbrCVRows-1-iRow];
      }
    }
  }

  // third row arrangement
  {
    faceRowArrang[2].resize(nbrCVRows);
    for (CFuint iRow = 0; iRow < nbrCVRows; ++iRow)
    {
      faceRowArrang[2][iRow].resize(iRow+1);
      for (CFuint iFace = 0; iFace < iRow+1; ++iFace)
      {
        faceRowArrang[2][iRow][iFace] = faceRowArrang[1][nbrCVRows-1-iRow+iFace][nbrCVRows-1-iRow];
      }
    }
  }

//   for (CFuint iRot = 0; iRot < 3; ++iRot)
//   {
// CF_DEBUG_OBJ(iRot);
//     for (CFuint iRow = 0; iRow < faceRowArrang[iRot].size(); ++iRow)
//     {
// CF_DEBUG_OBJ(iRow);
//       for (CFuint iFace = 0; iFace < faceRowArrang[iRot][iRow].size(); ++iFace)
//       {
// CF_DEBUG_OBJ(faceRowArrang[iRot][iRow][iFace]);
//       }
//     }
//   }

  // number of rotable CV face groups
  const CFuint nbrRotCVFaceGroups = nbLocFacesPerCellFace/3;

  // storage of face groups
  vector< vector< CFuint > > rotCVFaceIdxs(nbrRotCVFaceGroups);
  for (CFuint iRotGroup = 0; iRotGroup < nbrRotCVFaceGroups; ++iRotGroup)
  {
    rotCVFaceIdxs[iRotGroup].resize(3);
  }

  // assign face groups
  CFuint faceIdx = 0;
  for (CFuint iRow = 0; faceIdx < nbrRotCVFaceGroups; ++iRow)
  {
    for (CFuint iFace = 0; iFace < iRow+1 && faceIdx < nbrRotCVFaceGroups; ++iFace, ++faceIdx)
    {
      for (CFuint iRot = 0; iRot < 3; ++iRot)
      {
        rotCVFaceIdxs[faceIdx][iRot] = faceRowArrang[iRot][iRow][iFace];
      }
    }
  }

  // number of CV connectivities through SV faces
  const CFuint nbrOrients = 30;

  // fill m_extSVFaceCVConnPerOrient and m_eSVFaceLocESVFaceConnPerOrient
  CFuint iOrient = 0;
  m_extSVFaceCVConnPerOrient.resize(nbrOrients);
  m_eSVFaceLocESVFaceConnPerOrient.resize(nbrOrients);
  for (CFuint iSVFaceL = 0; iSVFaceL < nbrCellFaces; ++iSVFaceL)
  {
    vector< CFuint > extSVFaceCVConnL = m_extSVFaceCVConn[iSVFaceL];
    for (CFuint iSVFaceR = iSVFaceL; iSVFaceR < nbrCellFaces; ++iSVFaceR)
    {
      vector< CFuint > extSVFaceCVConnR = m_extSVFaceCVConn[iSVFaceR];
      for (CFuint iRot = 0; iRot < 3; ++iRot, ++iOrient)
      {
        // resize m_eSVFaceLocESVFaceConnPerOrient[iOrient]
        m_eSVFaceLocESVFaceConnPerOrient[iOrient].resize(2);

        // set local SV face idxs
        m_eSVFaceLocESVFaceConnPerOrient[iOrient][LEFT ] = iSVFaceL;
        m_eSVFaceLocESVFaceConnPerOrient[iOrient][RIGHT] = iSVFaceR;

        // resize m_extSVFaceCVConnPerOrient[iOrient]
        m_extSVFaceCVConnPerOrient[iOrient].resize(nbLocFacesPerCellFace);

        // create external face CV connectivity
        for (CFuint iExtFace = 0; iExtFace < nbLocFacesPerCellFace; ++iExtFace)
        {
          // resize the connectivity
          m_extSVFaceCVConnPerOrient[iOrient][iExtFace].resize(2);

          // put neighbouring cv IDs in connectivity
          m_extSVFaceCVConnPerOrient[iOrient][iExtFace][LEFT ] = extSVFaceCVConnL[iExtFace];
          m_extSVFaceCVConnPerOrient[iOrient][iExtFace][RIGHT] = extSVFaceCVConnR[invFaceIdxes[iExtFace]];
        }

        // rotate the right face
        for (CFuint iRotGroup = 0; iRotGroup < nbrRotCVFaceGroups; ++iRotGroup)
        {
          // indexes of faces to be rotated
          const CFuint face0Idx = rotCVFaceIdxs[iRotGroup][0];
          const CFuint face1Idx = rotCVFaceIdxs[iRotGroup][1];
          const CFuint face2Idx = rotCVFaceIdxs[iRotGroup][2];

          // rotate faces
          const CFuint swap = invFaceIdxes[face2Idx];
          invFaceIdxes[face2Idx] = invFaceIdxes[face1Idx];
          invFaceIdxes[face1Idx] = invFaceIdxes[face0Idx];
          invFaceIdxes[face0Idx] = swap;
        }
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////

void TetraSpectralFVElementData::createSVFaceNodeConnectivity()
{
  CFAUTOTRACE;

  // number of SV faces
  const CFuint nbrSVFaces = 4;

  // resize the variable
  m_svFaceNodeConn.resize(nbrSVFaces);
  for (CFuint iSVFace = 0; iSVFace < nbrSVFaces; ++iSVFace)
  {
    m_svFaceNodeConn[iSVFace].resize(3);
  }

  // fill the variable
  // first SV face
  m_svFaceNodeConn[0][0] = 0;
  m_svFaceNodeConn[0][1] = 2;
  m_svFaceNodeConn[0][2] = 1;

  // second SV face
  m_svFaceNodeConn[1][0] = 0;
  m_svFaceNodeConn[1][1] = 1;
  m_svFaceNodeConn[1][2] = 3;

  // third SV face
  m_svFaceNodeConn[2][0] = 1;
  m_svFaceNodeConn[2][1] = 2;
  m_svFaceNodeConn[2][2] = 3;

  // fourth SV face
  m_svFaceNodeConn[3][0] = 0;
  m_svFaceNodeConn[3][1] = 3;
  m_svFaceNodeConn[3][2] = 2;
}

//////////////////////////////////////////////////////////////////////

void TetraSpectralFVElementData::createSVFaceNodeConnectivityPerOrient()
{
  CFAUTOTRACE;

  // number of SV faces
  const CFuint nbrSVFaces = m_svFaceNodeConn.size();

  // number of possible orientations
  const CFuint nbrOrient = 30;

  // resize the variable
  m_svFaceNodeConnPerOrientation.resize(nbrOrient);
  for (CFuint iOrient = 0; iOrient < nbrOrient; ++iOrient)
  {
    m_svFaceNodeConnPerOrientation[iOrient].resize(2);
    for (CFuint iCell = 0; iCell < 2; ++iCell)
      m_svFaceNodeConnPerOrientation[iOrient][iCell].resize(3);
  }

  // fill the variable
  CFuint iOrient = 0;
  for (CFuint iSVFaceL = 0; iSVFaceL < nbrSVFaces; ++iSVFaceL)
  {
    vector< CFuint > faceNodesL = m_svFaceNodeConn[iSVFaceL];
    for (CFuint iSVFaceR = iSVFaceL; iSVFaceR < nbrSVFaces; ++iSVFaceR)
    {
      vector< CFuint > faceNodesR(3);
      faceNodesR[0] = m_svFaceNodeConn[iSVFaceR][0];
      faceNodesR[1] = m_svFaceNodeConn[iSVFaceR][2];
      faceNodesR[2] = m_svFaceNodeConn[iSVFaceR][1];

      for (CFuint iRot = 0; iRot < 3; ++iRot, ++iOrient)
      {
        for (CFuint iNode = 0; iNode < 3; ++iNode)
        {
          m_svFaceNodeConnPerOrientation[iOrient][LEFT ][iNode] = faceNodesL[iNode];
          m_svFaceNodeConnPerOrientation[iOrient][RIGHT][iNode] = faceNodesR[iNode];
        }

// CF_DEBUG_OBJ(iOrient);
//         for (CFuint iNode = 0; iNode < 3; ++iNode)
//         {
// CF_DEBUG_OBJ(m_svFaceNodeConnPerOrientation[iOrient][LEFT ][iNode]);
// CF_DEBUG_OBJ(m_extSVFaceCVConnPerOrient[iOrient][iNode][LEFT ]);
//         }
//
//         for (CFuint iNode = 0; iNode < 3; ++iNode)
//         {
// CF_DEBUG_OBJ(m_svFaceNodeConnPerOrientation[iOrient][RIGHT][iNode]);
// CF_DEBUG_OBJ(m_extSVFaceCVConnPerOrient[iOrient][iNode][RIGHT]);
//         }

        // rotate nodes of right face to new orientation
        CFuint swap = faceNodesR[2];
        faceNodesR[2] = faceNodesR[1];
        faceNodesR[1] = faceNodesR[0];
        faceNodesR[0] = swap;
      }

      cf_assert(faceNodesR[0] == m_svFaceNodeConn[iSVFaceR][0]);
      cf_assert(faceNodesR[1] == m_svFaceNodeConn[iSVFaceR][2]);
      cf_assert(faceNodesR[2] == m_svFaceNodeConn[iSVFaceR][1]);
    }
  }
}

//////////////////////////////////////////////////////////////////////

void TetraSpectralFVElementData::createSVFaceNodeConnectivityPerOrientNoSymm()
{
  CFAUTOTRACE;

  // number of SV faces
  const CFuint nbrSVFaces = m_svFaceNodeConn.size();

  // number of possible orientations
  const CFuint nbrOrient = 48;

  // resize the variable
  m_svFaceNodeConnPerOrientationNoSymm.resize(nbrOrient);
  for (CFuint iOrient = 0; iOrient < nbrOrient; ++iOrient)
  {
    m_svFaceNodeConnPerOrientationNoSymm[iOrient].resize(2);
    for (CFuint iCell = 0; iCell < 2; ++iCell)
      m_svFaceNodeConnPerOrientationNoSymm[iOrient][iCell].resize(3);
  }

  // fill the variable
  CFuint iOrient = 0;
  for (CFuint iSVFaceL = 0; iSVFaceL < nbrSVFaces; ++iSVFaceL)
  {
    vector< CFuint > faceNodesL = m_svFaceNodeConn[iSVFaceL];
    for (CFuint iSVFaceR = 0; iSVFaceR < nbrSVFaces; ++iSVFaceR)
    {
      vector< CFuint > faceNodesR(3);
      faceNodesR[0] = m_svFaceNodeConn[iSVFaceR][0];
      faceNodesR[1] = m_svFaceNodeConn[iSVFaceR][2];
      faceNodesR[2] = m_svFaceNodeConn[iSVFaceR][1];

      for (CFuint iRot = 0; iRot < 3; ++iRot, ++iOrient)
      {
        for (CFuint iNode = 0; iNode < 3; ++iNode)
        {
          m_svFaceNodeConnPerOrientationNoSymm[iOrient][LEFT ][iNode] = faceNodesL[iNode];
          m_svFaceNodeConnPerOrientationNoSymm[iOrient][RIGHT][iNode] = faceNodesR[iNode];
        }

// CF_DEBUG_OBJ(iOrient);
//         for (CFuint iNode = 0; iNode < 3; ++iNode)
//         {
// CF_DEBUG_OBJ(m_svFaceNodeConnPerOrientation[iOrient][LEFT ][iNode]);
// CF_DEBUG_OBJ(m_extSVFaceCVConnPerOrient[iOrient][iNode][LEFT ]);
//         }
        //
//         for (CFuint iNode = 0; iNode < 3; ++iNode)
//         {
// CF_DEBUG_OBJ(m_svFaceNodeConnPerOrientation[iOrient][RIGHT][iNode]);
// CF_DEBUG_OBJ(m_extSVFaceCVConnPerOrient[iOrient][iNode][RIGHT]);
//         }

        // rotate nodes of right face to new orientation
        CFuint swap = faceNodesR[2];
        faceNodesR[2] = faceNodesR[1];
        faceNodesR[1] = faceNodesR[0];
        faceNodesR[0] = swap;
      }

      cf_assert(faceNodesR[0] == m_svFaceNodeConn[iSVFaceR][0]);
      cf_assert(faceNodesR[1] == m_svFaceNodeConn[iSVFaceR][2]);
      cf_assert(faceNodesR[2] == m_svFaceNodeConn[iSVFaceR][1]);
    }
  }
}

//////////////////////////////////////////////////////////////////////

void TetraSpectralFVElementData::computeFaceFractionsOfCVs()
{
  CFAUTOTRACE;

  // number of SV faces
  const CFuint nbrSVFaces = 4;

  // number of CVs at a SV face
  const CFuint nbrCVsAtSVFace = (m_polyOrder + 1)*(m_polyOrder + 2)/2;

  // resize the vector containing the CV face fractions
  m_faceFracCV.resize(nbrSVFaces);
  for (CFuint iFace = 0; iFace < nbrSVFaces; ++iFace)
  {
    m_faceFracCV[iFace].resize(nbrCVsAtSVFace);
  }

  // fill in CV fractions first SV face (first nbrCVsAtSVFace faces in m_localFaceNodeConn make up one SV face)
  for (CFuint iCV = 0; iCV < nbrCVsAtSVFace; ++iCV)
  {
    // number of nodes in CV face
    const CFuint nbrNodes = m_localFaceNodeConn[iCV].size();

    // create vector containing the CV face node coordinates
    vector< RealVector > nodeCoords(nbrNodes);
    for (CFuint iNode = 0; iNode < nbrNodes; ++iNode)
    {
      // node ID
      const CFuint nodeID = m_localFaceNodeConn[iCV][iNode];

      // set node coordinates
      nodeCoords[iNode].resize(m_dimensionality);
      nodeCoords[iNode] = m_localNodeCoord[nodeID];
    }

    // compute face surface fraction
    m_faceFracCV[0][iCV] = computePolygonSurface(nodeCoords)/0.5;
//     CF_DEBUG_OBJ(m_faceFracCV[0][iCV]);
  }

  // fill in CV fractions of other SV faces
  for (CFuint iFace = 1; iFace < nbrSVFaces; ++iFace)
  {
    for (CFuint iCV = 0; iCV < nbrCVsAtSVFace; ++iCV)
    {
      m_faceFracCV[iFace][iCV] = m_faceFracCV[0][iCV];
    }
  }
//   CF_DEBUG_EXIT;
}

//////////////////////////////////////////////////////////////////////

void TetraSpectralFVElementData::createPolyExponents()
{
  CFAUTOTRACE;

  // number of polynomial terms
  const CFuint nbrPolyTerms = (m_polyOrder + 1)*(m_polyOrder + 2)*(m_polyOrder + 3)/6;

  // resize the variable
  m_polyExponents.resize(nbrPolyTerms);
  for (CFuint iTerm = 0; iTerm < nbrPolyTerms; ++iTerm)
  {
    m_polyExponents[iTerm].resize(3);
  }

  // define exponents
  CFuint iTerm = 0;
  const CFuint polyOrderP1 = m_polyOrder+1;
  for (CFuint iXYZ = 0; iXYZ < polyOrderP1; ++iXYZ)
  {
    for (CFuint iX = 0; iX < iXYZ+1; ++iX)
    {
      for (CFuint iY = 0; iY < iXYZ-iX+1; ++iY,++iTerm)
      {
        m_polyExponents[iTerm][KSI] = iX;
        m_polyExponents[iTerm][ETA] = iY;
        m_polyExponents[iTerm][ZTA] = iXYZ-iX-iY;
/*        CF_DEBUG_OBJ(iTerm);
        CF_DEBUG_OBJ(m_polyExponents[iTerm][KSI]);
        CF_DEBUG_OBJ(m_polyExponents[iTerm][ETA]);
        CF_DEBUG_OBJ(m_polyExponents[iTerm][ZTA]);*/
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////

void TetraSpectralFVElementData::computePolyCoefs()
{
  CFAUTOTRACE;

//   // set the dimensionality and order of the simplex integrator
//   m_sIntegrator.setDimensionality(m_dimensionality);
//   m_sIntegrator.setIntegratorOrder(m_polyOrder);

  // set the dimensionality and order of the tensor product integrator
  m_tpIntegrator.setDimensionality(m_dimensionality);
  m_tpIntegrator.setIntegratorOrder(m_polyOrder);

  // number of control volumes in a SV
  const CFuint nbrOfCVs = m_localCVNodeConn.size();

  // matrix for linear system
  RealMatrix polyCoefSystem(nbrOfCVs,nbrOfCVs);

  // fill the linear system matrix
  for (CFuint iCV = 0; iCV < nbrOfCVs; ++iCV)
  {
    // number of control volume nodes
    const CFuint nbrCVNodes = m_localCVNodeConn[iCV].size();

//     // CV decompositions into tetrahedrons
//     vector< vector< CFuint > > tetraDecomp;
//     setPolyHedronTetraDecomposition(nbrCVNodes,tetraDecomp);
//
//     // number of tetrahedrons in polyhedron
//     const CFuint nbrTetras = tetraDecomp.size();
//
//     // get CV node coordinates
//     vector< RealVector > cvNodeCoord(nbrCVNodes);
//     for (CFuint iNode = 0; iNode < nbrCVNodes; ++iNode)
//     {
//       cvNodeCoord[iNode].resize(m_dimensionality);
//       const CFuint nodeID = m_localCVNodeConn[iCV][iNode];
//       cvNodeCoord[iNode] = m_localNodeCoord[nodeID];
//     }

    // CV decompositions into boxes
    // first set the CV nodes
    vector< RealVector > polyHedronNodes(nbrCVNodes,RealVector(m_dimensionality));
    for (CFuint iNode = 0; iNode < nbrCVNodes; ++iNode)
    {
      const CFuint nodeID = m_localCVNodeConn[iCV][iNode];
      polyHedronNodes[iNode] = m_localNodeCoord[nodeID];
    }
    // get the coordinates of the boxes composing the CV
    vector< vector< RealVector > > boxDecomp;
    setPolyHedronBoxDecomposition(polyHedronNodes,boxDecomp);
    const CFuint nbrBoxes = boxDecomp.size();

//     // get quadrature nodes and wheights (using tetras)
//     vector< RealVector > qNodeCoord;
//     vector< CFreal > qWheights;
//     for (CFuint iTetra = 0; iTetra < nbrTetras; ++iTetra)
//     {
//       // get node coordinates of tetrahedron
//       vector< RealVector > tetraNodeCoord(4);
//       for (CFuint iNode = 0; iNode < 4; ++iNode)
//       {
//         tetraNodeCoord[iNode].resize(m_dimensionality);
//         tetraNodeCoord[iNode] = cvNodeCoord[tetraDecomp[iTetra][iNode]];
//       }
//
//       // get tetra quadrature nodes and wheights
//       vector< RealVector > qNodeCoordTetra = m_sIntegrator.getQuadPntsCoords  (tetraNodeCoord);
//       vector< CFreal >     qWheightsTetra  = m_sIntegrator.getQuadPntsWheights(tetraNodeCoord);
//
//       // add tetra quadrature nodes and wheights to global list
//       qNodeCoord.insert(qNodeCoord.end(),qNodeCoordTetra.begin(),qNodeCoordTetra.end());
//       qWheights .insert(qWheights .end(),qWheightsTetra .begin(),qWheightsTetra .end());
//     }

    // get quadrature nodes and wheights (using boxes)
    vector< RealVector > qNodeCoord;
    vector< CFreal > qWheights;
    for (CFuint iBox = 0; iBox < nbrBoxes; ++iBox)
    {
      // get box quadrature nodes and wheights
      vector< RealVector > qNodeCoordBox = m_tpIntegrator.getQuadPntsCoords  (boxDecomp[iBox]);
      vector< CFreal >     qWheightsBox  = m_tpIntegrator.getQuadPntsWheights(boxDecomp[iBox]);

      // add box quadrature nodes and wheights to global list
      qNodeCoord.insert(qNodeCoord.end(),qNodeCoordBox.begin(),qNodeCoordBox.end());
      qWheights .insert(qWheights .end(),qWheightsBox .begin(),qWheightsBox .end());
    }

    // compute left hand side of linear system
    const CFuint nbrQPnts = qWheights.size();
    for (CFuint iTerm = 0; iTerm < nbrOfCVs; ++iTerm)
    {
      polyCoefSystem(iCV,iTerm) = 0;

      for (CFuint iQPnt = 0; iQPnt < nbrQPnts; ++iQPnt)
      {
        polyCoefSystem(iCV,iTerm) += qWheights[iQPnt]
                                      *pow(qNodeCoord[iQPnt][KSI],m_polyExponents[iTerm][KSI])
                                      *pow(qNodeCoord[iQPnt][ETA],m_polyExponents[iTerm][ETA])
                                      *pow(qNodeCoord[iQPnt][ZTA],m_polyExponents[iTerm][ZTA]);
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
    rhsLinSys(iCV,iCV) = m_volFracCV[iCV]/6.0;
  }

  // multiply inverted linear system matrix with rhs
  RealMatrix polyCoefMatrix(nbrOfCVs,nbrOfCVs);
  polyCoefMatrix = invLinSysMatrix*rhsLinSys;

  // store polynomial coefficients in m_polyCoefSFV
  m_polyCoefSFV.resize(nbrOfCVs);
  for (CFuint iPoly = 0; iPoly < nbrOfCVs; ++iPoly)
  {
//     CF_DEBUG_OBJ(iPoly);
    m_polyCoefSFV[iPoly].resize(nbrOfCVs);
    for (CFuint iTerm = 0; iTerm < nbrOfCVs; ++iTerm)
    {
      m_polyCoefSFV[iPoly][iTerm] = polyCoefMatrix(iTerm,iPoly);
//       CF_DEBUG_OBJ(m_polyCoefSFV[iPoly][iTerm]);
    }
  }
//   CF_DEBUG_EXIT;
}

//////////////////////////////////////////////////////////////////////

void TetraSpectralFVElementData::computeIntFaceQuadPntsData()
{
  CFAUTOTRACE;

  // set the dimensionality and order of the tensor product integrator
  m_tpIntegrator.setDimensionality(static_cast<CFDim>(m_dimensionality-1));
  m_tpIntegrator.setIntegratorOrder(m_polyOrder);

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

    // number of quadrilaterals in CV face
    vector< vector< RealVector > > quadDecomp;
    setPolyGonQuadrilateralDecomposition(nodeCoord,quadDecomp);
    const CFuint nbrQuads = quadDecomp.size();

    // get face quadrature data
    for (CFuint iQuad = 0; iQuad < nbrQuads; ++iQuad)
    {
      // get quadrilateral quadrature nodes and wheights
      vector< RealVector > qNodeCoordQuad = m_tpIntegrator.getQuadPntsCoordsPlus1D  (quadDecomp[iQuad]);
      vector< CFreal >     qWheightsQuad  = m_tpIntegrator.getQuadPntsWheightsPlus1D(quadDecomp[iQuad]);

      // add quadrilateral quadrature nodes and wheights to global list
      m_intFaceQuadPntCoords[iIntFace]
          .insert(m_intFaceQuadPntCoords[iIntFace].end(),qNodeCoordQuad.begin(),qNodeCoordQuad.end());
      m_intFaceQuadWheights [iIntFace]
          .insert(m_intFaceQuadWheights [iIntFace].end(),qWheightsQuad .begin(),qWheightsQuad .end());
    }

    // The sum of these quadrature wheights should be one,
    // because the face surface is included in the normal vector!!!
//     const CFreal invQWheightsSum = 1.0/accumulate(m_intFaceQuadWheights[iIntFace].begin(),
//         m_intFaceQuadWheights[iIntFace].end()  ,0.0);
//     const CFuint nbrQPntsFace = m_intFaceQuadWheights[iIntFace].size();
//     for (CFuint iQPnt = 0; iQPnt < nbrQPntsFace; ++iQPnt)
//     {
//       m_intFaceQuadWheights[iIntFace][iQPnt] *= invQWheightsSum;
//     }

    // compute quadrature point polynomial values
    const CFuint nbrQPntsFace = m_intFaceQuadWheights[iIntFace].size();
    for (CFuint iQPnt = 0; iQPnt < nbrQPntsFace; ++iQPnt)
    {
      m_intFaceQuadPntPolyVals[iIntFace] = getSVPolyValsAtNode(m_intFaceQuadPntCoords[iIntFace]);
    }
  }
}

//////////////////////////////////////////////////////////////////////

void TetraSpectralFVElementData::createFaceFluxPntsConn()
{
 CFAUTOTRACE;

  // number of flux points at a face
  const CFuint nbrFaceFlxPnts = (m_polyOrder + 2)*(m_polyOrder + 3)/2;

  // resize m_faceFlxPntsConn (nbr of faces)
  m_faceFlxPntsConn.resize(4);

  // create connectivity
  // face 0
  m_faceFlxPntsConn[0].resize(nbrFaceFlxPnts);
  CFuint iFlx = 0;
  for (CFuint iLayer = 0; iFlx < nbrFaceFlxPnts; ++iLayer)
  {
    // number of flux points at face in this layer
    const CFuint nbrFlxPntsInLayer = iLayer+1;

    // set flux points in face
    CFuint flxID = iLayer*(iLayer+1)*(iLayer+2)/6;
    for (CFuint flxIdx = 0; flxIdx < nbrFlxPntsInLayer; ++flxIdx, ++iFlx)
    {
      m_faceFlxPntsConn[0][iFlx] = flxID; ++flxID;
    }
  }

  // face 1
  m_faceFlxPntsConn[1].resize(nbrFaceFlxPnts);
  iFlx = 0;
  for (CFuint iLayer = 0; iFlx < nbrFaceFlxPnts; ++iLayer)
  {
    // number of flux points at face in this layer
    const CFuint nbrFlxPntsInLayer = iLayer+1;

    // set flux points in face
    CFuint flxID = (iLayer+1)*(iLayer+2)*(iLayer+3)/6 - 1;
    for (CFuint flxIdx = 0; flxIdx < nbrFlxPntsInLayer; ++flxIdx, ++iFlx)
    {
      m_faceFlxPntsConn[1][iFlx] = flxID; flxID -= flxIdx + 2;
    }
  }

  // face 2
  // flux point rows in other direction
  const CFuint polyOrderP2 = m_polyOrder+2;
  vector< CFuint > elemRow2(polyOrderP2);
  CFuint flxID = (m_polyOrder+1)*(polyOrderP2)*(m_polyOrder+3)/6;
  for (CFuint iRow2 = 0; iRow2 < polyOrderP2; ++iRow2)
  {
      elemRow2[iRow2] = flxID; flxID += m_polyOrder+2-iRow2;
  }

  m_faceFlxPntsConn[2].resize(nbrFaceFlxPnts);
  iFlx = 0;
  for (CFuint iRow = 0; iRow < polyOrderP2; ++iRow)
  {
    for (CFuint flxIdx = 0; flxIdx < iRow+1; ++flxIdx, ++iFlx)
    {
      m_faceFlxPntsConn[2][iFlx] = elemRow2[iRow-flxIdx]; ++elemRow2[iRow-flxIdx];
    }
  }

  // face 3
  m_faceFlxPntsConn[3].resize(nbrFaceFlxPnts);
  iFlx = 0;
  for (CFuint iLayer = 0; iFlx < nbrFaceFlxPnts; ++iLayer)
  {
    // number of flux points at face in this layer
    const CFuint nbrFlxPntsInLayer = iLayer+1;

    // set flux points in face
    CFuint flxID = iLayer*(iLayer+1)*(iLayer+2)/6 + nbrFlxPntsInLayer - 1;
    for (CFuint flxIdx = 0; flxIdx < nbrFlxPntsInLayer; ++flxIdx, ++iFlx)
    {
      m_faceFlxPntsConn[3][iFlx] = flxID; flxID += nbrFlxPntsInLayer - 1 - flxIdx;
    }
  }

/*  for (CFuint iFace = 0; iFace < 4; ++iFace)
  {
    CF_DEBUG_OBJ(iFace);
    for (CFuint iFlx = 0; iFlx < nbrFaceFlxPnts; ++iFlx)
    {
      CF_DEBUG_OBJ(m_faceFlxPntsConn[iFace][iFlx]);
    }
  }*/
}

//////////////////////////////////////////////////////////////////////

void TetraSpectralFVElementData::createFaceFluxPolyNodeWheightCoord()
{
  CFAUTOTRACE;

  // get the number of face flux points
  const CFuint nbrFaceFlxPnts = getNbrOfFaceFlxPnts();

  // resize m_faceFluxPolyNodeWheightCoord
  m_faceFluxPolyNodeWheightCoord.resize(nbrFaceFlxPnts);

  // use the first SV face (for which the last local coordinate should be zero)
  for (CFuint iFlx = 0; iFlx < nbrFaceFlxPnts; ++iFlx)
  {
    // get node ID
    const CFuint nodeID = m_faceFlxPntsConn[0][iFlx];

    // resize
    m_faceFluxPolyNodeWheightCoord[iFlx].resize(m_dimensionality);

    // compute wheight coordinates
    m_faceFluxPolyNodeWheightCoord[iFlx][0] = 1.0 - m_fluxPolyNodeCoord[nodeID].sum();
    // !!! the KSI and ETA axis are inverted for the SV face with respect to the SV axis
    m_faceFluxPolyNodeWheightCoord[iFlx][1] = m_fluxPolyNodeCoord[nodeID][1];
    m_faceFluxPolyNodeWheightCoord[iFlx][2] = m_fluxPolyNodeCoord[nodeID][0];
  }
}

//////////////////////////////////////////////////////////////////////

void TetraSpectralFVElementData::createFluxPolyExponents()
{
  CFAUTOTRACE;

  // number of polynomial terms
  const CFuint nbrPolyTerms = (m_polyOrder + 2)*(m_polyOrder + 3)*(m_polyOrder + 4)/6;

  // resize the variable
  m_fluxPolyExponents.resize(nbrPolyTerms);
  for (CFuint iTerm = 0; iTerm < nbrPolyTerms; ++iTerm)
  {
    m_fluxPolyExponents[iTerm].resize(3);
  }

  // define exponents
  const CFuint polyOrderP2 = m_polyOrder+2;
  CFuint iTerm = 0;
  for (CFuint iXYZ = 0; iXYZ < polyOrderP2; ++iXYZ)
  {
    for (CFuint iX = 0; iX < iXYZ+1; ++iX)
    {
      for (CFuint iY = 0; iY < iXYZ-iX+1; ++iY,++iTerm)
      {
        m_fluxPolyExponents[iTerm][KSI] = iX;
        m_fluxPolyExponents[iTerm][ETA] = iY;
        m_fluxPolyExponents[iTerm][ZTA] = iXYZ-iX-iY;
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////

void TetraSpectralFVElementData::createVolumeTermsTensor()
{
  CFAUTOTRACE;

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
  const CFuint dimUns = static_cast<CFuint>(m_dimensionality);
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

  // set the dimensionality and order of the tensor product integrator
  const CFDim dimM1 = static_cast<CFDim>(m_dimensionality-1);
  m_tpIntegrator.setDimensionality(dimM1);
  const CFPolyOrder::Type polyOrderP1 = static_cast<CFPolyOrder::Type>(m_polyOrder+1);
  m_tpIntegrator.setIntegratorOrder(polyOrderP1);

  // loop over internal faces
  for (CFuint iFace = 0; iFace < nbrLocalIntFaces; ++iFace)
  {
    // local face ID
    const CFuint faceID = iFace + nbrLocalExtFaces;

    // number of face nodes
    const CFuint nbrFaceNodes = m_localFaceNodeConn[faceID].size();

    // get node coordinates
    vector< RealVector > faceNodeCoord(nbrFaceNodes,RealVector(m_dimensionality));
    for (CFuint iNode = 0; iNode < nbrFaceNodes; ++iNode)
    {
      const CFuint nodeID = m_localFaceNodeConn[faceID][iNode];
      faceNodeCoord[iNode] = m_localNodeCoord[nodeID];
    }

    // get decomposition into quadrilaterals
    vector< vector< RealVector > > quadDecomp;
    setPolyGonQuadrilateralDecomposition(faceNodeCoord,quadDecomp);

    // number of quadrilaterals in CV face
    const CFuint nbrQuads = quadDecomp.size();

    // get quadrature nodes and wheights
    vector< RealVector > qNodeCoord;
    vector< RealVector > qNodeNormals;
    vector< CFreal >     qWheights;
    for (CFuint iQuad = 0; iQuad < nbrQuads; ++iQuad)
    {
      // get quadrilateral quadrature nodes and wheights
      vector< RealVector > qNodeCoordQuad   = m_tpIntegrator.getQuadPntsCoordsPlus1D     (quadDecomp[iQuad]);
      vector< RealVector > qNodeNormalsQuad = m_tpIntegrator.getQuadPntsUnitNormalsPlus1D(quadDecomp[iQuad]);
      vector< CFreal >     qWheightsQuad    = m_tpIntegrator.getQuadPntsWheightsPlus1D   (quadDecomp[iQuad]);

      // add quadrilateral quadrature nodes and wheights to global list
      qNodeCoord  .insert(qNodeCoord.end()  ,qNodeCoordQuad  .begin(),qNodeCoordQuad  .end());
      qNodeNormals.insert(qNodeNormals.end(),qNodeNormalsQuad.begin(),qNodeNormalsQuad.end());
      qWheights   .insert(qWheights .end()  ,qWheightsQuad   .begin(),qWheightsQuad   .end());
    }

    // compute integral of flux polynomial over face
    const CFuint nbrQPnts = qWheights.size();
    for (CFuint iFlx = 0; iFlx < nbrFlxPnts; ++iFlx)
    {
      RealVector integral(m_dimensionality);
      for (CFuint iQPnt = 0; iQPnt < nbrQPnts; ++iQPnt)
      {
        CFreal term = 0.0;
        for (CFuint iTerm = 0; iTerm < nbrFlxPnts; ++iTerm)
        {
          CFreal term2 = m_fluxPolyCoefs[iFlx][iTerm];
          for (CFuint iCoor = 0; iCoor < dimUns; ++iCoor)
          {
            term2 *= pow(qNodeCoord[iQPnt][iCoor],m_fluxPolyExponents[iTerm][iCoor]);
          }
          term += term2;
        }
        integral += qWheights[iQPnt]*term*qNodeNormals[iQPnt];
      }

      // add contribution to volume term tensor
      for (CFuint iCoor = 0; iCoor < dimUns; ++iCoor)
      {
        m_volTermTensor[m_localIntFaceCVConn[iFace][LEFT ]][iFlx][iCoor] -= integral[iCoor];
        m_volTermTensor[m_localIntFaceCVConn[iFace][RIGHT]][iFlx][iCoor] += integral[iCoor];
      }
    }
  }

/*  for (CFuint iCV = 0; iCV < nbrCVs; ++iCV)
  {
  CF_DEBUG_OBJ(iCV);
  for (CFuint iFlx = 0; iFlx < nbrFlxPnts; ++iFlx)
  {
  for (CFuint iCoor = 0; iCoor < dimUns; ++iCoor)
  {
  CF_DEBUG_OBJ(m_volTermTensor[iCV][iFlx][iCoor]);
}
}
}*/
}

//////////////////////////////////////////////////////////////////////

void TetraSpectralFVElementData::setInterpolationNodeSet(const CFPolyOrder::Type order,
                                                         vector< RealVector >& nodalSet)
{
  // number of nodes
  const CFuint nbrNodes = (order+1)*(order+2)*(order+3)/6;

  // return variable
  nodalSet.resize(nbrNodes);
  for (CFuint iNode = 0; iNode < nbrNodes; ++iNode)
  {
    nodalSet[iNode].resize(3);
  }

  // fill the vector containing the flux polynomial node coordinates
  // nodes from Hesthaven and Teng, STABLE SPECTRAL METHODS ON TETRAHEDRAL ELEMENTS.
  // SIAM J. Sci. Comput. (2000); 21(6):2352-2380.
  switch (order)
  {
    case CFPolyOrder::ORDER0:
    {
      nodalSet[0][KSI] = 1.0/4.0;
      nodalSet[0][ETA] = 1.0/4.0;
      nodalSet[0][ZTA] = 1.0/4.0;
    } break;
    case CFPolyOrder::ORDER1:
    {
      nodalSet[0][KSI] = 0.0;
      nodalSet[0][ETA] = 0.0;
      nodalSet[0][ZTA] = 0.0;

      nodalSet[1][KSI] = 1.0;
      nodalSet[1][ETA] = 0.0;
      nodalSet[1][ZTA] = 0.0;

      nodalSet[2][KSI] = 0.0;
      nodalSet[2][ETA] = 1.0;
      nodalSet[2][ZTA] = 0.0;

      nodalSet[3][KSI] = 0.0;
      nodalSet[3][ETA] = 0.0;
      nodalSet[3][ZTA] = 1.0;
    } break;
    case CFPolyOrder::ORDER2:
    {
      nodalSet[0][KSI] = 0.0;
      nodalSet[0][ETA] = 0.0;
      nodalSet[0][ZTA] = 0.0;

      nodalSet[1][KSI] = 0.5;
      nodalSet[1][ETA] = 0.0;
      nodalSet[1][ZTA] = 0.0;

      nodalSet[2][KSI] = 0.0;
      nodalSet[2][ETA] = 0.5;
      nodalSet[2][ZTA] = 0.0;

      nodalSet[3][KSI] = 0.0;
      nodalSet[3][ETA] = 0.0;
      nodalSet[3][ZTA] = 0.5;

      nodalSet[4][KSI] = 1.0;
      nodalSet[4][ETA] = 0.0;
      nodalSet[4][ZTA] = 0.0;

      nodalSet[5][KSI] = 0.5;
      nodalSet[5][ETA] = 0.5;
      nodalSet[5][ZTA] = 0.0;

      nodalSet[6][KSI] = 0.0;
      nodalSet[6][ETA] = 1.0;
      nodalSet[6][ZTA] = 0.0;

      nodalSet[7][KSI] = 0.5;
      nodalSet[7][ETA] = 0.0;
      nodalSet[7][ZTA] = 0.5;

      nodalSet[8][KSI] = 0.0;
      nodalSet[8][ETA] = 0.5;
      nodalSet[8][ZTA] = 0.5;

      nodalSet[9][KSI] = 0.0;
      nodalSet[9][ETA] = 0.0;
      nodalSet[9][ZTA] = 1.0;
    } break;
    case CFPolyOrder::ORDER3:
    {
      const CFreal oEminusA = 0.7236067977;
      const CFreal a = 1.0 - oEminusA;

      nodalSet[0][KSI] = 0.0;
      nodalSet[0][ETA] = 0.0;
      nodalSet[0][ZTA] = 0.0;

      nodalSet[1][KSI] = a;
      nodalSet[1][ETA] = 0.0;
      nodalSet[1][ZTA] = 0.0;

      nodalSet[2][KSI] = 0.0;
      nodalSet[2][ETA] = a;
      nodalSet[2][ZTA] = 0.0;

      nodalSet[3][KSI] = 0.0;
      nodalSet[3][ETA] = 0.0;
      nodalSet[3][ZTA] = a;

      nodalSet[4][KSI] = oEminusA;
      nodalSet[4][ETA] = 0.0;
      nodalSet[4][ZTA] = 0.0;

      nodalSet[5][KSI] = 1.0/3.0;
      nodalSet[5][ETA] = 1.0/3.0;
      nodalSet[5][ZTA] = 0.0;

      nodalSet[6][KSI] = 0.0;
      nodalSet[6][ETA] = oEminusA;
      nodalSet[6][ZTA] = 0.0;

      nodalSet[7][KSI] = 1.0/3.0;
      nodalSet[7][ETA] = 0.0;
      nodalSet[7][ZTA] = 1.0/3.0;

      nodalSet[8][KSI] = 0.0;
      nodalSet[8][ETA] = 1.0/3.0;
      nodalSet[8][ZTA] = 1.0/3.0;

      nodalSet[9][KSI] = 0.0;
      nodalSet[9][ETA] = 0.0;
      nodalSet[9][ZTA] = oEminusA;

      nodalSet[10][KSI] = 1.0;
      nodalSet[10][ETA] = 0.0;
      nodalSet[10][ZTA] = 0.0;

      nodalSet[11][KSI] = oEminusA;
      nodalSet[11][ETA] = a;
      nodalSet[11][ZTA] = 0.0;

      nodalSet[12][KSI] = a;
      nodalSet[12][ETA] = oEminusA;
      nodalSet[12][ZTA] = 0.0;

      nodalSet[13][KSI] = 0.0;
      nodalSet[13][ETA] = 1.0;
      nodalSet[13][ZTA] = 0.0;

      nodalSet[14][KSI] = oEminusA;
      nodalSet[14][ETA] = 0.0;
      nodalSet[14][ZTA] = a;

      nodalSet[15][KSI] = 1.0/3.0;
      nodalSet[15][ETA] = 1.0/3.0;
      nodalSet[15][ZTA] = 1.0/3.0;

      nodalSet[16][KSI] = 0.0;
      nodalSet[16][ETA] = oEminusA;
      nodalSet[16][ZTA] = a;

      nodalSet[17][KSI] = a;
      nodalSet[17][ETA] = 0.0;
      nodalSet[17][ZTA] = oEminusA;

      nodalSet[18][KSI] = 0.0;
      nodalSet[18][ETA] = a;
      nodalSet[18][ZTA] = oEminusA;

      nodalSet[19][KSI] = 0.0;
      nodalSet[19][ETA] = 0.0;
      nodalSet[19][ZTA] = 1.0;
    } break;
    case CFPolyOrder::ORDER4:
    {
      const CFreal oEminusA = 0.8273268354;
      const CFreal a = 1.0 - oEminusA;

      const CFreal oEminusB = 0.5257599664;
      const CFreal bDivTwo = 0.5*(1.0 - oEminusB);

      nodalSet[0][KSI] = 0.0;
      nodalSet[0][ETA] = 0.0;
      nodalSet[0][ZTA] = 0.0;

      nodalSet[1][KSI] = a;
      nodalSet[1][ETA] = 0.0;
      nodalSet[1][ZTA] = 0.0;

      nodalSet[2][KSI] = 0.0;
      nodalSet[2][ETA] = a;
      nodalSet[2][ZTA] = 0.0;

      nodalSet[3][KSI] = 0.0;
      nodalSet[3][ETA] = 0.0;
      nodalSet[3][ZTA] = a;

      nodalSet[4][KSI] = 0.5;
      nodalSet[4][ETA] = 0.0;
      nodalSet[4][ZTA] = 0.0;

      nodalSet[5][KSI] = bDivTwo;
      nodalSet[5][ETA] = bDivTwo;
      nodalSet[5][ZTA] = 0.0;

      nodalSet[6][KSI] = 0.0;
      nodalSet[6][ETA] = 0.5;
      nodalSet[6][ZTA] = 0.0;

      nodalSet[7][KSI] = bDivTwo;
      nodalSet[7][ETA] = 0.0;
      nodalSet[7][ZTA] = bDivTwo;

      nodalSet[8][KSI] = 0.0;
      nodalSet[8][ETA] = bDivTwo;
      nodalSet[8][ZTA] = bDivTwo;

      nodalSet[9][KSI] = 0.0;
      nodalSet[9][ETA] = 0.0;
      nodalSet[9][ZTA] = 0.5;

      nodalSet[10][KSI] = oEminusA;
      nodalSet[10][ETA] = 0.0;
      nodalSet[10][ZTA] = 0.0;

      nodalSet[11][KSI] = oEminusB;
      nodalSet[11][ETA] = bDivTwo;
      nodalSet[11][ZTA] = 0.0;

      nodalSet[12][KSI] = bDivTwo;
      nodalSet[12][ETA] = oEminusB;
      nodalSet[12][ZTA] = 0.0;

      nodalSet[13][KSI] = 0.0;
      nodalSet[13][ETA] = oEminusA;
      nodalSet[13][ZTA] = 0.0;

      nodalSet[14][KSI] = oEminusB;
      nodalSet[14][ETA] = 0.0;
      nodalSet[14][ZTA] = bDivTwo;

      nodalSet[15][KSI] = 1.0/4.0;
      nodalSet[15][ETA] = 1.0/4.0;
      nodalSet[15][ZTA] = 1.0/4.0;

      nodalSet[16][KSI] = 0.0;
      nodalSet[16][ETA] = oEminusB;
      nodalSet[16][ZTA] = bDivTwo;

      nodalSet[17][KSI] = bDivTwo;
      nodalSet[17][ETA] = 0.0;
      nodalSet[17][ZTA] = oEminusB;

      nodalSet[18][KSI] = 0.0;
      nodalSet[18][ETA] = bDivTwo;
      nodalSet[18][ZTA] = oEminusB;

      nodalSet[19][KSI] = 0.0;
      nodalSet[19][ETA] = 0.0;
      nodalSet[19][ZTA] = oEminusA;

      nodalSet[20][KSI] = 1.0;
      nodalSet[20][ETA] = 0.0;
      nodalSet[20][ZTA] = 0.0;

      nodalSet[21][KSI] = oEminusA;
      nodalSet[21][ETA] = a;
      nodalSet[21][ZTA] = 0.0;

      nodalSet[22][KSI] = 0.5;
      nodalSet[22][ETA] = 0.5;
      nodalSet[22][ZTA] = 0.0;

      nodalSet[23][KSI] = a;
      nodalSet[23][ETA] = oEminusA;
      nodalSet[23][ZTA] = 0.0;

      nodalSet[24][KSI] = 0.0;
      nodalSet[24][ETA] = 1.0;
      nodalSet[24][ZTA] = 0.0;

      nodalSet[25][KSI] = oEminusA;
      nodalSet[25][ETA] = 0.0;
      nodalSet[25][ZTA] = a;

      nodalSet[26][KSI] = oEminusB;
      nodalSet[26][ETA] = bDivTwo;
      nodalSet[26][ZTA] = bDivTwo;

      nodalSet[27][KSI] = bDivTwo;
      nodalSet[27][ETA] = oEminusB;
      nodalSet[27][ZTA] = bDivTwo;

      nodalSet[28][KSI] = 0.0;
      nodalSet[28][ETA] = oEminusA;
      nodalSet[28][ZTA] = a;

      nodalSet[29][KSI] = 0.5;
      nodalSet[29][ETA] = 0.0;
      nodalSet[29][ZTA] = 0.5;

      nodalSet[30][KSI] = bDivTwo;
      nodalSet[30][ETA] = bDivTwo;
      nodalSet[30][ZTA] = oEminusB;

      nodalSet[31][KSI] = 0.0;
      nodalSet[31][ETA] = 0.5;
      nodalSet[31][ZTA] = 0.5;

      nodalSet[32][KSI] = a;
      nodalSet[32][ETA] = 0.0;
      nodalSet[32][ZTA] = oEminusA;

      nodalSet[33][KSI] = 0.0;
      nodalSet[33][ETA] = a;
      nodalSet[33][ZTA] = oEminusA;

      nodalSet[34][KSI] = 0.0;
      nodalSet[34][ETA] = 0.0;
      nodalSet[34][ZTA] = 1.0;
    } break;
    default:
    {
      throw Common::NotImplementedException (FromHere(),"Higher-order nodal set not defined!");
    }
  }
}

//////////////////////////////////////////////////////////////////////

void TetraSpectralFVElementData::
      setPolyHedronTetraDecomposition(const CFuint nbrNodes,vector< vector< CFuint > >& tetraDecomp)
{
  // set tetra decomposition
  switch (nbrNodes)
  {
    case 4:
    {
      tetraDecomp.resize(1);

      tetraDecomp[0].resize(4);
      tetraDecomp[0][0] = 0;
      tetraDecomp[0][1] = 1;
      tetraDecomp[0][2] = 2;
      tetraDecomp[0][3] = 3;
    } break;
    case 8:
    {
      tetraDecomp.resize(5);

      // for SV3_3D schemes of first family
      tetraDecomp[0].resize(4);
      tetraDecomp[0][0] = 0;
      tetraDecomp[0][1] = 1;
      tetraDecomp[0][2] = 2;
      tetraDecomp[0][3] = 5;

      tetraDecomp[1].resize(4);
      tetraDecomp[1][0] = 5;
      tetraDecomp[1][1] = 2;
      tetraDecomp[1][2] = 7;
      tetraDecomp[1][3] = 6;

      tetraDecomp[2].resize(4);
      tetraDecomp[2][0] = 0;
      tetraDecomp[2][1] = 5;
      tetraDecomp[2][2] = 7;
      tetraDecomp[2][3] = 4;

      tetraDecomp[3].resize(4);
      tetraDecomp[3][0] = 0;
      tetraDecomp[3][1] = 2;
      tetraDecomp[3][2] = 3;
      tetraDecomp[3][3] = 7;

      tetraDecomp[4].resize(4);
      tetraDecomp[4][0] = 0;
      tetraDecomp[4][1] = 5;
      tetraDecomp[4][2] = 2;
      tetraDecomp[4][3] = 7;


/*      // for SV3_3D schemes of second family
      tetraDecomp[0].resize(4);
      tetraDecomp[0][0] = 0;
      tetraDecomp[0][1] = 1;
      tetraDecomp[0][2] = 3;
      tetraDecomp[0][3] = 4;

      tetraDecomp[1].resize(4);
      tetraDecomp[1][0] = 1;
      tetraDecomp[1][1] = 2;
      tetraDecomp[1][2] = 3;
      tetraDecomp[1][3] = 6;

      tetraDecomp[2].resize(4);
      tetraDecomp[2][0] = 1;
      tetraDecomp[2][1] = 4;
      tetraDecomp[2][2] = 5;
      tetraDecomp[2][3] = 6;

      tetraDecomp[3].resize(4);
      tetraDecomp[3][0] = 3;
      tetraDecomp[3][1] = 4;
      tetraDecomp[3][2] = 6;
      tetraDecomp[3][3] = 7;

      tetraDecomp[4].resize(4);
      tetraDecomp[4][0] = 1;
      tetraDecomp[4][1] = 3;
      tetraDecomp[4][2] = 4;
      tetraDecomp[4][3] = 6;*/
    } break;
    case 11:
    {
      tetraDecomp.resize(9);

      // for SV3_3D schemes of first family
      tetraDecomp[0].resize(4);
      tetraDecomp[0][0] = 0;
      tetraDecomp[0][1] = 1;
      tetraDecomp[0][2] = 2;
      tetraDecomp[0][3] = 5;

      tetraDecomp[1].resize(4);
      tetraDecomp[1][0] = 5;
      tetraDecomp[1][1] = 2;
      tetraDecomp[1][2] = 7;
      tetraDecomp[1][3] = 6;

      tetraDecomp[2].resize(4);
      tetraDecomp[2][0] = 0;
      tetraDecomp[2][1] = 5;
      tetraDecomp[2][2] = 7;
      tetraDecomp[2][3] = 4;

      tetraDecomp[3].resize(4);
      tetraDecomp[3][0] = 0;
      tetraDecomp[3][1] = 2;
      tetraDecomp[3][2] = 3;
      tetraDecomp[3][3] = 7;

      tetraDecomp[4].resize(4);
      tetraDecomp[4][0] = 0;
      tetraDecomp[4][1] = 5;
      tetraDecomp[4][2] = 2;
      tetraDecomp[4][3] = 7;

      tetraDecomp[5].resize(4);
      tetraDecomp[5][0] = 10;
      tetraDecomp[5][1] = 9;
      tetraDecomp[5][2] = 8;
      tetraDecomp[5][3] = 5;

      tetraDecomp[6].resize(4);
      tetraDecomp[6][0] = 7;
      tetraDecomp[6][1] = 10;
      tetraDecomp[6][2] = 8;
      tetraDecomp[6][3] = 5;

      tetraDecomp[7].resize(4);
      tetraDecomp[7][0] = 6;
      tetraDecomp[7][1] = 10;
      tetraDecomp[7][2] = 7;
      tetraDecomp[7][3] = 5;

      tetraDecomp[8].resize(4);
      tetraDecomp[8][0] = 4;
      tetraDecomp[8][1] = 7;
      tetraDecomp[8][2] = 8;
      tetraDecomp[8][3] = 5;


/*      // for SV3_3D schemes of second family
      tetraDecomp[0].resize(4);
      tetraDecomp[0][0] = 0;
      tetraDecomp[0][1] = 1;
      tetraDecomp[0][2] = 3;
      tetraDecomp[0][3] = 4;

      tetraDecomp[1].resize(4);
      tetraDecomp[1][0] = 1;
      tetraDecomp[1][1] = 2;
      tetraDecomp[1][2] = 3;
      tetraDecomp[1][3] = 6;

      tetraDecomp[2].resize(4);
      tetraDecomp[2][0] = 1;
      tetraDecomp[2][1] = 4;
      tetraDecomp[2][2] = 5;
      tetraDecomp[2][3] = 6;

      tetraDecomp[3].resize(4);
      tetraDecomp[3][0] = 3;
      tetraDecomp[3][1] = 4;
      tetraDecomp[3][2] = 6;
      tetraDecomp[3][3] = 7;

      tetraDecomp[4].resize(4);
      tetraDecomp[4][0] = 1;
      tetraDecomp[4][1] = 3;
      tetraDecomp[4][2] = 4;
      tetraDecomp[4][3] = 6;

      tetraDecomp[5].resize(4);
      tetraDecomp[5][0] = 4;
      tetraDecomp[5][1] = 5;
      tetraDecomp[5][2] = 6;
      tetraDecomp[5][3] = 9;

      tetraDecomp[6].resize(4);
      tetraDecomp[6][0] = 4;
      tetraDecomp[6][1] = 6;
      tetraDecomp[6][2] = 7;
      tetraDecomp[6][3] = 9;

      tetraDecomp[7].resize(4);
      tetraDecomp[7][0] = 4;
      tetraDecomp[7][1] = 7;
      tetraDecomp[7][2] = 8;
      tetraDecomp[7][3] = 9;

      tetraDecomp[8].resize(4);
      tetraDecomp[8][0] = 6;
      tetraDecomp[8][1] = 7;
      tetraDecomp[8][2] = 9;
      tetraDecomp[8][3] = 10;*/
    } break;
    default:
    {
      throw Common::NotImplementedException (FromHere(),"Unsupported polyhedron form...");
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void TetraSpectralFVElementData::setPolyHedronEdgeNodeConn(const CFuint nbrNodes,vector< vector< CFuint > >& edgeToNodes)
{
  // set tetra decomposition
  switch (nbrNodes)
  {
    case 4:
    {
      edgeToNodes.resize(6);

      CFuint iEdge = 0;

      edgeToNodes[iEdge].resize(2);
      edgeToNodes[iEdge][0] = 0;
      edgeToNodes[iEdge][1] = 1;
      ++iEdge;

      edgeToNodes[iEdge].resize(2);
      edgeToNodes[iEdge][0] = 1;
      edgeToNodes[iEdge][1] = 2;
      ++iEdge;

      edgeToNodes[iEdge].resize(2);
      edgeToNodes[iEdge][0] = 2;
      edgeToNodes[iEdge][1] = 0;
      ++iEdge;

      edgeToNodes[iEdge].resize(2);
      edgeToNodes[iEdge][0] = 3;
      edgeToNodes[iEdge][1] = 0;
      ++iEdge;

      edgeToNodes[iEdge].resize(2);
      edgeToNodes[iEdge][0] = 3;
      edgeToNodes[iEdge][1] = 1;
      ++iEdge;

      edgeToNodes[iEdge].resize(2);
      edgeToNodes[iEdge][0] = 3;
      edgeToNodes[iEdge][1] = 2;
      ++iEdge;

    } break;
    case 8:
    {
      edgeToNodes.resize(15);

      CFuint iEdge = 0;

      // for SV3_3D schemes of first family
      edgeToNodes[iEdge].resize(2);
      edgeToNodes[iEdge][0] = 0;
      edgeToNodes[iEdge][1] = 1;
      ++iEdge;

      edgeToNodes[iEdge].resize(2);
      edgeToNodes[iEdge][0] = 1;
      edgeToNodes[iEdge][1] = 2;
      ++iEdge;

      edgeToNodes[iEdge].resize(2);
      edgeToNodes[iEdge][0] = 2;
      edgeToNodes[iEdge][1] = 3;
      ++iEdge;

      edgeToNodes[iEdge].resize(2);
      edgeToNodes[iEdge][0] = 3;
      edgeToNodes[iEdge][1] = 0;
      ++iEdge;

      edgeToNodes[iEdge].resize(2);
      edgeToNodes[iEdge][0] = 4;
      edgeToNodes[iEdge][1] = 5;
      ++iEdge;

      edgeToNodes[iEdge].resize(2);
      edgeToNodes[iEdge][0] = 5;
      edgeToNodes[iEdge][1] = 6;
      ++iEdge;

      edgeToNodes[iEdge].resize(2);
      edgeToNodes[iEdge][0] = 6;
      edgeToNodes[iEdge][1] = 7;
      ++iEdge;

      edgeToNodes[iEdge].resize(2);
      edgeToNodes[iEdge][0] = 7;
      edgeToNodes[iEdge][1] = 4;
      ++iEdge;

      edgeToNodes[iEdge].resize(2);
      edgeToNodes[iEdge][0] = 0;
      edgeToNodes[iEdge][1] = 4;
      ++iEdge;

      edgeToNodes[iEdge].resize(2);
      edgeToNodes[iEdge][0] = 1;
      edgeToNodes[iEdge][1] = 5;
      ++iEdge;

      edgeToNodes[iEdge].resize(2);
      edgeToNodes[iEdge][0] = 2;
      edgeToNodes[iEdge][1] = 6;
      ++iEdge;

      edgeToNodes[iEdge].resize(2);
      edgeToNodes[iEdge][0] = 3;
      edgeToNodes[iEdge][1] = 7;
      ++iEdge;

      edgeToNodes[iEdge].resize(2);
      edgeToNodes[iEdge][0] = 2;
      edgeToNodes[iEdge][1] = 5;
      ++iEdge;

      edgeToNodes[iEdge].resize(2);
      edgeToNodes[iEdge][0] = 2;
      edgeToNodes[iEdge][1] = 7;
      ++iEdge;

      edgeToNodes[iEdge].resize(2);
      edgeToNodes[iEdge][0] = 5;
      edgeToNodes[iEdge][1] = 7;
      ++iEdge;

/*      // for SV3_3D schemes of second family
      edgeToNodes[iEdge].resize(2);
      edgeToNodes[iEdge][0] = 0;
      edgeToNodes[iEdge][1] = 1;
      ++iEdge;

      edgeToNodes[iEdge].resize(2);
      edgeToNodes[iEdge][0] = 1;
      edgeToNodes[iEdge][1] = 2;
      ++iEdge;

      edgeToNodes[iEdge].resize(2);
      edgeToNodes[iEdge][0] = 2;
      edgeToNodes[iEdge][1] = 3;
      ++iEdge;

      edgeToNodes[iEdge].resize(2);
      edgeToNodes[iEdge][0] = 3;
      edgeToNodes[iEdge][1] = 0;
      ++iEdge;

      edgeToNodes[iEdge].resize(2);
      edgeToNodes[iEdge][0] = 4;
      edgeToNodes[iEdge][1] = 5;
      ++iEdge;

      edgeToNodes[iEdge].resize(2);
      edgeToNodes[iEdge][0] = 5;
      edgeToNodes[iEdge][1] = 6;
      ++iEdge;

      edgeToNodes[iEdge].resize(2);
      edgeToNodes[iEdge][0] = 6;
      edgeToNodes[iEdge][1] = 7;
      ++iEdge;

      edgeToNodes[iEdge].resize(2);
      edgeToNodes[iEdge][0] = 7;
      edgeToNodes[iEdge][1] = 4;
      ++iEdge;

      edgeToNodes[iEdge].resize(2);
      edgeToNodes[iEdge][0] = 0;
      edgeToNodes[iEdge][1] = 4;
      ++iEdge;

      edgeToNodes[iEdge].resize(2);
      edgeToNodes[iEdge][0] = 1;
      edgeToNodes[iEdge][1] = 5;
      ++iEdge;

      edgeToNodes[iEdge].resize(2);
      edgeToNodes[iEdge][0] = 2;
      edgeToNodes[iEdge][1] = 6;
      ++iEdge;

      edgeToNodes[iEdge].resize(2);
      edgeToNodes[iEdge][0] = 3;
      edgeToNodes[iEdge][1] = 7;
      ++iEdge;

      edgeToNodes[iEdge].resize(2);
      edgeToNodes[iEdge][0] = 1;
      edgeToNodes[iEdge][1] = 6;
      ++iEdge;

      edgeToNodes[iEdge].resize(2);
      edgeToNodes[iEdge][0] = 4;
      edgeToNodes[iEdge][1] = 6;
      ++iEdge;

      edgeToNodes[iEdge].resize(2);
      edgeToNodes[iEdge][0] = 3;
      edgeToNodes[iEdge][1] = 6;
      ++iEdge;*/
    } break;
    case 11:
    {
      edgeToNodes.resize(19);

      CFuint iEdge = 0;

      // for SV3_3D schemes of first family
      edgeToNodes[iEdge].resize(2);
      edgeToNodes[iEdge][0] = 0;
      edgeToNodes[iEdge][1] = 1;
      ++iEdge;

      edgeToNodes[iEdge].resize(2);
      edgeToNodes[iEdge][0] = 1;
      edgeToNodes[iEdge][1] = 2;
      ++iEdge;

      edgeToNodes[iEdge].resize(2);
      edgeToNodes[iEdge][0] = 2;
      edgeToNodes[iEdge][1] = 3;
      ++iEdge;

      edgeToNodes[iEdge].resize(2);
      edgeToNodes[iEdge][0] = 3;
      edgeToNodes[iEdge][1] = 0;
      ++iEdge;

      edgeToNodes[iEdge].resize(2);
      edgeToNodes[iEdge][0] = 7;
      edgeToNodes[iEdge][1] = 8;
      ++iEdge;

      edgeToNodes[iEdge].resize(2);
      edgeToNodes[iEdge][0] = 8;
      edgeToNodes[iEdge][1] = 9;
      ++iEdge;

      edgeToNodes[iEdge].resize(2);
      edgeToNodes[iEdge][0] = 9;
      edgeToNodes[iEdge][1] = 10;
      ++iEdge;

      edgeToNodes[iEdge].resize(2);
      edgeToNodes[iEdge][0] = 10;
      edgeToNodes[iEdge][1] = 7;
      ++iEdge;

      edgeToNodes[iEdge].resize(2);
      edgeToNodes[iEdge][0] = 3;
      edgeToNodes[iEdge][1] = 7;
      ++iEdge;

      edgeToNodes[iEdge].resize(2);
      edgeToNodes[iEdge][0] = 0;
      edgeToNodes[iEdge][1] = 4;
      ++iEdge;

      edgeToNodes[iEdge].resize(2);
      edgeToNodes[iEdge][0] = 4;
      edgeToNodes[iEdge][1] = 8;
      ++iEdge;

      edgeToNodes[iEdge].resize(2);
      edgeToNodes[iEdge][0] = 1;
      edgeToNodes[iEdge][1] = 5;
      ++iEdge;

      edgeToNodes[iEdge].resize(2);
      edgeToNodes[iEdge][0] = 5;
      edgeToNodes[iEdge][1] = 9;
      ++iEdge;

      edgeToNodes[iEdge].resize(2);
      edgeToNodes[iEdge][0] = 2;
      edgeToNodes[iEdge][1] = 6;
      ++iEdge;

      edgeToNodes[iEdge].resize(2);
      edgeToNodes[iEdge][0] = 6;
      edgeToNodes[iEdge][1] = 10;
      ++iEdge;

      edgeToNodes[iEdge].resize(2);
      edgeToNodes[iEdge][0] = 4;
      edgeToNodes[iEdge][1] = 5;
      ++iEdge;

      edgeToNodes[iEdge].resize(2);
      edgeToNodes[iEdge][0] = 5;
      edgeToNodes[iEdge][1] = 6;
      ++iEdge;

      edgeToNodes[iEdge].resize(2);
      edgeToNodes[iEdge][0] = 2;
      edgeToNodes[iEdge][1] = 0;
      ++iEdge;

      edgeToNodes[iEdge].resize(2);
      edgeToNodes[iEdge][0] = 10;
      edgeToNodes[iEdge][1] = 8;
      ++iEdge;


//       // for SV3_3D schemes of second family
//       edgeToNodes[iEdge].resize(2);
//       edgeToNodes[iEdge][0] = 0;
//       edgeToNodes[iEdge][1] = 1;
//       ++iEdge;
//
//       edgeToNodes[iEdge].resize(2);
//       edgeToNodes[iEdge][0] = 1;
//       edgeToNodes[iEdge][1] = 2;
//       ++iEdge;
//
//       edgeToNodes[iEdge].resize(2);
//       edgeToNodes[iEdge][0] = 2;
//       edgeToNodes[iEdge][1] = 3;
//       ++iEdge;
//
//       edgeToNodes[iEdge].resize(2);
//       edgeToNodes[iEdge][0] = 3;
//       edgeToNodes[iEdge][1] = 0;
//       ++iEdge;
//
//       edgeToNodes[iEdge].resize(2);
//       edgeToNodes[iEdge][0] = 7;
//       edgeToNodes[iEdge][1] = 8;
//       ++iEdge;
//
//       edgeToNodes[iEdge].resize(2);
//       edgeToNodes[iEdge][0] = 8;
//       edgeToNodes[iEdge][1] = 9;
//       ++iEdge;
//
//       edgeToNodes[iEdge].resize(2);
//       edgeToNodes[iEdge][0] = 9;
//       edgeToNodes[iEdge][1] = 10;
//       ++iEdge;
//
//       edgeToNodes[iEdge].resize(2);
//       edgeToNodes[iEdge][0] = 10;
//       edgeToNodes[iEdge][1] = 7;
//       ++iEdge;
//
//       edgeToNodes[iEdge].resize(2);
//       edgeToNodes[iEdge][0] = 3;
//       edgeToNodes[iEdge][1] = 7;
//       ++iEdge;
//
//       edgeToNodes[iEdge].resize(2);
//       edgeToNodes[iEdge][0] = 0;
//       edgeToNodes[iEdge][1] = 4;
//       ++iEdge;
//
//       edgeToNodes[iEdge].resize(2);
//       edgeToNodes[iEdge][0] = 4;
//       edgeToNodes[iEdge][1] = 8;
//       ++iEdge;
//
//       edgeToNodes[iEdge].resize(2);
//       edgeToNodes[iEdge][0] = 1;
//       edgeToNodes[iEdge][1] = 5;
//       ++iEdge;
//
//       edgeToNodes[iEdge].resize(2);
//       edgeToNodes[iEdge][0] = 5;
//       edgeToNodes[iEdge][1] = 9;
//       ++iEdge;
//
//       edgeToNodes[iEdge].resize(2);
//       edgeToNodes[iEdge][0] = 2;
//       edgeToNodes[iEdge][1] = 6;
//       ++iEdge;
//
//       edgeToNodes[iEdge].resize(2);
//       edgeToNodes[iEdge][0] = 6;
//       edgeToNodes[iEdge][1] = 10;
//       ++iEdge;
//
//       edgeToNodes[iEdge].resize(2);
//       edgeToNodes[iEdge][0] = 4;
//       edgeToNodes[iEdge][1] = 5;
//       ++iEdge;
//
//       edgeToNodes[iEdge].resize(2);
//       edgeToNodes[iEdge][0] = 5;
//       edgeToNodes[iEdge][1] = 6;
//       ++iEdge;
//
//       edgeToNodes[iEdge].resize(2);
//       edgeToNodes[iEdge][0] = 1;
//       edgeToNodes[iEdge][1] = 3;
//       ++iEdge;
//
//       edgeToNodes[iEdge].resize(2);
//       edgeToNodes[iEdge][0] = 9;
//       edgeToNodes[iEdge][1] = 7;
//       ++iEdge;

    } break;
    default:
    {
      throw Common::NotImplementedException (FromHere(),"Unsupported polyhedron form...");
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void TetraSpectralFVElementData::setCFLConvDiffRatio()
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
      throw Common::NotImplementedException (FromHere(),"Higher-order tetrahedral SV element not defined!");
    }
  }
}

//////////////////////////////////////////////////////////////////////

void TetraSpectralFVElementData::setPolyHedronBoxDecomposition(const vector< RealVector > polyHedronNodes,
                                                               vector< vector< RealVector > >& boxDecomp)
{
  // number of nodes
  const CFuint nbrNodes = polyHedronNodes.size();

  // set coordinates of box nodes in decomposition
  switch(nbrNodes)
  {
    case 4:
    {
      boxDecomp.resize(1);
      boxDecomp[0].resize(8,RealVector(3));

      // set box coordinates
      boxDecomp[0][0] = polyHedronNodes[0];
      boxDecomp[0][1] = polyHedronNodes[1];
      boxDecomp[0][2] = polyHedronNodes[2];
      boxDecomp[0][3] = polyHedronNodes[2];
      boxDecomp[0][4] = polyHedronNodes[3];
      boxDecomp[0][5] = polyHedronNodes[3];
      boxDecomp[0][6] = polyHedronNodes[3];
      boxDecomp[0][7] = polyHedronNodes[3];
    } break;
    case 8:
    {
      boxDecomp.resize(1);
      boxDecomp[0].resize(8,RealVector(3));

      // set box coordinates
      boxDecomp[0][0] = polyHedronNodes[0];
      boxDecomp[0][1] = polyHedronNodes[1];
      boxDecomp[0][2] = polyHedronNodes[2];
      boxDecomp[0][3] = polyHedronNodes[3];
      boxDecomp[0][4] = polyHedronNodes[4];
      boxDecomp[0][5] = polyHedronNodes[5];
      boxDecomp[0][6] = polyHedronNodes[6];
      boxDecomp[0][7] = polyHedronNodes[7];
    } break;
    case 11:
    {
      boxDecomp.resize(2);
      boxDecomp[0].resize(8,RealVector(3));
      boxDecomp[1].resize(8,RealVector(3));

      // additional node
      RealVector addNode = 0.5*(polyHedronNodes[3] + polyHedronNodes[7]);

      // set box coordinates
      boxDecomp[0][0] = polyHedronNodes[0];
      boxDecomp[0][1] = polyHedronNodes[1];
      boxDecomp[0][2] = polyHedronNodes[2];
      boxDecomp[0][3] = polyHedronNodes[3];
      boxDecomp[0][4] = polyHedronNodes[4];
      boxDecomp[0][5] = polyHedronNodes[5];
      boxDecomp[0][6] = polyHedronNodes[6];
      boxDecomp[0][7] = addNode;

      boxDecomp[1][0] = polyHedronNodes[4];
      boxDecomp[1][1] = polyHedronNodes[5];
      boxDecomp[1][2] = polyHedronNodes[6];
      boxDecomp[1][3] = addNode;
      boxDecomp[1][4] = polyHedronNodes[8];
      boxDecomp[1][5] = polyHedronNodes[9];
      boxDecomp[1][6] = polyHedronNodes[10];
      boxDecomp[1][7] = polyHedronNodes[7];
    } break;
  }
}

//////////////////////////////////////////////////////////////////////////////

void TetraSpectralFVElementData::setPolyGonQuadrilateralDecomposition(const vector< RealVector > polygonNodes,
                                                                      vector< vector< RealVector > >& quadDecomp)
{
  // number of nodes
  const CFuint nbrNodes = polygonNodes.size();

  // set coordinates of quadrilateral nodes in decomposition
  switch(nbrNodes)
  {
    case 3:
    {
      quadDecomp.resize(1);
      quadDecomp[0].resize(4,RealVector(3));

      // set quad coordinates
      quadDecomp[0][0] = polygonNodes[0];
      quadDecomp[0][1] = polygonNodes[1];
      quadDecomp[0][2] = polygonNodes[2];
      quadDecomp[0][3] = polygonNodes[2];
    } break;
    case 4:
    {
      quadDecomp.resize(1);
      quadDecomp[0].resize(4,RealVector(3));

      // set quad coordinates
      quadDecomp[0][0] = polygonNodes[0];
      quadDecomp[0][1] = polygonNodes[1];
      quadDecomp[0][2] = polygonNodes[2];
      quadDecomp[0][3] = polygonNodes[3];
    } break;
    case 5:
    {
      quadDecomp.resize(2);
      quadDecomp[0].resize(4,RealVector(3));
      quadDecomp[1].resize(4,RealVector(3));

      // additional node
      RealVector addNode = 0.5*(polygonNodes[2] + polygonNodes[3]);

      // set quad coordinates
      quadDecomp[0][0] = polygonNodes[0];
      quadDecomp[0][1] = polygonNodes[1];
      quadDecomp[0][2] = polygonNodes[2];
      quadDecomp[0][3] = addNode;

      quadDecomp[1][0] = polygonNodes[0];
      quadDecomp[1][1] = addNode;
      quadDecomp[1][2] = polygonNodes[3];
      quadDecomp[1][3] = polygonNodes[4];
    } break;
  }
}

//////////////////////////////////////////////////////////////////////////////

void TetraSpectralFVElementData::createFaceOutputPntCellMappedCoords()
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
  m_faceOutputPntFaceMappedCoords.resize(0);
  CFreal ksi = 0.0;
  for (CFuint iKsi = 0; iKsi < m_polyOrder+1; ++iKsi)
  {
    CFreal eta = 0.0;
    for (CFuint iEta = 0; iEta < m_polyOrder+1-iKsi; ++iEta)
    {
      RealVector mapCoord(2);
      mapCoord[KSI] = ksi;
      mapCoord[ETA] = eta;
      m_faceOutputPntFaceMappedCoords.push_back(mapCoord);
      eta += dKsi;
    }
    ksi += dKsi;
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
      const CFreal fun0 = m_faceOutputPntFaceMappedCoords[iPnt][KSI];
      const CFreal fun1 = m_faceOutputPntFaceMappedCoords[iPnt][ETA];
      const CFreal fun2 = 1 - fun0 - fun1;
      m_faceOutputPntCellMappedCoords[iFace].push_back(fun0*faceNodeCoords[0]+fun1*faceNodeCoords[1]+fun2*faceNodeCoords[2]);
    }
  }
}

//////////////////////////////////////////////////////////////////////

void TetraSpectralFVElementData::createFaceOutputPntConn()
{
  CFLog(WARN,"Output on boundary faces not implemented for 3D, connectivity is not created\n");
/*  CFAUTOTRACE;

  m_faceOutputPntConn.resize(0);
  for (CFuint iCell = 0; iCell < m_polyOrder; ++iCell)
  {
    vector<CFuint> cellNode(2);
    cellNode[0] = iCell;
    cellNode[1] = iCell+1;
    m_faceOutputPntConn.push_back(cellNode);
  }*/
}

//////////////////////////////////////////////////////////////////////

  } // namespace SpectralFV

} // namespace COOLFluiD

