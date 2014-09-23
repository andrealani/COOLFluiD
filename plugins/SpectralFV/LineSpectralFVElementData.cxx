#include <fstream>


#include "Common/CFLog.hh"
#include "Environment/DirPaths.hh"
#include "Common/NotImplementedException.hh"
#include "MathTools/RealMatrix.hh"
#include "SpectralFV/LineSpectralFVElementData.hh"
#include "MathTools/MathConsts.hh"

//////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Common;

namespace COOLFluiD {

  namespace SpectralFV {

//////////////////////////////////////////////////////////////////////

LineSpectralFVElementData::LineSpectralFVElementData() :
  SpectralFVElementData()
{
  m_shape = CFGeoShape::LINE;
  m_dimensionality = DIM_1D;
}

//////////////////////////////////////////////////////////////////////

LineSpectralFVElementData::LineSpectralFVElementData(CFPolyOrder::Type polyOrder)
{
  m_shape = CFGeoShape::LINE;
  m_dimensionality = DIM_1D;
  m_polyOrder = polyOrder;

  resetSpectralFVElementData();
}

//////////////////////////////////////////////////////////////////////

LineSpectralFVElementData::~LineSpectralFVElementData()
{
}

//////////////////////////////////////////////////////////////////////

void LineSpectralFVElementData::computeSVNodeLocalCoords()
{
  m_svNodeCoords.resize(3);

  // first node
  m_svNodeCoords[0].resize(1);
  m_svNodeCoords[0][KSI] = -1.0;

  // second node
  m_svNodeCoords[1].resize(1);
  m_svNodeCoords[1][KSI] =  1.0;
}

//////////////////////////////////////////////////////////////////////

void LineSpectralFVElementData::createLocalNodeCoord()
{
  CFAUTOTRACE;

  // Path where to read files from
  std::string filename;
  const std::string append = Environment::DirPaths::getInstance().getBaseDir().string() + "/data/spectralfv/";

  // number of local nodes
  const CFuint nbrLocalNodes = m_polyOrder+2;

  // resize the vector containing the node coordinates
  m_localNodeCoord.resize(nbrLocalNodes);
  for (CFuint iNode = 0; iNode < nbrLocalNodes; ++iNode)
  {
    m_localNodeCoord[iNode].resize(1);
  }

  switch (m_polyOrder)
  {
    case CFPolyOrder::ORDER0:
    {
      // fill the matrix containing the node coordinates
      m_localNodeCoord[0][KSI] = -1.0;

      m_localNodeCoord[1][KSI] = 1.0;
    } break;
    case CFPolyOrder::ORDER1:
    {
      // fill the matrix containing the node coordinates
      m_localNodeCoord[0][KSI] = -1.0;

      m_localNodeCoord[1][KSI] =  0.0;

      m_localNodeCoord[2][KSI] =  1.0;
    } break;
    case CFPolyOrder::ORDER2:
    {
      // partition parameters
      CFreal alpha;

      // read defining parameters from file
      filename = append + "SV3LINEDEF.DAT";
      ifstream inputFile;
      inputFile.open(filename.c_str(), ios::in);

      if (inputFile.is_open())
      {
        inputFile >> alpha;
      }
      else
      {
        alpha = 0.58;
      }
      inputFile.close();

      // fill the matrix containing the node coordinates
      m_localNodeCoord[0][KSI] = -1.0;

      m_localNodeCoord[1][KSI] = -alpha;

      m_localNodeCoord[2][KSI] =  alpha;

      m_localNodeCoord[3][KSI] =  1.0;
    } break;
    case CFPolyOrder::ORDER3:
    {
      // partition parameters
      CFreal alpha;

      // read defining parameters from file
      filename = append + "SV4LINEDEF.DAT";
      ifstream inputFile;
      inputFile.open(filename.c_str(), ios::in);

      if (inputFile.is_open())
      {
        inputFile >> alpha;
      }
      else
      {
        alpha = 0.78;
      }
      inputFile.close();

      // fill the matrix containing the node coordinates
      m_localNodeCoord[0][KSI] = -1.0;

      m_localNodeCoord[1][KSI] = -alpha;

      m_localNodeCoord[2][KSI] =  0.0;

      m_localNodeCoord[3][KSI] =  alpha;

      m_localNodeCoord[4][KSI] =  1.0;
    } break;
    default:
      throw Common::NotImplementedException (FromHere(),"Only solution orders up to 3 have been implemented for spectral FV!");
  }
}

//////////////////////////////////////////////////////////////////////

void LineSpectralFVElementData::createLocalFaceNodeConn()
{
  CFAUTOTRACE;

  // number of local faces
  const CFuint nbrLocalFaces = m_polyOrder + 2;

  // resize the local face-node connectivity
  m_localFaceNodeConn.resize(nbrLocalFaces);
  for (CFuint iFace = 0; iFace < nbrLocalFaces; ++iFace)
  {
    m_localFaceNodeConn[iFace].resize(1);
  }

  // external faces
  m_localFaceNodeConn[0][0] = 0;
  m_localFaceNodeConn[1][0] = m_polyOrder + 1;

  // internal faces
  const CFuint polyOrderUns = static_cast<CFuint>(m_polyOrder);
  for (CFuint iIntFace = 0; iIntFace < polyOrderUns; ++iIntFace)
  {
    m_localFaceNodeConn[iIntFace+2][0] = iIntFace + 1;
  }
}

//////////////////////////////////////////////////////////////////////

void LineSpectralFVElementData::createLocalCVNodeConn()
{
  CFAUTOTRACE;

  // number of CVs
  const CFuint nbrCVs = m_polyOrder + 1;

  // resize the local CV-node connectivity
  m_localCVNodeConn.resize(nbrCVs);
  for (CFuint iCV = 0; iCV < nbrCVs; ++iCV)
  {
    m_localCVNodeConn[iCV].resize(2);
  }

  // create local connectivity
  for (CFuint iCV = 0; iCV < nbrCVs; ++iCV)
  {
    m_localCVNodeConn[iCV][0] = iCV;
    m_localCVNodeConn[iCV][1] = iCV + 1;
  }
}

//////////////////////////////////////////////////////////////////////

void LineSpectralFVElementData::computeLocalFaceNormals()
{
  CFAUTOTRACE;

  // get number of local internal faces
  const CFuint nbrLocalIntFaces = m_polyOrder;

  // compute local internal face normals
  m_intFaceQuadPntNorm.resize(nbrLocalIntFaces);
  for (CFuint iIntFace = 2; iIntFace < nbrLocalIntFaces; ++iIntFace)
  {
    m_intFaceQuadPntNorm[iIntFace].resize(1,RealVector(1));
    m_intFaceQuadPntNorm[iIntFace][0][KSI] = -1.0;
  }
}

//////////////////////////////////////////////////////////////////////

void LineSpectralFVElementData::createLocalIntFaceCVConn()
{
  CFAUTOTRACE;

  // resize the local internal face-CV connectivity
  const CFuint polyOrderUns = static_cast<CFuint>(m_polyOrder);
  m_localIntFaceCVConn.resize(polyOrderUns);
  for (CFuint iIntFace = 0; iIntFace < polyOrderUns; ++iIntFace)
  {
    m_localIntFaceCVConn[iIntFace].resize(2);
  }

  // create local internal face-CV connectivity
  for (CFuint iIntFace = 0; iIntFace < polyOrderUns; ++iIntFace)
  {
    m_localIntFaceCVConn[iIntFace][LEFT ] = iIntFace + 1;
    m_localIntFaceCVConn[iIntFace][RIGHT] = iIntFace;
  }
}

//////////////////////////////////////////////////////////////////////

void LineSpectralFVElementData::computeVolumeFractionsOfCVs()
{
  CFAUTOTRACE;

  // number of CVs in SV
  const CFuint nbrCVs = m_localCVNodeConn.size();

  // resize the vector containing the CV volume fractions
  m_volFracCV.resize(nbrCVs);

  // compute CV volume fractions
  for (CFuint iCV = 0; iCV < nbrCVs; ++iCV)
  {
    const CFuint node1 = m_localCVNodeConn[iCV][0];
    const CFuint node2 = m_localCVNodeConn[iCV][1];
    m_volFracCV[iCV] = (m_localNodeCoord[node2][KSI] - m_localNodeCoord[node1][KSI])/2.0;
  }
}

//////////////////////////////////////////////////////////////////////

void LineSpectralFVElementData::createExtFaceNodeLocalCoords()
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////

void LineSpectralFVElementData::createLocalExtFaceCVConn()
{
  CFAUTOTRACE;

  // number of external faces
  const CFuint nbrExtFaces = 2;

  // resize m_localExtFaceCVConn
  m_localExtFaceCVConn.resize(nbrExtFaces);

  // create local external face - CV connectivity
  // first SV face
  m_localExtFaceCVConn[0] = 0;

  // second SV face
  m_localExtFaceCVConn[1] = m_polyOrder;
}

//////////////////////////////////////////////////////////////////////

void LineSpectralFVElementData::createSVFaceLocalExtFaceConn()
{
  CFAUTOTRACE;

  // resize m_svFaceLocalExtFaceConn
  m_svFaceLocalExtFaceConn.resize(2);
  m_svFaceLocalExtFaceConn[0].resize(1);
  m_svFaceLocalExtFaceConn[1].resize(1);

  // create SV face - local external face connectivity
  // first SV face
  m_svFaceLocalExtFaceConn[0][0] = 0;

  // second SV face
  m_svFaceLocalExtFaceConn[1][0] = 1;
}

//////////////////////////////////////////////////////////////////////

void LineSpectralFVElementData::createSVFaceNodeConnectivity()
{
  CFAUTOTRACE;

  // number of SV faces
  const CFuint nbrSVFaces = 2;

  // resize the variable
  m_svFaceNodeConn.resize(nbrSVFaces);
  for (CFuint iSVFace = 0; iSVFace < nbrSVFaces; ++iSVFace)
  {
    m_svFaceNodeConn[iSVFace].resize(1);
  }

  // fill the variable
  // first SV face
  m_svFaceNodeConn[0][0] = 0;

  // second SV face
  m_svFaceNodeConn[1][0] = 1;
}

//////////////////////////////////////////////////////////////////////

void LineSpectralFVElementData::computeExtFaceLocalNormals()
{
  CFAUTOTRACE;

  // number of SV faces
  const CFuint nbrSVFaces = 2;

  // resize the variable
  m_extFaceLocalNorm.resize(nbrSVFaces);
  for (CFuint iSVFace = 0; iSVFace < nbrSVFaces; ++iSVFace)
  {
    m_extFaceLocalNorm[iSVFace].resize(1);
  }

  // fill the variable
  // first SV face
  m_extFaceLocalNorm[0][KSI] = -1.0;

  // second SV face
  m_extFaceLocalNorm[1][KSI] = +1.0;
}

//////////////////////////////////////////////////////////////////////

void LineSpectralFVElementData::createSVFaceNodeConnectivityPerOrient()
{
  CFAUTOTRACE;

  // number of SV faces
  const CFuint nbrSVFaces = m_svFaceNodeConn.size();

  // number of possible orientations
  const CFuint nbrOrient = 3;

  // resize the variable
  m_svFaceNodeConnPerOrientation.resize(nbrOrient);
  for (CFuint iOrient = 0; iOrient < nbrOrient; ++iOrient)
  {
    m_svFaceNodeConnPerOrientation[iOrient].resize(2);
  }

  // fill the variable
  CFuint iOrient = 0;
  for (CFuint iSVFaceL = 0; iSVFaceL < nbrSVFaces; ++iSVFaceL)
  {
    for (CFuint iSVFaceR = iSVFaceL; iSVFaceR < nbrSVFaces; ++iSVFaceR, ++iOrient)
    {
      m_svFaceNodeConnPerOrientation[iOrient][LEFT ] = m_svFaceNodeConn[iSVFaceL];
      m_svFaceNodeConnPerOrientation[iOrient][RIGHT] = m_svFaceNodeConn[iSVFaceR];
    }
  }
}

//////////////////////////////////////////////////////////////////////

void LineSpectralFVElementData::createSVFaceNodeConnectivityPerOrientNoSymm()
{
  CFAUTOTRACE;

  // number of SV faces
  const CFuint nbrSVFaces = m_svFaceNodeConn.size();

  // number of possible orientations
  const CFuint nbrOrient = 4;

  // resize the variable
  m_svFaceNodeConnPerOrientationNoSymm.resize(nbrOrient);
  for (CFuint iOrient = 0; iOrient < nbrOrient; ++iOrient)
  {
    m_svFaceNodeConnPerOrientationNoSymm[iOrient].resize(2);
  }

  // fill the variable
  CFuint iOrient = 0;
  for (CFuint iSVFaceL = 0; iSVFaceL < nbrSVFaces; ++iSVFaceL)
  {
    for (CFuint iSVFaceR = 0; iSVFaceR < nbrSVFaces; ++iSVFaceR, ++iOrient)
    {
      m_svFaceNodeConnPerOrientationNoSymm[iOrient][LEFT ] = m_svFaceNodeConn[iSVFaceL];
      m_svFaceNodeConnPerOrientationNoSymm[iOrient][RIGHT] = m_svFaceNodeConn[iSVFaceR];
    }
  }
}

//////////////////////////////////////////////////////////////////////

/*void LineSpectralFVElementData::createSVFaceNodeCoordsPerOrient()
{
  CFAUTOTRACE;

  // number of orientations
  const CFuint nbrOrients = 3;

  // number of cell faces
  const CFuint nbrCellFaces = m_svFaceNodeCoords.size();

  // number of face nodes
  cf_assert(nbrCellFaces > 0);
  const CFuint nbrFaceNodes = m_svFaceNodeCoords[0].size();

  // resize m_svFaceNodeCoordsPerOrient
  m_svFaceNodeCoordsPerOrient.resize(nbrOrients);
  for (CFuint iOrient = 0; iOrient < nbrOrients; ++iOrient)
  {
    m_svFaceNodeCoordsPerOrient[iOrient].resize(2);
    for (CFuint iSide = 0; iSide < 2; ++iSide)
    {
      m_svFaceNodeCoordsPerOrient[iOrient][iSide].resize(nbrFaceNodes);
      for (CFuint iNode = 0; iNode < nbrFaceNodes; ++iNode)
      {
        m_svFaceNodeCoordsPerOrient[iOrient][iSide][iNode].resize(m_dimensionality);
      }
    }
  }

  // create SV face-node coordinates per face connectivity orientation
  CFuint iOrient = 0;
  for (CFuint iCellFaceLeft = 0; iCellFaceLeft < nbrCellFaces; ++iCellFaceLeft)
  {
    for (CFuint iCellFaceRight = iCellFaceLeft; iCellFaceRight < nbrCellFaces; ++iCellFaceRight, ++iOrient)
    {
      m_svFaceNodeCoordsPerOrient[iOrient][LEFT ][0] = m_svFaceNodeCoords[iCellFaceLeft ][0];
      m_svFaceNodeCoordsPerOrient[iOrient][RIGHT][0] = m_svFaceNodeCoords[iCellFaceRight][0];
    }
  }
}

//////////////////////////////////////////////////////////////////////
*/
void LineSpectralFVElementData::computeFaceFractionsOfCVs()
{
  CFAUTOTRACE;

  // resize the vector containing the CV face fractions
  m_faceFracCV.resize(2);
  m_faceFracCV[0].resize(1);
  m_faceFracCV[1].resize(1);

  // fill in CV fractions
  m_faceFracCV[0][0] = 1.0;
  m_faceFracCV[1][0] = 1.0;
}

//////////////////////////////////////////////////////////////////////

void LineSpectralFVElementData::createPolyExponents()
{
  CFAUTOTRACE;

  // number of polynomial terms
  const CFuint nbrPolyTerms = m_polyOrder + 1;

  // resize the variable
  m_polyExponents.resize(nbrPolyTerms);
  for (CFuint iTerm = 0; iTerm < nbrPolyTerms; ++iTerm)
  {
    m_polyExponents[iTerm].resize(1);
  }

  // define exponents
  for (CFuint iP = 0; iP < nbrPolyTerms; ++iP)
  {
    m_polyExponents[iP][KSI] = iP;
  }
}

//////////////////////////////////////////////////////////////////////

void LineSpectralFVElementData::computePolyCoefs()
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

    // get node coordinates
    vector< RealVector > cvNodeCoord(nbrCVNodes);
    for (CFuint iNode = 0; iNode < nbrCVNodes; ++iNode)
    {
      cvNodeCoord[iNode].resize(1);
      const CFuint nodeID = m_localCVNodeConn[iCV][iNode];
      cvNodeCoord[iNode] = m_localNodeCoord[nodeID];
    }

    // get quadrature nodes and wheights
    vector< RealVector > qNodeCoord = m_sIntegrator.getQuadPntsCoords(cvNodeCoord);
    vector< CFreal > qWheights = m_sIntegrator.getQuadPntsWheights(cvNodeCoord);
    const CFuint nbrQNodes = qWheights.size();

    for (CFuint iTerm = 0; iTerm < nbrOfCVs; ++iTerm)
    {
      polyCoefSystem(iCV,iTerm) = 0;

      for (CFuint iQNode = 0; iQNode < nbrQNodes; ++iQNode)
      {
        polyCoefSystem(iCV,iTerm) += qWheights[iQNode]*pow(qNodeCoord[iQNode][KSI],m_polyExponents[iTerm][KSI]);
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
    rhsLinSys(iCV,iCV) = 2.0*m_volFracCV[iCV];
  }

  // multiply imverted linear system matrix with rhs
  RealMatrix polyCoefMatrix(nbrOfCVs,nbrOfCVs);
  polyCoefMatrix = invLinSysMatrix*rhsLinSys;

  // resize m_polyCoefSFV
  m_polyCoefSFV.resize(nbrOfCVs);
  for (CFuint iPoly = 0; iPoly < nbrOfCVs; ++iPoly)
  {
    m_polyCoefSFV[iPoly].resize(nbrOfCVs);
  }

  // store polynomial coefficients in m_polyCoefSFV
  for (CFuint iPoly = 0; iPoly < nbrOfCVs; ++iPoly)
  {
    for (CFuint iTerm = 0; iTerm < nbrOfCVs; ++iTerm)
    {
      m_polyCoefSFV[iPoly][iTerm] = polyCoefMatrix(iTerm,iPoly);
    }
  }
}

//////////////////////////////////////////////////////////////////////

void LineSpectralFVElementData::createFaceFluxPntsConn()
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////

void LineSpectralFVElementData::createFluxPolyExponents()
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////

void LineSpectralFVElementData::setInterpolationNodeSet(const CFPolyOrder::Type order,vector< RealVector >& nodalSet)
{
  // number of nodes
  const CFuint nbrNodes = order+1;

  // fill the vector containing the flux polynomial node coordinates
  // Gauss-Lobatto nodes
  nodalSet.resize(nbrNodes);
  for (CFuint iNode = 0; iNode < nbrNodes; ++iNode)
  {
    nodalSet[iNode].resize(1);
    nodalSet[iNode][KSI] = -cos(iNode*MathTools::MathConsts::CFrealPi()/nbrNodes);
  }
}

//////////////////////////////////////////////////////////////////////

void LineSpectralFVElementData::setCFLConvDiffRatio()
{
  CFAUTOTRACE;

CFLog(WARN,"Convective/diffusive CFL ratio not set for linear elements");
  m_cflConvDiffRatio = 0.0;
}

//////////////////////////////////////////////////////////////////////

void LineSpectralFVElementData::createFaceOutputPntCellMappedCoords()
{
  // number of points on a face
// unused //  const CFuint nbrFacePnts = 1;

  // face mapped coordinates of uniform distribution of points
  m_faceOutputPntFaceMappedCoords.resize(0);
  RealVector mapCoord(0);
  m_faceOutputPntFaceMappedCoords.push_back(mapCoord);

  // compute cell mapped coordinates for distribution on each face
  m_faceOutputPntCellMappedCoords.resize(2);
  m_faceOutputPntCellMappedCoords[0].resize(0);
  m_faceOutputPntCellMappedCoords[0].push_back(m_svFaceNodeCoords[0][0]);
  m_faceOutputPntCellMappedCoords[1].resize(0);
  m_faceOutputPntCellMappedCoords[1].push_back(m_svFaceNodeCoords[1][0]);
}

//////////////////////////////////////////////////////////////////////

void LineSpectralFVElementData::createFaceOutputPntConn()
{
  CFAUTOTRACE;

  m_faceOutputPntConn.resize(0);
  vector<CFuint> cellNode(1);
  cellNode[0] = 0;
  m_faceOutputPntConn.push_back(cellNode);
}

//////////////////////////////////////////////////////////////////////

  } // namespace SpectralFV

} // namespace COOLFluiD

