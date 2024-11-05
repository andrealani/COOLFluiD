#include <fstream>


#include "Common/CFLog.hh"
#include "Environment/DirPaths.hh"
#include "Common/NotImplementedException.hh"
#include "MathTools/RealMatrix.hh"
#include "FluxReconstructionMethod/HexaFluxReconstructionElementData.hh"
#include "FluxReconstructionMethod/TensorProductGaussIntegrator.hh"

//////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Common;

namespace COOLFluiD {

  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////

HexaFluxReconstructionElementData::HexaFluxReconstructionElementData() :
  FluxReconstructionElementData()
{
  m_shape = CFGeoShape::HEXA;
  m_dimensionality = DIM_3D;
}

//////////////////////////////////////////////////////////////////////

HexaFluxReconstructionElementData::HexaFluxReconstructionElementData(CFPolyOrder::Type polyOrder)
{
  m_shape = CFGeoShape::HEXA;
  m_dimensionality = DIM_3D;
  m_polyOrder = polyOrder;
  
  m_solPntsLocalCoord1D.resize(polyOrder+1);
  m_flxPntsLocalCoord1D.resize(polyOrder+1);
  // Use a default solution and flux point distribution: Gauss Legendre.  
  std::vector<CFreal> coords;
  coords.resize(polyOrder+1);
  switch(polyOrder)
    {
      case CFPolyOrder::ORDER0:
      {
	coords[0] = 0.0;
      } break;
      case CFPolyOrder::ORDER1:
      {
	coords[0] = -1./sqrt(3.);
	coords[1] = +1./sqrt(3.);
      } break;
      case CFPolyOrder::ORDER2:
      {
	coords[0] = -sqrt(3./5.);
	coords[1] = 0.0;
	coords[2] = +sqrt(3./5.);
      } break;
      case CFPolyOrder::ORDER3:
      {
	coords[0] = -sqrt((3.+2.*sqrt(6./5.))/7.);
	coords[1] = -sqrt((3.-2.*sqrt(6./5.))/7.);
	coords[2] = +sqrt((3.-2.*sqrt(6./5.))/7.);
	coords[3] = +sqrt((3.+2.*sqrt(6./5.))/7.);
      } break;
      case CFPolyOrder::ORDER4:
      {
	coords[0] = -sqrt(5.+2.*sqrt(10./7.))/3.;
	coords[1] = -sqrt(5.-2.*sqrt(10./7.))/3.;
	coords[2] = 0.0;
	coords[3] = +sqrt(5.-2.*sqrt(10./7.))/3.;
	coords[4] = +sqrt(5.+2.*sqrt(10./7.))/3.;
      } break;
      case CFPolyOrder::ORDER5:
      {
        coords[0] = -0.9324695142031521;
	coords[1] = -0.6612093864662645;
	coords[2] = -0.2386191860831969;
	coords[3] = 0.2386191860831969;
	coords[4] = 0.6612093864662645;
	coords[5] = 0.9324695142031521;
      } break;
      case CFPolyOrder::ORDER6:
      {
        coords[0] = -0.949107912342759;
	coords[1] = -0.741531185599394;
	coords[2] = -0.405845151377397;
	coords[3] = 0.0;
	coords[4] = 0.405845151377397;
	coords[5] = 0.741531185599394;
	coords[6] = 0.949107912342759;
      } break;
      case CFPolyOrder::ORDER7:
      {
        coords[0] = -0.960289856497536;
	coords[1] = -0.796666477413627;
	coords[2] = -0.525532409916329;
	coords[3] = -0.183434642495650;
	coords[4] = 0.183434642495650;
	coords[5] = 0.525532409916329;
	coords[6] = 0.796666477413627;
	coords[7] = 0.960289856497536;
      } break;
      case CFPolyOrder::ORDER8:
      {
        coords[0] = -0.968160239507626;
	coords[1] = -0.836031107326636;
	coords[2] = -0.613371432700590;
	coords[3] = -0.324253423403809;
	coords[4] = 0.0;
	coords[5] = 0.324253423403809;
	coords[6] = 0.613371432700590;
	coords[7] = 0.836031107326636;
	coords[8] = 0.968160239507626;
      } break;
      case CFPolyOrder::ORDER9:
      {
        coords[0] = -0.973906528517172;
	coords[1] = -0.865063366688985;
	coords[2] = -0.679409568299024;
	coords[3] = -0.433395394129247;
	coords[4] = -0.148874338981631;
	coords[5] = 0.148874338981631;
	coords[6] = 0.433395394129247;
	coords[7] = 0.679409568299024;
	coords[8] = 0.865063366688985;
	coords[9] = 0.973906528517172;
      } break;
      case CFPolyOrder::ORDER10:
      {
        coords[0] = -0.978228658146057;
	coords[1] = -0.887062599768095;
	coords[2] = -0.730152005574049;
	coords[3] = -0.519096129110681;
	coords[4] = -0.269543155952345;
	coords[5] = 0.0;
	coords[6] = 0.269543155952345;
	coords[7] = 0.519096129110681;
	coords[8] = 0.730152005574049;
	coords[9] = 0.887062599768095;
	coords[10] = 0.978228658146057;
      } break;
      default:
      {
        throw Common::NotImplementedException (FromHere(),"Gauss Legendre not implemented for order "
                                      + StringOps::to_str(polyOrder) + ".");
      }
    }
  m_solPntsLocalCoord1D = coords;
  m_flxPntsLocalCoord1D = coords;

  resetFluxReconstructionElementData();
}

//////////////////////////////////////////////////////////////////////

HexaFluxReconstructionElementData::HexaFluxReconstructionElementData(CFPolyOrder::Type polyOrder, 
								     Common::SafePtr< BasePointDistribution > solPntDist, 
								     Common::SafePtr< BasePointDistribution > flxPntDist)
{
  m_shape = CFGeoShape::HEXA;
  m_dimensionality = DIM_3D;
  m_polyOrder = polyOrder;
  m_solPntsLocalCoord1D = solPntDist->getLocalCoords1D(polyOrder);
  m_flxPntsLocalCoord1D = flxPntDist->getLocalCoords1D(polyOrder);
  m_solPntDistribution = solPntDist;
  m_flxPntDistribution = flxPntDist;

  resetFluxReconstructionElementData();
}

//////////////////////////////////////////////////////////////////////

HexaFluxReconstructionElementData::~HexaFluxReconstructionElementData()
{
}

//////////////////////////////////////////////////////////////////////

void HexaFluxReconstructionElementData::createFlxPntsLocalCoords()
{
  CFAUTOTRACE;

  // number of flux points in 1D
  const CFuint nbrFlxPnts1D = m_flxPntsLocalCoord1D.size();

  // set flux point local coordinates
  m_flxPntsLocalCoords.resize(0);
  
    RealVector flxCoords(3);
    //faces ZTA=-1 and 1
    flxCoords[ZTA] = -1;
    for (CFuint iKsi = 0; iKsi < nbrFlxPnts1D; ++iKsi)
    {
        for (CFuint iEta = 0; iEta < nbrFlxPnts1D; ++iEta)
        {
          flxCoords[KSI] = m_flxPntsLocalCoord1D[iKsi];
          flxCoords[ETA] = m_flxPntsLocalCoord1D[iEta];
          m_flxPntsLocalCoords.push_back(flxCoords);
        }
    }
    flxCoords[ZTA] = 1;
    for (CFuint iKsi = 0; iKsi < nbrFlxPnts1D; ++iKsi)
    {
        for (CFuint iEta = 0; iEta < nbrFlxPnts1D; ++iEta)
        {
            flxCoords[KSI] = m_flxPntsLocalCoord1D[iKsi];
            flxCoords[ETA] = m_flxPntsLocalCoord1D[iEta];
            m_flxPntsLocalCoords.push_back(flxCoords);
        }
    }
    //faces KSI=-1 and 1
    flxCoords[KSI] = -1;
    for (CFuint iEta = 0; iEta < nbrFlxPnts1D; ++iEta)
    {
        for (CFuint iZta = 0; iZta < nbrFlxPnts1D; ++iZta)
        {
            flxCoords[ZTA] = m_flxPntsLocalCoord1D[iZta];
            flxCoords[ETA] = m_flxPntsLocalCoord1D[iEta];
            m_flxPntsLocalCoords.push_back(flxCoords);
        }
    }
    flxCoords[KSI] = 1;
    for (CFuint iEta = 0; iEta < nbrFlxPnts1D; ++iEta)
    {
        for (CFuint iZta = 0; iZta < nbrFlxPnts1D; ++iZta)
        {
            flxCoords[ZTA] = m_flxPntsLocalCoord1D[iZta];
            flxCoords[ETA] = m_flxPntsLocalCoord1D[iEta];
            m_flxPntsLocalCoords.push_back(flxCoords);
        }
    }
    //faces ETA=-1 and 1
    flxCoords[ETA] = -1;
    for (CFuint iZta = 0; iZta < nbrFlxPnts1D; ++iZta)
    {
        for (CFuint iKsi = 0; iKsi < nbrFlxPnts1D; ++iKsi)
        {
            flxCoords[ZTA] = m_flxPntsLocalCoord1D[iZta];
            flxCoords[KSI] = m_flxPntsLocalCoord1D[iKsi];
            m_flxPntsLocalCoords.push_back(flxCoords);
        }
    }
    flxCoords[ETA] = 1;
    for (CFuint iZta = 0; iZta < nbrFlxPnts1D; ++iZta)
    {
        for (CFuint iKsi = 0; iKsi < nbrFlxPnts1D; ++iKsi)
        {
            flxCoords[ZTA] = m_flxPntsLocalCoord1D[iZta];
            flxCoords[KSI] = m_flxPntsLocalCoord1D[iKsi];
            m_flxPntsLocalCoords.push_back(flxCoords);
        }
    }
  cf_assert(m_flxPntsLocalCoords.size() == 6*nbrFlxPnts1D*nbrFlxPnts1D);
  
}

//////////////////////////////////////////////////////////////////////

void HexaFluxReconstructionElementData::createSolPntsLocalCoords()
{
  CFAUTOTRACE;

  // number of solution points in 1D
  const CFuint nbrSolPnts1D = m_solPntsLocalCoord1D.size();

  // set solution point local coordinates
  m_solPntsLocalCoords.resize(0);
  for (CFuint iKsi = 0; iKsi < nbrSolPnts1D; ++iKsi)
  {
    for (CFuint iEta = 0; iEta < nbrSolPnts1D; ++iEta)
    {
      for (CFuint iZta = 0; iZta < nbrSolPnts1D; ++iZta)
      {
        RealVector solCoords(3);
        solCoords[KSI] = m_solPntsLocalCoord1D[iKsi];
        solCoords[ETA] = m_solPntsLocalCoord1D[iEta];
        solCoords[ZTA] = m_solPntsLocalCoord1D[iZta];
        m_solPntsLocalCoords.push_back(solCoords);
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////

void HexaFluxReconstructionElementData::createFaceFlxPntsFaceLocalCoords()
{
  CFAUTOTRACE;

  // number of solution points in 1D
  const CFuint nbrFlxPnts1D = m_flxPntsLocalCoord1D.size();

  // set face flux point face local coordinates
  m_faceFlxPntsFaceLocalCoords.resize(0);
  for (CFuint iKsi = 0; iKsi < nbrFlxPnts1D; ++iKsi)
  {
    for (CFuint iEta = 0; iEta < nbrFlxPnts1D; ++iEta)
    {
      RealVector flxCoord(2);
      flxCoord[KSI] = m_flxPntsLocalCoord1D[iKsi];
      flxCoord[ETA] = m_flxPntsLocalCoord1D[iEta];
      m_faceFlxPntsFaceLocalCoords.push_back(flxCoord);
    }
  }
}

//////////////////////////////////////////////////////////////////////

void HexaFluxReconstructionElementData::createSolPolyExponents()
{
  CFAUTOTRACE;

  // number of solution points in 1D
  const CFuint nbrSolPnts1D = m_solPntsLocalCoord1D.size();

  // define exponents
  m_solPolyExponents.resize(0);
  for (CFuint iKsi = 0; iKsi < nbrSolPnts1D; ++iKsi)
  {
    for (CFuint iEta = 0; iEta < nbrSolPnts1D; ++iEta)
    {
      for (CFuint iZta = 0; iZta < nbrSolPnts1D; ++iZta)
      {
        vector< CFint > solPolyExps(3);
        solPolyExps[KSI] = iKsi;
        solPolyExps[ETA] = iEta;
        solPolyExps[ZTA] = iZta;
        m_solPolyExponents.push_back(solPolyExps);
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////

void HexaFluxReconstructionElementData::createNodePolyExponents()
{
  CFAUTOTRACE;

  // number of solution points in 1D
  const CFuint nbrNodes1D = 2;

  // define exponents
  m_nodePolyExponents.resize(0);
  for (CFuint iKsi = 0; iKsi < nbrNodes1D; ++iKsi)
  {
    for (CFuint iEta = 0; iEta < nbrNodes1D; ++iEta)
    {
      for (CFuint iZta = 0; iZta < nbrNodes1D; ++iZta)
      {
        vector< CFint > nodePolyExps(3);
        nodePolyExps[KSI] = iKsi;
        nodePolyExps[ETA] = iEta;
        nodePolyExps[ZTA] = iZta;
        m_nodePolyExponents.push_back(nodePolyExps);
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////

void HexaFluxReconstructionElementData::createFaceFluxPntsConn()
{
  CFAUTOTRACE;

  // number of flux points in 1D
  const CFuint nbrFlxPnts1D = m_flxPntsLocalCoord1D.size();

  // resize m_faceFlxPntConn
  m_faceFlxPntConn.resize(6);

  // variable holding the face index
  CFuint faceIdx = 0;

  // zeroth face
  for (CFuint iEta = 0; iEta < nbrFlxPnts1D; ++iEta)
  {
    for (CFuint iKsi = 0; iKsi < nbrFlxPnts1D; ++iKsi)
    {
      m_faceFlxPntConn[faceIdx].push_back(nbrFlxPnts1D*iKsi + iEta);
    }
  }
  ++faceIdx;

  // first face
  for (CFuint iKsi = 0; iKsi < nbrFlxPnts1D; ++iKsi)
  {
    for (CFuint iEta = 0; iEta < nbrFlxPnts1D; ++iEta)
    {
      m_faceFlxPntConn[faceIdx].push_back(nbrFlxPnts1D*nbrFlxPnts1D + nbrFlxPnts1D*iKsi + iEta);
    }
  }
  ++faceIdx;

  // second face
  for (CFuint iKsi = 0; iKsi < nbrFlxPnts1D; ++iKsi)
  {
    for (CFuint iZta = 0; iZta < nbrFlxPnts1D; ++iZta)
    {
      m_faceFlxPntConn[faceIdx].push_back(4*nbrFlxPnts1D*nbrFlxPnts1D + nbrFlxPnts1D*iZta + iKsi);
    }
  }
  ++faceIdx;

  // third face
  for (CFuint iEta = 0; iEta < nbrFlxPnts1D; ++iEta)
  {
    for (CFuint iZta = 0; iZta < nbrFlxPnts1D; ++iZta)
    {
      m_faceFlxPntConn[faceIdx].push_back(3*nbrFlxPnts1D*nbrFlxPnts1D + nbrFlxPnts1D*iEta + iZta); 
    }
  }
  ++faceIdx;
  
  // fourth face
  for (CFuint iKsi = 0; iKsi < nbrFlxPnts1D; ++iKsi)
  {
    const CFuint idxKsi = nbrFlxPnts1D-iKsi-1;
    for (CFuint iZta = 0; iZta < nbrFlxPnts1D; ++iZta)
    {
      m_faceFlxPntConn[faceIdx].push_back(5*nbrFlxPnts1D*nbrFlxPnts1D + nbrFlxPnts1D*iZta + idxKsi); 
    }
  }
  ++faceIdx;
  
  // fifth face
  for (CFuint iZta = 0; iZta < nbrFlxPnts1D; ++iZta)
  {
    for (CFuint iEta = 0; iEta < nbrFlxPnts1D; ++iEta)
    {
      m_faceFlxPntConn[faceIdx].push_back(2*nbrFlxPnts1D*nbrFlxPnts1D + nbrFlxPnts1D*iEta + iZta); 
    }
  }
}

//////////////////////////////////////////////////////////////////////

void HexaFluxReconstructionElementData::createFaceFluxPntsConnPerOrient()
{
  CFAUTOTRACE;
  // number of faces
  const CFuint nbrFaces = m_faceNodeConn.size();
  CFLog(VERBOSE,"nbrFaces:" << nbrFaces << "\n");
  // number of solution points in 1D
  const CFuint nbrSolPnts1D = m_solPntsLocalCoord1D.size();
  CFLog(VERBOSE,"nbrSolPnts1D:" << nbrSolPnts1D << "\n");
  // number of face flux points
  const CFuint nbrFaceFlxPnts = getNbrOfFaceFlxPnts();
  CFLog(VERBOSE,"nbrFaceFlxPnts:" << nbrFaceFlxPnts << "\n");
  // flux point indexes for inverted face
  vector< CFuint > invFlxIdxs;
  for (CFuint iKsi = 0; iKsi < nbrSolPnts1D; ++iKsi)
  {
    const CFuint idxKsi = nbrSolPnts1D - iKsi - 1;
    for (CFuint iEta = 0; iEta < nbrSolPnts1D; ++iEta)
    {
      invFlxIdxs.push_back(nbrSolPnts1D*idxKsi+iEta);
    }
  }
  // number of rotatable flux point groups
  const CFuint nbrRotFlxGroups = nbrFaceFlxPnts/4;
  CFLog(VERBOSE,"nbrRotFlxGroups:" << nbrRotFlxGroups << "\n");
  // maximum number of flux points in a line of flux points
  const CFuint maxNbrFlxPntsInLine = (nbrSolPnts1D+1)/2;
  CFLog(VERBOSE,"maxNbrFlxPntsInLine:" << maxNbrFlxPntsInLine << "\n");
  // storage of flux point rotatable groups
  vector< vector< CFuint > > rotFlxIdxs(nbrRotFlxGroups);
  CFuint iRotGroup = 0;
  for (CFuint iFlxLine = 0; iRotGroup < nbrRotFlxGroups; ++iFlxLine)
  {
    for (CFuint iFlx = 0;
         iFlx < maxNbrFlxPntsInLine && iRotGroup < nbrRotFlxGroups;
         ++iFlx, ++iRotGroup)
    {
      const CFuint idxFlxLine = nbrSolPnts1D - 1 - iFlxLine;
      const CFuint idxFlx     = nbrSolPnts1D - 1 - iFlx    ;

      rotFlxIdxs[iRotGroup].push_back(nbrSolPnts1D*iFlx      +iFlxLine  );
      rotFlxIdxs[iRotGroup].push_back(nbrSolPnts1D*idxFlxLine+iFlx      );
      rotFlxIdxs[iRotGroup].push_back(nbrSolPnts1D*idxFlx    +idxFlxLine);
      rotFlxIdxs[iRotGroup].push_back(nbrSolPnts1D*iFlxLine  +idxFlx    );
    }
  }

  // number of possible orientations
  const CFuint nbrOrient = 84;

  // create data structure
  m_faceFlxPntConnPerOrient.resize(nbrOrient);
  CFuint iOrient = 0;
  for (CFuint iFaceL = 0; iFaceL < nbrFaces; ++iFaceL)
  {
    vector< CFuint > faceFlxConnL = m_faceFlxPntConn[iFaceL];
    for (CFuint iFaceR = iFaceL; iFaceR < nbrFaces; ++iFaceR)
    {
      vector< CFuint > faceFlxConnR = m_faceFlxPntConn[iFaceR];
      for (CFuint iRot = 0; iRot < 4; ++iRot, ++iOrient)
      {
        m_faceFlxPntConnPerOrient[iOrient].resize(2);
        for (CFuint iFlx = 0; iFlx < nbrFaceFlxPnts; ++iFlx)
        {
          m_faceFlxPntConnPerOrient[iOrient][LEFT ]
              .push_back(faceFlxConnL[iFlx            ]);
          m_faceFlxPntConnPerOrient[iOrient][RIGHT]
              .push_back(faceFlxConnR[invFlxIdxs[iFlx]]);
        }
        // rotate the right face
        for (CFuint iRotGroup = 0; iRotGroup < nbrRotFlxGroups; ++iRotGroup)
        {
          // indexes of flux points to be rotated
          const CFuint flx0Idx = rotFlxIdxs[iRotGroup][0];
          const CFuint flx1Idx = rotFlxIdxs[iRotGroup][1];
          const CFuint flx2Idx = rotFlxIdxs[iRotGroup][2];
          const CFuint flx3Idx = rotFlxIdxs[iRotGroup][3];
          //CFLog(VERBOSE,"flxrotIdx" << flx0Idx << flx1Idx << flx2Idx << flx3Idx << "\n");

          // rotate flux points
          const CFuint swap = invFlxIdxs[flx3Idx];
          invFlxIdxs[flx3Idx] = invFlxIdxs[flx2Idx];
          invFlxIdxs[flx2Idx] = invFlxIdxs[flx1Idx];
          invFlxIdxs[flx1Idx] = invFlxIdxs[flx0Idx];
          invFlxIdxs[flx0Idx] = swap;
        }
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////

void HexaFluxReconstructionElementData::createCellNodeCoords()
{
  CFAUTOTRACE;

  m_cellNodeCoords.resize(8);

  // first node
  m_cellNodeCoords[0].resize(3);
  m_cellNodeCoords[0][KSI] = -1.0;
  m_cellNodeCoords[0][ETA] = -1.0;
  m_cellNodeCoords[0][ZTA] = -1.0;

  // second node
  m_cellNodeCoords[1].resize(3);
  m_cellNodeCoords[1][KSI] = +1.0;
  m_cellNodeCoords[1][ETA] = -1.0;
  m_cellNodeCoords[1][ZTA] = -1.0;

  // third node
  m_cellNodeCoords[2].resize(3);
  m_cellNodeCoords[2][KSI] = +1.0;
  m_cellNodeCoords[2][ETA] = +1.0;
  m_cellNodeCoords[2][ZTA] = -1.0;

  // fourth node
  m_cellNodeCoords[3].resize(3);
  m_cellNodeCoords[3][KSI] = -1.0;
  m_cellNodeCoords[3][ETA] = +1.0;
  m_cellNodeCoords[3][ZTA] = -1.0;

  // fifth node
  m_cellNodeCoords[4].resize(3);
  m_cellNodeCoords[4][KSI] = -1.0;
  m_cellNodeCoords[4][ETA] = -1.0;
  m_cellNodeCoords[4][ZTA] = +1.0;

  // sixth node
  m_cellNodeCoords[5].resize(3);
  m_cellNodeCoords[5][KSI] = +1.0;
  m_cellNodeCoords[5][ETA] = -1.0;
  m_cellNodeCoords[5][ZTA] = +1.0;

  // seventh node
  m_cellNodeCoords[6].resize(3);
  m_cellNodeCoords[6][KSI] = +1.0;
  m_cellNodeCoords[6][ETA] = +1.0;
  m_cellNodeCoords[6][ZTA] = +1.0;

  // eighth node
  m_cellNodeCoords[7].resize(3);
  m_cellNodeCoords[7][KSI] = -1.0;
  m_cellNodeCoords[7][ETA] = +1.0;
  m_cellNodeCoords[7][ZTA] = +1.0;
}

//////////////////////////////////////////////////////////////////////

void HexaFluxReconstructionElementData::createFaceNodeConnectivity()
{
  CFAUTOTRACE;

  m_faceNodeConn.resize(6);

  m_faceNodeConn[0].resize(4);
  m_faceNodeConn[0][0] = 0;
  m_faceNodeConn[0][1] = 3;
  m_faceNodeConn[0][2] = 2;
  m_faceNodeConn[0][3] = 1;

  m_faceNodeConn[1].resize(4);
  m_faceNodeConn[1][0] = 4;
  m_faceNodeConn[1][1] = 5;
  m_faceNodeConn[1][2] = 6;
  m_faceNodeConn[1][3] = 7;

  m_faceNodeConn[2].resize(4);
  m_faceNodeConn[2][0] = 0;
  m_faceNodeConn[2][1] = 1;
  m_faceNodeConn[2][2] = 5;
  m_faceNodeConn[2][3] = 4;

  m_faceNodeConn[3].resize(4);
  m_faceNodeConn[3][0] = 1;
  m_faceNodeConn[3][1] = 2;
  m_faceNodeConn[3][2] = 6;
  m_faceNodeConn[3][3] = 5;

  m_faceNodeConn[4].resize(4);
  m_faceNodeConn[4][0] = 2;
  m_faceNodeConn[4][1] = 3;
  m_faceNodeConn[4][2] = 7;
  m_faceNodeConn[4][3] = 6;

  m_faceNodeConn[5].resize(4);
  m_faceNodeConn[5][0] = 0;
  m_faceNodeConn[5][1] = 4;
  m_faceNodeConn[5][2] = 7;
  m_faceNodeConn[5][3] = 3;
}

//////////////////////////////////////////////////////////////////////

void HexaFluxReconstructionElementData::createFaceMappedCoordDir()
{
  CFAUTOTRACE;

  m_faceMappedCoordDir.resize(6);

  m_faceMappedCoordDir[0] = -1;
  m_faceMappedCoordDir[1] = +1;
  m_faceMappedCoordDir[2] = -1;
  m_faceMappedCoordDir[3] = +1;
  m_faceMappedCoordDir[4] = +1;
  m_faceMappedCoordDir[5] = -1;
}

//////////////////////////////////////////////////////////////////////

void HexaFluxReconstructionElementData::createFluxPntFluxDim()
{
  CFAUTOTRACE;

  m_flxPntFlxDim.resize(6*m_flxPntsLocalCoord1D.size()*m_flxPntsLocalCoord1D.size());
  
  for (CFuint iFace = 0; iFace < m_faceFlxPntConn.size(); ++iFace)
  {
    for (CFuint iFlx = 0; iFlx < m_faceFlxPntConn[iFace].size(); ++iFlx)
    {
      if (iFace == 3 || iFace == 5)
      {
        m_flxPntFlxDim[m_faceFlxPntConn[iFace][iFlx]] = 0;
      }
      else if (iFace == 2 || iFace == 4)
      {
	m_flxPntFlxDim[m_faceFlxPntConn[iFace][iFlx]] = 1;
      }
      else
      {
	m_flxPntFlxDim[m_faceFlxPntConn[iFace][iFlx]] = 2;
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////

void HexaFluxReconstructionElementData::createFaceNormals()
{
  CFAUTOTRACE;

  m_faceNormals.resize(6);
  for (CFuint iFace = 0; iFace < 6; ++iFace)
  {
    m_faceNormals[iFace].resize(3);
  }

  m_faceNormals[0][0] = 0.;
  m_faceNormals[0][1] = 0.;
  m_faceNormals[0][2] = -1.;
  m_faceNormals[1][0] = 0.;
  m_faceNormals[1][1] = 0.;
  m_faceNormals[1][2] = 1.;
  m_faceNormals[2][0] = 0.;
  m_faceNormals[2][1] = -1.;
  m_faceNormals[2][2] = 0.;
  m_faceNormals[3][0] = 1.;
  m_faceNormals[3][1] = 0.;
  m_faceNormals[3][2] = 0.;
  m_faceNormals[4][0] = 0.;
  m_faceNormals[4][1] = 1.;
  m_faceNormals[4][2] = 0.;
  m_faceNormals[5][0] = -1.;
  m_faceNormals[5][1] = 0.;
  m_faceNormals[5][2] = 0.;
}

//////////////////////////////////////////////////////////////////////

void HexaFluxReconstructionElementData::createFaceNodeConnectivityPerOrient()
{
  CFAUTOTRACE;

  // number of faces
  const CFuint nbrFaces = m_faceNodeConn.size();

  // number of possible orientations
  const CFuint nbrOrient = 84;

  // resize the variables
  m_faceNodeConnPerOrient.resize(nbrOrient);
  m_faceConnPerOrient.resize(nbrOrient);
  m_faceMappedCoordDirPerOrient.resize(nbrOrient);
  for (CFuint iOrient = 0; iOrient < nbrOrient; ++iOrient)
  {
    m_faceNodeConnPerOrient[iOrient].resize(2);
    m_faceConnPerOrient[iOrient].resize(2);
    m_faceMappedCoordDirPerOrient[iOrient].resize(2);
    for (CFuint iSide = 0; iSide < 2; ++iSide)
    {
      m_faceNodeConnPerOrient[iOrient][iSide].resize(4);
    }
  }

  // fill the variable
  CFuint iOrient = 0;
  for (CFuint iFaceL = 0; iFaceL < nbrFaces; ++iFaceL)
  {
    vector< CFuint > faceNodesL = m_faceNodeConn[iFaceL];
    for (CFuint iFaceR = iFaceL; iFaceR < nbrFaces; ++iFaceR)
    {
      vector< CFuint > faceNodesR(4);
      faceNodesR[0] = m_faceNodeConn[iFaceR][1];
      faceNodesR[1] = m_faceNodeConn[iFaceR][0];
      faceNodesR[2] = m_faceNodeConn[iFaceR][3];
      faceNodesR[3] = m_faceNodeConn[iFaceR][2];

      for (CFuint iRot = 0; iRot < 4; ++iRot, ++iOrient)
      {
        m_faceConnPerOrient[iOrient][LEFT ] = iFaceL;
        m_faceConnPerOrient[iOrient][RIGHT] = iFaceR;

        m_faceMappedCoordDirPerOrient[iOrient][LEFT ] = m_faceMappedCoordDir[iFaceL];
        m_faceMappedCoordDirPerOrient[iOrient][RIGHT] = -m_faceMappedCoordDir[iFaceR];

        for (CFuint iNode = 0; iNode < 4; ++iNode)
        {
          m_faceNodeConnPerOrient[iOrient][LEFT ][iNode] = faceNodesL[iNode];
          m_faceNodeConnPerOrient[iOrient][RIGHT][iNode] = faceNodesR[iNode];
        }

        // rotate nodes of right face to new orientation
        CFuint swap = faceNodesR[3];
        faceNodesR[3] = faceNodesR[2];
        faceNodesR[2] = faceNodesR[1];
        faceNodesR[1] = faceNodesR[0];
        faceNodesR[0] = swap;
      }

      cf_assert(faceNodesR[0] == m_faceNodeConn[iFaceR][1]);
      cf_assert(faceNodesR[1] == m_faceNodeConn[iFaceR][0]);
      cf_assert(faceNodesR[2] == m_faceNodeConn[iFaceR][3]);
      cf_assert(faceNodesR[3] == m_faceNodeConn[iFaceR][2]);
    }
  }
  cf_assert(iOrient == nbrOrient);
}

//////////////////////////////////////////////////////////////////////

void HexaFluxReconstructionElementData::createFaceIntegrationCoefs()
{
  CFAUTOTRACE;

  // number of solution points in 1D
  const CFuint nbrFlxPnts1D = m_flxPntsLocalCoord1D.size();

  // number of flux points on a face
  const CFuint nbrFlxPnts = nbrFlxPnts1D*nbrFlxPnts1D;

  // resize m_faceIntegrationCoefs
  m_faceIntegrationCoefs.resize(nbrFlxPnts);

  // create TensorProductGaussIntegrator
  TensorProductGaussIntegrator tpIntegrator(DIM_2D,m_polyOrder);

  // create face node local coordinates
  vector< RealVector > nodeCoord(4);
  nodeCoord[0].resize(2);
  nodeCoord[0][KSI] = -1.0;
  nodeCoord[0][ETA] = -1.0;
  nodeCoord[1].resize(2);
  nodeCoord[1][KSI] = +1.0;
  nodeCoord[1][ETA] = -1.0;
  nodeCoord[2].resize(2);
  nodeCoord[2][KSI] = +1.0;
  nodeCoord[2][ETA] = +1.0;
  nodeCoord[3].resize(2);
  nodeCoord[3][KSI] = -1.0;
  nodeCoord[3][ETA] = +1.0;

  // get quadrature point coordinates and wheights
  vector< RealVector > quadPntCoords   = tpIntegrator.getQuadPntsCoords  (nodeCoord);
  vector< CFreal     > quadPntWheights = tpIntegrator.getQuadPntsWheights(nodeCoord);
  const CFuint nbrQPnts = quadPntCoords.size();
  cf_assert(quadPntWheights.size() == nbrQPnts);

  // compute the coefficients for integration over a face
  // loop over flux points
  CFuint iFlx = 0;
  for (CFuint iFlxKsi = 0; iFlxKsi < nbrFlxPnts1D; ++iFlxKsi)
  {
    const CFreal ksiFlx = m_flxPntsLocalCoord1D[iFlxKsi];
    for (CFuint iFlxEta = 0; iFlxEta < nbrFlxPnts1D; ++iFlxEta, ++iFlx)
    {
      const CFreal etaFlx = m_flxPntsLocalCoord1D[iFlxEta];

      m_faceIntegrationCoefs[iFlx] = 0.0;
      for (CFuint iQPnt = 0; iQPnt < nbrQPnts; ++iQPnt)
      {
        // quadrature point local coordinate on the face
        const CFreal ksiQPnt = quadPntCoords[iQPnt][KSI];
        const CFreal etaQPnt = quadPntCoords[iQPnt][ETA];

        // evaluate polynomial value in quadrature point
        CFreal quadPntPolyVal = 1.;
        for (CFuint iFac = 0; iFac < nbrFlxPnts1D; ++iFac)
        {
          if (iFac != iFlxKsi)
          {
            const CFreal ksiFac = m_flxPntsLocalCoord1D[iFac];
            quadPntPolyVal *= (ksiQPnt-ksiFac)/(ksiFlx-ksiFac);
          }
        }
        for (CFuint iFac = 0; iFac < nbrFlxPnts1D; ++iFac)
        {
          if (iFac != iFlxEta)
          {
            const CFreal etaFac = m_flxPntsLocalCoord1D[iFac];
            quadPntPolyVal *= (etaQPnt-etaFac)/(etaFlx-etaFac);
          }
        }

        // add contribution of quadrature point to integration coefficient
        m_faceIntegrationCoefs[iFlx] += quadPntWheights[iQPnt]*quadPntPolyVal;
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////

void HexaFluxReconstructionElementData::createCellAvgSolCoefs()
{
  CFAUTOTRACE;

  // number of solution points
  const CFuint nbrSolPnts = getNbrOfSolPnts();

  // resize m_cellAvgSolCoefs
  m_cellAvgSolCoefs.resize(nbrSolPnts);

  // create TensorProductGaussIntegrator
  TensorProductGaussIntegrator tpIntegrator(DIM_3D,m_polyOrder);

  // create cell node local coordinates
  vector< RealVector > nodeCoord(8);
  nodeCoord[0].resize(3);
  nodeCoord[0][KSI] = -1.0;
  nodeCoord[0][ETA] = -1.0;
  nodeCoord[0][ZTA] = -1.0;
  nodeCoord[1].resize(3);
  nodeCoord[1][KSI] = +1.0;
  nodeCoord[1][ETA] = -1.0;
  nodeCoord[1][ZTA] = -1.0;
  nodeCoord[2].resize(3);
  nodeCoord[2][KSI] = +1.0;
  nodeCoord[2][ETA] = +1.0;
  nodeCoord[2][ZTA] = -1.0;
  nodeCoord[3].resize(3);
  nodeCoord[3][KSI] = -1.0;
  nodeCoord[3][ETA] = +1.0;
  nodeCoord[3][ZTA] = -1.0;
  nodeCoord[4].resize(3);
  nodeCoord[4][KSI] = -1.0;
  nodeCoord[4][ETA] = -1.0;
  nodeCoord[4][ZTA] = +1.0;
  nodeCoord[5].resize(3);
  nodeCoord[5][KSI] = +1.0;
  nodeCoord[5][ETA] = -1.0;
  nodeCoord[5][ZTA] = +1.0;
  nodeCoord[6].resize(3);
  nodeCoord[6][KSI] = +1.0;
  nodeCoord[6][ETA] = +1.0;
  nodeCoord[6][ZTA] = +1.0;
  nodeCoord[7].resize(3);
  nodeCoord[7][KSI] = -1.0;
  nodeCoord[7][ETA] = +1.0;
  nodeCoord[7][ZTA] = +1.0;

  // get quadrature point coordinates and wheights
  vector< RealVector > quadPntCoords   = tpIntegrator.getQuadPntsCoords  (nodeCoord);
  vector< CFreal     > quadPntWheights = tpIntegrator.getQuadPntsWheights(nodeCoord);
  const CFuint nbrQPnts = quadPntCoords.size();
  cf_assert(quadPntWheights.size() == nbrQPnts);

  // get the solution polynomial values at the quadrature points
  vector< vector< CFreal > > quadPntPolyVals = getSolPolyValsAtNode(quadPntCoords);

  // compute the coefficients for integration over a face
  // loop over solution points
  for (CFuint iSol = 0; iSol < nbrSolPnts; ++iSol)
  {
    m_cellAvgSolCoefs[iSol] = 0.0;
    for (CFuint iQPnt = 0; iQPnt < nbrQPnts; ++iQPnt)
    {
      m_cellAvgSolCoefs[iSol] += quadPntWheights[iQPnt]*quadPntPolyVals[iQPnt][iSol];
    }
    m_cellAvgSolCoefs[iSol] *= 0.125;
  }
}

//////////////////////////////////////////////////////////////////////

void HexaFluxReconstructionElementData::createCellCenterDerivCoefs()
{
  CFAUTOTRACE;

  // center coordinate
  vector< RealVector > centerCoord(1,RealVector(3));
  centerCoord[0][KSI] = 0.0;
  centerCoord[0][ETA] = 0.0;
  centerCoord[0][ZTA] = 0.0;

  vector< vector< vector< CFreal > > > polyDerivs =
                                            getSolPolyDerivsAtNode(centerCoord);

  // number of solution points
  const CFuint nbrSolPnts = getNbrOfSolPnts();

  // set polynomial derivatives
  m_cellCenterDerivCoefs.resize(3);
  m_cellCenterDerivCoefs[KSI].resize(nbrSolPnts);
  m_cellCenterDerivCoefs[ETA].resize(nbrSolPnts);
  m_cellCenterDerivCoefs[ZTA].resize(nbrSolPnts);
  for (CFuint iPoly = 0; iPoly < nbrSolPnts; ++iPoly)
  {
    m_cellCenterDerivCoefs[KSI][iPoly] = polyDerivs[0][KSI][iPoly];
    m_cellCenterDerivCoefs[ETA][iPoly] = polyDerivs[0][ETA][iPoly];
    m_cellCenterDerivCoefs[ZTA][iPoly] = polyDerivs[0][ZTA][iPoly];
  }
}

//////////////////////////////////////////////////////////////////////

void HexaFluxReconstructionElementData::setInterpolationNodeSet(const CFPolyOrder::Type order,
                                                        vector< RealVector >& nodalSet)
{
  std::vector<CFreal> coords;
  coords.resize(order+1);
   
  if(m_solPntDistribution.isNotNull()){
    coords = m_solPntDistribution->getLocalCoords1D(order);
  } else {
    // Use a default solution point distribution: Gauss Legendre.  
    switch(order)
      {
	case CFPolyOrder::ORDER0:
	{
	  coords[0] = 0.0;
	} break;
	case CFPolyOrder::ORDER1:
	{
	  coords[0] = -1./sqrt(3.);
	  coords[1] = +1./sqrt(3.);
	} break;
	case CFPolyOrder::ORDER2:
	{
	  coords[0] = -sqrt(3./5.);
	  coords[1] = 0.0;
	  coords[2] = +sqrt(3./5.);
	} break;
	case CFPolyOrder::ORDER3:
	{
	  coords[0] = -sqrt((3.+2.*sqrt(6./5.))/7.);
	  coords[1] = -sqrt((3.-2.*sqrt(6./5.))/7.);
	  coords[2] = +sqrt((3.-2.*sqrt(6./5.))/7.);
	  coords[3] = +sqrt((3.+2.*sqrt(6./5.))/7.);
	} break;
	case CFPolyOrder::ORDER4:
	{
	  coords[0] = -sqrt(5.+2.*sqrt(10./7.))/3.;
	  coords[1] = -sqrt(5.-2.*sqrt(10./7.))/3.;
	  coords[2] = 0.0;
	  coords[3] = +sqrt(5.-2.*sqrt(10./7.))/3.;
	  coords[4] = +sqrt(5.+2.*sqrt(10./7.))/3.;
	} break;
	case CFPolyOrder::ORDER5:
	{
	  coords[0] = -0.9324695142031521;
	  coords[1] = -0.6612093864662645;
	  coords[2] = -0.2386191860831969;
	  coords[3] = 0.2386191860831969;
	  coords[4] = 0.6612093864662645;
	  coords[5] = 0.9324695142031521;
	} break;
	case CFPolyOrder::ORDER6:
      {
        coords[0] = -0.949107912342759;
	coords[1] = -0.741531185599394;
	coords[2] = -0.405845151377397;
	coords[3] = 0.0;
	coords[4] = 0.405845151377397;
	coords[5] = 0.741531185599394;
	coords[6] = 0.949107912342759;
      } break;
      case CFPolyOrder::ORDER7:
      {
        coords[0] = -0.960289856497536;
	coords[1] = -0.796666477413627;
	coords[2] = -0.525532409916329;
	coords[3] = -0.183434642495650;
	coords[4] = 0.183434642495650;
	coords[5] = 0.525532409916329;
	coords[6] = 0.796666477413627;
	coords[7] = 0.960289856497536;
      } break;
      case CFPolyOrder::ORDER8:
      {
        coords[0] = -0.968160239507626;
	coords[1] = -0.836031107326636;
	coords[2] = -0.613371432700590;
	coords[3] = -0.324253423403809;
	coords[4] = 0.0;
	coords[5] = 0.324253423403809;
	coords[6] = 0.613371432700590;
	coords[7] = 0.836031107326636;
	coords[8] = 0.968160239507626;
      } break;
      case CFPolyOrder::ORDER9:
      {
        coords[0] = -0.973906528517172;
	coords[1] = -0.865063366688985;
	coords[2] = -0.679409568299024;
	coords[3] = -0.433395394129247;
	coords[4] = -0.148874338981631;
	coords[5] = 0.148874338981631;
	coords[6] = 0.433395394129247;
	coords[7] = 0.679409568299024;
	coords[8] = 0.865063366688985;
	coords[9] = 0.973906528517172;
      } break;
      case CFPolyOrder::ORDER10:
      {
        coords[0] = -0.978228658146057;
	coords[1] = -0.887062599768095;
	coords[2] = -0.730152005574049;
	coords[3] = -0.519096129110681;
	coords[4] = -0.269543155952345;
	coords[5] = 0.0;
	coords[6] = 0.269543155952345;
	coords[7] = 0.519096129110681;
	coords[8] = 0.730152005574049;
	coords[9] = 0.887062599768095;
	coords[10] = 0.978228658146057;
      } break;
	default:
	{
	  throw Common::NotImplementedException (FromHere(),"Gauss Legendre not implemented for order "
					+ StringOps::to_str(order) + ".");
	}
      }
  }
  
  const CFuint nbrPnts1D = order+1;

  // set solution point local coordinates
  nodalSet.resize(0);
  for (CFuint iKsi = 0; iKsi < nbrPnts1D; ++iKsi)
  {
    for (CFuint iEta = 0; iEta < nbrPnts1D; ++iEta)
    {
      for (CFuint iZta = 0; iZta < nbrPnts1D; ++iZta)
      {
	RealVector node(3);
	node[KSI] = coords[iKsi];
	node[ETA] = coords[iEta];
	node[ZTA] = coords[iZta];
	nodalSet.push_back(node);
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////

void HexaFluxReconstructionElementData::setCFLConvDiffRatio()
{
  CFAUTOTRACE;

  switch(m_polyOrder)
  {
    case CFPolyOrder::ORDER0:
    {
      m_cflConvDiffRatio = 2.0;//4.0; //2.0; // check this!  2/2
    } break;
    case CFPolyOrder::ORDER1:
    {
      m_cflConvDiffRatio = 4.0;//4.0;//6.5; //100.0; //200.0; // check this! 4/6
    } break;
    case CFPolyOrder::ORDER2:
    {
      m_cflConvDiffRatio = 6.0;//6.0;//17.0; //6.0; // check this! 6/10
    } break;
    case CFPolyOrder::ORDER3:
    {
      m_cflConvDiffRatio = 8.0;//8.0;//25.0; //500.0; //8.0; // check this! 8/14
    } break;
    case CFPolyOrder::ORDER4:
    {
      m_cflConvDiffRatio = 10.0;//10.0;//50.0; //10.0; // check this! 10/18
    } break;
    case CFPolyOrder::ORDER5:
    {
      m_cflConvDiffRatio = 12.0;//22.0;//12.0;//100.0; //1200.0; // check this! 12/22
    } break;
    case CFPolyOrder::ORDER6:
    {
      m_cflConvDiffRatio = 14.0;//14.0;//200.0; //14.0; // check this! 14/26
    } break;
    case CFPolyOrder::ORDER7:
    {
      m_cflConvDiffRatio = 16.0;//16.0;//2500.0; //16.0; // check this! 16/30
    } break;
    case CFPolyOrder::ORDER8:
    {
      m_cflConvDiffRatio = 18.0;//18.0;//800.0; //18.0; // check this! 18/34
    } break;
    case CFPolyOrder::ORDER9:
    {
      m_cflConvDiffRatio = 20.0;//20.0;//5500.0; //20.0; // check this! 20/38
    } break;
    case CFPolyOrder::ORDER10:
    {
      m_cflConvDiffRatio = 22.0;//22.0;//3200.0; //22.0; // check this! 22/42
    } break;
    default:
    {
      throw Common::NotImplementedException (FromHere(),"Higher-order quadrilateral FR cell not defined!");
    }
  }
}

//////////////////////////////////////////////////////////////////////

void HexaFluxReconstructionElementData::createFaceOutputPntCellMappedCoords()
{
  // number of points on a face
  const CFuint nbrFacePnts1D = m_polyOrder == 0 ? 2 :m_polyOrder + 1;

  // face mapped coordinates of uniform distribution of points
  vector<RealVector> faceMapCoords;
  const CFreal dKsiEta = 0 ? 2.0 : 2.0/m_polyOrder;
  CFreal ksi = -1.0;
  m_faceOutputPntFaceMappedCoords.resize(0);
  for (CFuint iKsi = 0; iKsi < nbrFacePnts1D; ++iKsi, ksi += dKsiEta)
  {
    CFreal eta = -1.0;
    for (CFuint iEta = 0; iEta < nbrFacePnts1D; ++iEta, eta += dKsiEta)
    {
      RealVector mapCoord(2);
      mapCoord[KSI] = ksi;
      mapCoord[ETA] = eta;
      m_faceOutputPntFaceMappedCoords.push_back(mapCoord);
    }
  }
  const CFuint nbrFacePnts = m_faceOutputPntFaceMappedCoords.size();
  cf_assert(nbrFacePnts == nbrFacePnts1D*nbrFacePnts1D);

  // compute cell mapped coordinates for distribution on each face
  const CFuint nbrCellFaces = getNbrCellFaces();
  m_faceOutputPntCellMappedCoords.resize(nbrCellFaces);
  for (CFuint iFace = 0; iFace < nbrCellFaces; ++iFace)
  {
    // current face node coordinates
    const vector<RealVector>& faceNodeCoords = m_faceNodeCoords[iFace];
    m_faceOutputPntCellMappedCoords[iFace].resize(0);
    for (CFuint iPnt = 0; iPnt < nbrFacePnts; ++iPnt)
    {
      const CFreal fun0 = 0.25*(1.0-m_faceOutputPntFaceMappedCoords[iPnt][KSI])
                              *(1.0-m_faceOutputPntFaceMappedCoords[iPnt][ETA]);
      const CFreal fun1 = 0.25*(1.0+m_faceOutputPntFaceMappedCoords[iPnt][KSI])
                              *(1.0-m_faceOutputPntFaceMappedCoords[iPnt][ETA]);
      const CFreal fun2 = 0.25*(1.0+m_faceOutputPntFaceMappedCoords[iPnt][KSI])
                              *(1.0+m_faceOutputPntFaceMappedCoords[iPnt][ETA]);
      const CFreal fun3 = 0.25*(1.0-m_faceOutputPntFaceMappedCoords[iPnt][KSI])
                              *(1.0+m_faceOutputPntFaceMappedCoords[iPnt][ETA]);
      m_faceOutputPntCellMappedCoords[iFace].push_back(fun0*faceNodeCoords[0]+
                                                       fun1*faceNodeCoords[1]+
                                                       fun2*faceNodeCoords[2]+
                                                       fun3*faceNodeCoords[3]);
    }
  }
}

//////////////////////////////////////////////////////////////////////

void HexaFluxReconstructionElementData::createFaceOutputPntConn()
{
  CFAUTOTRACE;

  // number of nodes 1D
  const CFuint nbrNodes1D = m_polyOrder + 1;

  m_faceOutputPntConn.resize(0);
  for (CFuint iKsi = 0; iKsi < static_cast<CFuint>(m_polyOrder); ++iKsi)
  {
    for (CFuint iEta = 0; iEta < static_cast<CFuint>(m_polyOrder); ++iEta)
    {
      vector< CFuint > cellNodesConn(4);
      cellNodesConn[0] = (iKsi  )*nbrNodes1D + iEta  ;
      cellNodesConn[1] = (iKsi+1)*nbrNodes1D + iEta  ;
      cellNodesConn[2] = (iKsi+1)*nbrNodes1D + iEta+1;
      cellNodesConn[3] = (iKsi  )*nbrNodes1D + iEta+1;
      m_faceOutputPntConn.push_back(cellNodesConn);
    }
  }
}

//////////////////////////////////////////////////////////////////////

void HexaFluxReconstructionElementData::createVandermondeMatrix()
{
  CFAUTOTRACE;
  
  const CFuint nbrSolPnts = m_solPntsLocalCoords.size();
  // number of solution points in 1D
  const CFuint nbrSolPnts1D = m_solPntsLocalCoord1D.size();
  
  CFuint order;
  CFuint idx;
  
  m_vandermonde.resize(nbrSolPnts,nbrSolPnts);
  m_vandermondeInv.resize(nbrSolPnts,nbrSolPnts);
  
  if(m_polyOrder != CFPolyOrder::ORDER0)
  {
    for (CFuint iSol = 0; iSol < nbrSolPnts; ++iSol)
    {
      CFuint modalDof = 0;
      RealVector orderIdx(nbrSolPnts1D);
      orderIdx = 0.0;
      
      for (CFuint iKsi = 0; iKsi < nbrSolPnts1D; ++iKsi)
      {
        for (CFuint iEta = 0; iEta < nbrSolPnts1D; ++iEta)
        {
          for (CFuint iZta = 0; iZta < nbrSolPnts1D; ++iZta, ++modalDof)
          {
            order = max(iKsi,iEta);
            order = max(order,iZta);
            idx = order*order*order+orderIdx[order];
            orderIdx[order] += 1;
              
            m_vandermonde(iSol,idx) = evaluateLegendre(m_solPntsLocalCoords[iSol][KSI],iKsi)*evaluateLegendre(m_solPntsLocalCoords[iSol][ETA],iEta)*evaluateLegendre(m_solPntsLocalCoords[iSol][ZTA],iZta);
          }
        }
      }
//       for (CFuint iOrder = 0; iOrder < nbrSolPnts1D; ++iOrder)
//       {
//         for (CFuint iOrderKsi = 0; iOrderKsi < iOrder; ++iOrderKsi, ++modalDof)
//         {
// 	  CFuint iOrderEta = iOrder;
// 	  CFuint iOrderZta = iOrder;
//           m_vandermonde(iSol,modalDof) = evaluateLegendre(m_solPntsLocalCoords[iSol][KSI],iOrderKsi)*evaluateLegendre(m_solPntsLocalCoords[iSol][ETA],iOrderEta)*evaluateLegendre(m_solPntsLocalCoords[iSol][ZTA],iOrderZta);
//         }
//         for (CFuint iOrderEta = 0; iOrderEta < iOrder+1; ++iOrderEta, ++modalDof)
//         {
// 	  CFuint iOrderKsi = iOrder;
// 	  CFuint iOrderZta = iOrder;
//           m_vandermonde(iSol,modalDof) = evaluateLegendre(m_solPntsLocalCoords[iSol][KSI],iOrderKsi)*evaluateLegendre(m_solPntsLocalCoords[iSol][ETA],iOrderEta)*evaluateLegendre(m_solPntsLocalCoords[iSol][ZTA],iOrderZta);
//         }
//         for (CFuint iOrderZta = 0; iOrderZta < iOrder+2; ++iOrderZta, ++modalDof)
//         {
// 	  CFuint iOrderKsi = iOrder;
// 	  CFuint iOrderEta = iOrder;
//           m_vandermonde(iSol,modalDof) = evaluateLegendre(m_solPntsLocalCoords[iSol][KSI],iOrderKsi)*evaluateLegendre(m_solPntsLocalCoords[iSol][ETA],iOrderEta)*evaluateLegendre(m_solPntsLocalCoords[iSol][ZTA],iOrderZta);
//         }
//       }
      cf_assert(modalDof == nbrSolPnts);
    }
    
    InvertMatrix(m_vandermonde,m_vandermondeInv);
  }
}

//////////////////////////////////////////////////////////////////////

void HexaFluxReconstructionElementData::createFlxSolDependencies()
{
  CFAUTOTRACE;
  
  const CFuint nbrSolPnts = m_solPntsLocalCoords.size();
  const CFuint nbrFlxPnts = m_flxPntsLocalCoords.size();
  // number of solution points in 1D
  const CFuint nbrSolPnts1D = m_solPntsLocalCoord1D.size();

  m_solFlxDep.resize(nbrSolPnts);
  m_solSolDep.resize(nbrSolPnts);
  m_flxSolDep.resize(nbrFlxPnts);

  CFuint iSol = 0;
  
  for (CFuint iKsi = 0; iKsi < nbrSolPnts1D; ++iKsi)
  {
    for (CFuint iEta = 0; iEta < nbrSolPnts1D; ++iEta)
    {
      for (CFuint iZta = 0; iZta < nbrSolPnts1D; ++iZta, ++iSol)
      {
        m_solFlxDep[iSol].push_back(iEta + iKsi*nbrSolPnts1D);
        m_solFlxDep[iSol].push_back(iEta + iKsi*nbrSolPnts1D + nbrSolPnts1D*nbrSolPnts1D);
        m_solFlxDep[iSol].push_back(iZta + iEta*nbrSolPnts1D + 2*nbrSolPnts1D*nbrSolPnts1D);
        m_solFlxDep[iSol].push_back(iZta + iEta*nbrSolPnts1D + 3*nbrSolPnts1D*nbrSolPnts1D);
        m_solFlxDep[iSol].push_back(iKsi + iZta*nbrSolPnts1D + 4*nbrSolPnts1D*nbrSolPnts1D);
        m_solFlxDep[iSol].push_back(iKsi + iZta*nbrSolPnts1D + 5*nbrSolPnts1D*nbrSolPnts1D);        

        m_flxSolDep[iEta + iKsi*nbrSolPnts1D].push_back(iSol);
        m_flxSolDep[iEta + iKsi*nbrSolPnts1D + nbrSolPnts1D*nbrSolPnts1D].push_back(iSol);
        m_flxSolDep[iZta + iEta*nbrSolPnts1D + 2*nbrSolPnts1D*nbrSolPnts1D].push_back(iSol);
        m_flxSolDep[iZta + iEta*nbrSolPnts1D + 3*nbrSolPnts1D*nbrSolPnts1D].push_back(iSol);
        m_flxSolDep[iKsi + iZta*nbrSolPnts1D + 4*nbrSolPnts1D*nbrSolPnts1D].push_back(iSol);
        m_flxSolDep[iKsi + iZta*nbrSolPnts1D + 5*nbrSolPnts1D*nbrSolPnts1D].push_back(iSol);

        for (CFuint jSol = 0; jSol < nbrSolPnts1D; ++jSol)
        {
          m_solSolDep[iSol].push_back(iKsi*nbrSolPnts1D*nbrSolPnts1D + iEta*nbrSolPnts1D + jSol);
          if (iSol != iZta + iKsi*nbrSolPnts1D*nbrSolPnts1D + jSol*nbrSolPnts1D) m_solSolDep[iSol].push_back(iZta + iKsi*nbrSolPnts1D*nbrSolPnts1D + jSol*nbrSolPnts1D);
          if (iSol != iZta + jSol*nbrSolPnts1D*nbrSolPnts1D + iEta*nbrSolPnts1D) m_solSolDep[iSol].push_back(iZta + jSol*nbrSolPnts1D*nbrSolPnts1D + iEta*nbrSolPnts1D);
        }
      }
    }
  }
}
//////////////////////////////////////////////////////////////////////

void HexaFluxReconstructionElementData::createFaceFlxPntsLocalCoordsPerType()
{
  //only needed for meshes with prisms or hybrid grids
  CFAUTOTRACE;
  // number of solution points in 1D
  const CFuint nbrFlxPnts1D = m_flxPntsLocalCoord1D.size();
  // number of solution points in triag
  const CFuint nbrFlxPntsTriag = (m_polyOrder+1)*(m_polyOrder+2)/2;

  // set face flux point face local coordinates on Triag Face (not applicable)
  m_faceFlxPntsLocalCoordsPerType.resize(2);
  m_faceFlxPntsLocalCoordsPerType[0].resize(0);
  

  // set face flux point face local coordinates on Quad Face (reference -1 1 quad)
  m_faceFlxPntsLocalCoordsPerType[1].resize(0);
  for (CFuint iKsi = 0; iKsi < nbrFlxPnts1D; ++iKsi)
  {
    for (CFuint iEta = 0; iEta < nbrFlxPnts1D; ++iEta)
    {
      RealVector flxCoord(2);
      flxCoord[KSI] = m_flxPntsLocalCoord1D[iKsi];
      flxCoord[ETA] = m_flxPntsLocalCoord1D[iEta];
      (m_faceFlxPntsLocalCoordsPerType[1]).push_back(flxCoord);
    }
  }

}
//////////////////////////////////////////////////////////////////////

void HexaFluxReconstructionElementData::createFaceIntegrationCoefsPerType()
{
  CFAUTOTRACE;

  // number of solution points in 1D
  const CFuint nbrFlxPnts1D = m_flxPntsLocalCoord1D.size();

  // set face integration coeffs coordinates on Triag Face
  m_faceIntegrationCoefsPerType.resize(2);

  m_faceIntegrationCoefsPerType[0].resize(0);

  // set face integration coeffs on Quad Face (reference -1 1 quad)

  // number of flux points on a face
  const CFuint nbrFlxPnts = nbrFlxPnts1D*nbrFlxPnts1D;

  // resize m_faceIntegrationCoefsPerType
  m_faceIntegrationCoefsPerType[1].resize(nbrFlxPnts);

  // create TensorProductGaussIntegrator
  TensorProductGaussIntegrator tpIntegrator(DIM_2D,m_polyOrder);

  // create face node local coordinates
  vector< RealVector > nodeCoord(4);
  nodeCoord[0].resize(2);
  nodeCoord[0][KSI] = -1.0;
  nodeCoord[0][ETA] = -1.0;
  nodeCoord[1].resize(2);
  nodeCoord[1][KSI] = +1.0;
  nodeCoord[1][ETA] = -1.0;
  nodeCoord[2].resize(2);
  nodeCoord[2][KSI] = +1.0;
  nodeCoord[2][ETA] = +1.0;
  nodeCoord[3].resize(2);
  nodeCoord[3][KSI] = -1.0;
  nodeCoord[3][ETA] = +1.0;

  // get quadrature point coordinates and wheights
  vector< RealVector > quadPntCoords   = tpIntegrator.getQuadPntsCoords  (nodeCoord);
  vector< CFreal     > quadPntWheights = tpIntegrator.getQuadPntsWheights(nodeCoord);
  const CFuint nbrQPnts = quadPntCoords.size();
  cf_assert(quadPntWheights.size() == nbrQPnts);

  // compute the coefficients for integration over a face
  // loop over flux points
  CFuint iFlx = 0;
  for (CFuint iFlxKsi = 0; iFlxKsi < nbrFlxPnts1D; ++iFlxKsi)
  {
    const CFreal ksiFlx = m_flxPntsLocalCoord1D[iFlxKsi];
    for (CFuint iFlxEta = 0; iFlxEta < nbrFlxPnts1D; ++iFlxEta, ++iFlx)
    {
      const CFreal etaFlx = m_flxPntsLocalCoord1D[iFlxEta];

      m_faceIntegrationCoefsPerType[1][iFlx] = 0.0;
      for (CFuint iQPnt = 0; iQPnt < nbrQPnts; ++iQPnt)
      {
        // quadrature point local coordinate on the face
        const CFreal ksiQPnt = quadPntCoords[iQPnt][KSI];
        const CFreal etaQPnt = quadPntCoords[iQPnt][ETA];

        // evaluate polynomial value in quadrature point
        CFreal quadPntPolyVal = 1.;
        for (CFuint iFac = 0; iFac < nbrFlxPnts1D; ++iFac)
        {
          if (iFac != iFlxKsi)
          {
            const CFreal ksiFac = m_flxPntsLocalCoord1D[iFac];
            quadPntPolyVal *= (ksiQPnt-ksiFac)/(ksiFlx-ksiFac);
          }
        }
        for (CFuint iFac = 0; iFac < nbrFlxPnts1D; ++iFac)
        {
          if (iFac != iFlxEta)
          {
            const CFreal etaFac = m_flxPntsLocalCoord1D[iFac];
            quadPntPolyVal *= (etaQPnt-etaFac)/(etaFlx-etaFac);
          }
        }

        // add contribution of quadrature point to integration coefficient
        m_faceIntegrationCoefsPerType[1][iFlx] += quadPntWheights[iQPnt]*quadPntPolyVal;
      }
    }
  }

}

//////////////////////////////////////////////////////////////////////

  } // namespace FluxReconstructionMethod

} // namespace COOLFluiD