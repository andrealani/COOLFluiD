#include <fstream>


#include "Common/CFLog.hh"
#include "Environment/DirPaths.hh"
#include "Common/NotImplementedException.hh"
#include "MathTools/RealMatrix.hh"
#include "SpectralFD/QuadSpectralFDElementData.hh"
#include "SpectralFD/TensorProductGaussIntegrator.hh"

//////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Common;

namespace COOLFluiD {

  namespace SpectralFD {

//////////////////////////////////////////////////////////////////////

QuadSpectralFDElementData::QuadSpectralFDElementData() :
  SpectralFDElementData()
{
  m_shape = CFGeoShape::QUAD;
  m_dimensionality = DIM_2D;
}

//////////////////////////////////////////////////////////////////////

QuadSpectralFDElementData::QuadSpectralFDElementData(CFPolyOrder::Type polyOrder)
{
  m_shape = CFGeoShape::QUAD;
  m_dimensionality = DIM_2D;
  m_polyOrder = polyOrder;

  resetSpectralFDElementData();
}

//////////////////////////////////////////////////////////////////////

QuadSpectralFDElementData::~QuadSpectralFDElementData()
{
}

//////////////////////////////////////////////////////////////////////

void QuadSpectralFDElementData::createFlxPntsLocalCoords()
{
  CFAUTOTRACE;

  // number of solution points in 1D
  const CFuint nbrSolPnts1D = m_solPntsLocalCoord1D.size();

  // number of flux points in 1D
  const CFuint nbrFlxPnts1D = m_flxPntsLocalCoord1D.size();

  // set flux point local coordinates
  m_flxPntsLocalCoords.resize(0);
  // ksi-flux points
  for (CFuint iEta = 0; iEta < nbrSolPnts1D; ++iEta)
  {
    for (CFuint iKsi = 0; iKsi < nbrFlxPnts1D; ++iKsi)
    {
      RealVector flxCoords(2);
      flxCoords[KSI] = m_flxPntsLocalCoord1D[iKsi];
      flxCoords[ETA] = m_solPntsLocalCoord1D[iEta];
      m_flxPntsLocalCoords.push_back(flxCoords);
    }
  }
  // eta-flux points
  for (CFuint iKsi = 0; iKsi < nbrSolPnts1D; ++iKsi)
  {
    for (CFuint iEta = 0; iEta < nbrFlxPnts1D; ++iEta)
    {
      RealVector flxCoords(2);
      flxCoords[KSI] = m_solPntsLocalCoord1D[iKsi];
      flxCoords[ETA] = m_flxPntsLocalCoord1D[iEta];
      m_flxPntsLocalCoords.push_back(flxCoords);
    }
  }
}

//////////////////////////////////////////////////////////////////////

void QuadSpectralFDElementData::createSolPntsLocalCoords()
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
      RealVector solCoords(2);
      solCoords[KSI] = m_solPntsLocalCoord1D[iKsi];
      solCoords[ETA] = m_solPntsLocalCoord1D[iEta];
      m_solPntsLocalCoords.push_back(solCoords);
    }
  }
}

//////////////////////////////////////////////////////////////////////

void QuadSpectralFDElementData::createFaceFlxPntsFaceLocalCoords()
{
  CFAUTOTRACE;

  // number of solution points in 1D
  const CFuint nbrSolPnts1D = m_solPntsLocalCoord1D.size();

  // set face flux point face local coordinates
  m_faceFlxPntsFaceLocalCoords.resize(0);
  for (CFuint iFlx = 0; iFlx < nbrSolPnts1D; ++iFlx)
  {
    RealVector flxCoord(1);
    flxCoord[KSI] = m_solPntsLocalCoord1D[iFlx];
    m_faceFlxPntsFaceLocalCoords.push_back(flxCoord);
  }
}

//////////////////////////////////////////////////////////////////////

void QuadSpectralFDElementData::createFlxPolyExponents()
{
  CFAUTOTRACE;

  // number of solution points in 1D
  const CFuint nbrSolPnts1D = m_solPntsLocalCoord1D.size();

  // number of flux points in 1D
  const CFuint nbrFlxPnts1D = m_flxPntsLocalCoord1D.size();

  // set flux point local coordinates
  m_flxPolyExponents.resize(0);
  m_flxPolyExponents.resize(2);
  // ksi-flux points
  for (CFuint iEta = 0; iEta < nbrSolPnts1D; ++iEta)
  {
    for (CFuint iKsi = 0; iKsi < nbrFlxPnts1D; ++iKsi)
    {
      vector< CFint > flxPolyExps(2);
      flxPolyExps[KSI] = iKsi;
      flxPolyExps[ETA] = iEta;
      m_flxPolyExponents[KSI].push_back(flxPolyExps);
    }
  }
  // eta-flux points
  for (CFuint iKsi = 0; iKsi < nbrSolPnts1D; ++iKsi)
  {
    for (CFuint iEta = 0; iEta < nbrFlxPnts1D; ++iEta)
    {
      vector< CFint > flxPolyExps(2);
      flxPolyExps[KSI] = iKsi;
      flxPolyExps[ETA] = iEta;
      m_flxPolyExponents[ETA].push_back(flxPolyExps);
    }
  }
}

//////////////////////////////////////////////////////////////////////

void QuadSpectralFDElementData::createSolPolyExponents()
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
      vector< CFint > solPolyExps(2);
      solPolyExps[KSI] = iKsi;
      solPolyExps[ETA] = iEta;
      m_solPolyExponents.push_back(solPolyExps);
    }
  }
}

//////////////////////////////////////////////////////////////////////

void QuadSpectralFDElementData::createFlxPntMatrixIdxForReconstruction()
{
  // number of solution points in 1D
  const CFuint nbrSolPnts1D = m_solPntsLocalCoord1D.size();

  // number of flux points in 1D
  const CFuint nbrFlxPnts1D = m_flxPntsLocalCoord1D.size();

  // set indices
  m_flxPntMatrixIdxForReconstruction.resize(0);
  // ksi-flux points
  for (CFuint iEta = 0; iEta < nbrSolPnts1D; ++iEta)
  {
    for (CFuint iKsi = 0; iKsi < nbrFlxPnts1D; ++iKsi)
    {
      m_flxPntMatrixIdxForReconstruction.push_back(iKsi);
    }
  }
  // eta-flux points
  for (CFuint iKsi = 0; iKsi < nbrSolPnts1D; ++iKsi)
  {
    for (CFuint iEta = 0; iEta < nbrFlxPnts1D; ++iEta)
    {
      m_flxPntMatrixIdxForReconstruction.push_back(iEta);
    }
  }
}

//////////////////////////////////////////////////////////////////////

void QuadSpectralFDElementData::createSolPntIdxsForReconstruction()
{
  // number of solution points in 1D
  const CFuint nbrSolPnts1D = m_solPntsLocalCoord1D.size();

  // number of flux points in 1D
  const CFuint nbrFlxPnts1D = m_flxPntsLocalCoord1D.size();

  // total number of flux polynomials
  const CFuint totNbrFlxPnts = getNbrOfFlxPnts();

  // set indices
  m_solPntIdxsForReconstruction.resize(0);
  m_solPntIdxsForReconstruction.resize(totNbrFlxPnts);
  CFuint flxIdx = 0;
  // ksi-flux points
  for (CFuint iEta = 0; iEta < nbrSolPnts1D; ++iEta)
  {
    for (CFuint iKsi = 0; iKsi < nbrFlxPnts1D; ++iKsi, ++flxIdx)
    {
      for (CFuint iPnt = 0; iPnt < nbrSolPnts1D; ++iPnt)
      {
        m_solPntIdxsForReconstruction[flxIdx].push_back(nbrSolPnts1D*iPnt+iEta);
      }
    }
  }
  // eta-flux points
  for (CFuint iKsi = 0; iKsi < nbrSolPnts1D; ++iKsi)
  {
    for (CFuint iEta = 0; iEta < nbrFlxPnts1D; ++iEta, ++flxIdx)
    {
      for (CFuint iPnt = 0; iPnt < nbrSolPnts1D; ++iPnt)
      {
        m_solPntIdxsForReconstruction[flxIdx].push_back(nbrSolPnts1D*iKsi+iPnt);
      }
    }
  }
  cf_assert(totNbrFlxPnts == flxIdx);

  // set indices for optimized reconstruction
  m_solPntIdxsForRecOptim.resize(0);
  m_solPntIdxsForRecOptim.resize(totNbrFlxPnts);
  flxIdx = 0;
  // ksi-flux points
  for (CFuint iEta = 0; iEta < nbrSolPnts1D; ++iEta)
  {
    for (CFuint iKsi = 0; iKsi < nbrFlxPnts1D; ++iKsi, ++flxIdx)
    {
      for (CFuint iPnt = 0; iPnt < m_solPntIdxsForRecFlxPnts1DOptim[iKsi].size(); ++iPnt)
      {
        const CFuint solIdx = m_solPntIdxsForRecFlxPnts1DOptim[iKsi][iPnt];
        m_solPntIdxsForRecOptim[flxIdx].push_back(nbrSolPnts1D*solIdx+iEta);
      }
    }
  }
  // eta-flux points
  for (CFuint iKsi = 0; iKsi < nbrSolPnts1D; ++iKsi)
  {
    for (CFuint iEta = 0; iEta < nbrFlxPnts1D; ++iEta, ++flxIdx)
    {
      for (CFuint iPnt = 0; iPnt < m_solPntIdxsForRecFlxPnts1DOptim[iEta].size(); ++iPnt)
      {
        const CFuint solIdx = m_solPntIdxsForRecFlxPnts1DOptim[iEta][iPnt];
        m_solPntIdxsForRecOptim[flxIdx].push_back(nbrSolPnts1D*iKsi+solIdx);
      }
    }
  }
  cf_assert(totNbrFlxPnts == flxIdx);
}

//////////////////////////////////////////////////////////////////////

void QuadSpectralFDElementData::createSolPntMatrixIdxForDerivation()
{
  CFAUTOTRACE;

  // number of solution points in 1D
  const CFuint nbrSolPnts1D = m_solPntsLocalCoord1D.size();

  // set indices
  m_solPntMatrixIdxForDerivation.resize(0);
  for (CFuint iKsi = 0; iKsi < nbrSolPnts1D; ++iKsi)
  {
    for (CFuint iEta = 0; iEta < nbrSolPnts1D; ++iEta)
    {
      vector< CFuint > solIdxs(2);
      solIdxs[KSI] = iKsi;
      solIdxs[ETA] = iEta;
      m_solPntMatrixIdxForDerivation.push_back(solIdxs);
    }
  }
}

//////////////////////////////////////////////////////////////////////

void QuadSpectralFDElementData::createFlxPntMatrixIdxForDerivation()
{
  CFAUTOTRACE;

  // number of solution points in 1D
  const CFuint nbrSolPnts1D = m_solPntsLocalCoord1D.size();

  // number of flux points in 1D
  const CFuint nbrFlxPnts1D = m_flxPntsLocalCoord1D.size();

  // set indices
  m_flxPntMatrixIdxForDerivation.resize(0);
  // ksi-flux points
  for (CFuint iEta = 0; iEta < nbrSolPnts1D; ++iEta)
  {
    for (CFuint iKsi = 0; iKsi < nbrFlxPnts1D; ++iKsi)
    {
      m_flxPntMatrixIdxForDerivation.push_back(iKsi);
    }
  }
  // eta-flux points
  for (CFuint iKsi = 0; iKsi < nbrSolPnts1D; ++iKsi)
  {
    for (CFuint iEta = 0; iEta < nbrFlxPnts1D; ++iEta)
    {
      m_flxPntMatrixIdxForDerivation.push_back(iEta);
    }
  }
/*  for (CFuint iFlx = 0; iFlx < m_flxPntMatrixIdxForDerivation.size(); ++iFlx)
  {
    CF_DEBUG_OBJ(m_flxPntMatrixIdxForDerivation[iFlx]);
  }*/
}

//////////////////////////////////////////////////////////////////////

void QuadSpectralFDElementData::createSolPntIdxsForDerivation()
{
  CFAUTOTRACE;

  // number of solution points in 1D
  const CFuint nbrSolPnts1D = m_solPntsLocalCoord1D.size();

  // number of flux points in 1D
  const CFuint nbrFlxPnts1D = m_flxPntsLocalCoord1D.size();

  // set indices
  m_solPntIdxsForDerivation.resize(0);
  // ksi-flux points
  for (CFuint iEta = 0; iEta < nbrSolPnts1D; ++iEta)
  {
    vector< CFuint > solPntIdxs;
    for (CFuint iSol = 0; iSol < nbrSolPnts1D; ++iSol)
    {
      solPntIdxs.push_back(nbrSolPnts1D*iSol + iEta);
    }
    for (CFuint iKsi = 0; iKsi < nbrFlxPnts1D; ++iKsi)
    {
      m_solPntIdxsForDerivation.push_back(solPntIdxs);
    }
  }
  // eta-flux points
  for (CFuint iKsi = 0; iKsi < nbrSolPnts1D; ++iKsi)
  {
    vector< CFuint > solPntIdxs;
    for (CFuint iSol = 0; iSol < nbrSolPnts1D; ++iSol)
    {
      solPntIdxs.push_back(nbrSolPnts1D*iKsi + iSol);
    }
    for (CFuint iEta = 0; iEta < nbrFlxPnts1D; ++iEta)
    {
      m_solPntIdxsForDerivation.push_back(solPntIdxs);
    }
  }
/*  for (CFuint iFlx = 0; iFlx < m_solPntIdxsForDerivation.size(); ++iFlx)
  {
    CF_DEBUG_OBJ(iFlx);
    for (CFuint iSol = 0; iSol < m_solPntIdxsForDerivation[iFlx].size(); ++iSol)
    {
      CF_DEBUG_OBJ(m_solPntIdxsForDerivation[iFlx][iSol]);
    }
  }*/
}

//////////////////////////////////////////////////////////////////////

void QuadSpectralFDElementData::createFlxPntIdxsForDerivation()
{
  CFAUTOTRACE;

  // number of solution points in 1D
  const CFuint nbrSolPnts1D = m_solPntsLocalCoord1D.size();

  // number of flux points in 1D
  const CFuint nbrFlxPnts1D = m_flxPntsLocalCoord1D.size();

  // number of solution points
  const CFuint nbrSolPnts = getNbrOfSolPnts();

  // number of flux points in one direction
  const CFuint nbrFlxPnts = nbrSolPnts1D*nbrFlxPnts1D;

  // set indices
  m_flxPntIdxsForDerivation.resize(0);
  m_flxPntIdxsForDerivation.resize(nbrSolPnts);
  CFuint iSol = 0;
  for (CFuint iKsi = 0; iKsi < nbrSolPnts1D; ++iKsi)
  {
    for (CFuint iEta = 0; iEta < nbrSolPnts1D; ++iEta, ++iSol)
    {
      for (CFuint iFlx = 0; iFlx < nbrFlxPnts1D; ++iFlx)
      {
        vector< CFuint > flxIdxs(2);
        flxIdxs[KSI] = iFlx + nbrFlxPnts1D*iEta;
        flxIdxs[ETA] = nbrFlxPnts + iFlx + nbrFlxPnts1D*iKsi;
        m_flxPntIdxsForDerivation[iSol].push_back(flxIdxs);
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////

void QuadSpectralFDElementData::createFlxPntDerivDir()
{
  CFAUTOTRACE;

  // number of solution points in 1D
  const CFuint nbrSolPnts1D = m_solPntsLocalCoord1D.size();

  // number of flux points in 1D
  const CFuint nbrFlxPnts1D = m_flxPntsLocalCoord1D.size();

  // set derivation directions
  m_flxPntDerivDir.resize(0);
  // ksi-flux points
  for (CFuint iEta = 0; iEta < nbrSolPnts1D; ++iEta)
  {
    for (CFuint iKsi = 0; iKsi < nbrFlxPnts1D; ++iKsi)
    {
      m_flxPntDerivDir.push_back(KSI);
    }
  }
  // eta-flux points
  for (CFuint iKsi = 0; iKsi < nbrSolPnts1D; ++iKsi)
  {
    for (CFuint iEta = 0; iEta < nbrFlxPnts1D; ++iEta)
    {
      m_flxPntDerivDir.push_back(ETA);
    }
  }
}

//////////////////////////////////////////////////////////////////////

void QuadSpectralFDElementData::createIntFlxPntIdxs()
{
  CFAUTOTRACE;

  // number of solution points in 1D
  const CFuint nbrSolPnts1D = m_solPntsLocalCoord1D.size();

  // number of flux points in 1D
  const CFuint nbrFlxPnts1D = m_flxPntsLocalCoord1D.size();

  // set internal flux points
  m_intFlxPntIdxs.resize(0);
  CFuint flxIdx = 0;
  // ksi-flux points
  for (CFuint iEta = 0; iEta < nbrSolPnts1D; ++iEta)
  {
    for (CFuint iKsi = 0; iKsi < nbrFlxPnts1D; ++iKsi, ++flxIdx)
    {
      if (iKsi != 0 && iKsi != nbrSolPnts1D)
      {
        m_intFlxPntIdxs.push_back(flxIdx);
      }
    }
  }
  // eta-flux points
  for (CFuint iKsi = 0; iKsi < nbrSolPnts1D; ++iKsi)
  {
    for (CFuint iEta = 0; iEta < nbrFlxPnts1D; ++iEta, ++flxIdx)
    {
      if (iEta != 0 && iEta != nbrSolPnts1D)
      {
        m_intFlxPntIdxs.push_back(flxIdx);
      }
    }
  }
  cf_assert(m_intFlxPntIdxs.size() == getNbrOfIntFlxPnts());
}

//////////////////////////////////////////////////////////////////////

void QuadSpectralFDElementData::createFaceFluxPntsConn()
{
  CFAUTOTRACE;

  // number of solution points in 1D
  const CFuint nbrSolPnts1D = m_solPntsLocalCoord1D.size();

  // number of flux points in 1D
  const CFuint nbrFlxPnts1D = m_flxPntsLocalCoord1D.size();

  // number of flux points per direction
  const CFuint nbrFlxPnts = nbrSolPnts1D*nbrFlxPnts1D;

  // resize m_faceFlxPntConn
  m_faceFlxPntConn.resize(4);

  // variable holding the face index
  CFuint faceIdx = 0;

  // first face
  for (CFuint iSol = 0; iSol < nbrSolPnts1D; ++iSol)
  {
    m_faceFlxPntConn[faceIdx].push_back(nbrFlxPnts+iSol*nbrFlxPnts1D);
  }
  ++faceIdx;


  // second face
  for (CFuint iSol = 0; iSol < nbrSolPnts1D; ++iSol)
  {
    m_faceFlxPntConn[faceIdx].push_back((iSol+1)*nbrFlxPnts1D-1);
  }
  ++faceIdx;

  // third face
  for (CFuint iSol = 0; iSol < nbrSolPnts1D; ++iSol)
  {
    const CFuint solIdxP1 = nbrSolPnts1D-iSol;
    m_faceFlxPntConn[faceIdx].push_back(nbrFlxPnts+solIdxP1*nbrFlxPnts1D-1);
  }
  ++faceIdx;

  // fourth face
  for (CFuint iSol = 0; iSol < nbrSolPnts1D; ++iSol)
  {
    const CFuint solIdx = nbrSolPnts1D-1-iSol;
    m_faceFlxPntConn[faceIdx].push_back(solIdx*nbrFlxPnts1D);
  }

/*  for (CFuint iFace = 0; iFace < m_faceFlxPntConn.size(); ++iFace)
  {
    for (CFuint iFlx = 0; iFlx < m_faceFlxPntConn[iFace].size(); ++iFlx)
    {
      CF_DEBUG_OBJ(m_faceFlxPntConn[iFace][iFlx]);
    }
  }*/
}

//////////////////////////////////////////////////////////////////////

void QuadSpectralFDElementData::createFaceFluxPntsConnPerOrient()
{
  CFAUTOTRACE;

  // number of orientations
  const CFuint nbrOrients = 10;

  // number of solution points in 1D
  const CFuint nbrSolPnts1D = m_solPntsLocalCoord1D.size();

  // create data structure
  m_faceFlxPntConnPerOrient.resize(nbrOrients);
  CFuint iOrient = 0;
  for (CFuint iFaceL = 0; iFaceL < 4; ++iFaceL)
  {
    for (CFuint iFaceR = iFaceL; iFaceR < 4; ++iFaceR, ++iOrient)
    {
      m_faceFlxPntConnPerOrient[iOrient].resize(2);
      for (CFuint iSol = 0; iSol < nbrSolPnts1D; ++iSol)
      {
        m_faceFlxPntConnPerOrient[iOrient][LEFT ]
            .push_back(m_faceFlxPntConn[iFaceL][iSol               ]);
        m_faceFlxPntConnPerOrient[iOrient][RIGHT]
            .push_back(m_faceFlxPntConn[iFaceR][nbrSolPnts1D-1-iSol]);
      }
    }
  }
/*  for (CFuint iOrient = 0; iOrient < m_faceFlxPntConnPerOrient.size(); ++iOrient)
  {
    CF_DEBUG_OBJ(iOrient);
    for (CFuint iFlx = 0; iFlx < m_faceFlxPntConnPerOrient[iOrient][LEFT ].size(); ++iFlx)
    {
      CF_DEBUG_OBJ(m_faceFlxPntConnPerOrient[iOrient][LEFT ][iFlx]);
    }
    for (CFuint iFlx = 0; iFlx < m_faceFlxPntConnPerOrient[iOrient][RIGHT].size(); ++iFlx)
    {
      CF_DEBUG_OBJ(m_faceFlxPntConnPerOrient[iOrient][RIGHT][iFlx]);
    }
  }*/
}

//////////////////////////////////////////////////////////////////////

void QuadSpectralFDElementData::createCellNodeCoords()
{
  CFAUTOTRACE;

  m_cellNodeCoords.resize(4);

  // first node
  m_cellNodeCoords[0].resize(2);
  m_cellNodeCoords[0][KSI] = -1.0;
  m_cellNodeCoords[0][ETA] = -1.0;

  // second node
  m_cellNodeCoords[1].resize(2);
  m_cellNodeCoords[1][KSI] = +1.0;
  m_cellNodeCoords[1][ETA] = -1.0;

  // third node
  m_cellNodeCoords[2].resize(2);
  m_cellNodeCoords[2][KSI] = +1.0;
  m_cellNodeCoords[2][ETA] = +1.0;

  // fourth node
  m_cellNodeCoords[3].resize(2);
  m_cellNodeCoords[3][KSI] = -1.0;
  m_cellNodeCoords[3][ETA] = +1.0;
}

//////////////////////////////////////////////////////////////////////

void QuadSpectralFDElementData::createFaceNodeConnectivity()
{
  CFAUTOTRACE;

  m_faceNodeConn.resize(4);

  m_faceNodeConn[0].resize(2);
  m_faceNodeConn[0][0] = 0;
  m_faceNodeConn[0][1] = 1;

  m_faceNodeConn[1].resize(2);
  m_faceNodeConn[1][0] = 1;
  m_faceNodeConn[1][1] = 2;

  m_faceNodeConn[2].resize(2);
  m_faceNodeConn[2][0] = 2;
  m_faceNodeConn[2][1] = 3;

  m_faceNodeConn[3].resize(2);
  m_faceNodeConn[3][0] = 3;
  m_faceNodeConn[3][1] = 0;
}

//////////////////////////////////////////////////////////////////////

void QuadSpectralFDElementData::createFaceMappedCoordDir()
{
  CFAUTOTRACE;

  m_faceMappedCoordDir.resize(4);

  m_faceMappedCoordDir[0] = -1;
  m_faceMappedCoordDir[1] = +1;
  m_faceMappedCoordDir[2] = +1;
  m_faceMappedCoordDir[3] = -1;
}

//////////////////////////////////////////////////////////////////////

void QuadSpectralFDElementData::createFaceNodeConnectivityPerOrient()
{
  CFAUTOTRACE;

  // number of faces
  const CFuint nbrFaces = m_faceNodeConn.size();

  // number of possible orientations
  const CFuint nbrOrient = 10;

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
      m_faceNodeConnPerOrient[iOrient][iSide].resize(2);
    }
  }

  // fill the variable
  CFuint iOrient = 0;
  for (CFuint iFaceL = 0; iFaceL < nbrFaces; ++iFaceL)
  {
    for (CFuint iFaceR = iFaceL; iFaceR < nbrFaces; ++iFaceR, ++iOrient)
    {
      m_faceConnPerOrient[iOrient][LEFT ] = iFaceL;
      m_faceConnPerOrient[iOrient][RIGHT] = iFaceR;

      m_faceMappedCoordDirPerOrient[iOrient][LEFT ] = m_faceMappedCoordDir[iFaceL];
      m_faceMappedCoordDirPerOrient[iOrient][RIGHT] = -m_faceMappedCoordDir[iFaceR];

      for (CFuint iNode = 0; iNode < 2; ++iNode)
      {
        m_faceNodeConnPerOrient[iOrient][LEFT ][iNode] = m_faceNodeConn[iFaceL][iNode  ];
        m_faceNodeConnPerOrient[iOrient][RIGHT][iNode] = m_faceNodeConn[iFaceR][1-iNode];
      }
    }
  }
  cf_assert(iOrient == nbrOrient);
}

//////////////////////////////////////////////////////////////////////

void QuadSpectralFDElementData::createFaceIntegrationCoefs()
{
  CFAUTOTRACE;

  // number of flux points on a face
  const CFuint nbrFlxPnts = m_solPntsLocalCoord1D.size();

  // resize m_faceIntegrationCoefs
  m_faceIntegrationCoefs.resize(nbrFlxPnts);

  // create TensorProductGaussIntegrator
  TensorProductGaussIntegrator tpIntegrator(DIM_1D,m_polyOrder);

  // create face node local coordinates
  vector< RealVector > nodeCoord(2);
  nodeCoord[0].resize(1);
  nodeCoord[0][KSI] = -1.0;
  nodeCoord[1].resize(1);
  nodeCoord[1][KSI] = +1.0;

  // get quadrature point coordinates and wheights
  vector< RealVector > quadPntCoords   = tpIntegrator.getQuadPntsCoords  (nodeCoord);
  vector< CFreal     > quadPntWheights = tpIntegrator.getQuadPntsWheights(nodeCoord);
  const CFuint nbrQPnts = quadPntCoords.size();
  cf_assert(quadPntWheights.size() == nbrQPnts);

  // compute the coefficients for integration over a face
  // loop over flux points
  for (CFuint iFlx = 0; iFlx < nbrFlxPnts; ++iFlx)
  {
    m_faceIntegrationCoefs[iFlx] = 0.0;

    const CFreal ksiFlx = m_solPntsLocalCoord1D[iFlx];
    for (CFuint iQPnt = 0; iQPnt < nbrQPnts; ++iQPnt)
    {
      // quadrature point local coordinate on the face
      const CFreal ksiQPnt = quadPntCoords[iQPnt][KSI];

      // evaluate polynomial value in quadrature point
      CFreal quadPntPolyVal = 1.;
      for (CFuint iFac = 0; iFac < nbrFlxPnts; ++iFac)
      {
        if (iFac != iFlx)
        {
          const CFreal ksiFac = m_solPntsLocalCoord1D[iFac];
          quadPntPolyVal *= (ksiQPnt-ksiFac)/(ksiFlx-ksiFac);
        }
      }

      // add contribution of quadrature point to integration coefficient
      m_faceIntegrationCoefs[iFlx] += quadPntWheights[iQPnt]*quadPntPolyVal;
    }
  }
}

//////////////////////////////////////////////////////////////////////

void QuadSpectralFDElementData::createCellAvgSolCoefs()
{
  CFAUTOTRACE;

  // number of solution points
  const CFuint nbrSolPnts = getNbrOfSolPnts();

  // resize m_cellAvgSolCoefs
  m_cellAvgSolCoefs.resize(nbrSolPnts);

  // create TensorProductGaussIntegrator
  TensorProductGaussIntegrator tpIntegrator(DIM_2D,m_polyOrder);

  // create cell node local coordinates
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
    m_cellAvgSolCoefs[iSol] *= 0.25;
  }
}

//////////////////////////////////////////////////////////////////////

void QuadSpectralFDElementData::createCellCenterDerivCoefs()
{
  CFAUTOTRACE;

  // center coordinate
  vector< RealVector > centerCoord(1,RealVector(2));
  centerCoord[0][KSI] = 0.0;
  centerCoord[0][ETA] = 0.0;

  vector< vector< vector< CFreal > > > polyDerivs =
      getSolPolyDerivsAtNode(centerCoord);

  // number of solution points
  const CFuint nbrSolPnts = getNbrOfSolPnts();

  // set polynomial derivatives
  m_cellCenterDerivCoefs.resize(2);
  m_cellCenterDerivCoefs[KSI].resize(nbrSolPnts);
  m_cellCenterDerivCoefs[ETA].resize(nbrSolPnts);
  for (CFuint iPoly = 0; iPoly < nbrSolPnts; ++iPoly)
  {
    m_cellCenterDerivCoefs[KSI][iPoly] = polyDerivs[0][KSI][iPoly];
    m_cellCenterDerivCoefs[ETA][iPoly] = polyDerivs[0][ETA][iPoly];
  }
}

//////////////////////////////////////////////////////////////////////

void QuadSpectralFDElementData::setInterpolationNodeSet(const CFPolyOrder::Type order,
                                                        vector< RealVector >& nodalSet)
{
  // number of points in one direction
  const CFuint nbrPnts1D = order+1;

  // set node coordinates
  nodalSet.resize(0);
  if (order == CFPolyOrder::ORDER0)
  {
    RealVector node(2);
    node[KSI] = 0.0;
    node[ETA] = 0.0;
    nodalSet.push_back(node);
  }
  else
  {
    for (CFuint iKsi = 0; iKsi < nbrPnts1D; ++iKsi)
    {
      for (CFuint iEta = 0; iEta < nbrPnts1D; ++iEta)
      {
        RealVector node(2);
        node[KSI] = -cos(iKsi*MathTools::MathConsts::CFrealPi()/order);
        node[ETA] = -cos(iEta*MathTools::MathConsts::CFrealPi()/order);
        nodalSet.push_back(node);
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////

void QuadSpectralFDElementData::setCFLConvDiffRatio()
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
    case CFPolyOrder::ORDER4:
    {
      m_cflConvDiffRatio = 50.0; // check this!
    } break;
    case CFPolyOrder::ORDER5:
    {
      m_cflConvDiffRatio = 50.0; // check this!
    } break;
    default:
    {
      throw Common::NotImplementedException (FromHere(),"Higher-order quadrilateral SD cell not defined!");
    }
  }
}

//////////////////////////////////////////////////////////////////////

void QuadSpectralFDElementData::createFaceOutputPntCellMappedCoords()
{
  // number of points on a face
  const CFuint nbrFacePnts = m_polyOrder == 0 ? 2 :m_polyOrder + 1;

  // face mapped coordinates of uniform distribution of points
  const CFreal dKsi = m_polyOrder == 0 ? 2.0 : 2.0/m_polyOrder;
  CFreal ksi = -1.0;
  m_faceOutputPntFaceMappedCoords.resize(0);
  for (CFuint iPnt = 0; iPnt < nbrFacePnts; ++iPnt, ksi += dKsi)
  {
    RealVector mapCoord(1);
    mapCoord[KSI] = ksi;
    m_faceOutputPntFaceMappedCoords.push_back(mapCoord);
  }

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
      const CFreal fun0 = 0.5*(1.0-m_faceOutputPntFaceMappedCoords[iPnt][KSI]);
      const CFreal fun1 = 0.5*(1.0+m_faceOutputPntFaceMappedCoords[iPnt][KSI]);
      m_faceOutputPntCellMappedCoords[iFace].push_back(fun0*faceNodeCoords[0]+fun1*faceNodeCoords[1]);
    }
  }
}

//////////////////////////////////////////////////////////////////////

void QuadSpectralFDElementData::createFaceOutputPntConn()
{
  CFAUTOTRACE;

  m_faceOutputPntConn.resize(0);
  for (CFuint iCell = 0; iCell < static_cast<CFuint>(m_polyOrder); ++iCell)
  {
    vector<CFuint> cellNode(2);
    cellNode[0] = iCell;
    cellNode[1] = iCell+1;
    m_faceOutputPntConn.push_back(cellNode);
  }
}

//////////////////////////////////////////////////////////////////////

  } // namespace SpectralFD

} // namespace COOLFluiD
