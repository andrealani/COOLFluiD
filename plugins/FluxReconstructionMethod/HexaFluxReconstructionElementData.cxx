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
//  CFLog(VERBOSE, m_flxPntsLocalCoords[0] << "\n");
//  CFLog(VERBOSE, m_flxPntsLocalCoords[1] << "\n");
//  CFLog(VERBOSE, m_flxPntsLocalCoords[2] << "\n");
//  CFLog(VERBOSE, m_flxPntsLocalCoords[3] << "\n");
//  CFLog(VERBOSE, m_flxPntsLocalCoords[4] << "\n");
//  CFLog(VERBOSE, m_flxPntsLocalCoords[5] << "\n");
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
/*  for (CFuint i = 0; i < m_solPntsLocalCoords.size(); ++i)
  {
    CF_DEBUG_OBJ(m_solPntsLocalCoords[i]);
  }*/
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
/*  for (CFuint i = 0; i < m_faceFlxPntsFaceLocalCoords.size(); ++i)
  {
    CF_DEBUG_OBJ(m_faceFlxPntsFaceLocalCoords[i]);
  }*/
}

//////////////////////////////////////////////////////////////////////

void HexaFluxReconstructionElementData::createFlxPolyExponents()
{
  CFAUTOTRACE;

//   // number of solution points in 1D
//   const CFuint nbrSolPnts1D = m_solPntsLocalCoord1D.size();
// 
//   // number of flux points in 1D
//   const CFuint nbrFlxPnts1D = m_flxPntsLocalCoord1D.size();
// 
//   // set flux point local coordinates
//   m_flxPolyExponents.resize(0);
//   m_flxPolyExponents.resize(3);
//   // ksi-flux points
//   for (CFuint iEta = 0; iEta < nbrSolPnts1D; ++iEta)
//   {
//     for (CFuint iZta = 0; iZta < nbrSolPnts1D; ++iZta)
//     {
//       for (CFuint iKsi = 0; iKsi < nbrFlxPnts1D; ++iKsi)
//       {
//         vector< CFint > flxPolyExps(3);
//         flxPolyExps[KSI] = iKsi;
//         flxPolyExps[ETA] = iEta;
//         flxPolyExps[ZTA] = iZta;
//         m_flxPolyExponents[KSI].push_back(flxPolyExps);
//       }
//     }
//   }
//   // eta-flux points
//   for (CFuint iZta = 0; iZta < nbrSolPnts1D; ++iZta)
//   {
//     for (CFuint iKsi = 0; iKsi < nbrSolPnts1D; ++iKsi)
//     {
//       for (CFuint iEta = 0; iEta < nbrFlxPnts1D; ++iEta)
//       {
//         vector< CFint > flxPolyExps(3);
//         flxPolyExps[KSI] = iKsi;
//         flxPolyExps[ETA] = iEta;
//         flxPolyExps[ZTA] = iZta;
//         m_flxPolyExponents[ETA].push_back(flxPolyExps);
//       }
//     }
//   }
//   // zta-flux points
//   for (CFuint iKsi = 0; iKsi < nbrSolPnts1D; ++iKsi)
//   {
//     for (CFuint iEta = 0; iEta < nbrSolPnts1D; ++iEta)
//     {
//       for (CFuint iZta = 0; iZta < nbrFlxPnts1D; ++iZta)
//       {
//         vector< CFint > flxPolyExps(3);
//         flxPolyExps[KSI] = iKsi;
//         flxPolyExps[ETA] = iEta;
//         flxPolyExps[ZTA] = iZta;
//         m_flxPolyExponents[ZTA].push_back(flxPolyExps);
//       }
//     }
//   }
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
/*  for (CFuint i = 0; i < m_solPolyExponents.size(); ++i)
  {
    CF_DEBUG_OBJ(m_solPolyExponents[i][KSI]);
    CF_DEBUG_OBJ(m_solPolyExponents[i][ETA]);
    CF_DEBUG_OBJ(m_solPolyExponents[i][ZTA]);
  }*/
}

//////////////////////////////////////////////////////////////////////

void HexaFluxReconstructionElementData::createFlxPntMatrixIdxForReconstruction()
{
//   // number of solution points in 1D
//   const CFuint nbrSolPnts1D = m_solPntsLocalCoord1D.size();
// 
//   // number of flux points in 1D
//   const CFuint nbrFlxPnts1D = m_flxPntsLocalCoord1D.size();
// 
//   // set indices
//   m_flxPntMatrixIdxForReconstruction.resize(0);
//   // ksi-flux points
//   for (CFuint iEta = 0; iEta < nbrSolPnts1D; ++iEta)
//   {
//     for (CFuint iZta = 0; iZta < nbrSolPnts1D; ++iZta)
//     {
//       for (CFuint iKsi = 0; iKsi < nbrFlxPnts1D; ++iKsi)
//       {
//         m_flxPntMatrixIdxForReconstruction.push_back(iKsi);
//       }
//     }
//   }
//   // eta-flux points
//   for (CFuint iZta = 0; iZta < nbrSolPnts1D; ++iZta)
//   {
//     for (CFuint iKsi = 0; iKsi < nbrSolPnts1D; ++iKsi)
//     {
//       for (CFuint iEta = 0; iEta < nbrFlxPnts1D; ++iEta)
//       {
//         m_flxPntMatrixIdxForReconstruction.push_back(iEta);
//       }
//     }
//   }
//   // zta-flux points
//   for (CFuint iKsi = 0; iKsi < nbrSolPnts1D; ++iKsi)
//   {
//     for (CFuint iEta = 0; iEta < nbrSolPnts1D; ++iEta)
//     {
//       for (CFuint iZta = 0; iZta < nbrFlxPnts1D; ++iZta)
//       {
//         m_flxPntMatrixIdxForReconstruction.push_back(iZta);
//       }
//     }
//   }
// /*  for (CFuint i = 0; i < m_flxPntMatrixIdxForReconstruction.size(); ++i)
//   {
//     CF_DEBUG_OBJ(m_flxPntMatrixIdxForReconstruction[i]);
//   }*/
}

//////////////////////////////////////////////////////////////////////

void HexaFluxReconstructionElementData::createSolPntIdxsForReconstruction()
{
//   // number of solution points in 1D
//   const CFuint nbrSolPnts1D = m_solPntsLocalCoord1D.size();
//   const CFuint nbrSolPnts1DSq = nbrSolPnts1D*nbrSolPnts1D;
// 
//   // number of flux points in 1D
//   const CFuint nbrFlxPnts1D = m_flxPntsLocalCoord1D.size();
// 
//   // total number of flux polynomials
//   const CFuint totNbrFlxPnts = getNbrOfFlxPnts();
// 
//   // set indices
//   m_solPntIdxsForReconstruction.resize(0);
//   m_solPntIdxsForReconstruction.resize(totNbrFlxPnts);
//   CFuint flxIdx = 0;
//   // ksi-flux points
//   for (CFuint iEta = 0; iEta < nbrSolPnts1D; ++iEta)
//   {
//     for (CFuint iZta = 0; iZta < nbrSolPnts1D; ++iZta)
//     {
//       for (CFuint iKsi = 0; iKsi < nbrFlxPnts1D; ++iKsi, ++flxIdx)
//       {
//         for (CFuint iPnt = 0; iPnt < nbrSolPnts1D; ++iPnt)
//         {
//           m_solPntIdxsForReconstruction[flxIdx]
//               .push_back(nbrSolPnts1DSq*iPnt+nbrSolPnts1D*iEta+iZta);
//         }
//       }
//     }
//   }
//   // eta-flux points
//   for (CFuint iZta = 0; iZta < nbrSolPnts1D; ++iZta)
//   {
//     for (CFuint iKsi = 0; iKsi < nbrSolPnts1D; ++iKsi)
//     {
//       for (CFuint iEta = 0; iEta < nbrFlxPnts1D; ++iEta, ++flxIdx)
//       {
//         for (CFuint iPnt = 0; iPnt < nbrSolPnts1D; ++iPnt)
//         {
//           m_solPntIdxsForReconstruction[flxIdx]
//               .push_back(nbrSolPnts1DSq*iKsi+nbrSolPnts1D*iPnt+iZta);
//         }
//       }
//     }
//   }
//   // zta-flux points
//   for (CFuint iKsi = 0; iKsi < nbrSolPnts1D; ++iKsi)
//   {
//     for (CFuint iEta = 0; iEta < nbrSolPnts1D; ++iEta)
//     {
//       for (CFuint iZta = 0; iZta < nbrFlxPnts1D; ++iZta, ++flxIdx)
//       {
//         for (CFuint iPnt = 0; iPnt < nbrSolPnts1D; ++iPnt)
//         {
//           m_solPntIdxsForReconstruction[flxIdx]
//               .push_back(nbrSolPnts1DSq*iKsi+nbrSolPnts1D*iEta+iPnt);
//         }
//       }
//     }
//   }
// //   for (CFuint i = 0; i < m_solPntIdxsForReconstruction.size(); ++i)
// //   {
// //     CF_DEBUG_OBJ(i);
// //     for (CFuint j = 0; j < m_solPntIdxsForReconstruction[i].size(); ++j)
// //     {
// //       CF_DEBUG_OBJ(m_solPntIdxsForReconstruction[i][j]);
// //     }
// //   }
// 
//   cf_assert(totNbrFlxPnts == flxIdx);
// 
//   // set indices for optimized reconstruction
//   m_solPntIdxsForRecOptim.resize(0);
//   m_solPntIdxsForRecOptim.resize(totNbrFlxPnts);
//   flxIdx = 0;
//   for (CFuint iEta = 0; iEta < nbrSolPnts1D; ++iEta)
//   {
//     for (CFuint iZta = 0; iZta < nbrSolPnts1D; ++iZta)
//     {
//       for (CFuint iKsi = 0; iKsi < nbrFlxPnts1D; ++iKsi, ++flxIdx)
//       {
//         for (CFuint iPnt = 0; iPnt < m_solPntIdxsForRecFlxPnts1DOptim[iKsi].size(); ++iPnt)
//         {
//           const CFuint solIdx = m_solPntIdxsForRecFlxPnts1DOptim[iKsi][iPnt];
//           m_solPntIdxsForRecOptim[flxIdx]
//               .push_back(nbrSolPnts1DSq*solIdx+nbrSolPnts1D*iEta+iZta);
//         }
//       }
//     }
//   }
//   // eta-flux points
//   for (CFuint iZta = 0; iZta < nbrSolPnts1D; ++iZta)
//   {
//     for (CFuint iKsi = 0; iKsi < nbrSolPnts1D; ++iKsi)
//     {
//       for (CFuint iEta = 0; iEta < nbrFlxPnts1D; ++iEta, ++flxIdx)
//       {
//         for (CFuint iPnt = 0; iPnt < m_solPntIdxsForRecFlxPnts1DOptim[iEta].size(); ++iPnt)
//         {
//           const CFuint solIdx = m_solPntIdxsForRecFlxPnts1DOptim[iEta][iPnt];
//           m_solPntIdxsForRecOptim[flxIdx]
//               .push_back(nbrSolPnts1DSq*iKsi+nbrSolPnts1D*solIdx+iZta);
//         }
//       }
//     }
//   }
//   // zta-flux points
//   for (CFuint iKsi = 0; iKsi < nbrSolPnts1D; ++iKsi)
//   {
//     for (CFuint iEta = 0; iEta < nbrSolPnts1D; ++iEta)
//     {
//       for (CFuint iZta = 0; iZta < nbrFlxPnts1D; ++iZta, ++flxIdx)
//       {
//         for (CFuint iPnt = 0; iPnt < m_solPntIdxsForRecFlxPnts1DOptim[iZta].size(); ++iPnt)
//         {
//           const CFuint solIdx = m_solPntIdxsForRecFlxPnts1DOptim[iZta][iPnt];
//           m_solPntIdxsForRecOptim[flxIdx]
//               .push_back(nbrSolPnts1DSq*iKsi+nbrSolPnts1D*iEta+solIdx);
//         }
//       }
//     }
//   }
//   cf_assert(totNbrFlxPnts == flxIdx);
}

//////////////////////////////////////////////////////////////////////

void HexaFluxReconstructionElementData::createSolPntMatrixIdxForDerivation()
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
      for (CFuint iZta = 0; iZta < nbrSolPnts1D; ++iZta)
      {
        vector< CFuint > solIdxs(3);
        solIdxs[KSI] = iKsi;
        solIdxs[ETA] = iEta;
        solIdxs[ZTA] = iZta;
        m_solPntMatrixIdxForDerivation.push_back(solIdxs);
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////

void HexaFluxReconstructionElementData::createFlxPntMatrixIdxForDerivation()
{
  CFAUTOTRACE;

//   // number of solution points in 1D
//   const CFuint nbrSolPnts1D = m_solPntsLocalCoord1D.size();
// 
//   // number of flux points in 1D
//   const CFuint nbrFlxPnts1D = m_flxPntsLocalCoord1D.size();
// 
//   // set indices
//   m_flxPntMatrixIdxForDerivation.resize(0);
//   // ksi-flux points
//   for (CFuint iEta = 0; iEta < nbrSolPnts1D; ++iEta)
//   {
//     for (CFuint iZta = 0; iZta < nbrSolPnts1D; ++iZta)
//     {
//       for (CFuint iKsi = 0; iKsi < nbrFlxPnts1D; ++iKsi)
//       {
//         m_flxPntMatrixIdxForDerivation.push_back(iKsi);
//       }
//     }
//   }
//   // eta-flux points
//   for (CFuint iZta = 0; iZta < nbrSolPnts1D; ++iZta)
//   {
//     for (CFuint iKsi = 0; iKsi < nbrSolPnts1D; ++iKsi)
//     {
//       for (CFuint iEta = 0; iEta < nbrFlxPnts1D; ++iEta)
//       {
//         m_flxPntMatrixIdxForDerivation.push_back(iEta);
//       }
//     }
//   }
//   // zta-flux points
//   for (CFuint iKsi = 0; iKsi < nbrSolPnts1D; ++iKsi)
//   {
//     for (CFuint iEta = 0; iEta < nbrSolPnts1D; ++iEta)
//     {
//       for (CFuint iZta = 0; iZta < nbrFlxPnts1D; ++iZta)
//       {
//         m_flxPntMatrixIdxForDerivation.push_back(iZta);
//       }
//     }
//   }
// /*  for (CFuint iFlx = 0; iFlx < m_flxPntMatrixIdxForDerivation.size(); ++iFlx)
//   {
//     CF_DEBUG_OBJ(m_flxPntMatrixIdxForDerivation[iFlx]);
//   }*/
}

//////////////////////////////////////////////////////////////////////

void HexaFluxReconstructionElementData::createSolPntIdxsForDerivation()
{
  CFAUTOTRACE;

//   // number of solution points in 1D
//   const CFuint nbrSolPnts1D = m_solPntsLocalCoord1D.size();
//   const CFuint nbrSolPnts1DSq = nbrSolPnts1D*nbrSolPnts1D;
// 
//   // number of flux points in 1D
//   const CFuint nbrFlxPnts1D = m_flxPntsLocalCoord1D.size();
// 
//   // set indices
//   m_solPntIdxsForDerivation.resize(0);
//   // ksi-flux points
//   for (CFuint iEta = 0; iEta < nbrSolPnts1D; ++iEta)
//   {
//     for (CFuint iZta = 0; iZta < nbrSolPnts1D; ++iZta)
//     {
//       vector< CFuint > solPntIdxs;
//       for (CFuint iSol = 0; iSol < nbrSolPnts1D; ++iSol)
//       {
//         solPntIdxs.push_back(nbrSolPnts1DSq*iSol+nbrSolPnts1D*iEta+iZta);
//       }
//       for (CFuint iKsi = 0; iKsi < nbrFlxPnts1D; ++iKsi)
//       {
//         m_solPntIdxsForDerivation.push_back(solPntIdxs);
//       }
//     }
//   }
//   // eta-flux points
//   for (CFuint iZta = 0; iZta < nbrSolPnts1D; ++iZta)
//   {
//     for (CFuint iKsi = 0; iKsi < nbrSolPnts1D; ++iKsi)
//     {
//       vector< CFuint > solPntIdxs;
//       for (CFuint iSol = 0; iSol < nbrSolPnts1D; ++iSol)
//       {
//         solPntIdxs.push_back(nbrSolPnts1DSq*iKsi+nbrSolPnts1D*iSol+iZta);
//       }
//       for (CFuint iEta = 0; iEta < nbrFlxPnts1D; ++iEta)
//       {
//         m_solPntIdxsForDerivation.push_back(solPntIdxs);
//       }
//     }
//   }
//   // zta-flux points
//   for (CFuint iKsi = 0; iKsi < nbrSolPnts1D; ++iKsi)
//   {
//     for (CFuint iEta = 0; iEta < nbrSolPnts1D; ++iEta)
//     {
//       vector< CFuint > solPntIdxs;
//       for (CFuint iSol = 0; iSol < nbrSolPnts1D; ++iSol)
//       {
//         solPntIdxs.push_back(nbrSolPnts1DSq*iKsi+nbrSolPnts1D*iEta+iSol);
//       }
//       for (CFuint iZta = 0; iZta < nbrFlxPnts1D; ++iZta)
//       {
//         m_solPntIdxsForDerivation.push_back(solPntIdxs);
//       }
//     }
//   }
// /*  for (CFuint iFlx = 0; iFlx < m_solPntIdxsForDerivation.size(); ++iFlx)
//   {
//     CF_DEBUG_OBJ(iFlx);
//     for (CFuint iSol = 0; iSol < m_solPntIdxsForDerivation[iFlx].size(); ++iSol)
//     {
//       CF_DEBUG_OBJ(m_solPntIdxsForDerivation[iFlx][iSol]);
//     }
//   }*/
}

//////////////////////////////////////////////////////////////////////

void HexaFluxReconstructionElementData::createFlxPntIdxsForDerivation()
{
  CFAUTOTRACE;

//   // number of solution points in 1D
//   const CFuint nbrSolPnts1D = m_solPntsLocalCoord1D.size();
// 
//   // number of flux points in 1D
//   const CFuint nbrFlxPnts1D = m_flxPntsLocalCoord1D.size();
// 
//   // number of solution points
//   const CFuint nbrSolPnts = getNbrOfSolPnts();
// 
//   // number of flux points in one direction
//   const CFuint nbrFlxPnts = nbrSolPnts1D*nbrFlxPnts1D;
// 
//   // set indices
//   m_flxPntIdxsForDerivation.resize(0);
//   m_flxPntIdxsForDerivation.resize(nbrSolPnts);
//   CFuint iSol = 0;
//   for (CFuint iKsi = 0; iKsi < nbrSolPnts1D; ++iKsi)
//   {
//     for (CFuint iEta = 0; iEta < nbrSolPnts1D; ++iEta)
//     {
//       for (CFuint iZta = 0; iZta < nbrSolPnts1D; ++iZta, ++iSol)
//       {
//         for (CFuint iFlx = 0; iFlx < nbrFlxPnts1D; ++iFlx)
//         {
//           vector< CFuint > flxIdxs(3);
//           flxIdxs[KSI] = iFlx + nbrFlxPnts1D*nbrSolPnts1D*iEta + nbrFlxPnts1D*iZta;
//           flxIdxs[ETA] = nbrFlxPnts + iFlx + nbrFlxPnts1D*nbrSolPnts1D*iZta + nbrFlxPnts1D*iKsi;
//           flxIdxs[ZTA] = 2*nbrFlxPnts + iFlx + nbrFlxPnts1D*nbrSolPnts1D*iKsi + nbrFlxPnts1D*iEta;
//           m_flxPntIdxsForDerivation[iSol].push_back(flxIdxs);
//         }
//       }
//     }
//   }
}

//////////////////////////////////////////////////////////////////////

void HexaFluxReconstructionElementData::createFlxPntDerivDir()
{
  CFAUTOTRACE;

//   // number of solution points in 1D
//   const CFuint nbrSolPnts1D = m_solPntsLocalCoord1D.size();
// 
//   // number of flux points in 1D
//   const CFuint nbrFlxPnts1D = m_flxPntsLocalCoord1D.size();
// 
//   // set derivation directions
//   m_flxPntDerivDir.resize(0);
//   // ksi-flux points
//   for (CFuint iEta = 0; iEta < nbrSolPnts1D; ++iEta)
//   {
//     for (CFuint iZta = 0; iZta < nbrSolPnts1D; ++iZta)
//     {
//       for (CFuint iKsi = 0; iKsi < nbrFlxPnts1D; ++iKsi)
//       {
//         m_flxPntDerivDir.push_back(KSI);
//       }
//     }
//   }
//   // eta-flux points
//   for (CFuint iZta = 0; iZta < nbrSolPnts1D; ++iZta)
//   {
//     for (CFuint iKsi = 0; iKsi < nbrSolPnts1D; ++iKsi)
//     {
//       for (CFuint iEta = 0; iEta < nbrFlxPnts1D; ++iEta)
//       {
//         m_flxPntDerivDir.push_back(ETA);
//       }
//     }
//   }
//   // zta-flux points
//   for (CFuint iKsi = 0; iKsi < nbrSolPnts1D; ++iKsi)
//   {
//     for (CFuint iEta = 0; iEta < nbrSolPnts1D; ++iEta)
//     {
//       for (CFuint iZta = 0; iZta < nbrFlxPnts1D; ++iZta)
//       {
//         m_flxPntDerivDir.push_back(ZTA);
//       }
//     }
//   }
}

//////////////////////////////////////////////////////////////////////

void HexaFluxReconstructionElementData::createIntFlxPntIdxs()
{
  CFAUTOTRACE;

//   // number of solution points in 1D
//   const CFuint nbrSolPnts1D = m_solPntsLocalCoord1D.size();
// 
//   // number of flux points in 1D
//   const CFuint nbrFlxPnts1D = m_flxPntsLocalCoord1D.size();
// 
//   // set internal flux points
//   m_intFlxPntIdxs.resize(0);
//   CFuint flxIdx = 0;
//   // ksi-flux points
//   for (CFuint iEta = 0; iEta < nbrSolPnts1D; ++iEta)
//   {
//     for (CFuint iZta = 0; iZta < nbrSolPnts1D; ++iZta)
//     {
//       for (CFuint iKsi = 0; iKsi < nbrFlxPnts1D; ++iKsi, ++flxIdx)
//       {
//         if (iKsi != 0 && iKsi != nbrSolPnts1D)
//         {
//           m_intFlxPntIdxs.push_back(flxIdx);
//         }
//       }
//     }
//   }
//   // eta-flux points
//   for (CFuint iZta = 0; iZta < nbrSolPnts1D; ++iZta)
//   {
//     for (CFuint iKsi = 0; iKsi < nbrSolPnts1D; ++iKsi)
//     {
//       for (CFuint iEta = 0; iEta < nbrFlxPnts1D; ++iEta, ++flxIdx)
//       {
//         if (iEta != 0 && iEta != nbrSolPnts1D)
//         {
//           m_intFlxPntIdxs.push_back(flxIdx);
//         }
//       }
//     }
//   }
//   // zta-flux points
//   for (CFuint iKsi = 0; iKsi < nbrSolPnts1D; ++iKsi)
//   {
//     for (CFuint iEta = 0; iEta < nbrSolPnts1D; ++iEta)
//     {
//       for (CFuint iZta = 0; iZta < nbrFlxPnts1D; ++iZta, ++flxIdx)
//       {
//         if (iZta != 0 && iZta != nbrSolPnts1D)
//         {xConnL
//           m_intFlxPntIdxs.push_back(flxIdx);
//         }
//       }
//     }
//   }
//   cf_assert(m_intFlxPntIdxs.size() == getNbrOfIntFlxPnts());
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


//   // number of solution points in 1D
//   const CFuint nbrSolPnts1D = m_solPntsLocalCoord1D.size();
// 
//   // number of flux points in 1D
//   const CFuint nbrFlxPnts1D = m_flxPntsLocalCoord1D.size();
// 
//   // number of flux points per direction
//   const CFuint nbrFlxPnts = nbrSolPnts1D*nbrSolPnts1D*nbrFlxPnts1D;
// 
//   // resize m_faceFlxPntConn
//   m_faceFlxPntConn.resize(6);
// 
//   // variable holding the face index
//   CFuint faceIdx = 0;
// 
//   // first face
//   for (CFuint iEta = 0; iEta < nbrSolPnts1D; ++iEta)
//   {
//     for (CFuint iKsi = 0; iKsi < nbrSolPnts1D; ++iKsi)
//     {
//       m_faceFlxPntConn[faceIdx]
//           .push_back(2*nbrFlxPnts+nbrSolPnts1D*nbrFlxPnts1D*iKsi+nbrFlxPnts1D*iEta);
//     }
//   }
//   ++faceIdx;
// 
//   // second face
//   for (CFuint iKsi = 0; iKsi < nbrSolPnts1D; ++iKsi)
//   {
//     for (CFuint iEta = 0; iEta < nbrSolPnts1D; ++iEta)
//     {
//       m_faceFlxPntConn[faceIdx]
//           .push_back(2*nbrFlxPnts+nbrSolPnts1D*nbrFlxPnts1D*iKsi+nbrFlxPnts1D*iEta+nbrSolPnts1D);
//     }
//   }
//   ++faceIdx;
// 
//   // third face
//   for (CFuint iKsi = 0; iKsi < nbrSolPnts1D; ++iKsi)
//   {
//     for (CFuint iZta = 0; iZta < nbrSolPnts1D; ++iZta)
//     {
//       m_faceFlxPntConn[faceIdx]
//           .push_back(nbrFlxPnts+nbrSolPnts1D*nbrFlxPnts1D*iZta+nbrFlxPnts1D*iKsi);
//     }
//   }
//   ++faceIdx;
// 
//   // fourth face
//   for (CFuint iEta = 0; iEta < nbrSolPnts1D; ++iEta)
//   {
//     for (CFuint iZta = 0; iZta < nbrSolPnts1D; ++iZta)
//     {
//       m_faceFlxPntConn[faceIdx]
//           .push_back(nbrSolPnts1D*nbrFlxPnts1D*iEta+nbrFlxPnts1D*iZta+nbrSolPnts1D);
//     }
//   }
//   ++faceIdx;
// 
//   // fifth face
//   for (CFuint iKsi = 0; iKsi < nbrSolPnts1D; ++iKsi)
//   {
//     const CFuint idxKsi = nbrSolPnts1D-iKsi-1;
//     for (CFuint iZta = 0; iZta < nbrSolPnts1D; ++iZta)
//     {
//       m_faceFlxPntConn[faceIdx]
//           .push_back(nbrFlxPnts+nbrSolPnts1D*nbrFlxPnts1D*iZta+nbrFlxPnts1D*idxKsi+nbrSolPnts1D);
//     }
//   }
//   ++faceIdx;
// 
//   // sixth face
//   for (CFuint iZta = 0; iZta < nbrSolPnts1D; ++iZta)
//   {
//     for (CFuint iEta = 0; iEta < nbrSolPnts1D; ++iEta)
//     {
//       m_faceFlxPntConn[faceIdx]
//         .push_back(nbrSolPnts1D*nbrFlxPnts1D*iEta+nbrFlxPnts1D*iZta);
//     }
//   }
//   ++faceIdx;
// 
// 
// /*  for (CFuint iFace = 0; iFace < m_faceFlxPntConn.size(); ++iFace)
//   {
//     for (CFuint iFlx = 0; iFlx < m_faceFlxPntConn[iFace].size(); ++iFlx)
//     {
//       CF_DEBUG_OBJ(m_faceFlxPntConn[iFace][iFlx]);
//     }
//   }*/
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
//      CFLog(VERBOSE,"rotFlxIdxs_0:" << rotFlxIdxs[0][0] << rotFlxIdxs[0][1] << rotFlxIdxs[0][2] << rotFlxIdxs[0][3] << "\n");
//      CFLog(VERBOSE,"rotFlxIdxs_1:" << rotFlxIdxs[1][0] << rotFlxIdxs[1][1] << rotFlxIdxs[1][2] << rotFlxIdxs[1][3] << "\n");
//      CFLog(VERBOSE,"iRotGroup:" << iRotGroup << "\n");
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
//          CFLog(VERBOSE,"m_faceFlxPntConnPerOrient[" <<iOrient << "][LEFT]" << m_faceFlxPntConnPerOrient[iOrient][LEFT][0] << m_faceFlxPntConnPerOrient[iOrient][LEFT][1] << m_faceFlxPntConnPerOrient[iOrient][LEFT][2] << m_faceFlxPntConnPerOrient[iOrient][LEFT][3] << m_faceFlxPntConnPerOrient[iOrient][LEFT][4] << m_faceFlxPntConnPerOrient[iOrient][LEFT][5] << m_faceFlxPntConnPerOrient[iOrient][LEFT][6] << m_faceFlxPntConnPerOrient[iOrient][LEFT][7] << "\n");
//          CFLog(VERBOSE,"m_faceFlxPntConnPerOrient[" <<iOrient << "][RIGHT]" << m_faceFlxPntConnPerOrient[iOrient][RIGHT][0] << m_faceFlxPntConnPerOrient[iOrient][RIGHT][1] << m_faceFlxPntConnPerOrient[iOrient][RIGHT][2] << m_faceFlxPntConnPerOrient[iOrient][RIGHT][3] << m_faceFlxPntConnPerOrient[iOrient][RIGHT][4] << m_faceFlxPntConnPerOrient[iOrient][RIGHT][5] << m_faceFlxPntConnPerOrient[iOrient][RIGHT][6] << m_faceFlxPntConnPerOrient[iOrient][RIGHT][7] << "\n");
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
      //CFLog(VERBOSE,"m_faceFlxPntConnPerOrient[0][LEFT]" << m_faceFlxPntConnPerOrient[0][LEFT][0] << m_faceFlxPntConnPerOrient[0][LEFT][1] << m_faceFlxPntConnPerOrient[0][LEFT][2] << m_faceFlxPntConnPerOrient[0][LEFT][3] << m_faceFlxPntConnPerOrient[0][LEFT][4] << m_faceFlxPntConnPerOrient[0][LEFT][5] << m_faceFlxPntConnPerOrient[0][LEFT][6] << m_faceFlxPntConnPerOrient[0][LEFT][7] << m_faceFlxPntConnPerOrient[0][LEFT][8] << "\n");
      //CFLog(VERBOSE,"m_faceFlxPntConnPerOrient[0][RIGHT]" << m_faceFlxPntConnPerOrient[0][RIGHT][0] << m_faceFlxPntConnPerOrient[0][RIGHT][1] << m_faceFlxPntConnPerOrient[0][RIGHT][2] << m_faceFlxPntConnPerOrient[0][RIGHT][3] << m_faceFlxPntConnPerOrient[0][RIGHT][4] << m_faceFlxPntConnPerOrient[0][RIGHT][5] << m_faceFlxPntConnPerOrient[0][RIGHT][6] << m_faceFlxPntConnPerOrient[0][RIGHT][7] << m_faceFlxPntConnPerOrient[0][RIGHT][8] << "\n");
    //CFLog(VERBOSE,iOrient << "\n");
  }
// /*  for (CFuint iOrient = 0; iOrient < m_faceFlxPntConnPerOrient.size(); ++iOrient)
//   {
//     CF_DEBUG_OBJ(iOrient);
//     for (CFuint iFlx = 0; iFlx < m_faceFlxPntConnPerOrient[iOrient][LEFT ].size(); ++iFlx)
//     {
//       CF_DEBUG_OBJ(m_faceFlxPntConnPerOrient[iOrient][LEFT ][iFlx]);
//     }
//     for (CFuint iFlx = 0; iFlx < m_faceFlxPntConnPerOrient[iOrient][RIGHT].size(); ++iFlx)
//     {
//       CF_DEBUG_OBJ(m_faceFlxPntConnPerOrient[iOrient][RIGHT][iFlx]);
//     }
//   }*/
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
  
  
//   // number of points in one direction
//   const CFuint nbrPnts1D = order+1;
// 
//   // set node coordinates
//   nodalSet.resize(0);
//   if (order == CFPolyOrder::ORDER0)
//   {
//     RealVector node(3);
//     node[KSI] = 0.0;
//     node[ETA] = 0.0;
//     node[ZTA] = 0.0;
//     nodalSet.push_back(node);
//   }
//   else
//   {
//     for (CFuint iKsi = 0; iKsi < nbrPnts1D; ++iKsi)
//     {
//       for (CFuint iEta = 0; iEta < nbrPnts1D; ++iEta)
//       {
//         for (CFuint iZta = 0; iZta < nbrPnts1D; ++iZta)
//         {
//           RealVector node(3);
//           node[KSI] = -cos(iKsi*MathTools::MathConsts::CFrealPi()/order);
//           node[ETA] = -cos(iEta*MathTools::MathConsts::CFrealPi()/order);
//           node[ZTA] = -cos(iZta*MathTools::MathConsts::CFrealPi()/order);
//           nodalSet.push_back(node);
//         }
//       }
//     }
//   }
}

//////////////////////////////////////////////////////////////////////

void HexaFluxReconstructionElementData::setCFLConvDiffRatio()
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

  } // namespace FluxReconstructionMethod

} // namespace COOLFluiD
