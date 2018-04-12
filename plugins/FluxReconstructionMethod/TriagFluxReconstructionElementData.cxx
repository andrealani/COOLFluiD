#include <fstream>


#include "Common/CFLog.hh"
#include "Environment/DirPaths.hh"
#include "Common/NotImplementedException.hh"
#include "MathTools/RealMatrix.hh"
#include "FluxReconstructionMethod/TriagFluxReconstructionElementData.hh"
#include "MathTools/MathConsts.hh"

//////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Common;

namespace COOLFluiD {

  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////

TriagFluxReconstructionElementData::TriagFluxReconstructionElementData() :
  FluxReconstructionElementData()
{
  m_shape = CFGeoShape::TRIAG;
  m_dimensionality = DIM_2D;
}

//////////////////////////////////////////////////////////////////////

TriagFluxReconstructionElementData::TriagFluxReconstructionElementData(CFPolyOrder::Type polyOrder)
{
  m_shape = CFGeoShape::TRIAG;
  m_dimensionality = DIM_2D;
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

TriagFluxReconstructionElementData::TriagFluxReconstructionElementData(CFPolyOrder::Type polyOrder, 
								       Common::SafePtr< BasePointDistribution > solPntDist, 
								       Common::SafePtr< BasePointDistribution > flxPntDist)
{
  m_shape = CFGeoShape::TRIAG;
  m_dimensionality = DIM_2D;
  m_polyOrder = polyOrder;
  m_solPntsLocalCoord1D = solPntDist->getLocalCoords1D(polyOrder);
  m_flxPntsLocalCoord1D = flxPntDist->getLocalCoords1D(polyOrder);
  m_solPntDistribution = solPntDist;
  m_flxPntDistribution = flxPntDist;

  resetFluxReconstructionElementData();
}

//////////////////////////////////////////////////////////////////////

TriagFluxReconstructionElementData::~TriagFluxReconstructionElementData()
{
}

//////////////////////////////////////////////////////////////////////

void TriagFluxReconstructionElementData::createFaceNormals()
{
  CFAUTOTRACE;
  
  // number of faces
  const CFuint nbrFaces = 3;

  // compute the normals
  m_faceNormals.resize(nbrFaces);
  for (CFuint iFace = 0; iFace < nbrFaces; ++iFace)
  {
    // resize the variable
    m_faceNormals[iFace].resize(m_dimensionality);

    // compute normal
    m_faceNormals[iFace][KSI] =  (m_faceNodeCoords[iFace][1][ETA] - m_faceNodeCoords[iFace][0][ETA]);
    m_faceNormals[iFace][ETA] = -(m_faceNodeCoords[iFace][1][KSI] - m_faceNodeCoords[iFace][0][KSI]);

  }

}

//////////////////////////////////////////////////////////////////////

void TriagFluxReconstructionElementData::createFaceNodeConnectivityPerOrient()
{
  CFAUTOTRACE;

  // number of faces
  const CFuint nbrFaces = m_faceNodeConn.size();

  // number of possible orientations
  const CFuint nbrOrient = 6;

  // resize the variable
  m_faceNodeConnPerOrient.resize(nbrOrient);
  m_faceConnPerOrient.resize(nbrOrient);
  m_faceMappedCoordDirPerOrient.resize(nbrOrient);
  for (CFuint iOrient = 0; iOrient < nbrOrient; ++iOrient)
  {
    m_faceNodeConnPerOrient[iOrient].resize(2);
    m_faceConnPerOrient[iOrient].resize(2);
    m_faceMappedCoordDirPerOrient[iOrient].resize(2);
    for (CFuint iCell = 0; iCell < 2; ++iCell)
      m_faceNodeConnPerOrient[iOrient][iCell].resize(2);
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
}

//////////////////////////////////////////////////////////////////////
  
void TriagFluxReconstructionElementData::createSolPolyExponents()
{
  CFAUTOTRACE;

  // define exponents
  m_solPolyExponents.resize(0);
  for (CFuint iKsi = 0; iKsi < static_cast< CFuint >(m_polyOrder)+1; ++iKsi)
  {
    for (CFuint iEta = 0; iEta < m_polyOrder+1-iKsi; ++iEta)
    {
      vector< CFint > solPolyExps(2);
      solPolyExps[KSI] = iKsi;
      solPolyExps[ETA] = iEta;
      m_solPolyExponents.push_back(solPolyExps);
    }
  }

  cf_assert(m_solPolyExponents.size() == static_cast<CFuint>((m_polyOrder+1)*(m_polyOrder+2)/2));
}

//////////////////////////////////////////////////////////////////////
  
void TriagFluxReconstructionElementData::createNodePolyExponents()
{
  CFAUTOTRACE;

  // define exponents
  m_nodePolyExponents.resize(0);
  for (CFuint iKsi = 0; iKsi < 2; ++iKsi)
  {
    for (CFuint iEta = 0; iEta < 2-iKsi; ++iEta)
    {
      vector< CFint > nodePolyExps(2);
      nodePolyExps[KSI] = iKsi;
      nodePolyExps[ETA] = iEta;
      m_nodePolyExponents.push_back(nodePolyExps);
    }
  }
}

//////////////////////////////////////////////////////////////////////

void TriagFluxReconstructionElementData::createFaceFluxPntsConn()
{
 CFAUTOTRACE;

//   // number of flux points at a face
//   const CFuint nbrFaceFlxPnts = m_polyOrder +2;
// 
//   // resize m_faceFlxPntsConn (nbr of faces)
//   m_faceFlxPntsConn.resize(3);
// 
//   // create connectivity
//   // face 0
//   m_faceFlxPntsConn[0].resize(nbrFaceFlxPnts);
//   CFuint flxPntID = 0;
//   for (CFuint iFlx = 0; iFlx < nbrFaceFlxPnts; ++iFlx)
//   {
//     m_faceFlxPntsConn[0][iFlx] = flxPntID;
//     flxPntID += iFlx+1;
//   }
// 
//   // face 1
//   m_faceFlxPntsConn[1].resize(nbrFaceFlxPnts);
//   flxPntID = m_faceFlxPntsConn[0][nbrFaceFlxPnts-1];
//   for (CFuint iFlx = 0; iFlx < nbrFaceFlxPnts; ++iFlx, ++flxPntID)
//   {
//     m_faceFlxPntsConn[1][iFlx] = flxPntID;
//   }
// 
//   // face 2
//   m_faceFlxPntsConn[2].resize(nbrFaceFlxPnts);
//   flxPntID = m_faceFlxPntsConn[1][nbrFaceFlxPnts-1];
//   for (CFuint iFlx = 0; iFlx < nbrFaceFlxPnts; ++iFlx)
//   {
//     m_faceFlxPntsConn[2][iFlx] = flxPntID;
//     flxPntID -= nbrFaceFlxPnts-iFlx;
//   }
// 
// /*  for (CFuint iFace = 0; iFace < 3; ++iFace)
//   {
//     CF_DEBUG_OBJ(iFace);
//     for (CFuint iFlx = 0; iFlx < nbrFaceFlxPnts; ++iFlx)
//     {
//       CF_DEBUG_OBJ(m_faceFlxPntsConn[iFace][iFlx]);
//     }
//   }*/
}

//////////////////////////////////////////////////////////////////////

void TriagFluxReconstructionElementData::createFluxPolyExponents()
{
  CFAUTOTRACE;

//   // helper variable
//   const CFuint polyOrderP2 = m_polyOrder + 2;
// 
//   // number of polynomial terms
//   const CFuint nbrPolyTerms = (m_polyOrder + 2)*(m_polyOrder + 3)/2;
// 
//   // resize the variable
//   m_fluxPolyExponents.resize(nbrPolyTerms);
//   for (CFuint iTerm = 0; iTerm < nbrPolyTerms; ++iTerm)
//   {
//     m_fluxPolyExponents[iTerm].resize(2);
//   }
// 
//   // define exponents
//   CFuint iTerm = 0;
//   for (CFuint iP = 0; iP < polyOrderP2; ++iP)
//   {
//     for (CFuint iY = 0; iY < iP+1; ++iY, ++iTerm)
//     {
//       m_fluxPolyExponents[iTerm][KSI] = iP-iY;
//       m_fluxPolyExponents[iTerm][ETA] = iY;
//     }
//   }
}

//////////////////////////////////////////////////////////////////////

void TriagFluxReconstructionElementData::setInterpolationNodeSet(const CFPolyOrder::Type order,
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
    for (CFuint iEta = 0; iEta < nbrPnts1D-iKsi; ++iEta)
    {
      RealVector node(2);
      node[KSI] = (coords[iKsi]+1.)/2.; //Map [-1;1] onto [0;1]
      node[ETA] = (coords[iEta]+1.)/2.*(1.-coords[iKsi])/2.; //Map for eta [0;1] onto [0;1-ksi'] with ksi' on [-;1]
      nodalSet.push_back(node);
    }
  }
}

//////////////////////////////////////////////////////////////////////

void TriagFluxReconstructionElementData::setCFLConvDiffRatio()
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

void TriagFluxReconstructionElementData::createFaceOutputPntCellMappedCoords()
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
  const CFuint nbrCellFaces = 3;
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

void TriagFluxReconstructionElementData::createFaceOutputPntConn()
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

void TriagFluxReconstructionElementData::createFlxPntsLocalCoords()
{
  // number of flux points in 1D
  const CFuint nbrFlxPnts1D = m_flxPntsLocalCoord1D.size();

  // set flux point local coordinates
  m_flxPntsLocalCoords.resize(0);
  // loop over 1D Flx pnts
  for (CFuint i = 0; i < nbrFlxPnts1D; ++i)
  {
    RealVector flxCoords(2);
    flxCoords[KSI] = -1;
    flxCoords[ETA] = m_flxPntsLocalCoord1D[i];
    m_flxPntsLocalCoords.push_back(flxCoords);
    
    flxCoords[KSI] = m_flxPntsLocalCoord1D[i];
    flxCoords[ETA] = -1;
    m_flxPntsLocalCoords.push_back(flxCoords);
    
    flxCoords[KSI] = m_flxPntsLocalCoord1D[i];
    flxCoords[ETA] = -m_flxPntsLocalCoord1D[i];
    m_flxPntsLocalCoords.push_back(flxCoords);
  }
}

//////////////////////////////////////////////////////////////////////
  
void TriagFluxReconstructionElementData::createSolPntsLocalCoords()
{   
  // number of solution points in 1D
  const CFuint nbrSolPnts1D = m_solPntsLocalCoord1D.size();

  // set solution point local coordinates
  m_solPntsLocalCoords.resize(0);
  for (CFuint iKsi = 0; iKsi < nbrSolPnts1D; ++iKsi)
  { 
    for (CFuint iEta = 0; iEta < nbrSolPnts1D-iKsi; ++iEta) //as ksi grows, create pnts on vertical closer to node (1;-1)
    {
      RealVector solCoords(2);
      solCoords[KSI] = (m_solPntsLocalCoord1D[iKsi]+1.)/2.; //Map [-1,1] to [0,1]
      solCoords[ETA] = (m_solPntsLocalCoord1D[iEta]+1.)/2.*(1-m_solPntsLocalCoord1D[iKsi])/2.; //Map for eta [0;1] onto [0;1-ksi'] with ksi' on [-;1]
      m_solPntsLocalCoords.push_back(solCoords);
    }
  }
    
}

//////////////////////////////////////////////////////////////////////
  
void TriagFluxReconstructionElementData::createFaceFlxPntsFaceLocalCoords()
{
  CFAUTOTRACE;

  // number of solution points in 1D
  const CFuint nbrFlxPnts1D = m_flxPntsLocalCoord1D.size();

  // set face flux point face local coordinates
  m_faceFlxPntsFaceLocalCoords.resize(0);
  for (CFuint iFlx = 0; iFlx < nbrFlxPnts1D; ++iFlx)
  {
    RealVector flxCoord(1);
    flxCoord[KSI] = m_flxPntsLocalCoord1D[iFlx]/2.+1.;
    m_faceFlxPntsFaceLocalCoords.push_back(flxCoord);
  }
}


//////////////////////////////////////////////////////////////////////
  
void TriagFluxReconstructionElementData::createFaceFluxPntsConnPerOrient()
{
  CFAUTOTRACE;

//   // number of orientations
//   const CFuint nbrOrients = 10;
// 
//   // number of solution points in 1D
//   const CFuint nbrSolPnts1D = m_solPntsLocalCoord1D.size();
// 
//   // create data structure
//   m_faceFlxPntConnPerOrient.resize(nbrOrients);
//   CFuint iOrient = 0;
//   for (CFuint iFaceL = 0; iFaceL < 4; ++iFaceL)
//   {
//     for (CFuint iFaceR = iFaceL; iFaceR < 4; ++iFaceR, ++iOrient)
//     {
//       m_faceFlxPntConnPerOrient[iOrient].resize(2);
//       for (CFuint iSol = 0; iSol < nbrSolPnts1D; ++iSol)
//       {
//         m_faceFlxPntConnPerOrient[iOrient][LEFT ]
//             .push_back(m_faceFlxPntConn[iFaceL][iSol               ]);
//         m_faceFlxPntConnPerOrient[iOrient][RIGHT]
//             .push_back(m_faceFlxPntConn[iFaceR][nbrSolPnts1D-1-iSol]);
//       }
//     }
//   }
}

//////////////////////////////////////////////////////////////////////
  
void TriagFluxReconstructionElementData::createCellNodeCoords()
{
  CFAUTOTRACE;

  m_cellNodeCoords.resize(3);

  // first node
  m_cellNodeCoords[0].resize(2);
  m_cellNodeCoords[0][KSI] = 0.;
  m_cellNodeCoords[0][ETA] = 0.;

  // second node
  m_cellNodeCoords[1].resize(2);
  m_cellNodeCoords[1][KSI] = 0.;
  m_cellNodeCoords[1][ETA] = +1.0;

  // third node
  m_cellNodeCoords[2].resize(2);
  m_cellNodeCoords[2][KSI] = +1.0;
  m_cellNodeCoords[2][ETA] = 0.;
}

//////////////////////////////////////////////////////////////////////
  
void TriagFluxReconstructionElementData::createFaceNodeConnectivity()
{
  CFAUTOTRACE;

  m_faceNodeConn.resize(3);

  m_faceNodeConn[0].resize(2);
  m_faceNodeConn[0][0] = 0;
  m_faceNodeConn[0][1] = 1;

  m_faceNodeConn[1].resize(2);
  m_faceNodeConn[1][0] = 1;
  m_faceNodeConn[1][1] = 2;

  m_faceNodeConn[2].resize(2);
  m_faceNodeConn[2][0] = 2;
  m_faceNodeConn[2][1] = 0;
}

//////////////////////////////////////////////////////////////////////
  
void TriagFluxReconstructionElementData::createFaceMappedCoordDir()
{
  CFAUTOTRACE;

  m_faceMappedCoordDir.resize(3);

  m_faceMappedCoordDir[0] = -1;
  m_faceMappedCoordDir[1] = 0.5;
  m_faceMappedCoordDir[2] = -1;
}

//////////////////////////////////////////////////////////////////////
  
  void TriagFluxReconstructionElementData::createFaceIntegrationCoefs(){};

//////////////////////////////////////////////////////////////////////
  
  void TriagFluxReconstructionElementData::createCellAvgSolCoefs(){};

//////////////////////////////////////////////////////////////////////
  
  void TriagFluxReconstructionElementData::createCellCenterDerivCoefs(){};


  } // namespace FluxReconstructionMethod

} // namespace COOLFluiD

