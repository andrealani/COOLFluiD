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
  //CFLog(VERBOSE,"computeLocalFaceNormals\n");
  
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
  //CFLog(VERBOSE,"createFaceNodeConnectivityPerOrient\n");

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

void TriagFluxReconstructionElementData::createSVFaceNodeConnectivityPerOrientNoSymm()
{
  CFAUTOTRACE;

//   // number of SV faces
//   const CFuint nbrSVFaces = m_svFaceNodeConn.size();
// 
//   // number of possible orientations
//   const CFuint nbrOrient = 9;
// 
//   // resize the variable
//   m_svFaceNodeConnPerOrientationNoSymm.resize(nbrOrient);
//   for (CFuint iOrient = 0; iOrient < nbrOrient; ++iOrient)
//   {
//     m_svFaceNodeConnPerOrientationNoSymm[iOrient].resize(2);
//     for (CFuint iCell = 0; iCell < 2; ++iCell)
//       m_svFaceNodeConnPerOrientationNoSymm[iOrient][iCell].resize(2);
//   }
// 
//   // fill the variable
//   CFuint iOrient = 0;
//   for (CFuint iSVFaceL = 0; iSVFaceL < nbrSVFaces; ++iSVFaceL)
//   {
//     for (CFuint iSVFaceR = 0; iSVFaceR < nbrSVFaces; ++iSVFaceR, ++iOrient)
//     {
//       for (CFuint iNode = 0; iNode < 2; ++iNode)
//       {
//         m_svFaceNodeConnPerOrientationNoSymm[iOrient][LEFT ][iNode] = m_svFaceNodeConn[iSVFaceL][iNode  ];
//         m_svFaceNodeConnPerOrientationNoSymm[iOrient][RIGHT][iNode] = m_svFaceNodeConn[iSVFaceR][1-iNode];
//       }
//     }
//   }
}

//////////////////////////////////////////////////////////////////////
  
void TriagFluxReconstructionElementData::createSolPolyExponents()
{
  CFAUTOTRACE;
  //CFLog(VERBOSE,"createSolPolyExponents\n");

  // define exponents
  m_solPolyExponents.resize(0);
  for (CFuint iKsi = 0; iKsi < m_polyOrder+1; ++iKsi)
  {
    for (CFuint iEta = 0; iEta < m_polyOrder+1-iKsi; ++iEta)
    {
      vector< CFint > solPolyExps(2);
      solPolyExps[KSI] = iKsi;
      solPolyExps[ETA] = iEta;
      m_solPolyExponents.push_back(solPolyExps);
    }
  }
  cf_assert(m_solPolyExponents.size() == (m_polyOrder+1)*(m_polyOrder+2)/2);
}

//////////////////////////////////////////////////////////////////////

// void TriagFluxReconstructionElementData::createPolyExponents()
// {
//   CFAUTOTRACE;
// 
//   // helper variable
//   const CFuint polyOrderP1 = m_polyOrder + 1;
// 
//   // number of polynomial terms
//   const CFuint nbrPolyTerms = (m_polyOrder + 1)*(m_polyOrder + 2)/2;
// 
//   // resize the variable
//   m_solPolyExponents.resize(nbrPolyTerms);
//   for (CFuint iTerm = 0; iTerm < nbrPolyTerms; ++iTerm)
//   {
//     m_solPolyExponents[iTerm].resize(2);
//   }
// 
//   // define exponents
//   CFuint iTerm = 0;
//   for (CFuint iP = 0; iP < polyOrderP1; ++iP)
//   {
//     for (CFuint iY = 0; iY < iP+1; ++iY, ++iTerm)
//     {
//       m_solPolyExponents[iTerm][KSI] = iP-iY;
//       m_solPolyExponents[iTerm][ETA] = iY;
//     }
//   }
// }

//////////////////////////////////////////////////////////////////////

void TriagFluxReconstructionElementData::computePolyCoefs()
{
  CFAUTOTRACE;

//   // set the dimensionality and order of the simplex integrator
//   m_sIntegrator.setDimensionality(m_dimensionality);
//   m_sIntegrator.setIntegratorOrder(m_polyOrder);
// 
//   // number of control volumes in a SV
//   const CFuint nbrOfCVs = m_localCVNodeConn.size();
// 
//   // matrix for linear system
//   RealMatrix polyCoefSystem(nbrOfCVs,nbrOfCVs);
// 
//   // fill the linear system matrix
//   for (CFuint iCV = 0; iCV < nbrOfCVs; ++iCV)
//   {
//     // number of control volume nodes
//     const CFuint nbrCVNodes = m_localCVNodeConn[iCV].size();
// 
//     // number of triangles in polygon
//     const CFuint nbrTriangles = nbrCVNodes - 2;
// 
//     // get node coordinates
//     vector< RealVector > cvNodeCoord(nbrCVNodes);
//     for (CFuint iNode = 0; iNode < nbrCVNodes; ++iNode)
//     {
//       cvNodeCoord[iNode].resize(m_dimensionality);
//       const CFuint nodeID = m_localCVNodeConn[iCV][iNode];
//       cvNodeCoord[iNode] = m_localNodeCoord[nodeID];
//     }
// 
//     // get quadrature nodes and wheights
//     vector< RealVector > qNodeCoord;
//     vector< CFreal > qWheights;
//     for (CFuint iTriangle = 0; iTriangle < nbrTriangles; ++iTriangle)
//     {
//       // get node coordinates of triangle
//       vector< RealVector > triagNodeCoord(3);
//       for (CFuint iNode = 0; iNode < 3; ++iNode)
//       {
//         triagNodeCoord[iNode].resize(m_dimensionality);
//       }
//       triagNodeCoord[0] = cvNodeCoord[0];
//       triagNodeCoord[1] = cvNodeCoord[iTriangle+1];
//       triagNodeCoord[2] = cvNodeCoord[iTriangle+2];
// 
//       // get triangle quadrature nodes and wheights
//       vector< RealVector > qNodeCoordTriangle = m_sIntegrator.getQuadPntsCoords  (triagNodeCoord);
//       vector< CFreal >     qWheightsTriangle  = m_sIntegrator.getQuadPntsWheights(triagNodeCoord);
// 
//       // add triangle quadrature nodes and wheights to global list
//       qNodeCoord.insert(qNodeCoord.end(),qNodeCoordTriangle.begin(),qNodeCoordTriangle.end());
//       qWheights .insert(qWheights .end(),qWheightsTriangle .begin(),qWheightsTriangle .end());
//     }
// 
//     // compute left hand side of linear system
//     const CFuint nbrQNodes = qWheights.size();
//     for (CFuint iTerm = 0; iTerm < nbrOfCVs; ++iTerm)
//     {
//       polyCoefSystem(iCV,iTerm) = 0;
// 
//       for (CFuint iQNode = 0; iQNode < nbrQNodes; ++iQNode)
//       {
//         polyCoefSystem(iCV,iTerm) += qWheights[iQNode]
//                                         *pow(qNodeCoord[iQNode][KSI],m_polyExponents[iTerm][KSI])
//                                         *pow(qNodeCoord[iQNode][ETA],m_polyExponents[iTerm][ETA]);
//       }
//     }
//   }
// 
//   // invert the linear system matrix
//   RealMatrix invLinSysMatrix(nbrOfCVs,nbrOfCVs);
//   InvertMatrix(polyCoefSystem,invLinSysMatrix);
// 
//   // right hand side of linear system
//   RealMatrix rhsLinSys(nbrOfCVs,nbrOfCVs,0.0);
//   for (CFuint iCV = 0; iCV < nbrOfCVs; ++iCV)
//   {
//     rhsLinSys(iCV,iCV) = 0.5*m_volFracCV[iCV];
//   }
// 
//   // multiply imverted linear system matrix with rhs
//   RealMatrix polyCoefMatrix(nbrOfCVs,nbrOfCVs);
//   polyCoefMatrix = invLinSysMatrix*rhsLinSys;
// 
//   // store polynomial coefficients in m_polyCoefSFV
//   m_polyCoefSFV.resize(nbrOfCVs);
//   for (CFuint iPoly = 0; iPoly < nbrOfCVs; ++iPoly)
//   {
// //    CF_DEBUG_OBJ(iPoly);
//     m_polyCoefSFV[iPoly].resize(nbrOfCVs);
//     for (CFuint iTerm = 0; iTerm < nbrOfCVs; ++iTerm)
//     {
//       m_polyCoefSFV[iPoly][iTerm] = polyCoefMatrix(iTerm,iPoly);
// //      cout << m_polyCoefSFV[iPoly][iTerm] << "*ksi" << m_polyExponents[iTerm][KSI] << "eta" << m_polyExponents[iTerm][ETA] << " + ";
//     }
// //    cout << endl;
//   }
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

void TriagFluxReconstructionElementData::createFaceFluxPolyNodeWheightCoord()
{
   CFAUTOTRACE;
// 
//   // get the number of face flux points
//   const CFuint nbrFaceFlxPnts = getNbrOfFaceFlxPnts();
// 
//   // resize m_faceFluxPolyNodeWheightCoord
//   m_faceFluxPolyNodeWheightCoord.resize(nbrFaceFlxPnts);
// 
//   // dimensionality -1
//   const CFuint dimM1 = m_dimensionality - 1;
// 
//   // use the first SV face (for which the last local coordinate should be zero)
//   for (CFuint iFlx = 0; iFlx < nbrFaceFlxPnts; ++iFlx)
//   {
//     // get node ID
//     const CFuint nodeID = m_faceFlxPntsConn[0][iFlx];
// 
//     // resize
//     m_faceFluxPolyNodeWheightCoord[iFlx].resize(m_dimensionality);
// 
//     // compute wheight coordinates
//     m_faceFluxPolyNodeWheightCoord[iFlx][KSI] = 1.0 - m_fluxPolyNodeCoord[nodeID].sum();
//     for (CFuint iCoor = 0; iCoor < dimM1; ++iCoor)
//     {
//       m_faceFluxPolyNodeWheightCoord[iFlx][iCoor+1] = m_fluxPolyNodeCoord[nodeID][iCoor];
//     }
//   }
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
  //CFLog(VERBOSE,"setInterpolationNodeSet\n");

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
  
  
//   // number of nodes
//   const CFuint nbrNodes = (order+1)*(order+2)/2;
// 
//   // return variable
//   nodalSet.resize(nbrNodes);
//   for (CFuint iNode = 0; iNode < nbrNodes; ++iNode)
//   {
//     nodalSet[iNode].resize(2);
//   }
// 
//   // fill the vector containing the flux polynomial node coordinates
//   // Legendre-Gauss-Lobatto nodes from Hesthaven, FROM ELECTROSTATICS TO ALMOST OPTIMAL NODAL SETS
//   // FOR POLYNOMIAL INTERPOLATION IN A SIMPLEX. SIAM J. Numer. Anal. (1998); 35(2):655-676.
//   switch (order)
//   {
//     case CFPolyOrder::ORDER0:
//     {
//       nodalSet[0][KSI] = 1.0/3.0;
//       nodalSet[0][ETA] = 1.0/3.0;
//     } break;
//     case CFPolyOrder::ORDER1:
//     {
//       nodalSet[0][KSI] = 0.0;
//       nodalSet[0][ETA] = 0.0;
// 
//       nodalSet[1][KSI] = 1.0;
//       nodalSet[1][ETA] = 0.0;
// 
//       nodalSet[2][KSI] = 0.0;
//       nodalSet[2][ETA] = 1.0;
//     } break;
//     case CFPolyOrder::ORDER2:
//     {
//       nodalSet[0][KSI] = 0.0;
//       nodalSet[0][ETA] = 0.0;
// 
//       nodalSet[1][KSI] = 0.5;
//       nodalSet[1][ETA] = 0.0;
// 
//       nodalSet[2][KSI] = 0.0;
//       nodalSet[2][ETA] = 0.5;
// 
//       nodalSet[3][KSI] = 1.0;
//       nodalSet[3][ETA] = 0.0;
// 
//       nodalSet[4][KSI] = 0.5;
//       nodalSet[4][ETA] = 0.5;
// 
//       nodalSet[5][KSI] = 0.0;
//       nodalSet[5][ETA] = 1.0;
//     } break;
//     case CFPolyOrder::ORDER3:
//     {
//       nodalSet[0][KSI] = 0.0;
//       nodalSet[0][ETA] = 0.0;
// 
//       nodalSet[1][KSI] = 0.25;
//       nodalSet[1][ETA] = 0.0;
// 
//       nodalSet[2][KSI] = 0.0;
//       nodalSet[2][ETA] = 0.25;
// 
//       nodalSet[3][KSI] = 0.75;
//       nodalSet[3][ETA] = 0.0;
// 
//       nodalSet[4][KSI] = 1.0/3.0;
//       nodalSet[4][ETA] = 1.0/3.0;
// 
//       nodalSet[5][KSI] = 0.0;
//       nodalSet[5][ETA] = 0.75;
// 
//       nodalSet[6][KSI] = 1.0;
//       nodalSet[6][ETA] = 0.0;
// 
//       nodalSet[7][KSI] = 0.75;
//       nodalSet[7][ETA] = 0.25;
// 
//       nodalSet[8][KSI] = 0.25;
//       nodalSet[8][ETA] = 0.75;
// 
//       nodalSet[9][KSI] = 0.0;
//       nodalSet[9][ETA] = 1.0;
//     } break;
//     case CFPolyOrder::ORDER4:
//     {
//       nodalSet[0][KSI] = 0.0;
//       nodalSet[0][ETA] = 0.0;
// 
//       nodalSet[1][KSI] = 0.5*(1-cos(MathTools::MathConsts::CFrealPi()/4.0));
//       nodalSet[1][ETA] = 0.0;
// 
//       nodalSet[2][KSI] = 0.0;
//       nodalSet[2][ETA] = 0.5*(1-cos(MathTools::MathConsts::CFrealPi()/4.0));
// 
//       nodalSet[3][KSI] = 0.5;
//       nodalSet[3][ETA] = 0.0;
// 
//       nodalSet[4][KSI] = 0.2371200168;
//       nodalSet[4][ETA] = 0.2371200168;
// 
//       nodalSet[5][KSI] = 0.0;
//       nodalSet[5][ETA] = 0.5;
// 
//       nodalSet[6][KSI] = 0.5*(1+cos(MathTools::MathConsts::CFrealPi()/4.0));
//       nodalSet[6][ETA] = 0.0;
// 
//       nodalSet[7][KSI] = 0.5257599664;
//       nodalSet[7][ETA] = 0.2371200168;
// 
//       nodalSet[8][KSI] = 0.2371200168;
//       nodalSet[8][ETA] = 0.5257599664;
// 
//       nodalSet[9][KSI] = 0.0;
//       nodalSet[9][ETA] = 0.5*(1+cos(MathTools::MathConsts::CFrealPi()/4.0));
// 
//       nodalSet[10][KSI] = 1.0;
//       nodalSet[10][ETA] = 0.0;
// 
//       nodalSet[11][KSI] = 0.5*(1+cos(MathTools::MathConsts::CFrealPi()/4.0));
//       nodalSet[11][ETA] = 0.5*(1-cos(MathTools::MathConsts::CFrealPi()/4.0));
// 
//       nodalSet[12][KSI] = 0.5;
//       nodalSet[12][ETA] = 0.5;
// 
//       nodalSet[13][KSI] = 0.5*(1-cos(MathTools::MathConsts::CFrealPi()/4.0));
//       nodalSet[13][ETA] = 0.5*(1+cos(MathTools::MathConsts::CFrealPi()/4.0));
// 
//       nodalSet[14][KSI] = 0.0;
//       nodalSet[14][ETA] = 1.0;
//     } break;
//     default:
//     {
//       throw Common::NotImplementedException (FromHere(),"Higher-order nodal set not defined!");
//     }
//   }
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
  //CFLog(VERBOSE,"createFaceOutputPntCellMappedCoords\n");
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
      //RealVector temp = fun0*faceNodeCoords[0]+fun1*faceNodeCoords[1];
      //CFLog(VERBOSE,"Face output Pnt coord: (" << temp[0] << " , " << temp[1] << ")\n");
    }
  }
}

//////////////////////////////////////////////////////////////////////

void TriagFluxReconstructionElementData::createFaceOutputPntConn()
{
  CFAUTOTRACE;
  //CFLog(VERBOSE,"createFaceOutputPntConn\n");

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
  //CFLog(VERBOSE,"createFlxPntsLocalCoords\n");
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
  //CFLog(VERBOSE,"createSolPntsLocalCoords\n");
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
  //CFLog(VERBOSE,"createFaceFlxPntsFaceLocalCoords\n");

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
   
void TriagFluxReconstructionElementData::createFlxPntDerivDir()
{
//   CFLog(VERBOSE,"createFlxPntDerivDir\n");
//   // number of flux points in 1D
//   const CFuint nbrFlxPnts1D = m_flxPntsLocalCoord1D.size();
// 
//   // set derivation directions
//   m_flxPntDerivDir.resize(0);
//   // loop over 1D Flx pnts
//   for (CFuint i = 0; i < nbrFlxPnts1D; ++i)
//   {
//     // Side at ksi = -1
//     m_flxPntDerivDir.push_back(KSI);
//     
//     // Side at eta = -1
//     m_flxPntDerivDir.push_back(ETA);
//     
//     // Side at 45deg
//     m_flxPntDerivDir.push_back(ETA);
//   }
}

//////////////////////////////////////////////////////////////////////
  
  void TriagFluxReconstructionElementData::createIntFlxPntIdxs(){}//internal Flx pnts not necessary


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
  
void TriagFluxReconstructionElementData::createFlxPolyExponents()//Shouldn't be necessary
{
}

//////////////////////////////////////////////////////////////////////
  
  void TriagFluxReconstructionElementData::createFlxPntMatrixIdxForReconstruction(){}

//////////////////////////////////////////////////////////////////////
  
  void TriagFluxReconstructionElementData::createSolPntIdxsForReconstruction(){}

//////////////////////////////////////////////////////////////////////
  
void TriagFluxReconstructionElementData::createSolPntMatrixIdxForDerivation()
{
  //CFLog(VERBOSE,"createSolPntMatrixIdxForDerivation\n");

  // number of solution points in 1D
  const CFuint nbrSolPnts1D = m_solPntsLocalCoord1D.size();

  // set indices
  m_solPntMatrixIdxForDerivation.resize(0);
  for (CFuint iKsi = 0; iKsi < nbrSolPnts1D; ++iKsi)
  {
    for (CFuint iEta = 0; iEta < nbrSolPnts1D-iKsi; ++iEta)
    {
      vector< CFuint > solIdxs(2);
      solIdxs[KSI] = iKsi;
      solIdxs[ETA] = iEta;
      m_solPntMatrixIdxForDerivation.push_back(solIdxs);
    }
  }
}

//////////////////////////////////////////////////////////////////////
  
  void TriagFluxReconstructionElementData::createFlxPntMatrixIdxForDerivation(){}

//////////////////////////////////////////////////////////////////////
  
  void TriagFluxReconstructionElementData::createSolPntIdxsForDerivation(){}

//////////////////////////////////////////////////////////////////////
  
  void TriagFluxReconstructionElementData::createFlxPntIdxsForDerivation(){}

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

