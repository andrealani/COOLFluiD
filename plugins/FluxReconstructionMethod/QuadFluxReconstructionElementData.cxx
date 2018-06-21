#include <fstream>


#include "Common/CFLog.hh"
#include "Environment/DirPaths.hh"
#include "Common/NotImplementedException.hh"
#include "MathTools/RealMatrix.hh"
#include "FluxReconstructionMethod/QuadFluxReconstructionElementData.hh"
#include "FluxReconstructionMethod/TensorProductGaussIntegrator.hh"

//////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Framework;

namespace COOLFluiD {

  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////

QuadFluxReconstructionElementData::QuadFluxReconstructionElementData() :
  FluxReconstructionElementData()
{
  m_shape = CFGeoShape::QUAD;
  m_dimensionality = DIM_2D;
  createTensorProductIdx();
}

//////////////////////////////////////////////////////////////////////

QuadFluxReconstructionElementData::QuadFluxReconstructionElementData(CFPolyOrder::Type polyOrder)
{
  m_shape = CFGeoShape::QUAD;
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
  createTensorProductIdx();
}

//////////////////////////////////////////////////////////////////////

QuadFluxReconstructionElementData::QuadFluxReconstructionElementData(CFPolyOrder::Type polyOrder, 
								     Common::SafePtr< BasePointDistribution > solPntDist, 
								     Common::SafePtr< BasePointDistribution > flxPntDist)
{
  m_shape = CFGeoShape::QUAD;
  m_dimensionality = DIM_2D;
  m_polyOrder = polyOrder;
  m_solPntsLocalCoord1D = solPntDist->getLocalCoords1D(polyOrder);
  m_flxPntsLocalCoord1D = flxPntDist->getLocalCoords1D(polyOrder);
  m_solPntDistribution = solPntDist;
  m_flxPntDistribution = flxPntDist;

  resetFluxReconstructionElementData();
  createTensorProductIdx();
}

//////////////////////////////////////////////////////////////////////

QuadFluxReconstructionElementData::~QuadFluxReconstructionElementData()
{
}

//////////////////////////////////////////////////////////////////////

void QuadFluxReconstructionElementData::createFlxPntsLocalCoords()
{
  CFAUTOTRACE;
  
  // number of flux points in 1D
  const CFuint nbrFlxPnts1D = m_flxPntsLocalCoord1D.size();

  // set flux point local coordinates
  m_flxPntsLocalCoords.resize(0); 
  
  
  for (CFuint i = 0; i < nbrFlxPnts1D; ++i)
  {
    RealVector flxCoords(2);
    flxCoords[KSI] = -1;
    flxCoords[ETA] = m_flxPntsLocalCoord1D[i];
    m_flxPntsLocalCoords.push_back(flxCoords);
    flxCoords[KSI] = 1;
    m_flxPntsLocalCoords.push_back(flxCoords);
    flxCoords[ETA] = -1;
    flxCoords[KSI] = m_flxPntsLocalCoord1D[i];
    m_flxPntsLocalCoords.push_back(flxCoords);
    flxCoords[ETA] = 1;
    m_flxPntsLocalCoords.push_back(flxCoords);
  }
  cf_assert(m_flxPntsLocalCoords.size() == 4*nbrFlxPnts1D);
}

//////////////////////////////////////////////////////////////////////

void QuadFluxReconstructionElementData::createSolPntsLocalCoords()
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
  cf_assert(m_solPntsLocalCoords.size() == nbrSolPnts1D*nbrSolPnts1D);
}

//////////////////////////////////////////////////////////////////////

void QuadFluxReconstructionElementData::createFlxSolDependencies()
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
    for (CFuint iEta = 0; iEta < nbrSolPnts1D; ++iEta, ++iSol)
    {
      m_solFlxDep[iSol].push_back(iEta*4);
      m_solFlxDep[iSol].push_back(iEta*4 + 1);
      m_solFlxDep[iSol].push_back(iKsi*4 + 2);
      m_solFlxDep[iSol].push_back(iKsi*4 + 3);

      m_flxSolDep[iEta*4].push_back(iSol);
      m_flxSolDep[iEta*4 + 1].push_back(iSol);
      m_flxSolDep[iKsi*4 + 2].push_back(iSol);
      m_flxSolDep[iKsi*4 + 3].push_back(iSol);

      for (CFuint jSol = 0; jSol < nbrSolPnts1D; ++jSol)
      {
        m_solSolDep[iSol].push_back(iKsi*nbrSolPnts1D + jSol);
        if (iSol != iEta + jSol*nbrSolPnts1D) m_solSolDep[iSol].push_back(iEta + jSol*nbrSolPnts1D);
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////

void QuadFluxReconstructionElementData::createFaceFlxPntsFaceLocalCoords()
{
  CFAUTOTRACE;

  // number of flux points in 1D
  const CFuint nbrFlxPnts1D = m_flxPntsLocalCoord1D.size();

  // set face flux point face local coordinates
  m_faceFlxPntsFaceLocalCoords.resize(0);
  for (CFuint iFlx = 0; iFlx < nbrFlxPnts1D; ++iFlx)
  {
    RealVector flxCoord(1);
    flxCoord[KSI] = m_flxPntsLocalCoord1D[iFlx];
    m_faceFlxPntsFaceLocalCoords.push_back(flxCoord);
  }
}

//////////////////////////////////////////////////////////////////////

void QuadFluxReconstructionElementData::createSolPolyExponents()
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

void QuadFluxReconstructionElementData::createNodePolyExponents()
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
      vector< CFint > nodePolyExps(2);
      nodePolyExps[KSI] = iKsi;
      nodePolyExps[ETA] = iEta;
      m_nodePolyExponents.push_back(nodePolyExps);
    }
  }
}

//////////////////////////////////////////////////////////////////////

void QuadFluxReconstructionElementData::createFaceFluxPntsConn()
{
  CFAUTOTRACE;

  // number of flux points in 1D
  const CFuint nbrFlxPnts1D = m_flxPntsLocalCoord1D.size();

  // resize m_faceFlxPntConn
  m_faceFlxPntConn.resize(4);

  // variable holding the face index
  CFuint faceIdx = 0;

  // zeroth face
  for (CFuint iSol = 0; iSol < nbrFlxPnts1D; ++iSol)
  {
    m_faceFlxPntConn[faceIdx].push_back(2+iSol*4);// 4*nbrFlxPnts1D-2-iSol
  }
  ++faceIdx;


  // first face
  for (CFuint iSol = 0; iSol < nbrFlxPnts1D; ++iSol)
  {
    m_faceFlxPntConn[faceIdx].push_back(1+iSol*4);// 4*nbrFlxPnts1D-3-iSol
  }
  ++faceIdx;

  // second face
  for (CFuint iSol = 0; iSol < nbrFlxPnts1D; ++iSol)
  {
    m_faceFlxPntConn[faceIdx].push_back(4*nbrFlxPnts1D-1-4*iSol);// 3+iSol*4
  }
  ++faceIdx;

  // third face
  for (CFuint iSol = 0; iSol < nbrFlxPnts1D; ++iSol)
  {
    m_faceFlxPntConn[faceIdx].push_back(4*nbrFlxPnts1D-4-4*iSol);//iSol*4 
  }
}

//////////////////////////////////////////////////////////////////////

void QuadFluxReconstructionElementData::createFaceFluxPntsConnPerOrient()
{
  CFAUTOTRACE;

  // number of orientations
  const CFuint nbrOrients = 10;

  // number of solution points in 1D
  const CFuint nbrFlxPnts1D = m_flxPntsLocalCoord1D.size();

  // create data structure
  m_faceFlxPntConnPerOrient.resize(nbrOrients);
  CFuint iOrient = 0;
  for (CFuint iFaceL = 0; iFaceL < 4; ++iFaceL)
  {
    for (CFuint iFaceR = iFaceL; iFaceR < 4; ++iFaceR, ++iOrient)
    {
      m_faceFlxPntConnPerOrient[iOrient].resize(2);
      for (CFuint iSol = 0; iSol < nbrFlxPnts1D; ++iSol)
      {
        m_faceFlxPntConnPerOrient[iOrient][LEFT ]
            .push_back(m_faceFlxPntConn[iFaceL][iSol               ]);
        m_faceFlxPntConnPerOrient[iOrient][RIGHT]
            .push_back(m_faceFlxPntConn[iFaceR][nbrFlxPnts1D-1-iSol]);
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////

void QuadFluxReconstructionElementData::createCellNodeCoords()
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

void QuadFluxReconstructionElementData::createFaceNodeConnectivity()
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

void QuadFluxReconstructionElementData::createFluxPntFluxDim()
{
  CFAUTOTRACE;

  m_flxPntFlxDim.resize(4*m_flxPntsLocalCoord1D.size());
  
  for (CFuint iFace = 0; iFace < m_faceFlxPntConn.size(); ++iFace)
  {
    for (CFuint iFlx = 0; iFlx < m_flxPntsLocalCoord1D.size(); ++iFlx)
    {
      if (iFace == 1 || iFace == 3)
      {
        m_flxPntFlxDim[m_faceFlxPntConn[iFace][iFlx]] = 0;
      }
      else
      {
	m_flxPntFlxDim[m_faceFlxPntConn[iFace][iFlx]] = 1;
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////

void QuadFluxReconstructionElementData::createFaceMappedCoordDir()
{
  CFAUTOTRACE;

  m_faceMappedCoordDir.resize(4);

  m_faceMappedCoordDir[0] = -1;
  m_faceMappedCoordDir[1] = 1;
  m_faceMappedCoordDir[2] = 1;
  m_faceMappedCoordDir[3] = -1;
}

//////////////////////////////////////////////////////////////////////

void QuadFluxReconstructionElementData::createFaceNormals()
{
  CFAUTOTRACE;

  m_faceNormals.resize(4);
  
  for (CFuint iFace = 0; iFace < 4; ++iFace)
  {
    m_faceNormals[iFace].resize(2);
  }

  m_faceNormals[0][0] = 0.;
  m_faceNormals[0][1] = -1.;
  m_faceNormals[1][0] = 1.;
  m_faceNormals[1][1] = 0.;
  m_faceNormals[2][0] = 0.;
  m_faceNormals[2][1] = 1.;
  m_faceNormals[3][0] = -1.;
  m_faceNormals[3][1] = 0.;
}

//////////////////////////////////////////////////////////////////////

void QuadFluxReconstructionElementData::createFaceNodeConnectivityPerOrient()
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

void QuadFluxReconstructionElementData::createFaceIntegrationCoefs()
{
  CFAUTOTRACE;

  // number of flux points on a face
  const CFuint nbrFlxPnts = m_flxPntsLocalCoord1D.size();

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

    const CFreal ksiFlx = m_flxPntsLocalCoord1D[iFlx];
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
          const CFreal ksiFac = m_flxPntsLocalCoord1D[iFac];
          quadPntPolyVal *= (ksiQPnt-ksiFac)/(ksiFlx-ksiFac);
        }
      }

      // add contribution of quadrature point to integration coefficient
      m_faceIntegrationCoefs[iFlx] += quadPntWheights[iQPnt]*quadPntPolyVal;
    }
  }
}

//////////////////////////////////////////////////////////////////////

void QuadFluxReconstructionElementData::createCellAvgSolCoefs()
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

void QuadFluxReconstructionElementData::createCellCenterDerivCoefs()
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

void QuadFluxReconstructionElementData::setInterpolationNodeSet(const CFPolyOrder::Type order,
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
      RealVector node(2);
      node[KSI] = coords[iKsi];
      node[ETA] = coords[iEta];
      nodalSet.push_back(node);
    }
  }
}

//////////////////////////////////////////////////////////////////////

void QuadFluxReconstructionElementData::setCFLConvDiffRatio()
{
  CFAUTOTRACE;

  switch(m_polyOrder)
  {
    case CFPolyOrder::ORDER0:
    {
      m_cflConvDiffRatio = 2.0; //4.0; // check this!
    } break;
    case CFPolyOrder::ORDER1:
    {
      m_cflConvDiffRatio = 4.0; //6.5; // check this!
    } break;
    case CFPolyOrder::ORDER2:
    {
      m_cflConvDiffRatio = 6.0; //17.0; // check this!
    } break;
    case CFPolyOrder::ORDER3:
    {
      m_cflConvDiffRatio = 8.0; //25.0; // check this!
    } break;
    case CFPolyOrder::ORDER4:
    {
      m_cflConvDiffRatio = 10.0; //50.0; // check this!
    } break;
    case CFPolyOrder::ORDER5:
    {
      m_cflConvDiffRatio = 12.0; //50.0; // check this!
    } break;
    case CFPolyOrder::ORDER6:
    {
      m_cflConvDiffRatio = 14.0; //50.0; // check this!
    } break;
    case CFPolyOrder::ORDER7:
    {
      m_cflConvDiffRatio = 16.0; //50.0; // check this!
    } break;
    case CFPolyOrder::ORDER8:
    {
      m_cflConvDiffRatio = 18.0; //50.0; // check this!
    } break;
    case CFPolyOrder::ORDER9:
    {
      m_cflConvDiffRatio = 20.0; //50.0; // check this!
    } break;
    case CFPolyOrder::ORDER10:
    {
      m_cflConvDiffRatio = 22.0; //50.0; // check this!
    } break;
    default:
    {
      throw Common::NotImplementedException (FromHere(),"Higher-order quadrilateral FR cell not defined!");
    }
  }
}

//////////////////////////////////////////////////////////////////////

void QuadFluxReconstructionElementData::createFaceOutputPntCellMappedCoords()
{
  //CFLog(VERBOSE,"createFaceOutputPntCellMappedCoords\n");
  // number of points on a face
  const CFuint nbrFacePnts = m_polyOrder == 0 ? 1 :m_polyOrder;

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

void QuadFluxReconstructionElementData::createFaceOutputPntConn()
{
  CFAUTOTRACE;
  //CFLog(VERBOSE,"createFaceOutputPntConn\n");

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

void QuadFluxReconstructionElementData::createTensorProductIdx()
{
  CFAUTOTRACE;
  
  const CFuint nbrSolPnts = m_solPntsLocalCoords.size();
  // number of solution points in 1D
  const CFuint nbrSolPnts1D = m_solPntsLocalCoord1D.size();
  
  m_tensorProdIdx.resize(nbrSolPnts);
  for (CFuint iSol = 0; iSol < nbrSolPnts; ++iSol)
  {
    m_tensorProdIdx[iSol].resize(2);
  }

  CFuint iSol = 0;
  for (CFuint iKsi = 0; iKsi < nbrSolPnts1D; ++iKsi)
  {
    for (CFuint iEta = 0; iEta < nbrSolPnts1D; ++iEta, ++iSol)
    {
      for (CFuint jKsi = 0; jKsi < nbrSolPnts1D-1; ++jKsi)
      {
	if (iKsi == jKsi)
	{
	  jKsi += 1;
	}
        m_tensorProdIdx[iSol][KSI].push_back(iEta+nbrSolPnts1D*jKsi);
      }
      for (CFuint jEta = 0; jEta < nbrSolPnts1D-1; ++jEta)
      {
	if (iEta == jEta)
	{
	  jEta += 1;
	}
        m_tensorProdIdx[iSol][ETA].push_back(jEta+nbrSolPnts1D*iKsi);
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////

void QuadFluxReconstructionElementData::createVandermondeMatrix()
{
  CFAUTOTRACE;
  
  const CFuint nbrSolPnts = m_solPntsLocalCoords.size();
  // number of solution points in 1D
  const CFuint nbrSolPnts1D = m_solPntsLocalCoord1D.size();
  
  m_vandermonde.resize(nbrSolPnts,nbrSolPnts);
  m_vandermondeInv.resize(nbrSolPnts,nbrSolPnts);
  
  if(m_polyOrder != CFPolyOrder::ORDER0 && m_polyOrder != CFPolyOrder::ORDER1)
  {
    for (CFuint iSol = 0; iSol < nbrSolPnts; ++iSol)
    {
      CFuint modalDof = 0;
      for (CFuint iOrder = 0; iOrder < nbrSolPnts1D; ++iOrder)
      {
        for (CFuint iOrderKsi = 0; iOrderKsi < iOrder; ++iOrderKsi, ++modalDof)
        {
	  CFuint iOrderEta = iOrder;
          m_vandermonde(iSol,modalDof) = evaluateLegendre(m_solPntsLocalCoords[iSol][KSI],iOrderKsi)*evaluateLegendre(m_solPntsLocalCoords[iSol][ETA],iOrderEta);
        }
        for (CFuint iOrderEta = 0; iOrderEta < iOrder+1; ++iOrderEta, ++modalDof)
        {
	  CFuint iOrderKsi = iOrder;
          m_vandermonde(iSol,modalDof) = evaluateLegendre(m_solPntsLocalCoords[iSol][KSI],iOrderKsi)*evaluateLegendre(m_solPntsLocalCoords[iSol][ETA],iOrderEta);
        }
      }
      cf_assert(modalDof == nbrSolPnts);
    }
    
    InvertMatrix(m_vandermonde,m_vandermondeInv);
  }
}

//////////////////////////////////////////////////////////////////////

  } // namespace FluxReconstructionMethod

} // namespace COOLFluiD
