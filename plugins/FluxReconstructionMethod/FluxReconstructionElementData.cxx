#include "Common/StringOps.hh"
#include "Common/NotImplementedException.hh"
#include "Common/CFLog.hh"
#include "FluxReconstructionMethod/FluxReconstructionElementData.hh"
#include "FluxReconstructionMethod/FluxReconstructionSolver.hh"

#include "FluxReconstructionMethod/FluxReconstruction.hh"


//////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Framework;

namespace COOLFluiD {

  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////

FluxReconstructionElementData::FluxReconstructionElementData() :
  m_dimensionality(),
  m_shape(),
  m_polyOrder(),
  m_solPntsLocalCoord1D(),
  m_flxPntsLocalCoord1D(),
  m_derivCoefsSolPnts1D(),
  m_solPolyDerivCoefsFlxPnts1D(),
  m_solPolyDerivCoefsSolPnts1D(),
  m_solPntsLocalCoords(),
  m_flxPntsLocalCoords(),
  m_faceFlxPntsFaceLocalCoords(),
  m_flxPntDerivDir(),
  m_allSolPntIdxs(),
  m_allFlxPntIdxs(),
  m_intFlxPntIdxs(),
  m_faceFlxPntConn(),
  m_faceFlxPntConnPerOrient(),
  m_faceFlxPntCellMappedCoords(),
  m_faceFlxPntCellMappedCoordsPerOrient(),
  m_cellNodeCoords(),
  m_faceNodeCoords(),
  m_faceNodeConn(),
  m_faceMappedCoordDir(),
  m_faceNodeConnPerOrient(),
  m_faceConnPerOrient(),
  m_faceMappedCoordDirPerOrient(),
  m_faceNodeCoordsPerOrient(),
  m_solPolyExponents(),
  m_solPolyCoefs(),
  m_nodePolyExponents(),
  m_nodePolyCoefs(),
  m_initPntsCoords(),
  m_initTransfMatrix(),
  m_faceIntegrationCoefs(),
  m_cellAvgSolCoefs(),
  m_cellCenterDerivCoefs(),
  m_flxPolyExponents(),
  m_flxPolyCoefs(),
  m_coefSolPolyDerivInSolPnts(),
  m_coefSolPolyDerivInNodes(),
  m_coefSolPolyDerivInFlxPnts(),
  m_cflConvDiffRatio(),
  m_faceOutputPntCellMappedCoords(),
  m_faceOutputPntSolPolyCoef(),
  m_faceOutputPntSolDerivCoef(),
  m_faceOutputPntConn(),
  m_solPntDistribution(),
  m_flxPntDistribution(),
  m_faceNormals(),
  m_coefSolPolyInFlxPnts(),
  m_flxPntFlxDim(),
  m_vandermonde(),
  m_vandermondeInv(),
  m_coefSolPolyInNodes()
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////

FluxReconstructionElementData::~FluxReconstructionElementData()
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////

void FluxReconstructionElementData::setPolyOrder(CFPolyOrder::Type polyOrder)
{
  m_polyOrder = polyOrder;
  resetFluxReconstructionElementData();
}

//////////////////////////////////////////////////////////////////////

void FluxReconstructionElementData::setFlxPntDistribution(Common::SafePtr< BasePointDistribution > flxPntDist)
{
  m_flxPntsLocalCoord1D = flxPntDist->getLocalCoords1D(m_polyOrder);
  resetFluxReconstructionElementData();
}

//////////////////////////////////////////////////////////////////////

void FluxReconstructionElementData::setSolPntDistribution(Common::SafePtr< BasePointDistribution > solPntDist)
{
  m_solPntsLocalCoord1D = solPntDist->getLocalCoords1D(m_polyOrder);
  resetFluxReconstructionElementData();
}

//////////////////////////////////////////////////////////////////////

void FluxReconstructionElementData::resetFluxReconstructionElementData()
{
  CFAUTOTRACE;

  createFlxPntsLocalCoords();
  createSolPntsLocalCoords();
  createFaceFlxPntsFaceLocalCoords();
  createSolPolyExponents();
  computeSolPolyCoefs();
  computeInitPntsCoords();
  computeInitTransfMatrix();
  createAllSolPntIdxs();
  createAllFlxPntIdxs();
  createFaceNodeConnectivity();
  createFaceMappedCoordDir();
  createFaceNodeConnectivityPerOrient();
  createFaceFluxPntsConn();
  createFaceFluxPntsConnPerOrient();
  createCellNodeCoords();
  createFaceNodeCoords();
  createFaceNodeCoordsPerOrient();
  createNodePolyExponents();
  computeNodePolyCoefs();
  createFaceIntegrationCoefs();
  createCellAvgSolCoefs();
  createCellCenterDerivCoefs();
  setCFLConvDiffRatio();
  createCoefSolPolyDerivInSolPnts();
  createCoefSolPolyInFlxPnts();
  createFaceFlxPntsCellLocalCoords();
  createFaceOutputPntCellMappedCoords();
  createFaceOutputPntSolPolyAndDerivCoef();
  createFaceOutputPntConn();
  createFaceNormals();
  createFluxPntFluxDim();
  createVandermondeMatrix();
  createCoefSolPolyInNodes();
  createCoefSolPolyDerivInNodes();
}

//////////////////////////////////////////////////////////////////////

void FluxReconstructionElementData::computeDerivCoefsSolPnts1D()
{
  CFAUTOTRACE;

  // number of solution points
  const CFuint nbrSolPnts = m_polyOrder+1;

  // number of flux points
  const CFuint nbrFlxPnts = 2;

  // resize m_recCoefsFlxPnts1D
  m_derivCoefsSolPnts1D.resize(nbrSolPnts,nbrFlxPnts);

  // compute derivation coefficients
  for (CFuint iSol = 0; iSol < nbrSolPnts; ++iSol)
  {
    const CFreal ksiSol = m_solPntsLocalCoord1D[iSol];
    for (CFuint iFlx = 0; iFlx < nbrFlxPnts; ++iFlx)
    {
      const CFreal ksiFlx = 2.*iFlx-1.;
      m_derivCoefsSolPnts1D(iSol,iFlx) = 0.;
      for (CFuint iTerm = 0; iTerm < nbrFlxPnts; ++iTerm)
      {
        if (iTerm != iFlx)
        {
          const CFreal ksiTerm = m_flxPntsLocalCoord1D[iTerm];
          CFreal term = 1./(ksiFlx-ksiTerm);
          for (CFuint iFac = 0; iFac < nbrFlxPnts; ++iFac)
          {
            if (iFac != iFlx && iFac != iTerm)
            {
              const CFreal ksiFac = 2.*iFac-1.;
              term *= (ksiSol-ksiFac)/(ksiFlx-ksiFac);
            }
          }
          m_derivCoefsSolPnts1D(iSol,iFlx) += term;
        }
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////

void FluxReconstructionElementData::computeSolPolyDerivCoefsSolPnts1D()
{
  CFAUTOTRACE;

  // number of solution points
  const CFuint nbrSolPnts = m_polyOrder+1;

  // number of flux points
  const CFuint nbrFlxPnts = 2;

  // resize m_solPolyDerivCoefsFlxPnts1D
  m_solPolyDerivCoefsFlxPnts1D.resize(nbrFlxPnts,nbrSolPnts);

  // compute derivation coefficients
  for (CFuint iFlx = 0; iFlx < nbrFlxPnts; ++iFlx)
  {
    const CFreal ksiFlx = 2.*iFlx-1.;
    for (CFuint iSol = 0; iSol < nbrSolPnts; ++iSol)
    {
      const CFreal ksiSol = m_solPntsLocalCoord1D[iSol];
      m_solPolyDerivCoefsFlxPnts1D(iFlx,iSol) = 0.;
      for (CFuint iTerm = 0; iTerm < nbrSolPnts; ++iTerm)
      {
        if (iTerm != iSol)
        {
          const CFreal ksiTerm = m_solPntsLocalCoord1D[iTerm];
          CFreal term = 1./(ksiSol-ksiTerm);
          for (CFuint iFac = 0; iFac < nbrSolPnts; ++iFac)
          {
            if (iFac != iSol && iFac != iTerm)
            {
              const CFreal ksiFac = m_solPntsLocalCoord1D[iFac];
              term *= (ksiFlx-ksiFac)/(ksiSol-ksiFac);
            }
          }
          m_solPolyDerivCoefsFlxPnts1D(iFlx,iSol) += term;
        }
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////

void FluxReconstructionElementData::computeSolPolyDerivCoefsFlxPnts1D()
{
  CFAUTOTRACE;

  // number of solution points
  const CFuint nbrSolPnts = m_polyOrder+1;

  // resize m_solPolyDerivCoefsSolPnts1D
  m_solPolyDerivCoefsSolPnts1D.resize(nbrSolPnts,nbrSolPnts);

  // compute derivation coefficients
  for (CFuint iSol = 0; iSol < nbrSolPnts; ++iSol)
  {
    const CFreal ksiSol = m_solPntsLocalCoord1D[iSol];
    for (CFuint iSol2 = 0; iSol2 < nbrSolPnts; ++iSol2)
    {
      const CFreal ksiSol2 = m_solPntsLocalCoord1D[iSol2];
      m_solPolyDerivCoefsSolPnts1D(iSol,iSol2) = 0.;
      for (CFuint iTerm = 0; iTerm < nbrSolPnts; ++iTerm)
      {
        if (iTerm != iSol2)
        {
          const CFreal ksiTerm = m_solPntsLocalCoord1D[iTerm];
          CFreal term = 1./(ksiSol2-ksiTerm);
          for (CFuint iFac = 0; iFac < nbrSolPnts; ++iFac)
          {
            if (iFac != iSol2 && iFac != iTerm)
            {
              const CFreal ksiFac = m_solPntsLocalCoord1D[iFac];
              term *= (ksiSol-ksiFac)/(ksiSol2-ksiFac);
            }
          }
          m_solPolyDerivCoefsSolPnts1D(iSol,iSol2) += term;
        }
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////

vector< vector< CFreal > >
FluxReconstructionElementData::getSolPolyValsAtNode(vector< RealVector > nodeLocalCoords)
{
  cf_assert(nodeLocalCoords.size() > 0);
  cf_assert(nodeLocalCoords[0].size() == static_cast<CFuint>(m_dimensionality));

  // number of nodes passed
  const CFuint nbrNodes = nodeLocalCoords.size();

  // number of polynomials
  const CFuint nbrPolys = getNbrOfSolPnts();

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
        CFreal term = m_solPolyCoefs[iPoly][iTerm];

        // loop over coordinates x-y-z
        for (CFuint iCoor = 0; iCoor < static_cast<CFuint>(m_dimensionality); ++iCoor)
        {
          term *= pow(nodeLocalCoords[iNode][iCoor],m_solPolyExponents[iTerm][iCoor]);
        }

        // add term to polynomial value
        polyValsAtNodes[iNode][iPoly] += term;
      }
    }
  }
  return polyValsAtNodes;
}

//////////////////////////////////////////////////////////////////////

vector< vector< CFreal > >
FluxReconstructionElementData::getNodePolyValsAtPnt(vector< RealVector > pntLocalCoords)
{
  cf_assert(pntLocalCoords.size() > 0);
  cf_assert(pntLocalCoords[0].size() == static_cast<CFuint>(m_dimensionality));

  // number of nodes passed
  const CFuint nbrNodes = pntLocalCoords.size();

  // number of polynomials
  const CFuint nbrPolys = getNbrCornerNodes();

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
        CFreal term = m_nodePolyCoefs[iPoly][iTerm];

        // loop over coordinates x-y-z
        for (CFuint iCoor = 0; iCoor < static_cast<CFuint>(m_dimensionality); ++iCoor)
        {
          term *= pow(pntLocalCoords[iNode][iCoor],m_nodePolyExponents[iTerm][iCoor]);
        }

        // add term to polynomial value
        polyValsAtNodes[iNode][iPoly] += term;
      }
    }
  }
  return polyValsAtNodes;
}

//////////////////////////////////////////////////////////////////////

vector< vector< vector< CFreal > > >
    FluxReconstructionElementData::getSolPolyDerivsAtNode(vector< RealVector > nodeLocalCoords)
{
  //CFLog(VERBOSE,"getSolPolyDerivsAtNode\n");
  cf_assert(nodeLocalCoords.size() > 0);
  cf_assert(nodeLocalCoords[0].size() == static_cast<CFuint>(m_dimensionality));

  // number of nodes passed
  const CFuint nbrNodes = nodeLocalCoords.size();

  // number of polynomials
  const CFuint nbrPolys = getNbrOfSolPnts();

  // return variable
  vector< vector< vector< CFreal > > > polyDerivsAtNodes(nbrNodes);

  // compute polynomial values at given nodes
  for (CFuint iNode = 0; iNode < nbrNodes; ++iNode)
  {
    // resize
    polyDerivsAtNodes[iNode].resize(m_dimensionality);

    // loop over derivatives
    for (CFuint iDeriv = 0; iDeriv < static_cast<CFuint>(m_dimensionality); ++iDeriv)
    {
      // resize
      polyDerivsAtNodes[iNode][iDeriv].resize(nbrPolys);

      // loop over polynomials
      for (CFuint iPoly = 0; iPoly < nbrPolys; ++iPoly)
      {
        // loop over terms
        for (CFuint iTerm = 0; iTerm < nbrPolys; ++iTerm)
        {
          CFreal term = m_solPolyCoefs[iPoly][iTerm];

          // loop over coordinates x-y-z
          for (CFuint iCoor = 0; iCoor < static_cast<CFuint>(m_dimensionality); ++iCoor)
          {
            if (iCoor != iDeriv)
            {
              term *= pow(nodeLocalCoords[iNode][iCoor],m_solPolyExponents[iTerm][iCoor]);
            }
            else
            {
              if (m_solPolyExponents[iTerm][iCoor] != 0)
              {
                term *= m_solPolyExponents[iTerm][iCoor]
                        * pow(nodeLocalCoords[iNode][iCoor],m_solPolyExponents[iTerm][iCoor]-1.0);
              }
              else
              {
                term = 0.0;
              }
            }
          }

          // add term to polynomial value
          polyDerivsAtNodes[iNode][iDeriv][iPoly] += term;
        }
      }
    }
  }

  return polyDerivsAtNodes;
}

//////////////////////////////////////////////////////////////////////

void FluxReconstructionElementData::createAllSolPntIdxs()
{
  CFAUTOTRACE;

  // number of solution points
  const CFuint nbrSolPnts = getNbrOfSolPnts();

  // put indexes for all solution points in list
  m_allSolPntIdxs.resize(0);
  for (CFuint iSol = 0; iSol < nbrSolPnts; ++iSol)
  {
    m_allSolPntIdxs.push_back(iSol);
  }
}

//////////////////////////////////////////////////////////////////////

void FluxReconstructionElementData::createAllFlxPntIdxs()
{
  CFAUTOTRACE;

  // number of flux points
  const CFuint nbrFlxPnts = getNbrOfFlxPnts();

  // put indexes for all flux points in list
  m_allFlxPntIdxs.resize(0);
  for (CFuint iFlx = 0; iFlx < nbrFlxPnts; ++iFlx)
  {
    m_allFlxPntIdxs.push_back(iFlx);
  }
}

//////////////////////////////////////////////////////////////////////

void FluxReconstructionElementData::computeInitPntsCoords()
{
  CFAUTOTRACE;

  // set nodal set
  setInterpolationNodeSet(m_polyOrder,m_initPntsCoords);
}

//////////////////////////////////////////////////////////////////////

void FluxReconstructionElementData::computeInitTransfMatrix()
{
  CFAUTOTRACE;

  /// @note this could be done by just evaluating the lagrangian polynomials
  /// associated with the initialization points. Done here like this because
  /// this implementation is independent of cell type.

  // number of solution points
  const CFuint nbrSolPnts = getNbrOfSolPnts();

  // evaluate polynomials at initialization points
  vector< vector< CFreal > > lhs = getSolPolyValsAtNode(m_initPntsCoords);

  // fill in left hand side matrix
  RealMatrix lhsRM(nbrSolPnts,nbrSolPnts);
  for (CFuint iInit = 0; iInit < nbrSolPnts; ++iInit)
  {
    for (CFuint iSol = 0; iSol < nbrSolPnts; ++iSol)
    {
      lhsRM(iInit,iSol) = lhs[iInit][iSol];
    }
  }

  // compute values in transformation matrix
  m_initTransfMatrix.resize(nbrSolPnts,nbrSolPnts);
  InvertMatrix(lhsRM,m_initTransfMatrix);
}

//////////////////////////////////////////////////////////////////////

void FluxReconstructionElementData::createFaceNodeCoords()
{
  CFAUTOTRACE;

  // number of cell faces
  const CFuint nbrCellFaces = getNbrCellFaces();

  // create face-node coordinates
  m_faceNodeCoords.resize(nbrCellFaces);
  for (CFuint iFace = 0; iFace < nbrCellFaces; ++iFace)
  {
    // number of face nodes
    const CFuint nbrFaceNodes = m_faceNodeConn[iFace].size();

    m_faceNodeCoords[iFace].resize(nbrFaceNodes);
    for (CFuint iNode = 0; iNode < nbrFaceNodes; ++iNode)
    {
      m_faceNodeCoords[iFace][iNode].resize(m_dimensionality);

      // node ID
      const CFuint nodeID = m_faceNodeConn[iFace][iNode];

      // set coordinate
      m_faceNodeCoords[iFace][iNode] = m_cellNodeCoords[nodeID];
    }
  }
}

//////////////////////////////////////////////////////////////////////

void FluxReconstructionElementData::createFaceNodeCoordsPerOrient()
{
  CFAUTOTRACE;

  // number of orientations
  const CFuint nbrOrients = m_faceNodeConnPerOrient.size();

  // number of face nodes
  cf_assert(getNbrCellFaces() > 0);
  const CFuint nbrFaceNodes = m_faceNodeCoords[0].size();

  // create face-node coordinates per face connectivity orientation
  m_faceNodeCoordsPerOrient.resize(nbrOrients);
  for (CFuint iOrient = 0; iOrient < nbrOrients; ++iOrient)
  {
    // resize
    m_faceNodeCoordsPerOrient[iOrient].resize(2);
    for (CFuint iSide = 0; iSide < 2; ++iSide)
    {
      m_faceNodeCoordsPerOrient[iOrient][iSide].resize(nbrFaceNodes);
      for (CFuint iNode = 0; iNode < nbrFaceNodes; ++iNode)
      {
        m_faceNodeCoordsPerOrient[iOrient][iSide][iNode].resize(m_dimensionality);
      }
    }

    // fill variable
    for (CFuint iNode = 0; iNode < nbrFaceNodes; ++iNode)
    {
      // node IDs
      const CFuint nodeIDL = m_faceNodeConnPerOrient[iOrient][LEFT ][iNode];
      const CFuint nodeIDR = m_faceNodeConnPerOrient[iOrient][RIGHT][iNode];

      // set coordinates
      m_faceNodeCoordsPerOrient[iOrient][LEFT ][iNode] = m_cellNodeCoords[nodeIDL];
      m_faceNodeCoordsPerOrient[iOrient][RIGHT][iNode] = m_cellNodeCoords[nodeIDR];
    }
  }
}

//////////////////////////////////////////////////////////////////////

void FluxReconstructionElementData::computeSolPolyCoefs()
{
  CFAUTOTRACE;

  // number of solution polynomials
  const CFuint nbrSolPolys = getNbrOfSolPnts();

  // resize m_solPolyCoefs
  m_solPolyCoefs.resize(nbrSolPolys);
  for (CFuint iPoly = 0; iPoly < nbrSolPolys; ++iPoly)
  {
    m_solPolyCoefs[iPoly].resize(nbrSolPolys);
  }

  // variable for LHS of linear system
  RealMatrix lhs   (nbrSolPolys,nbrSolPolys);
  RealMatrix lhsInv(nbrSolPolys,nbrSolPolys);

  // create LHS matrix
  for (CFuint iPoly = 0; iPoly < nbrSolPolys; ++iPoly)
  {
    for (CFuint iTerm = 0; iTerm < nbrSolPolys; ++iTerm)
    {
      lhs(iPoly,iTerm) = 1.0;
      for (CFuint iCoor = 0; iCoor < static_cast<CFuint>(m_dimensionality); ++ iCoor)
      {
        lhs(iPoly,iTerm) *= pow(m_solPntsLocalCoords[iPoly][iCoor],m_solPolyExponents[iTerm][iCoor]);
      }
    }
  }

  // invert the LHS matrix
  InvertMatrix(lhs,lhsInv);

  // store solution polynomial coefficients
  for (CFuint iPoly = 0; iPoly < nbrSolPolys; ++iPoly)
  {
    for (CFuint iTerm = 0; iTerm < nbrSolPolys; ++iTerm)
    {
      m_solPolyCoefs[iPoly][iTerm] = lhsInv(iTerm,iPoly);
    }
  }
}

//////////////////////////////////////////////////////////////////////

void FluxReconstructionElementData::computeNodePolyCoefs()
{
  CFAUTOTRACE;

  // number of solution polynomials
  const CFuint nbrNodePolys = getNbrCornerNodes();

  // resize m_solPolyCoefs
  m_nodePolyCoefs.resize(nbrNodePolys);
  for (CFuint iPoly = 0; iPoly < nbrNodePolys; ++iPoly)
  {
    m_nodePolyCoefs[iPoly].resize(nbrNodePolys);
  }

  // variable for LHS of linear system
  RealMatrix lhs   (nbrNodePolys,nbrNodePolys);
  RealMatrix lhsInv(nbrNodePolys,nbrNodePolys);

  // create LHS matrix
  for (CFuint iPoly = 0; iPoly < nbrNodePolys; ++iPoly)
  {
    for (CFuint iTerm = 0; iTerm < nbrNodePolys; ++iTerm)
    {
      lhs(iPoly,iTerm) = 1.0;
      for (CFuint iCoor = 0; iCoor < static_cast<CFuint>(m_dimensionality); ++ iCoor)
      {
        lhs(iPoly,iTerm) *= pow(m_cellNodeCoords[iPoly][iCoor],m_nodePolyExponents[iTerm][iCoor]);
      }
    }
  }

  // invert the LHS matrix
  InvertMatrix(lhs,lhsInv);

  // store solution polynomial coefficients
  for (CFuint iPoly = 0; iPoly < nbrNodePolys; ++iPoly)
  {
    for (CFuint iTerm = 0; iTerm < nbrNodePolys; ++iTerm)
    {
      m_nodePolyCoefs[iPoly][iTerm] = lhsInv(iTerm,iPoly);
    }
  }
}

//////////////////////////////////////////////////////////////////////

void FluxReconstructionElementData::createCoefSolPolyDerivInSolPnts()
{
  //CFLog(VERBOSE,"createCoefSolPolyDerivInFlxPnts\n");
  m_coefSolPolyDerivInSolPnts = getSolPolyDerivsAtNode(m_solPntsLocalCoords);
}

//////////////////////////////////////////////////////////////////////

void FluxReconstructionElementData::createCoefSolPolyDerivInNodes()
{
  //CFLog(VERBOSE,"createCoefSolPolyDerivInFlxPnts\n");
  m_coefSolPolyDerivInNodes = getSolPolyDerivsAtNode(m_cellNodeCoords);
}

//////////////////////////////////////////////////////////////////////

void FluxReconstructionElementData::createCoefSolPolyInFlxPnts()
{
  //CFLog(VERBOSE,"createCoefSolPolyDerivInFlxPnts\n");
  m_coefSolPolyInFlxPnts = getSolPolyValsAtNode(m_flxPntsLocalCoords);
}

//////////////////////////////////////////////////////////////////////

void FluxReconstructionElementData::createCoefSolPolyInNodes()
{
  //CFLog(VERBOSE,"createCoefSolPolyDerivInFlxPnts\n");
  m_coefSolPolyInNodes = getSolPolyValsAtNode(m_cellNodeCoords);
}

//////////////////////////////////////////////////////////////////////

void FluxReconstructionElementData::createFaceFlxPntsCellLocalCoords()
{
  CFAUTOTRACE;
//   CFLog(VERBOSE,"createFaceFlxPntsCellLocalCoords\n");
  
  // number of face flux points
  const CFuint nbrFaceFlxPnts = getNbrOfFaceFlxPnts();

  // connectivity for all faces
  const CFuint nbrCellFaces = getNbrCellFaces();
  m_faceFlxPntCellMappedCoords.resize(nbrCellFaces);
  for (CFuint iFace = 0; iFace < nbrCellFaces; ++iFace)
  {
    for (CFuint iFlx = 0; iFlx < nbrFaceFlxPnts; ++iFlx)
    {
      const CFuint flxIdx = m_faceFlxPntConn[iFace][iFlx];
      m_faceFlxPntCellMappedCoords[iFace].push_back(m_flxPntsLocalCoords[flxIdx]);
    }
  }
  // connectivity for all orientations
  const CFuint nbrOrients = m_faceFlxPntConnPerOrient.size();
  m_faceFlxPntCellMappedCoordsPerOrient.resize(nbrOrients);
  for (CFuint iOrient = 0; iOrient < nbrOrients; ++iOrient)
  {
    m_faceFlxPntCellMappedCoordsPerOrient[iOrient].resize(2);
    for (CFuint iFlx = 0; iFlx < nbrFaceFlxPnts; ++iFlx)
    {
      const CFuint flxIdxL = m_faceFlxPntConnPerOrient[iOrient][LEFT ][iFlx];
      m_faceFlxPntCellMappedCoordsPerOrient[iOrient][LEFT ].push_back(m_flxPntsLocalCoords[flxIdxL]);
      const CFuint flxIdxR = m_faceFlxPntConnPerOrient[iOrient][RIGHT][iFlx];
      m_faceFlxPntCellMappedCoordsPerOrient[iOrient][RIGHT].push_back(m_flxPntsLocalCoords[flxIdxR]);
    }
  }
}

//////////////////////////////////////////////////////////////////////

void FluxReconstructionElementData::createFaceOutputPntSolPolyAndDerivCoef()
{
  CFAUTOTRACE;
  //CFLog(VERBOSE,"createFaceOutputPntSolPolyAndDerivCoef\n");

  // get number of cell faces
  const CFuint nbrCellFaces = getNbrCellFaces();
  cf_assert(m_faceOutputPntCellMappedCoords.size() == nbrCellFaces);

  // compute polynomial values in cell face output points
  m_faceOutputPntSolPolyCoef .resize(nbrCellFaces);
  m_faceOutputPntSolDerivCoef.resize(nbrCellFaces);
  for (CFuint iFace = 0; iFace < nbrCellFaces; ++iFace)
  {
    m_faceOutputPntSolPolyCoef [iFace] = getSolPolyValsAtNode  (m_faceOutputPntCellMappedCoords[iFace]);
    m_faceOutputPntSolDerivCoef[iFace] = getSolPolyDerivsAtNode(m_faceOutputPntCellMappedCoords[iFace]);
  }
}

//////////////////////////////////////////////////////////////////////

void FluxReconstructionElementData::InvertMatrix(RealMatrix A, RealMatrix& AI)
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

void FluxReconstructionElementData::SwapRows(RealMatrix& A, CFuint row1, CFuint row2)
{
  RealVector swapRow = A.getRow<RealVector>(row1);
  A.setRow(A.getRow<RealVector>(row2),row1);
  A.setRow(swapRow,row2);
}

//////////////////////////////////////////////////////////////////////

CFreal FluxReconstructionElementData::evaluateLegendre(CFreal coord, CFuint order)
{
  CFreal result = 0.0;
  switch(order)
    {
      case 0:
      {
	result = 1;

      } break;
      case 1:
      {
	result = coord;

      } break;
      case 2:
      {
	result = 0.5*(3.0*pow(coord,2)-1.0);

      } break;
      case 3:
      {
	result = 0.5*(5.0*pow(coord,3)-3.0*coord);

      } break;
      case 4:
      {
	result = 0.125*(35.0*pow(coord,4)-30.0*pow(coord,2)+3.0);

      } break;
      case 5:
      {
	result = 0.125*(63.0*pow(coord,5)-70.0*pow(coord,3)+15.0*coord);

      } break;
      default:
      {
        throw Common::NotImplementedException (FromHere(),"Legendre not implemented for order "
                                      + StringOps::to_str(order) + ".");
      }
    }
  return result;
}

//////////////////////////////////////////////////////////////////////

  } // namespace FluxReconstructionMethod

} // namespace COOLFluiD
