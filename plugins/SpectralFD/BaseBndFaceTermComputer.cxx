#include "Common/BadValueException.hh"
#include "Framework/BaseTerm.hh"
#include "Framework/MethodStrategyProvider.hh"
#include "Framework/ConsistencyException.hh"

#include "SpectralFD/BaseBndFaceTermComputer.hh"
#include "SpectralFD/SpectralFD.hh"
#include "SpectralFD/SpectralFDElementData.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace SpectralFD {

//////////////////////////////////////////////////////////////////////////////

Framework::MethodStrategyProvider<
    BaseBndFaceTermComputer,SpectralFDMethodData,BaseBndFaceTermComputer,SpectralFDModule >
    BaseBndFaceTermComputerProvider("BaseBndFaceTermComputer");

//////////////////////////////////////////////////////////////////////////////

BaseBndFaceTermComputer::BaseBndFaceTermComputer(const std::string& name) :
  SpectralFDMethodStrategy(name),
  socket_extraVars("meanflow",false),/// @todo get the name of this socket in another way, not hardcoded, to make it more general
  socket_faceJacobVecSizeFaceFlxPnts("faceJacobVecSizeFaceFlxPnts"),
  m_updateVarSet(CFNULL),
  m_diffusiveVarSet(CFNULL),
  m_statesReconstr(CFNULL),
  m_bcStateComputer(CFNULL),
  m_riemannFluxComputer(CFNULL),
  m_faceDiffFluxComputer(CFNULL),
  m_cflConvDiffRatio(),
  m_flxPntsRecCoefs(CFNULL),
  m_solPntsDerivCoefs(CFNULL),
  m_faceFlxPntConn(CFNULL),
  m_faceMappedCoordDir(CFNULL),
  m_flxPntMatrixIdxForReconstruction(CFNULL),
  m_solPntIdxsForReconstruction(CFNULL),
  m_flxPntMatrixIdxForDerivation(CFNULL),
  m_solPntIdxsForDerivation(CFNULL),
  m_faceIntegrationCoefs(CFNULL),
  m_faceFlxPntsFaceLocalCoords(CFNULL),
  m_faceFlxPntCellMappedCoords(CFNULL),
  m_face(),
  m_orient(),
  m_cellVolume(),
  m_unitNormalFlxPnts(),
  m_faceJacobVecSizeFlxPnts(),
  m_faceJacobVecAbsSizeFlxPnts(),
  m_faceInvCharLengths(),
  m_cellExtraVars(),
  m_flxPntIntSol(),
  m_flxPntGhostSol(),
  m_flxPntIntExtraVars(),
  m_flxPntGhostExtraVars(),
  m_flxPntCoords(),
  m_allSol(),
  m_allExtraVars(),
  m_flxPntIntRVSol(),
  m_flxPntGhostRVSol(),
  m_flxPntRiemannFlux(),
  m_flxPntIntGrads(),
  m_flxPntGhostGrads(),
  m_flxPntIntGradPtrs(),
  m_flxPntGhostGradPtrs(),
  m_backupPhysVar(),
  m_backupGhostStates(),
  m_nbrFlxPnts(),
  m_gradTerm(),
  m_nbrEqs(),
  m_dim(),
  m_nbrExtraVars(),
  m_faceMappedCoordsPntSet(),
  m_cellMappedCoordsPntSet(),
  m_solRecCoefsPntSet(),
  m_statesPntSet(),
  m_intSolPntSet(),
  m_ghostSolPntSet(),
  m_intSolPntSetPtrs(),
  m_ghostSolPntSetPtrs(),
  m_unitNormalPntSet(),
  m_faceJacobPntSet(),
  m_coordsPntSet(),
  m_gradsPntSet(),
  m_intGradsPntSet(),
  m_ghostGradsPntSet()
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

BaseBndFaceTermComputer::~BaseBndFaceTermComputer()
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

void BaseBndFaceTermComputer::configure ( Config::ConfigArgs& args )
{
  SpectralFDMethodStrategy::configure(args);

  // initialize socket_extraVars
}

//////////////////////////////////////////////////////////////////////////////

void BaseBndFaceTermComputer::setFaceTermData()
{
  // get the SpectralFDElementData
  vector< SpectralFDElementData* >& sdLocalData = getMethodData().getSDLocalData();

  // get derivation coefficients for the solution points
  m_solPntsDerivCoefs = sdLocalData[0]->getDerivCoefsSolPnts1D();

  // get indexes of internal flux points
  m_faceFlxPntConn = sdLocalData[0]->getFaceFlxPntConn();

  // get flux point mapped coordinate directions
  m_faceMappedCoordDir = sdLocalData[0]->getFaceMappedCoordDir();

  // get flux point index (in the matrix flxPntRecCoefs) for reconstruction
  m_flxPntMatrixIdxForReconstruction = sdLocalData[0]->getFlxPntMatrixIdxForReconstruction();

  // get flux point index (in the matrix m_solPntsDerivCoefs) for derivation
  m_flxPntMatrixIdxForDerivation = sdLocalData[0]->getFlxPntMatrixIdxForDerivation();

  // get solution point index (in the cell) for derivation
  m_solPntIdxsForDerivation = sdLocalData[0]->getSolPntIdxsForDerivation();

  // get coefficients for integration over a face
  m_faceIntegrationCoefs = sdLocalData[0]->getFaceIntegrationCoefs();

  // get face local coordinates of face flux points
  m_faceFlxPntsFaceLocalCoords = sdLocalData[0]->getFaceFlxPntsFaceLocalCoords();

  // get convective/diffusive CFL ratio
  m_cflConvDiffRatio = sdLocalData[0]->getCFLConvDiffRatio();

  // set number of flux points
  const CFuint nbrOrients = m_faceFlxPntConn->size();
  m_nbrFlxPnts.resize(nbrOrients);
  for (CFuint iOrient = 0; iOrient < nbrOrients; ++iOrient)
  {
    m_nbrFlxPnts[iOrient] = (*m_faceFlxPntConn)[iOrient].size();
  }

  // set interpolation type dependent data
  const std::string interpolationType = getMethodData().getInterpolationType();
  if (interpolationType == "standard")
  {
    // get reconstruction coefficients for the flux points
    m_flxPntsRecCoefs = sdLocalData[0]->getRecCoefsFlxPnts1D();

    // get solution point index (in the cell) for reconstruction
    m_solPntIdxsForReconstruction = sdLocalData[0]->getSolPntIdxsForReconstruction();
  }
  else if (interpolationType == "optimized")
  {
    // get reconstruction coefficients for the flux points
    m_flxPntsRecCoefs = sdLocalData[0]->getRecCoefsFlxPnts1DOptim();

    // get solution point index (in the cell) for reconstruction
    m_solPntIdxsForReconstruction = sdLocalData[0]->getSolPntIdxsForRecOptim();
  }
  else
  {
    throw BadValueException (FromHere(),"BaseVolTermComputer::setFaceTermData --> Interpolation type should be standard or optimized");
  }

  // get face flux point cell mapped coordinates
  m_faceFlxPntCellMappedCoords = sdLocalData[0]->getFaceFlxPntCellMappedCoords();
}

//////////////////////////////////////////////////////////////////////////////

void BaseBndFaceTermComputer::computeFaceData()
{
  // face ID
  const CFuint faceID = m_face->getID();

  // compute face Jacobian vectors
  vector< RealVector > faceJacobVecs = m_face->computeFaceJacobDetVectorAtMappedCoords(*m_faceFlxPntsFaceLocalCoords);

  // get face Jacobian vector sizes in the flux points
  DataHandle< vector< CFreal > >
      faceJacobVecSizeFaceFlxPnts = socket_faceJacobVecSizeFaceFlxPnts.getDataHandle();

  // set face Jacobian vector sizes and unit normals in the flux points
  for (CFuint iFlx = 0; iFlx < m_nbrFlxPnts[m_orient]; ++iFlx)
  {
    // get face Jacobian vector size
    m_faceJacobVecAbsSizeFlxPnts[iFlx] = faceJacobVecSizeFaceFlxPnts[faceID][iFlx];

    // set face Jacobian vector size with sign depending on mapped coordinate direction
    m_faceJacobVecSizeFlxPnts[iFlx] =
        m_faceJacobVecAbsSizeFlxPnts[iFlx]*(*m_faceMappedCoordDir)[m_orient];

    // set unit normal vector
    m_unitNormalFlxPnts[iFlx] = faceJacobVecs[iFlx]/m_faceJacobVecAbsSizeFlxPnts[iFlx];
  }

  // set the face ID in the BCStateComputer
  m_bcStateComputer->setFaceID(m_face->getID());

  // if necessary, compute the flux point coordinates
  if (m_bcStateComputer->needsSpatialCoordinates())
  {
    computeFlxPntCoords();
  }
}

//////////////////////////////////////////////////////////////////////////////

void BaseBndFaceTermComputer::computeNeighbourCellData()
{
  // get neighbouring cell
  GeometricEntity* cell = m_face->getNeighborGeo(0);

  // compute volume
  m_cellVolume = cell->computeVolume();

  // compute Jacobian determinants
  std::valarray<CFreal> jacobDets(m_nbrFlxPnts[m_orient]);
  jacobDets = cell->computeGeometricShapeFunctionJacobianDeterminant((*m_faceFlxPntCellMappedCoords)[m_orient]);

  // compute inverse characteristic lengths
  for (CFuint iFlx = 0; iFlx < m_nbrFlxPnts[m_orient]; ++iFlx)
  {
    m_faceInvCharLengths[iFlx] = 0.5*m_faceJacobVecAbsSizeFlxPnts[iFlx]/jacobDets[iFlx];
  }
}

//////////////////////////////////////////////////////////////////////////////

void BaseBndFaceTermComputer::reconstructFluxPntsStates(const vector< State* >& cellIntStates, bool onlyExtraVars)
{
  // compute internal states
  if (!onlyExtraVars)
  {
    m_statesReconstr->reconstructStates(cellIntStates,m_flxPntIntSol,*m_flxPntsRecCoefs,
                                        (*m_faceFlxPntConn)[m_orient],
                                        *m_flxPntMatrixIdxForReconstruction,
                                        *m_solPntIdxsForReconstruction);
  }

  // if needed, reconstruct the extra variables
  if (m_nbrExtraVars > 0)
  {
    cf_assert(socket_extraVars.isConnected());
    DataHandle<RealVector> extraVars = socket_extraVars.getDataHandle();

    // get extra vars at the solution points in the internal cell
    const CFuint nbrSolPnts = cellIntStates.size();
    for (CFuint iSol = 0; iSol < nbrSolPnts; ++iSol)
    {
      const CFuint stateID = cellIntStates[iSol]->getLocalID();
      m_cellExtraVars[iSol] = &extraVars[stateID];
    }

    // reconstruct extra variables in internal cell
    m_statesReconstr->reconstructExtraVars(m_cellExtraVars,m_flxPntIntExtraVars,*m_flxPntsRecCoefs,
                                           (*m_faceFlxPntConn)[m_orient],
                                           *m_flxPntMatrixIdxForReconstruction,
                                           *m_solPntIdxsForReconstruction);

    // copy values to ghost cell extra variables
    const CFuint nbrFlxPnts = m_flxPntIntExtraVars.size();
    for (CFuint iFlx = 0; iFlx < nbrFlxPnts; ++iFlx)
    {
      *m_flxPntGhostExtraVars[iFlx] = *m_flxPntIntExtraVars[iFlx];
    }
  }

  // compute ghost states
  if (!onlyExtraVars)
  {
    computeGhostStates();
  }
}

//////////////////////////////////////////////////////////////////////////////

void BaseBndFaceTermComputer::reconstructFluxPntsGradients(const vector< vector< RealVector >* >& cellIntGrads)
{
  // reconstruct internal gradients
  m_statesReconstr->reconstructGradients(cellIntGrads,m_flxPntIntGrads,*m_flxPntsRecCoefs,
                                         (*m_faceFlxPntConn)[m_orient],
                                         *m_flxPntMatrixIdxForReconstruction,
                                         *m_solPntIdxsForReconstruction);

  // compute the ghost gradients
  computeGhostGradients();
}

//////////////////////////////////////////////////////////////////////////////

void BaseBndFaceTermComputer::backupAndReconstructPhysVar
    (const CFuint iVar, const vector< State* >& cellIntStates)
{
  // backup
  backupPhysVar(iVar);

  // reconstruct internal states
  m_statesReconstr
      ->reconstructPhysVar(iVar,cellIntStates,m_flxPntIntSol,*m_flxPntsRecCoefs,
                           (*m_faceFlxPntConn)[m_orient],
                           *m_flxPntMatrixIdxForReconstruction,
                           *m_solPntIdxsForReconstruction);

  // compute ghost states
  computePerturbedGhostStates();
}

//////////////////////////////////////////////////////////////////////////////

void BaseBndFaceTermComputer::backupPhysVar(const CFuint iVar)
{
  cf_assert(m_nbrFlxPnts[m_orient] <= m_backupPhysVar.size());
  cf_assert(m_nbrFlxPnts[m_orient] <= m_backupGhostStates.size());
  for (CFuint iFlx = 0; iFlx < m_nbrFlxPnts[m_orient]; ++iFlx)
  {
    m_backupPhysVar[iFlx] = (*m_flxPntIntSol[iFlx])[iVar];
    m_backupGhostStates[iFlx] = *m_flxPntGhostSol[iFlx];
  }
}

//////////////////////////////////////////////////////////////////////////////

void BaseBndFaceTermComputer::restorePhysVar(const CFuint iVar)
{
  cf_assert(m_nbrFlxPnts[m_orient] <= m_backupPhysVar.size());
  cf_assert(m_nbrFlxPnts[m_orient] <= m_backupGhostStates.size());
  for (CFuint iFlx = 0; iFlx < m_nbrFlxPnts[m_orient]; ++iFlx)
  {
    (*m_flxPntIntSol[iFlx])[iVar] = m_backupPhysVar[iFlx];
    *m_flxPntGhostSol[iFlx] = m_backupGhostStates[iFlx];
  }
}

//////////////////////////////////////////////////////////////////////////////

void BaseBndFaceTermComputer::computeGhostStates()
{
  cf_assert(m_bcStateComputer.isNotNull());

  // set extra variables
  m_bcStateComputer->setExtraVars(&m_flxPntIntExtraVars);

  // compute ghost states
  m_bcStateComputer->computeGhostStates(m_flxPntIntSol,m_flxPntGhostSol,m_unitNormalFlxPnts,m_flxPntCoords);
}

//////////////////////////////////////////////////////////////////////////////

void BaseBndFaceTermComputer::computePerturbedGhostStates()
{
  // compute ghost states
  cf_assert(m_bcStateComputer.isNotNull());
  m_bcStateComputer->computeGhostStates(m_flxPntIntSol,m_flxPntGhostSol,m_unitNormalFlxPnts,m_flxPntCoords);
}

//////////////////////////////////////////////////////////////////////////////

void BaseBndFaceTermComputer::computeGhostGradients()
{
  cf_assert(m_bcStateComputer.isNotNull());

  // compute ghost gradients
  m_bcStateComputer->
      computeGhostGradients(m_flxPntIntGrads,m_flxPntGhostGrads,m_unitNormalFlxPnts,m_flxPntCoords);
}

//////////////////////////////////////////////////////////////////////////////

void BaseBndFaceTermComputer::computeFlxPntCoords()
{
  // compute flux point coordinates
  for (CFuint iFlx = 0; iFlx < m_nbrFlxPnts[m_orient]; ++iFlx)
  {
    m_flxPntCoords[iFlx] =
        m_face->computeCoordFromMappedCoord((*m_faceFlxPntsFaceLocalCoords)[iFlx]);
  }
}

//////////////////////////////////////////////////////////////////////////////

void BaseBndFaceTermComputer::computeConvFaceTerm(RealVector& resUpdates)
{

  // evaluate the riemann fluxes in all the flux points
  m_flxPntRiemannFlux = m_riemannFluxComputer->computeFlux(m_flxPntIntSol,m_flxPntGhostSol,
                                                           m_unitNormalFlxPnts,
                                                           m_nbrFlxPnts[m_orient]);

  // compute the actual face term
  computeFaceTermFromFlxPntFluxes(resUpdates);

  // change the sign of the updates to the residuals
  resUpdates *= -1.0;
}

//////////////////////////////////////////////////////////////////////////////

void BaseBndFaceTermComputer::computeConvFaceTermAndWaveSpeedUpdates(RealVector& resUpdates,
                                                                     CFreal& waveSpeedUpd)
{
// compute the face term
  computeConvFaceTerm(resUpdates);

  // compute the wave speed updates for the neighbouring cell
  waveSpeedUpd = 0.0;
  for (CFuint iFlx = 0; iFlx < m_nbrFlxPnts[m_orient]; ++iFlx)
  {
    const CFreal jacobXIntCoef = m_faceJacobVecAbsSizeFlxPnts[iFlx]*
                                 (*m_faceIntegrationCoefs)[iFlx];
    
    // transform update states to physical data to calculate eigenvalues
    m_updateVarSet->computePhysicalData(*m_flxPntIntSol[iFlx], m_pData);
    waveSpeedUpd += jacobXIntCoef*
        m_updateVarSet->getMaxAbsEigenValue(m_pData,m_unitNormalFlxPnts[iFlx]);
  }
}

//////////////////////////////////////////////////////////////////////////////

void BaseBndFaceTermComputer::computeGradientFaceTerm(vector< vector< RealVector > >& gradUpdates)
{
  // compute the averaged gradient variables in the flux points (stored in m_flxPntRiemannFlux)
  // (m_leftRVSol and m_rightRVSol are RealVector pointers to the left and right states)
  m_flxPntRiemannFlux = m_faceDiffFluxComputer->computeAvgGradVars(m_flxPntIntRVSol,m_flxPntGhostRVSol,
                                                                  m_nbrFlxPnts[m_orient]);

  // compute the actual volume term contribution to the gradient
  computeGradFaceTermFromFlxPntSol(gradUpdates);
}

//////////////////////////////////////////////////////////////////////////////

void BaseBndFaceTermComputer::computeGradientExtraVarsFaceTerm(vector< vector< RealVector > >& gradUpdates)
{
  // compute the averaged extra variables in the flux points (stored in m_flxPntRiemannFlux)
  for (CFuint iFlx = 0; iFlx < m_nbrFlxPnts[m_orient]; ++iFlx)
  {
    m_flxPntRiemannFlux[iFlx] = *m_flxPntIntExtraVars[iFlx];
    // ghost extra vars are the same as internal extra vars
  }

  // compute the actual volume term contribution to the gradient
  computeGradFaceTermFromFlxPntSol(gradUpdates);
}

//////////////////////////////////////////////////////////////////////////////

void BaseBndFaceTermComputer::computeDiffFaceTerm(RealVector& resUpdates)
{
  // compute the diffusive fluxes in the face flux points
  m_flxPntRiemannFlux = m_faceDiffFluxComputer->computeDiffFlux(m_flxPntIntGradPtrs,m_flxPntGhostGradPtrs,
                                                                m_flxPntIntRVSol,m_flxPntGhostRVSol,
                                                                m_faceInvCharLengths,m_unitNormalFlxPnts,
                                                                m_nbrFlxPnts[m_orient]);

  // compute the actual face term
  computeFaceTermFromFlxPntFluxes(resUpdates);
}

//////////////////////////////////////////////////////////////////////////////

void BaseBndFaceTermComputer::computeDiffFaceTermAndUpdateCoefContributions(RealVector& resUpdates,
                                                                            CFreal& updateCoefContr)
{
// compute the face term
  computeDiffFaceTerm(resUpdates);

  // compute the update coefficient contributions for the neighbouring cells
  const CFreal visc = 1.0;
  updateCoefContr = 0.0;
  for (CFuint iFlx = 0; iFlx < m_nbrFlxPnts[m_orient]; ++iFlx)
  {
    const CFreal jacobXJacobXIntCoef = m_faceJacobVecAbsSizeFlxPnts[iFlx]*
                                       m_faceJacobVecAbsSizeFlxPnts[iFlx]*
                                       (*m_faceIntegrationCoefs)[iFlx]*
                                       m_cflConvDiffRatio;
    updateCoefContr += visc*jacobXJacobXIntCoef/m_cellVolume;
  }
}

//////////////////////////////////////////////////////////////////////////////

void BaseBndFaceTermComputer::setPointSet(const std::vector< RealVector >& faceMappedCoordPntSet)
{
  // get the SpectralFDElementData
  vector< SpectralFDElementData* >& sdLocalData = getMethodData().getSDLocalData();

  // get face node cell-mapped coordinates
  SafePtr< vector< vector< RealVector > > > faceNodeCoords = sdLocalData[0]->getFaceNodeCoords();

  // get current face nodes
  const vector< RealVector >& currFaceNodes = (*faceNodeCoords)[m_orient];

  // number of points in point set
  const CFuint nbrPnts = faceMappedCoordPntSet.size();

  // compute cell mapped coordinates of point set
  m_faceMappedCoordsPntSet.resize(0);
  m_cellMappedCoordsPntSet.resize(0);
  for (CFuint iPnt = 0; iPnt < nbrPnts; ++iPnt)
  {
    // get point face-mapped coordinates
    const RealVector& faceMappedCoord = faceMappedCoordPntSet[iPnt];

    // store face-mapped coordinates
    m_faceMappedCoordsPntSet.push_back(faceMappedCoord);

    // compute face node wheights
    vector< CFreal > nodeWheights;
    if (m_dim == 2)
    {
      cf_assert(faceMappedCoord.size() == 1);
      nodeWheights.push_back(0.5*(1.0-faceMappedCoord[KSI]));
      nodeWheights.push_back(0.5*(1.0+faceMappedCoord[KSI]));
    }
    else if (m_dim == 3)
    {
      cf_assert(faceMappedCoord.size() == 2);
      nodeWheights.push_back(0.25*(1.0-faceMappedCoord[KSI])*(1.0-faceMappedCoord[ETA]));
      nodeWheights.push_back(0.25*(1.0+faceMappedCoord[KSI])*(1.0-faceMappedCoord[ETA]));
      nodeWheights.push_back(0.25*(1.0+faceMappedCoord[KSI])*(1.0+faceMappedCoord[ETA]));
      nodeWheights.push_back(0.25*(1.0-faceMappedCoord[KSI])*(1.0+faceMappedCoord[ETA]));
    }

    // compute cell mapped coordinates of given point
    const CFuint nbrFaceNodes = nodeWheights.size();
    cf_assert(currFaceNodes.size() == nbrFaceNodes);
    RealVector cellMappedCoord(m_dim);
    for (CFuint iNode = 0; iNode < nbrFaceNodes; ++iNode)
    {
      cellMappedCoord += nodeWheights[iNode]*currFaceNodes[iNode];
    }
    m_cellMappedCoordsPntSet.push_back(cellMappedCoord);
  }

  // compute solution reconstruction coefficients
  m_solRecCoefsPntSet = sdLocalData[0]->getSolPolyValsAtNode(m_cellMappedCoordsPntSet);

  // resize some variables
  if (nbrPnts > m_statesPntSet.size())
  {
    m_statesPntSet.resize(nbrPnts,RealVector(m_nbrEqs));
    m_unitNormalPntSet.resize(nbrPnts,RealVector(m_dim));
    m_faceJacobPntSet .resize(nbrPnts,RealVector(m_dim));
    m_coordsPntSet    .resize(nbrPnts,RealVector(m_dim));
    m_gradsPntSet.resize(nbrPnts,vector<RealVector>(m_nbrEqs,RealVector(m_dim)));

    CFuint iPnt = m_intSolPntSet.size();
    for (; iPnt < nbrPnts; ++iPnt)
    {
      m_intSolPntSet  .push_back(new State());
      m_ghostSolPntSet.push_back(new State());

      // set pointers to internal and ghost states
      m_intSolPntSetPtrs  .push_back(m_intSolPntSet  [iPnt]);
      m_ghostSolPntSetPtrs.push_back(m_ghostSolPntSet[iPnt]);

      vector<RealVector*> intGrads  ;
      vector<RealVector*> ghostGrads;
      for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
      {
        intGrads  .push_back(new RealVector(m_dim));
        ghostGrads.push_back(new RealVector(m_dim));
      }
      m_intGradsPntSet  .push_back(intGrads  );
      m_ghostGradsPntSet.push_back(ghostGrads);
    }
  }
  else if (nbrPnts < m_statesPntSet.size())
  {
    m_statesPntSet.resize(nbrPnts,RealVector(m_nbrEqs));
    m_unitNormalPntSet.resize(nbrPnts,RealVector(m_dim));
    m_faceJacobPntSet .resize(nbrPnts,RealVector(m_dim));
    m_coordsPntSet    .resize(nbrPnts,RealVector(m_dim));
    m_gradsPntSet.resize(nbrPnts,vector<RealVector>(m_nbrEqs,RealVector(m_dim)));

    CFuint iPnt = m_intSolPntSet.size()-1;
    for (; iPnt >= nbrPnts; --iPnt)
    {
      deletePtr(m_intSolPntSet  [iPnt]);
      deletePtr(m_ghostSolPntSet[iPnt]);

      for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
      {
        deletePtr(m_intGradsPntSet  [iPnt][iEq]);
        deletePtr(m_ghostGradsPntSet[iPnt][iEq]);
      }

      // remove last elements
      m_intSolPntSet      .pop_back();
      m_ghostSolPntSet    .pop_back();
      m_intSolPntSetPtrs  .pop_back();
      m_ghostSolPntSetPtrs.pop_back();
      m_intGradsPntSet    .pop_back();
      m_ghostGradsPntSet  .pop_back();
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void BaseBndFaceTermComputer::computeFacePntSetData()
{
  // compute face Jacobian vectors in the points set
  m_faceJacobPntSet = m_face->
      computeFaceJacobDetVectorAtMappedCoords(m_faceMappedCoordsPntSet);

  // set face Jacobian vector sizes and unit normals in the flux points
  const CFuint nbrPnts = m_faceJacobPntSet.size();
  cf_assert(nbrPnts == m_unitNormalPntSet.size());
  for (CFuint iPnt = 0; iPnt < nbrPnts; ++iPnt)
  {
    // set unit normal vector
    m_unitNormalPntSet[iPnt] = m_faceJacobPntSet[iPnt]/m_faceJacobPntSet[iPnt].norm2();
  }

  // set the face ID in the BCStateComputer
  m_bcStateComputer->setFaceID(m_face->getID());

  // if necessary, compute the flux point coordinates
  // the output point coordinates are stored in m_faceFlxPntsFaceLocalCoords
  // thus the coordinates are computed correctly
  if (m_bcStateComputer->needsSpatialCoordinates())
  {
    for (CFuint iPnt = 0; iPnt < nbrPnts; ++iPnt)
    {
      m_coordsPntSet[iPnt] =
          m_face->computeCoordFromMappedCoord(m_faceMappedCoordsPntSet[iPnt]);
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

vector< RealVector >& BaseBndFaceTermComputer::reconstructGivenPntsStates(const std::vector< Framework::State* >& cellIntStates)
{
  // number of output points
  const CFuint nbrPnts = m_solRecCoefsPntSet.size();
  cf_assert(nbrPnts == m_statesPntSet.size());
//   cf_assert(nbrPnts == m_intSolPntSetPtrs.size());
//   cf_assert(nbrPnts == m_ghostSolPntSetPtrs.size());

  // number of solution points/basis polynomials
  const CFuint nbrPolys = cellIntStates.size();
  for (CFuint iPnt = 0; iPnt < nbrPnts; ++iPnt)
  {
    *m_intSolPntSet[iPnt] = 0.0;
    cf_assert(nbrPolys == m_solRecCoefsPntSet[iPnt].size());
    for (CFuint iPoly = 0; iPoly < nbrPolys; ++iPoly)
    {
      *m_intSolPntSet[iPnt] += m_solRecCoefsPntSet[iPnt][iPoly]*(*cellIntStates[iPoly]);
    }
  }

  // compute ghost states
  m_bcStateComputer->computeGhostStates(m_intSolPntSet,m_ghostSolPntSet,m_unitNormalPntSet,m_coordsPntSet);

  // average internal and ghost states
  for (CFuint iPnt = 0; iPnt < nbrPnts; ++iPnt)
  {
    m_statesPntSet[iPnt] = 0.5*(*m_intSolPntSet[iPnt] + *m_ghostSolPntSet[iPnt]);
  }

  return m_statesPntSet;
}

//////////////////////////////////////////////////////////////////////////////

vector< vector< RealVector > >&
    BaseBndFaceTermComputer::reconstructGivenPntsGrads(const vector< vector< RealVector >* >& cellIntGrads)
{
  // number of points
  const CFuint nbrPnts = m_solRecCoefsPntSet.size();
  cf_assert(nbrPnts == m_gradsPntSet.size());

  // number of solution points/basis polynomials
  const CFuint nbrPolys = cellIntGrads.size();
  for (CFuint iPnt = 0; iPnt < nbrPnts; ++iPnt)
  {
    for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
    {
      // dereference gradient
      RealVector& grad = *m_intGradsPntSet[iPnt][iEq];

      // reconstruct gradients
      grad = 0.0;
      cf_assert(nbrPolys == m_solRecCoefsPntSet[iPnt].size());
      for (CFuint iPoly = 0; iPoly < nbrPolys; ++iPoly)
      {
        grad += m_solRecCoefsPntSet[iPnt][iPoly]*(*cellIntGrads[iPoly])[iEq];
      }
    }
  }

  // compute ghost gradients
  m_bcStateComputer->computeGhostGradients(m_intGradsPntSet,m_ghostGradsPntSet,m_unitNormalPntSet,m_coordsPntSet);

  // average internal and ghost gradients
  for (CFuint iPnt = 0; iPnt < nbrPnts; ++iPnt)
  {
    for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
    {
      m_gradsPntSet[iPnt][iEq] = 0.5*(*m_intGradsPntSet[iPnt][iEq] + *m_ghostGradsPntSet[iPnt][iEq]);
    }
  }

  return m_gradsPntSet;
}

//////////////////////////////////////////////////////////////////////////////

void BaseBndFaceTermComputer::computeFaceTermFromFlxPntFluxes(RealVector& resUpdates)
{
  // number of solution points one flux point contributes to
  cf_assert(m_solPntIdxsForDerivation->size() > 0);
  const CFuint nbrSolsPerFlx = (*m_solPntIdxsForDerivation)[0].size();

  // set updates to zero
  resUpdates = 0.0;

  // loop over flux points
  for (CFuint iFlx = 0; iFlx < m_nbrFlxPnts[m_orient]; ++iFlx)
  {
    // flux point index
    const CFuint flxIdx = (*m_faceFlxPntConn)[m_orient][iFlx];

    // flux point index in the matrix m_solPntsDerivCoefs
    const CFuint flxPntMatrixIdx = (*m_flxPntMatrixIdxForDerivation)[flxIdx];

    // evaluate flux projected on the projection vector
    const RealVector fluxXProjVect = m_flxPntRiemannFlux[iFlx]*
                                     m_faceJacobVecSizeFlxPnts[iFlx];

    // add contribution of this flux point to the solution points
    for (CFuint iSol = 0; iSol < nbrSolsPerFlx; ++iSol)
    {
      // first residual in solution point index
      CFuint resIdx = m_nbrEqs*(*m_solPntIdxsForDerivation)[flxIdx][iSol];
      for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar, ++resIdx)
      {
        resUpdates[resIdx] +=
            (*m_solPntsDerivCoefs)(iSol,flxPntMatrixIdx)*fluxXProjVect[iVar];
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void BaseBndFaceTermComputer::computeGradFaceTermFromFlxPntSol(vector< vector< RealVector > >& gradUpdates)
{
  // number of solution points one flux point contributes to
  cf_assert(m_solPntIdxsForDerivation->size() > 0);
  const CFuint nbrSolsPerFlx = (*m_solPntIdxsForDerivation)[0].size();

  // set updates to zero
  const CFuint nbrSolPnts = gradUpdates.size();
  for (CFuint iSol = 0; iSol < nbrSolPnts; ++iSol)
  {
    for (CFuint iGrad = 0; iGrad < m_nbrEqs; ++iGrad)
    {
      gradUpdates[iSol][iGrad] = 0.0;
    }
  }

  // loop over flux points
  for (CFuint iFlx = 0; iFlx < m_nbrFlxPnts[m_orient]; ++iFlx)
  {
    // flux point index
    const CFuint flxIdx = (*m_faceFlxPntConn)[m_orient][iFlx];

    // flux point index in the matrix m_solPntsDerivCoefs
    const CFuint flxPntMatrixIdx = (*m_flxPntMatrixIdxForDerivation)[flxIdx];

    for (CFuint iGrad = 0; iGrad < m_nbrEqs; ++iGrad)
    {
      // compute gradient term
      m_gradTerm = m_flxPntRiemannFlux[iFlx][iGrad]*
                   m_faceJacobVecSizeFlxPnts[iFlx]*
                   m_unitNormalFlxPnts[iFlx];

      // add contribution of this flux point to the solution points
      for (CFuint iSol = 0; iSol < nbrSolsPerFlx; ++iSol)
      {
        // solution point index
        const CFuint solIdx = (*m_solPntIdxsForDerivation)[flxIdx][iSol];
        gradUpdates[solIdx][iGrad] +=
            (*m_solPntsDerivCoefs)(iSol,flxPntMatrixIdx)*m_gradTerm;
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void BaseBndFaceTermComputer::setup()
{
  CFAUTOTRACE;

  // m_dimensionality and number of equations
  m_dim    = PhysicalModelStack::getActive()->getDim();
  m_nbrEqs = PhysicalModelStack::getActive()->getNbEq();

  // get the update varset
  m_updateVarSet = getMethodData().getUpdateVar();

  // get number of extra variables from the update variable set
  m_nbrExtraVars = m_updateVarSet->getExtraPhysicalVarsSize();
  if (m_nbrExtraVars > 0 && !socket_extraVars.isConnected())
  {
    throw ConsistencyException(FromHere(), "Extra variables socket not connected but required...");
  }

  // get the states reconstructor
  m_statesReconstr = getMethodData().getStatesReconstructor();

  // get the Riemann flux
  m_riemannFluxComputer = getMethodData().getRiemannFlux();

  // get the face diffusive flux
  m_faceDiffFluxComputer = getMethodData().getFaceDiffusiveFlux();

  // get the local spectral FD data
  vector< SpectralFDElementData* >& sdLocalData = getMethodData().getSDLocalData();
  cf_assert(sdLocalData.size() > 0);

  // number of solution points
  const CFuint nbrSolPnts = sdLocalData[0]->getNbrOfSolPnts();

  // resize m_cellExtraVars
  m_cellExtraVars.resize(nbrSolPnts);

  // number of flux points
  const CFuint nbrFlxPnts = sdLocalData[0]->getNbrOfFaceFlxPnts();

  m_flxPntCoords.resize(nbrFlxPnts);
  for (CFuint iFlx = 0; iFlx < nbrFlxPnts; ++iFlx)
  {
    m_flxPntCoords[iFlx].resize(m_dim);
  }

  // create internal and ghost states and extra variables
  for (CFuint iFlx = 0; iFlx < nbrFlxPnts; ++iFlx)
  {
    m_flxPntIntSol        .push_back(new State()                   );
    m_flxPntGhostSol      .push_back(new State()                   );
    m_flxPntIntExtraVars  .push_back(new RealVector(m_nbrExtraVars));
    m_flxPntGhostExtraVars.push_back(new RealVector(m_nbrExtraVars));
  }

  for (CFuint iFlx = 0; iFlx < nbrFlxPnts; ++iFlx)
  {
  m_flxPntIntSol[iFlx]->setLocalID(iFlx);
  m_flxPntGhostSol[iFlx]->setLocalID(iFlx);
  }

  // resize m_backupPhysVar
  m_backupPhysVar.resize(nbrFlxPnts);

  // resize m_faceJacobVecAbsSizeFlxPnts
  m_faceJacobVecAbsSizeFlxPnts.resize(nbrFlxPnts);

  // resize m_faceJacobVecSizeFlxPnts
  m_faceJacobVecSizeFlxPnts.resize(nbrFlxPnts);

  // resize m_faceInvCharLengths
  m_faceInvCharLengths.resize(nbrFlxPnts);

  // resize m_unitNormalFlxPnts
  m_unitNormalFlxPnts.resize(nbrFlxPnts,RealVector(m_dim));

  // resize m_flxPntRiemannFlux
  m_flxPntRiemannFlux.resize(nbrFlxPnts,RealVector(m_nbrEqs));

  // resize m_backupGhostStates
  m_backupGhostStates.resize(nbrFlxPnts);
  for (CFuint iFlx = 0; iFlx < nbrFlxPnts; ++iFlx)
  {
    m_backupGhostStates[iFlx].resize(m_nbrEqs);
  }

  // set m_allSol and m_allExtraVars pointers
  for (CFuint iFlx = 0; iFlx < nbrFlxPnts; ++iFlx)
  {
    m_allSol      .push_back(m_flxPntIntSol        [iFlx]);
    m_allSol      .push_back(m_flxPntGhostSol      [iFlx]);
    m_allExtraVars.push_back(m_flxPntIntExtraVars  [iFlx]);
    m_allExtraVars.push_back(m_flxPntGhostExtraVars[iFlx]);
  }

  // set RealVector pointers to states
  for (CFuint iFlx = 0; iFlx < nbrFlxPnts; ++iFlx)
  {
    m_flxPntIntRVSol  .push_back(m_flxPntIntSol  [iFlx]);
    m_flxPntGhostRVSol.push_back(m_flxPntGhostSol[iFlx]);
  }

  // setup variables for gradient computation
  if (getMethodData().hasDiffTerm())
  {
    // get the diffusive varset
    m_diffusiveVarSet = getMethodData().getDiffusiveVar();
  }

  // create gradients for flux points
  m_flxPntIntGrads  .resize(nbrFlxPnts);
  m_flxPntGhostGrads.resize(nbrFlxPnts);
  for (CFuint iFlx = 0; iFlx < nbrFlxPnts; ++iFlx)
  {
    for (CFuint iGrad = 0; iGrad < m_nbrEqs; ++iGrad)
    {
      m_flxPntIntGrads  [iFlx].push_back(new RealVector(m_dim));
      m_flxPntGhostGrads[iFlx].push_back(new RealVector(m_dim));
    }
  }

  // create pointers to gradients for quadrature points
  for (CFuint iFlx = 0; iFlx < nbrFlxPnts; ++iFlx)
  {
    m_flxPntIntGradPtrs  .push_back(&m_flxPntIntGrads  [iFlx]);
    m_flxPntGhostGradPtrs.push_back(&m_flxPntGhostGrads[iFlx]);
  }

  // resize m_gradTerm
  m_gradTerm.resize(m_dim);

  // set maximum number of flux points that will be passed at one time to the physical model
  const CFuint maxNbrFlx = 2*nbrFlxPnts;
  const CFuint prevMaxNbrStatesData = getMethodData().getMaxNbrStatesData();
  getMethodData().
      setMaxNbrStatesData(prevMaxNbrStatesData > maxNbrFlx ? prevMaxNbrStatesData : maxNbrFlx);

  // set maximum number of points in which the Riemann flux has to be evaluated at the same time
  const CFuint prevMaxNbrRFluxPnts = getMethodData().getMaxNbrRFluxPnts();
  getMethodData().
      setMaxNbrRFluxPnts
      (
       prevMaxNbrRFluxPnts > nbrFlxPnts ? prevMaxNbrRFluxPnts : nbrFlxPnts
      );
  

  // resize the physical data temporary vector
  SafePtr<BaseTerm> convTerm = PhysicalModelStack::getActive()->getImplementor()->getConvectiveTerm(); 
  convTerm->resizePhysicalData(m_pData);
  
}

//////////////////////////////////////////////////////////////////////////////

void BaseBndFaceTermComputer::unsetup()
{
  CFAUTOTRACE;

  for (CFuint iFlx = 0; iFlx < m_flxPntIntSol.size(); ++iFlx)
  {
    deletePtr(m_flxPntIntSol  [iFlx]);
    deletePtr(m_flxPntGhostSol[iFlx]);
    deletePtr(m_flxPntIntExtraVars  [iFlx]);
    deletePtr(m_flxPntGhostExtraVars[iFlx]);
  }
  m_flxPntIntSol  .resize(0);
  m_flxPntGhostSol.resize(0);
  m_flxPntIntExtraVars  .resize(0);
  m_flxPntGhostExtraVars.resize(0);

  for (CFuint iFlx = 0; iFlx < m_flxPntIntGrads.size(); ++iFlx)
  {
    for (CFuint iGrad = 0; iGrad < m_flxPntIntGrads[iFlx].size(); ++iGrad)
    {
      deletePtr(m_flxPntIntGrads  [iFlx][iGrad]);
      deletePtr(m_flxPntGhostGrads[iFlx][iGrad]);
    }
    m_flxPntIntGrads  [iFlx].resize(0);
    m_flxPntGhostGrads[iFlx].resize(0);
  }
  m_flxPntIntGrads  .resize(0);
  m_flxPntGhostGrads.resize(0);

  for (CFuint iPnt = 0; iPnt < m_intSolPntSet.size(); ++iPnt)
  {
    deletePtr(m_intSolPntSet  [iPnt]);
    deletePtr(m_ghostSolPntSet[iPnt]);
  }
  m_intSolPntSet  .resize(0);
  m_ghostSolPntSet.resize(0);
}

//////////////////////////////////////////////////////////////////////////////

vector<SafePtr<BaseDataSocketSink> >
    BaseBndFaceTermComputer::needsSockets()
{
  vector<SafePtr<BaseDataSocketSink> > result;

  result.push_back(&socket_extraVars                  );
  result.push_back(&socket_faceJacobVecSizeFaceFlxPnts);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

  }  // namespace SpectralFD

}  // namespace COOLFluiD
