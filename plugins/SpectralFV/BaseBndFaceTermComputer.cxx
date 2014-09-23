#include "Framework/MethodStrategyProvider.hh"

#include "SpectralFV/SpectralFV.hh"
#include "SpectralFV/BaseBndFaceTermComputer.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace SpectralFV {

//////////////////////////////////////////////////////////////////////////////

BaseBndFaceTermComputer::BaseBndFaceTermComputer(const std::string& name) :
  SpectralFVMethodStrategy(name),
  socket_faceSurf("faceSurf"),
  m_updateVarSet(CFNULL),
  m_diffusiveVarSet(CFNULL),
  m_statesReconstr(CFNULL),
  m_bcStateComputer(CFNULL),
  m_riemannFluxComputer(CFNULL),
  m_faceDiffFluxComputer(CFNULL),
  m_cflConvDiffRatio(),
  m_svFaceflxPntsRecCoefs(CFNULL),
  m_flxPntWheightCoordsSVFaces(CFNULL),
  m_face(),
  m_orient(),
  m_cellVolume(),
  m_surf(),
  m_avgNormal(),
  m_avgUnitNormal(),
  m_unitNormalFlxPnt(),
  m_flxPntIntSol(),
  m_flxPntGhostSol(),
  m_flxPntIntExtraVars(),
  m_flxPntGhostExtraVars(),
  m_flxPntCoords(),
  m_faceAvgSolInt(CFNULL),
  m_faceAvgExtraVarInt(CFNULL),
  m_allSol(),
  m_allExtraVars(),
  m_flxPntIntRVSol(),
  m_flxPntGhostRVSol(),
  m_flxPntRiemannFlx(),
  m_flxPntIntGradVars(),
  m_flxPntGhostGradVars(),
  m_flxPntIntGrads(),
  m_flxPntGhostGrads(),
  m_flxPntIntGradPtrs(),
  m_flxPntGhostGradPtrs(),
  m_backupPhysVar(),
  m_backupGhostStates(),
  m_nbrFlxPnts(),
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
//  m_faceJacobPntSet(),
  m_coordsPntSet(),
  m_gradsPntSet()
//  m_intGradsPntSet(),
//  m_ghostGradsPntSet()
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

BaseBndFaceTermComputer::~BaseBndFaceTermComputer()
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

void BaseBndFaceTermComputer::setFaceTermData()
{
  CFAUTOTRACE;

  // get the SpectralFVElementData
  vector< SpectralFVElementData* >& svLocalData = getMethodData().getSVLocalData();

  // get convective/diffusive CFL ratio
  m_cflConvDiffRatio = svLocalData[0]->getCFLConvDiffRatio();
}

//////////////////////////////////////////////////////////////////////////////

void BaseBndFaceTermComputer::computeFaceData()
{
  // compute the face surface
  DataHandle< CFreal > faceSurf = socket_faceSurf.getDataHandle();
  m_surf = faceSurf[m_face->getID()];

  // compute the face average normal vector
  m_avgNormal = m_face->computeAvgCellNormal();

  // compute the face average unit normal vector
  m_avgUnitNormal = m_avgNormal/m_surf;

  // set the face ID in the BCStateComputer
  m_bcStateComputer->setFaceID(m_face->getID());
}

//////////////////////////////////////////////////////////////////////////////

void BaseBndFaceTermComputer::computeNeighbourCellData()
{
  GeometricEntity* cell = m_face->getNeighborGeo(0);
  m_cellVolume = cell->computeVolume();
}

//////////////////////////////////////////////////////////////////////////////

void BaseBndFaceTermComputer::reconstructFluxPntsStates(const vector< State* >& cellIntStates)
{
  // compute internal states
  m_statesReconstr->reconstructStates(cellIntStates,m_flxPntIntSol,
                                      (*m_svFaceflxPntsRecCoefs)[m_orient],
                                      cellIntStates.size());

  // compute ghost states
  computeGhostStates();
}

//////////////////////////////////////////////////////////////////////////////

void BaseBndFaceTermComputer::reconstructFluxPntsGradients(const vector< vector< RealVector >* >& cellIntGrads,
                                                           const CFuint nbrCellGrads)
{
  // reconstruct internal gradients
  m_statesReconstr->reconstructGradients(cellIntGrads,m_flxPntIntGrads,
                                         (*m_svFaceflxPntsRecCoefs)[m_orient],
                                           nbrCellGrads);

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
      ->reconstructPhysVar(iVar,cellIntStates,m_flxPntIntSol,
                           (*m_svFaceflxPntsRecCoefs)[m_orient],
                             cellIntStates.size());

  // compute ghost states
  computePerturbedGhostStates();
}

//////////////////////////////////////////////////////////////////////////////

void BaseBndFaceTermComputer::backupPhysVar(const CFuint iVar)
{
  const CFuint nbrFlxPnts = m_flxPntIntSol.size();
  cf_assert(nbrFlxPnts <= m_backupPhysVar.size());
  cf_assert(nbrFlxPnts <= m_backupGhostStates.size());
  for (CFuint iFlx = 0; iFlx < nbrFlxPnts; ++iFlx)
  {
    m_backupPhysVar[iFlx] = (*m_flxPntIntSol[iFlx])[iVar];
    m_backupGhostStates[iFlx] = *m_flxPntGhostSol[iFlx];
  }
}

//////////////////////////////////////////////////////////////////////////////

void BaseBndFaceTermComputer::restorePhysVar(const CFuint iVar)
{
  const CFuint nbrFlxPnts = m_flxPntIntSol.size();
  cf_assert(nbrFlxPnts <= m_backupPhysVar.size());
  cf_assert(nbrFlxPnts <= m_backupGhostStates.size());
  for (CFuint iFlx = 0; iFlx < nbrFlxPnts; ++iFlx)
  {
    (*m_flxPntIntSol[iFlx])[iVar] = m_backupPhysVar[iFlx];
    *m_flxPntGhostSol[iFlx] = m_backupGhostStates[iFlx];
  }
}

//////////////////////////////////////////////////////////////////////////////

void BaseBndFaceTermComputer::computeGhostStates()
{
  cf_assert(m_bcStateComputer.isNotNull());

  // if necessary, compute the flux point coordinates
  if (m_bcStateComputer->needsSpatialCoordinates())
  {
    computeFlxPntCoords();
  }

  // compute ghost states
  const CFuint nbrFlxPnts = m_unitNormalFlxPnt.size();
  for (CFuint iFlx = 0; iFlx < nbrFlxPnts; ++iFlx)
  {
    m_unitNormalFlxPnt[iFlx] = m_avgUnitNormal;
  }
  m_bcStateComputer->computeGhostStates(m_flxPntIntSol,m_flxPntGhostSol,m_unitNormalFlxPnt,m_flxPntCoords);
}

//////////////////////////////////////////////////////////////////////////////

void BaseBndFaceTermComputer::computePerturbedGhostStates()
{
  // the flux point coordinates should have been computed before (unperturbed)
  // the unit normals should have been computed before (unperturbed)

  // compute ghost states
  cf_assert(m_bcStateComputer.isNotNull());
  m_bcStateComputer->computeGhostStates(m_flxPntIntSol,m_flxPntGhostSol,m_unitNormalFlxPnt,m_flxPntCoords);
}

//////////////////////////////////////////////////////////////////////////////

void BaseBndFaceTermComputer::computeGhostGradients()
{
  cf_assert(m_bcStateComputer.isNotNull());

  // if necessary, compute the flux point coordinates
  if (m_bcStateComputer->needsSpatialCoordinates())
  {
    computeFlxPntCoords();
  }

  // compute ghost gradients
  const CFuint nbrFlxPnts = m_unitNormalFlxPnt.size();
  for (CFuint iFlx = 0; iFlx < nbrFlxPnts; ++iFlx)
  {
    m_unitNormalFlxPnt[iFlx] = m_avgUnitNormal;
  }
  m_bcStateComputer->
      computeGhostGradients(m_flxPntIntGrads,m_flxPntGhostGrads,m_unitNormalFlxPnt,m_flxPntCoords);
}

//////////////////////////////////////////////////////////////////////////////

void BaseBndFaceTermComputer::computeFlxPntCoords()
{
  // number of flux points
  const CFuint nbrFlxPnts = m_flxPntCoords.size();

  // get the face nodes
  vector< Node* >* nodes = m_face->getNodes();

  // number of nodes in the face
  const CFuint nbrFaceNodes = nodes->size();

  // compute flux point coordinates
  for (CFuint iFlx = 0; iFlx < nbrFlxPnts; ++iFlx)
  {
    cf_assert(nbrFaceNodes == (*m_flxPntWheightCoordsSVFaces)[iFlx].size());

    // dereference flux point coordinates
    RealVector& bCoord = m_flxPntCoords[iFlx];

    bCoord = 0.0;
    for (CFuint iNode = 0; iNode < nbrFaceNodes; ++iNode)
    {
      bCoord += (*m_flxPntWheightCoordsSVFaces)[iFlx][iNode]*(*(*nodes)[iNode]);
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void BaseBndFaceTermComputer::computeConvFaceTerm(RealVector& resUpdates)
{
  // compute the states data in the flux points (and not for the face averaged states --> -1)
 /// @todo broken after release 2009.3
//   m_updateVarSet->computeStatesData(m_allSol,m_allExtraVars,m_allSol.size()-1);

  // evaluate the riemann fluxes in all the flux points
  m_flxPntRiemannFlx = m_riemannFluxComputer->computeFlux(m_flxPntIntSol,m_flxPntGhostSol,m_unitNormalFlxPnt,m_nbrFlxPnts[m_orient]);

  // compute the actual face term
  computeFaceFluxIntegralFromFlxPntFluxes(resUpdates);
}

//////////////////////////////////////////////////////////////////////////////

void BaseBndFaceTermComputer::computeConvFaceTermAndWaveSpeedUpdates(RealVector& resUpdates,
                                                                     CFreal& waveSpeedUpd)
{
  // compute the states data in the flux points and for the face averaged states
 /// @todo broken after release 2009.3
//   m_updateVarSet->computeStatesData(m_allSol,m_allExtraVars,m_allSol.size());

  // evaluate the riemann fluxes in all the flux points
  m_flxPntRiemannFlx = m_riemannFluxComputer->computeFlux(m_flxPntIntSol,m_flxPntGhostSol,m_unitNormalFlxPnt,m_nbrFlxPnts[m_orient]);

  // compute the actual face term
  computeFaceFluxIntegralFromFlxPntFluxes(resUpdates);

  // compute the wave speed updates for the neighbouring cells
  waveSpeedUpd = m_updateVarSet->getMaxAbsEigenValue(*m_faceAvgSolInt,m_avgUnitNormal)*m_surf;
}

//////////////////////////////////////////////////////////////////////////////

void BaseBndFaceTermComputer::computeGradientFaceTerm(vector< vector< RealVector > >& gradUpdates)
{
  // compute the averaged gradient variables in the flux points (stored in m_flxPntRiemannFlx)
  // (m_leftRVSol and m_rightRVSol are RealVector pointers to the left and right states)
  m_flxPntRiemannFlx = m_faceDiffFluxComputer->computeAvgGradVars(m_flxPntIntRVSol,m_flxPntGhostRVSol,
                                                                  m_nbrFlxPnts[m_orient]);

  // compute the actual volume term contribution to the gradient
  computeGradFaceTermFromFlxPntSol(gradUpdates);
}

//////////////////////////////////////////////////////////////////////////////

void BaseBndFaceTermComputer::computeDiffFaceTerm(RealVector& resUpdates)
{
  // compute the diffusive fluxes in the face flux points
  m_flxPntRiemannFlx = m_faceDiffFluxComputer->computeDiffFlux(m_flxPntIntGradPtrs,m_flxPntGhostGradPtrs,
                                                               m_flxPntIntRVSol,m_flxPntGhostRVSol,
                                                               m_cellVolume,m_cellVolume,
                                                               m_surf,m_unitNormalFlxPnt,
                                                               m_nbrFlxPnts[m_orient]);

  // compute the actual face term
  computeFaceFluxIntegralFromFlxPntFluxes(resUpdates);
}

//////////////////////////////////////////////////////////////////////////////

void BaseBndFaceTermComputer::computeDiffFaceTermAndUpdateCoefContributions(RealVector& resUpdates,
                                                                            CFreal& updateCoefContr)
{
  // compute the diffusive fluxes in the face flux points
  m_flxPntRiemannFlx = m_faceDiffFluxComputer->computeDiffFlux(m_flxPntIntGradPtrs,m_flxPntGhostGradPtrs,
                                                               m_flxPntIntRVSol,m_flxPntGhostRVSol,
                                                               m_cellVolume,m_cellVolume,
                                                               m_surf,m_unitNormalFlxPnt,
                                                               m_nbrFlxPnts[m_orient]);

  // compute the actual face term
  computeFaceFluxIntegralFromFlxPntFluxes(resUpdates);

  // compute the update coefficient contributions for the neighbouring cells
  CFLog(WARN,"Contribution of the diffusive terms to the update coefficients is not added...");
  updateCoefContr = 0.0;
}

//////////////////////////////////////////////////////////////////////////////

void BaseBndFaceTermComputer::setPointSet(const std::vector< RealVector >& faceMappedCoordPntSet)
{
  // get the SpectralFVElementData
  vector< SpectralFVElementData* >& svLocalData = getMethodData().getSVLocalData();

  // get face node cell-mapped coordinates
  SafePtr< vector< vector< RealVector > > > faceNodeCoords = svLocalData[0]->getSVFaceNodeCoords();

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
      throw Common::NotImplementedException(FromHere(),"Output on boundary faces not implemented for 3D SV");
//      cf_assert(faceMappedCoord.size() == 2);
//      nodeWheights.push_back(0.25*(1.0-faceMappedCoord[KSI])*(1.0-faceMappedCoord[ETA]));
//      nodeWheights.push_back(0.25*(1.0+faceMappedCoord[KSI])*(1.0-faceMappedCoord[ETA]));
//      nodeWheights.push_back(0.25*(1.0+faceMappedCoord[KSI])*(1.0+faceMappedCoord[ETA]));
//      nodeWheights.push_back(0.25*(1.0-faceMappedCoord[KSI])*(1.0+faceMappedCoord[ETA]));
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
  m_solRecCoefsPntSet = svLocalData[0]->getSVPolyValsAtNode(m_cellMappedCoordsPntSet);

  // resize some variables
  if (nbrPnts > m_statesPntSet.size())
  {
    m_statesPntSet.resize(nbrPnts,RealVector(m_nbrEqs));
    m_unitNormalPntSet.resize(nbrPnts,RealVector(m_dim));
//    m_faceJacobPntSet .resize(nbrPnts,RealVector(m_dim));
    m_coordsPntSet    .resize(nbrPnts,RealVector(m_dim));
//    m_gradsPntSet.resize(nbrPnts,vector<RealVector>(m_nbrEqs,RealVector(m_dim)));

    CFuint iPnt = m_intSolPntSet.size();
    for (; iPnt < nbrPnts; ++iPnt)
    {
      m_intSolPntSet  .push_back(new State());
      m_ghostSolPntSet.push_back(new State());

      // set pointers to internal and ghost states
      m_intSolPntSetPtrs  .push_back(m_intSolPntSet  [iPnt]);
      m_ghostSolPntSetPtrs.push_back(m_ghostSolPntSet[iPnt]);

//      vector<RealVector*> intGrads  ;
//      vector<RealVector*> ghostGrads;
//      for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
//      {
//        intGrads  .push_back(new RealVector(m_dim));
//        ghostGrads.push_back(new RealVector(m_dim));
//      }
//      m_intGradsPntSet  .push_back(intGrads  );
//      m_ghostGradsPntSet.push_back(ghostGrads);
    }
  }
  else if (nbrPnts < m_statesPntSet.size())
  {
    m_statesPntSet.resize(nbrPnts,RealVector(m_nbrEqs));
    m_unitNormalPntSet.resize(nbrPnts,RealVector(m_dim));
//    m_faceJacobPntSet .resize(nbrPnts,RealVector(m_dim));
    m_coordsPntSet    .resize(nbrPnts,RealVector(m_dim));
//    m_gradsPntSet.resize(nbrPnts,vector<RealVector>(m_nbrEqs,RealVector(m_dim)));

    CFuint iPnt = m_intSolPntSet.size()-1;
    for (; iPnt >= nbrPnts; --iPnt)
    {
      deletePtr(m_intSolPntSet  [iPnt]);
      deletePtr(m_ghostSolPntSet[iPnt]);

//      for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
//      {
//        deletePtr(m_intGradsPntSet  [iPnt][iEq]);
//        deletePtr(m_ghostGradsPntSet[iPnt][iEq]);
//      }

      // remove last elements
      m_intSolPntSet      .pop_back();
      m_ghostSolPntSet    .pop_back();
      m_intSolPntSetPtrs  .pop_back();
      m_ghostSolPntSetPtrs.pop_back();
//      m_intGradsPntSet    .pop_back();
//      m_ghostGradsPntSet  .pop_back();
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void BaseBndFaceTermComputer::computeFacePntSetData()
{
  // compute face Jacobian vectors in the points set
//  m_faceJacobPntSet = m_face->
//      computeFaceJacobDetVectorAtMappedCoords(m_faceMappedCoordsPntSet);

  // compute the face surface
  DataHandle< CFreal > faceSurf = socket_faceSurf.getDataHandle();
  m_surf = faceSurf[m_face->getID()];

  // compute the face average normal vector
  m_avgNormal = m_face->computeAvgCellNormal();

  // compute the face average unit normal vector
  m_avgUnitNormal = m_avgNormal/m_surf;

  // set unit normals in the flux points
  const CFuint nbrPnts = m_faceMappedCoordsPntSet.size();
  cf_assert(nbrPnts == m_unitNormalPntSet.size());
  for (CFuint iPnt = 0; iPnt < nbrPnts; ++iPnt)
  {
    // set unit normal vector
    m_unitNormalPntSet[iPnt] = m_avgUnitNormal;
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
/*  // number of points
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
*/
  return m_gradsPntSet;
}

//////////////////////////////////////////////////////////////////////////////

void BaseBndFaceTermComputer::setup()
{
  CFAUTOTRACE;

  // dimensionality and number of equations
  const CFuint dim   = PhysicalModelStack::getActive()->getDim();
  m_nbrEqs = PhysicalModelStack::getActive()->getNbEq();

  // get the update varset
  m_updateVarSet = getMethodData().getUpdateVar();

  // get number of extra variables from the update variable set
  m_nbrExtraVars = m_updateVarSet->getExtraPhysicalVarsSize();

  // get the diffusive varset
  if (getMethodData().hasDiffTerm())
  {
    m_diffusiveVarSet = getMethodData().getDiffusiveVar();
  }

  // get the states reconstructor
  m_statesReconstr = getMethodData().getStatesReconstructor();

  // get the Riemann flux
  m_riemannFluxComputer = getMethodData().getRiemannFlux();

  // get the face diffusive flux
  m_faceDiffFluxComputer = getMethodData().getFaceDiffusiveFlux();

  // create state and extra variables for face averaged solution
  m_faceAvgSolInt      = new State();
  m_faceAvgExtraVarInt = new RealVector(m_nbrExtraVars);

  // resize variables
  m_avgNormal.resize(dim);
  m_avgUnitNormal.resize(dim);
}

//////////////////////////////////////////////////////////////////////////////

void BaseBndFaceTermComputer::unsetup()
{
  CFAUTOTRACE;

  deletePtr(m_faceAvgSolInt     );
  deletePtr(m_faceAvgExtraVarInt);
}

//////////////////////////////////////////////////////////////////////////////

vector<SafePtr<BaseDataSocketSink> >
    BaseBndFaceTermComputer::needsSockets()
{
  vector<SafePtr<BaseDataSocketSink> > result;

  result.push_back(&socket_faceSurf);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

  }  // namespace SpectralFV

}  // namespace COOLFluiD
