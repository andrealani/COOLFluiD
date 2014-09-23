#include "Framework/MethodStrategyProvider.hh"

#include "SpectralFV/SpectralFV.hh"
#include "SpectralFV/BaseFaceTermComputer.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace SpectralFV {

//////////////////////////////////////////////////////////////////////////////

BaseFaceTermComputer::BaseFaceTermComputer(const std::string& name) :
  SpectralFVMethodStrategy(name),
  socket_faceSurf("faceSurf"),
  m_updateVarSet(CFNULL),
  m_diffusiveVarSet(CFNULL),
  m_statesReconstr(CFNULL),
  m_riemannFluxComputer(CFNULL),
  m_faceDiffFluxComputer(CFNULL),
  m_cflConvDiffRatio(),
  m_svFaceflxPntsRecCoefs(CFNULL),
  m_orient(),
  m_surf(),
  m_cellVolumes(),
  m_avgNormal(),
  m_avgUnitNormal(),
  m_unitNormalFlxPnt(),
  m_flxPntSol(),
  m_flxPntExtraVars(),
  m_faceAvgSol(),
  m_faceAvgExtraVars(),
  m_allSol(),
  m_allExtraVars(),
  m_flxPntRVSol(),
  m_flxPntRiemannFlx(),
  m_flxPntGradVars(),
  m_flxPntGrads(),
  m_flxPntGradPtrs(),
  m_backupPhysVar(),
  m_nbrFlxPnts(),
  m_nbrEqs(),
  m_nbrExtraVars()
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

BaseFaceTermComputer::~BaseFaceTermComputer()
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

void BaseFaceTermComputer::setFaceTermData()
{
  CFAUTOTRACE;

  // get the SpectralFVElementData
  vector< SpectralFVElementData* >& svLocalData = getMethodData().getSVLocalData();

  // get convective/diffusive CFL ratio
  m_cflConvDiffRatio = svLocalData[0]->getCFLConvDiffRatio();
}

//////////////////////////////////////////////////////////////////////////////

void BaseFaceTermComputer::computeFaceData()
{
  // compute the face surface
  DataHandle< CFreal > faceSurf = socket_faceSurf.getDataHandle();
  m_surf = faceSurf[m_face->getID()];

  // compute the face average normal vector
  m_avgNormal = m_face->computeAvgCellNormal();

  // compute the face average unit normal vector
  m_avgUnitNormal = m_avgNormal/m_surf;
}

//////////////////////////////////////////////////////////////////////////////

void BaseFaceTermComputer::computeNeighbourCellData()
{
  for (CFuint iSide = 0; iSide < 2; ++iSide)
  {
    GeometricEntity* cell = m_face->getNeighborGeo(iSide);
    m_cellVolumes[iSide] = cell->computeVolume();
  }
}

//////////////////////////////////////////////////////////////////////////////

void BaseFaceTermComputer::reconstructFluxPntsStates(const vector< State* >& cellLStates,
                                                     const vector< State* >& cellRStates)
{
  m_statesReconstr->reconstructStates(cellLStates,m_flxPntSol[LEFT ],
                                      (*m_svFaceflxPntsRecCoefs)[m_orient][LEFT ],
                                        cellLStates.size());
  m_statesReconstr->reconstructStates(cellRStates,m_flxPntSol[RIGHT],
                                      (*m_svFaceflxPntsRecCoefs)[m_orient][RIGHT],
                                        cellRStates.size());
}

//////////////////////////////////////////////////////////////////////////////

void BaseFaceTermComputer::reconstructFluxPntsGradients(const vector< vector< RealVector >* >& cellLGrads,
                                                        const vector< vector< RealVector >* >& cellRGrads,
                                                        const CFuint nbrCellLGrads, const CFuint nbrCellRGrads)
{
  m_statesReconstr->reconstructGradients(cellLGrads,m_flxPntGrads[LEFT ],
                                         (*m_svFaceflxPntsRecCoefs)[m_orient][LEFT ],
                                           nbrCellLGrads);
  m_statesReconstr->reconstructGradients(cellRGrads,m_flxPntGrads[RIGHT],
                                         (*m_svFaceflxPntsRecCoefs)[m_orient][RIGHT],
                                           nbrCellRGrads);
}

//////////////////////////////////////////////////////////////////////////////

void BaseFaceTermComputer::reconstructFluxPntsGradients(const CFuint side,
                                                        const vector< vector< RealVector >* >& cellGrads,
                                                        const CFuint nbrCellGrads)
{
  m_statesReconstr->reconstructGradients(cellGrads,m_flxPntGrads[side],
                                         (*m_svFaceflxPntsRecCoefs)[m_orient][side],
                                           nbrCellGrads);
}

//////////////////////////////////////////////////////////////////////////////

void BaseFaceTermComputer::backupAndReconstructPhysVar(const CFuint side,
                                                       const CFuint iVar,
                                                       const vector< State* >& cellStates)
{
  // backup
  backupPhysVar(side,iVar);

  // reconstruct
  m_statesReconstr
      ->reconstructPhysVar(iVar,cellStates,m_flxPntSol[side],
                           (*m_svFaceflxPntsRecCoefs)[m_orient][side],
                             cellStates.size());
}

//////////////////////////////////////////////////////////////////////////////

void BaseFaceTermComputer::backupPhysVar(const CFuint side, const CFuint iVar)
{
  const CFuint nbrFlxPnts = m_flxPntSol[side].size();
  cf_assert(nbrFlxPnts <= m_backupPhysVar.size());
  for (CFuint iFlx = 0; iFlx < nbrFlxPnts; ++iFlx)
  {
    m_backupPhysVar[iFlx] = (*m_flxPntSol[side][iFlx])[iVar];
  }
}

//////////////////////////////////////////////////////////////////////////////

void BaseFaceTermComputer::restorePhysVar(const CFuint side, const CFuint iVar)
{
  const CFuint nbrFlxPnts = m_flxPntSol[side].size();
  cf_assert(nbrFlxPnts <= m_backupPhysVar.size());
  for (CFuint iFlx = 0; iFlx < nbrFlxPnts; ++iFlx)
  {
    (*m_flxPntSol[side][iFlx])[iVar] = m_backupPhysVar[iFlx];
  }
}

//////////////////////////////////////////////////////////////////////////////

void BaseFaceTermComputer::computeConvFaceTerm(RealVector& resUpdates)
{
  // compute the states data in the flux points (and not for the face averaged states --> -2)
 /// @todo broken after release 2009.3
//   m_updateVarSet->computeStatesData(m_allSol,m_allExtraVars,m_allSol.size()-2);

  // evaluate the riemann fluxes in all the flux points
  const CFuint nbrFlxPnts = m_unitNormalFlxPnt.size();
  for (CFuint iFlx = 0; iFlx < nbrFlxPnts; ++iFlx)
  {
    m_unitNormalFlxPnt[iFlx] = m_avgUnitNormal;
  }
  m_flxPntRiemannFlx = m_riemannFluxComputer->computeFlux(m_flxPntSol[LEFT],m_flxPntSol[RIGHT],
                                                          m_unitNormalFlxPnt,m_nbrFlxPnts[m_orient]);

  // compute the actual face term
  computeFaceFluxIntegralFromFlxPntFluxes(resUpdates);
}

//////////////////////////////////////////////////////////////////////////////

void BaseFaceTermComputer::computeConvFaceTermAndWaveSpeedUpdates(RealVector& resUpdates,
                                                                  CFreal& waveSpeedUpdL,
                                                                  CFreal& waveSpeedUpdR)
{
  // compute the states data in the flux points and for the face averaged states
 /// @todo broken after release 2009.3
//  m_updateVarSet->computeStatesData(m_allSol,m_allExtraVars,m_allSol.size());

  // evaluate the riemann fluxes in all the flux points
  const CFuint nbrFlxPnts = m_unitNormalFlxPnt.size();
  for (CFuint iFlx = 0; iFlx < nbrFlxPnts; ++iFlx)
  {
    m_unitNormalFlxPnt[iFlx] = m_avgUnitNormal;
  }
  m_flxPntRiemannFlx = m_riemannFluxComputer->computeFlux(m_flxPntSol[LEFT],m_flxPntSol[RIGHT],
                                                          m_unitNormalFlxPnt,m_nbrFlxPnts[m_orient]);

  // compute the actual face term
  computeFaceFluxIntegralFromFlxPntFluxes(resUpdates);

  // compute the wave speed updates for the neighbouring cells
  waveSpeedUpdL = m_updateVarSet->getMaxAbsEigenValue(*m_faceAvgSol[LEFT ],m_avgUnitNormal)*m_surf;
  waveSpeedUpdR = m_updateVarSet->getMaxAbsEigenValue(*m_faceAvgSol[RIGHT],m_avgUnitNormal)*m_surf;
}

//////////////////////////////////////////////////////////////////////////////

void BaseFaceTermComputer::computeGradientFaceTerm(vector< vector< RealVector > >& gradUpdates)
{
  // compute the averaged gradient variables in the flux points (stored in m_flxPntRiemannFlx)
  // (m_leftRVSol and m_rightRVSol are RealVector pointers to the left and right states)
  const CFuint nbrFlxPnts = m_unitNormalFlxPnt.size();
  for (CFuint iFlx = 0; iFlx < nbrFlxPnts; ++iFlx)
  {
    m_unitNormalFlxPnt[iFlx] = m_avgUnitNormal;
  }
  m_flxPntRiemannFlx = m_faceDiffFluxComputer->computeAvgGradVars(m_flxPntRVSol[LEFT],m_flxPntRVSol[RIGHT],
                                                                  m_nbrFlxPnts[m_orient]);

  // compute the actual volume term contribution to the gradient
  computeGradFaceTermFromFlxPntSol(gradUpdates);
}

//////////////////////////////////////////////////////////////////////////////

void BaseFaceTermComputer::computeDiffFaceTerm(RealVector& resUpdates)
{
  // compute the diffusive fluxes in the face flux points
  const CFuint nbrFlxPnts = m_unitNormalFlxPnt.size();
  for (CFuint iFlx = 0; iFlx < nbrFlxPnts; ++iFlx)
  {
    m_unitNormalFlxPnt[iFlx] = m_avgUnitNormal;
  }
  m_flxPntRiemannFlx = m_faceDiffFluxComputer->computeDiffFlux(m_flxPntGradPtrs[LEFT],m_flxPntGradPtrs[RIGHT],
                                                               m_flxPntRVSol[LEFT],m_flxPntRVSol[RIGHT],
                                                               m_cellVolumes[LEFT],m_cellVolumes[RIGHT],
                                                               m_surf,m_unitNormalFlxPnt,
                                                               m_nbrFlxPnts[m_orient]);

  // compute the actual face term
  computeFaceFluxIntegralFromFlxPntFluxes(resUpdates);
}

//////////////////////////////////////////////////////////////////////////////

void BaseFaceTermComputer::computeDiffFaceTermAndUpdateCoefContributions(RealVector& resUpdates,
                                                                         CFreal& updateCoefContrL,
                                                                         CFreal& updateCoefContrR)
{
  // compute the diffusive fluxes in the face flux points
  const CFuint nbrFlxPnts = m_unitNormalFlxPnt.size();
  for (CFuint iFlx = 0; iFlx < nbrFlxPnts; ++iFlx)
  {
    m_unitNormalFlxPnt[iFlx] = m_avgUnitNormal;
  }
  m_flxPntRiemannFlx = m_faceDiffFluxComputer->computeDiffFlux(m_flxPntGradPtrs[LEFT],m_flxPntGradPtrs[RIGHT],
                                                               m_flxPntRVSol[LEFT],m_flxPntRVSol[RIGHT],
                                                               m_cellVolumes[LEFT],m_cellVolumes[RIGHT],
                                                               m_surf,m_unitNormalFlxPnt,
                                                               m_nbrFlxPnts[m_orient]);

  // compute the actual face term
  computeFaceFluxIntegralFromFlxPntFluxes(resUpdates);

  // compute the update coefficient contributions for the neighbouring cells
  CFLog(WARN,"Contribution of the diffusive terms to the update coefficients is not added...");
  updateCoefContrL = 0.0;
  updateCoefContrR = 0.0;
}

//////////////////////////////////////////////////////////////////////////////

void BaseFaceTermComputer::setup()
{
  CFAUTOTRACE;

  // dimensionality and number of variables
  const CFuint dim    = PhysicalModelStack::getActive()->getDim();
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

  // create states and extra variables for face averaged solutions
  m_faceAvgSol      .resize(2);
  m_faceAvgExtraVars.resize(2);
  m_faceAvgSol      [LEFT ] = new State();
  m_faceAvgSol      [RIGHT] = new State();
  m_faceAvgExtraVars[LEFT ] = new RealVector(m_nbrExtraVars);
  m_faceAvgExtraVars[RIGHT] = new RealVector(m_nbrExtraVars);

  // resize variables
  m_cellVolumes.resize(2);
  m_avgNormal.resize(dim);
  m_avgUnitNormal.resize(dim);
}

//////////////////////////////////////////////////////////////////////////////

void BaseFaceTermComputer::unsetup()
{
  CFAUTOTRACE;

  deletePtr(m_faceAvgSol      [LEFT ]);
  deletePtr(m_faceAvgSol      [RIGHT]);
  deletePtr(m_faceAvgExtraVars[LEFT ]);
  deletePtr(m_faceAvgExtraVars[RIGHT]);
  m_faceAvgSol      .resize(0);
  m_faceAvgExtraVars.resize(0);
}

//////////////////////////////////////////////////////////////////////////////

vector<SafePtr<BaseDataSocketSink> >
    BaseFaceTermComputer::needsSockets()
{
  vector<SafePtr<BaseDataSocketSink> > result;

  result.push_back(&socket_faceSurf);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

  }  // namespace SpectralFV

}  // namespace COOLFluiD
