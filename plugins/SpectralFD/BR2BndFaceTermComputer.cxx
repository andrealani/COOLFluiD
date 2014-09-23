#include "Common/BadValueException.hh"
#include "Framework/MethodStrategyProvider.hh"

#include "SpectralFD/BR2BndFaceTermComputer.hh"
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
    BR2BndFaceTermComputer,SpectralFDMethodData,BaseBndFaceTermComputer,SpectralFDModule >
    BR2BndFaceTermComputerProvider("BR2BndFaceTermComputer");

//////////////////////////////////////////////////////////////////////////////

void BR2BndFaceTermComputer::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< CFreal >("Alpha","Parameter of Bassi-Rebay 2 approach determining the amount of damping.");
}

//////////////////////////////////////////////////////////////////////////////

BR2BndFaceTermComputer::BR2BndFaceTermComputer(const std::string& name) :
  CompactBndFaceTermComputer(name),
  m_solPntsLocalCoords(CFNULL),
  m_backupIntStateGrad(),
  m_backupGhostStateGrad(),
  m_solJacobDet(),
  m_alpha(),
  m_liftingOperators(),
  m_liftingOperatorsPtrs(),
  m_liftOperatorTerms(),
  m_liftingOperator(),
  m_backupLiftingOperatorTerms()
{
  CFAUTOTRACE;
  addConfigOptionsTo(this);

  m_alpha = 0.0;
  setParameter("Alpha",&m_alpha);
}

//////////////////////////////////////////////////////////////////////////////

BR2BndFaceTermComputer::~BR2BndFaceTermComputer()
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

void BR2BndFaceTermComputer::configure ( Config::ConfigArgs& args )
{
  CompactBndFaceTermComputer::configure(args);
}

//////////////////////////////////////////////////////////////////////////////

void BR2BndFaceTermComputer::setFaceTermData()
{
  // call setFaceTermData of parent class
  CompactBndFaceTermComputer::setFaceTermData();

  // get the SpectralFDElementData
  vector< SpectralFDElementData* >& sdLocalData = getMethodData().getSDLocalData();

  // get solution point mapped coordinates for Jacobian determinant computation
  m_solPntsLocalCoords = sdLocalData[0]->getSolPntsLocalCoords();
}

//////////////////////////////////////////////////////////////////////////////

void BR2BndFaceTermComputer::computeNeighbourCellData()
{
  // call computeNeighbourCellData() of parent class
  CompactBndFaceTermComputer::computeNeighbourCellData();

  // get neighbouring cell
  GeometricEntity* cell = m_face->getNeighborGeo(0);

  // compute solution point Jacobian determinants
  m_solJacobDet = cell->computeGeometricShapeFunctionJacobianDeterminant(*m_solPntsLocalCoords);
}

//////////////////////////////////////////////////////////////////////////////

void BR2BndFaceTermComputer::reconstructFluxPntsSolPolyGrads(const vector< State* >& cellIntStates)
{
  // compute lifting operators
  for (CFuint iFlx = 0; iFlx < m_nbrFlxPnts[m_orient]; ++iFlx)
  {
    m_flxPntRiemannFlux[iFlx] = m_alpha*(*m_flxPntIntRVSol[iFlx] - *m_flxPntGhostRVSol[iFlx]);
  }
  computeGradFaceTermFromFlxPntSol(m_liftingOperators);

  // compute the internal gradients
  m_statesReconstr->computePolyGradients(cellIntStates,m_flxPntIntGrads,
                                         *m_coefSolPolyDerivInFlxPnts,
                                         (*m_faceFlxPntConn)[m_orient],
                                         *m_solPntIdxsSolPolyDerivInFlxPnts,
                                         m_cellInvJacobMatr);

  // compute the ghost gradients (this has to be done before adding the lifting operators, otherwise their contribution vanishes)
  computeGhostGradients();

  // subtract lifting operator terms
  const CFuint nbrSol = m_liftingOperators.size();
  for (CFuint iSol = 0; iSol < nbrSol; ++iSol)
  {
    for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
    {
      m_liftingOperators[iSol][iEq] /= m_solJacobDet[iSol];
    }
  }
  m_statesReconstr->reconstructGradients(m_liftingOperatorsPtrs,m_liftOperatorTerms,*m_flxPntsRecCoefs,
                                         (*m_faceFlxPntConn)[m_orient],
                                         *m_flxPntMatrixIdxForReconstruction,
                                         *m_solPntIdxsForReconstruction);
  for (CFuint iFlx = 0; iFlx < m_nbrFlxPnts[m_orient]; ++iFlx)
  {
    for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
    {
      *m_flxPntIntGrads[iFlx][iEq] -= *m_liftOperatorTerms[iFlx][iEq]*2.0;// factor 2 because of averaging in facediffusiveflux computer
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void BR2BndFaceTermComputer::backupAndReconstrFluxPntsSolPolyGrad(const CFuint iVar,
                                                                  const vector< State* >& cellIntStates)
{
  // backup
  backupGradPhysVar(iVar);

  // remove lifting operators
  for (CFuint iFlx = 0; iFlx < m_nbrFlxPnts[m_orient]; ++iFlx)
  {
    for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
    {
      *m_flxPntIntGrads[iFlx][iEq] += *m_liftOperatorTerms[iFlx][iEq]*2.0;
    }
  }

  // reconstruct
  m_statesReconstr->computePolyGradient(iVar,cellIntStates,m_flxPntIntGrads,
                                        *m_coefSolPolyDerivInFlxPnts,
                                        (*m_faceFlxPntConn)[m_orient],
                                        *m_solPntIdxsSolPolyDerivInFlxPnts,
                                        m_cellInvJacobMatr);

  // compute the ghost gradients
  computeGhostGradients();

  // compute lifting operators
  for (CFuint iFlx = 0; iFlx < m_nbrFlxPnts[m_orient]; ++iFlx)
  {
    m_flxPntRiemannFlux[iFlx] = m_alpha*(*m_flxPntIntRVSol[iFlx] - *m_flxPntGhostRVSol[iFlx]);
  }
  computeGradFaceTermFromFlxPntSol(m_liftingOperators);

  // subtract lifting operator terms
  const CFuint nbrSol = m_liftingOperators.size();
  for (CFuint iSol = 0; iSol < nbrSol; ++iSol)
  {
    for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
    {
      m_liftingOperators[iSol][iEq] /= m_solJacobDet[iSol];
    }
  }
  m_statesReconstr->reconstructGradients(m_liftingOperatorsPtrs,m_liftOperatorTerms,*m_flxPntsRecCoefs,
                                         (*m_faceFlxPntConn)[m_orient],
                                         *m_flxPntMatrixIdxForReconstruction,
                                         *m_solPntIdxsForReconstruction);
  for (CFuint iFlx = 0; iFlx < m_nbrFlxPnts[m_orient]; ++iFlx)
  {
    for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
    {
      *m_flxPntIntGrads[iFlx][iEq] -= *m_liftOperatorTerms[iFlx][iEq]*2.0;// factor 2 because of averaging in facediffusiveflux computer
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void BR2BndFaceTermComputer::backupGradPhysVar(const CFuint iVar)
{
  cf_assert(m_nbrFlxPnts[m_orient] <= m_backupIntStateGrad.size());
  cf_assert(m_nbrFlxPnts[m_orient] <= m_backupGhostStateGrad.size());
  cf_assert(m_nbrFlxPnts[m_orient] <= m_backupLiftingOperatorTerms.size());
  for (CFuint iFlx = 0; iFlx < m_nbrFlxPnts[m_orient]; ++iFlx)
  {
    for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
    {
      m_backupIntStateGrad[iFlx][iEq] = *m_flxPntIntGrads[iFlx][iEq];
      m_backupGhostStateGrad[iFlx][iEq] = *m_flxPntGhostGrads[iFlx][iEq];
      m_backupLiftingOperatorTerms[iFlx][iEq] = *m_liftOperatorTerms[iFlx][iEq];
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void BR2BndFaceTermComputer::restorePhysVarGrad(const CFuint iVar)
{
  cf_assert(m_nbrFlxPnts[m_orient] <= m_backupIntStateGrad.size());
  cf_assert(m_nbrFlxPnts[m_orient] <= m_backupGhostStateGrad.size());
  cf_assert(m_nbrFlxPnts[m_orient] <= m_backupLiftingOperatorTerms.size());
  for (CFuint iFlx = 0; iFlx < m_nbrFlxPnts[m_orient]; ++iFlx)
  {
    for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
    {
      *m_flxPntIntGrads[iFlx][iEq] = m_backupIntStateGrad[iFlx][iEq];
      *m_flxPntGhostGrads[iFlx][iEq] = m_backupGhostStateGrad[iFlx][iEq];
      *m_liftOperatorTerms[iFlx][iEq] = m_backupLiftingOperatorTerms[iFlx][iEq];
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void BR2BndFaceTermComputer::setup()
{
  CFAUTOTRACE;

  // call setup of parent class
  CompactBndFaceTermComputer::setup();

  // get the local spectral FD data
  vector< SpectralFDElementData* >& sdLocalData = getMethodData().getSDLocalData();
  cf_assert(sdLocalData.size() > 0);

  // number of solution points
  const CFuint nbrSolPnts = sdLocalData[0]->getNbrOfSolPnts();

  // number of flux points
  const CFuint nbrFlxPnts = sdLocalData[0]->getNbrOfFaceFlxPnts();

  // resize m_backupIntStateGrad and m_backupGhostStateGrad
  m_backupIntStateGrad  .resize(nbrFlxPnts);
  m_backupGhostStateGrad.resize(nbrFlxPnts);
  for (CFuint iFlx = 0; iFlx < nbrFlxPnts; ++iFlx)
  {
    m_backupIntStateGrad  [iFlx].resize(m_nbrEqs,RealVector(m_dim));
    m_backupGhostStateGrad[iFlx].resize(m_nbrEqs,RealVector(m_dim));
  }

  // resize variables
  m_solJacobDet         .resize(nbrSolPnts);
  m_liftingOperator     .resize(nbrSolPnts,RealVector(m_dim));
  m_liftingOperators    .resize(nbrSolPnts);
  m_liftingOperatorsPtrs.resize(nbrSolPnts);
  for (CFuint iSol = 0; iSol < nbrSolPnts; ++iSol)
  {
    m_liftingOperators    [iSol].resize(m_nbrEqs,RealVector(m_dim));
    m_liftingOperatorsPtrs[iSol] = &m_liftingOperators[iSol];
  }
  m_liftOperatorTerms         .resize(nbrFlxPnts);
  m_backupLiftingOperatorTerms.resize(nbrFlxPnts);
  for (CFuint iFlx = 0; iFlx < nbrFlxPnts; ++iFlx)
  {
    m_liftOperatorTerms[iFlx]         .resize(m_nbrEqs);
    m_backupLiftingOperatorTerms[iFlx].resize(m_nbrEqs,RealVector(m_dim));
    for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
    {
      m_liftOperatorTerms[iFlx][iEq] = new RealVector(m_dim);
    }
  }

  // check if damping value is sufficiently large
//   const CFuint nbrCellFaces = sdLocalData[0]->getNbrCellFaces();
//   m_alpha = nbrCellFaces > m_alpha ? nbrCellFaces : m_alpha;
//   if (sdLocalData[0]->getPolyOrder() == 0)
//   {
//     m_alpha = 1.0;
//   }
  m_alpha = 1.0;
}

//////////////////////////////////////////////////////////////////////////////

void BR2BndFaceTermComputer::unsetup()
{
  CFAUTOTRACE;

  deletePtr(m_matrInverter);

  for (CFuint iFlx = 0; iFlx < m_liftOperatorTerms.size(); ++iFlx)
  {
    for (CFuint iEq = 0; iEq < m_liftOperatorTerms[iFlx].size(); ++iEq)
    {
      deletePtr(m_liftOperatorTerms[iFlx][iEq]);
    }
  }

  // call unsetup of parent class
  CompactBndFaceTermComputer::unsetup();
}

//////////////////////////////////////////////////////////////////////////////

  }  // namespace SpectralFD

}  // namespace COOLFluiD
