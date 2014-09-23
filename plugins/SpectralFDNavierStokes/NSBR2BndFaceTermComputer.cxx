#include "Framework/MethodStrategyProvider.hh"

#include "SpectralFD/SpectralFDElementData.hh"

#include "SpectralFDNavierStokes/SpectralFDNavierStokes.hh"
#include "SpectralFDNavierStokes/NSBR2BndFaceTermComputer.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Physics::NavierStokes;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace SpectralFD {

//////////////////////////////////////////////////////////////////////////////

Framework::MethodStrategyProvider<
    NSBR2BndFaceTermComputer,SpectralFDMethodData,BaseBndFaceTermComputer,SpectralFDNavierStokesModule >
NSBR2BndFaceTermComputerProvider("NSBR2BndFaceTermComputer");

//////////////////////////////////////////////////////////////////////////////

NSBR2BndFaceTermComputer::NSBR2BndFaceTermComputer(const std::string& name) :
  BR2BndFaceTermComputer(name),
  m_navierStokesVarSet(CFNULL),
//   m_nsFaceDiffusiveFluxBR2(CFNULL),
  m_flxPntIntGradVarGrads(),
  m_flxPntGhostGradVarGrads(),
  m_intGradVarGradsPntSet()
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

NSBR2BndFaceTermComputer::~NSBR2BndFaceTermComputer()
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

void NSBR2BndFaceTermComputer::computeGhostGradients()
{
  cf_assert(m_bcStateComputer.isNotNull());

  // compute internal gradient variable gradients
  m_navierStokesVarSet->setGradientVarGradients(m_flxPntIntRVSol,
                                                m_flxPntIntGrads,
                                                m_flxPntIntGradVarGrads,
                                                m_nbrFlxPnts[m_orient]);

  // compute ghost gradient variables gradients
  m_bcStateComputer->computeGhostGradients(m_flxPntIntGradVarGrads,
                                           m_flxPntGhostGradVarGrads,
                                           m_unitNormalFlxPnts,
                                           m_flxPntCoords);

  // compute ghost state gradients
  m_navierStokesVarSet->setStateGradients(m_flxPntGhostRVSol,// m_flxPntIntRVSol,
                                          m_flxPntGhostGradVarGrads,
                                          m_flxPntGhostGrads,
                                          m_nbrFlxPnts[m_orient]);
}

//////////////////////////////////////////////////////////////////////////////

// void NSBR2BndFaceTermComputer::computeDiffFaceTerm(RealVector& resUpdates)
// {
//   // compute the diffusive fluxes in the face flux points
//   m_nsFaceDiffusiveFluxBR2->setBCStateComputer(m_bcStateComputer);
//   m_flxPntRiemannFlux = m_nsFaceDiffusiveFluxBR2->computeBndDiffFlux(m_flxPntIntGrads,
//                                                                     m_flxPntIntRVSol,m_flxPntGhostRVSol,
//                                                                     m_faceInvCharLengths,m_unitNormalFlxPnts,
//                                                                     m_flxPntCoords,m_nbrFlxPnts[m_orient]);
//
//   // compute the actual face term
//   computeFaceTermFromFlxPntFluxes(resUpdates);
// }

//////////////////////////////////////////////////////////////////////////////

void NSBR2BndFaceTermComputer::computeDiffFaceTermAndUpdateCoefContributions(RealVector& resUpdates,
                                                                             CFreal& updateCoefContr)
{
  // compute the face term
  computeDiffFaceTerm(resUpdates);

  // compute the update coefficient contributions for the neighbouring cells
  const CFreal dynVisc = m_navierStokesVarSet->getCurrDynViscosity();
  updateCoefContr = 0.0;
  for (CFuint iFlx = 0; iFlx < m_nbrFlxPnts[m_orient]; ++iFlx)
  {
    const CFreal jacobXJacobXIntCoef = m_faceJacobVecAbsSizeFlxPnts[iFlx]*
                                       m_faceJacobVecAbsSizeFlxPnts[iFlx]*
                                       (*m_faceIntegrationCoefs)[iFlx]*
                                       m_cflConvDiffRatio;
    const CFreal rho = m_navierStokesVarSet->getDensity(*m_flxPntIntSol[iFlx]);
    updateCoefContr += dynVisc*jacobXJacobXIntCoef/rho/m_cellVolume;
  }
}

//////////////////////////////////////////////////////////////////////////////

void NSBR2BndFaceTermComputer::setPointSet(const std::vector< RealVector >& faceMappedCoordPntSet)
{
  // call setPointSet of parent class
  CompactBndFaceTermComputer::setPointSet(faceMappedCoordPntSet);

  // number of points in point set
  const CFuint nbrPnts = faceMappedCoordPntSet.size();

  // resize some variables
  if (nbrPnts > m_intGradVarGradsPntSet.size())
  {
    CFuint iPnt = m_intGradVarGradsPntSet.size();
    for (; iPnt < nbrPnts; ++iPnt)
    {
      vector<RealVector*> intGrads  ;
      for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
      {
        intGrads.push_back(new RealVector(m_dim));
      }
      m_intGradVarGradsPntSet.push_back(intGrads  );
    }
  }
  else if (nbrPnts < m_intGradVarGradsPntSet.size())
  {
    CFuint iPnt = m_intGradVarGradsPntSet.size()-1;
    for (; iPnt >= nbrPnts; --iPnt)
    {
      for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
      {
        deletePtr(m_intGradVarGradsPntSet[iPnt][iEq]);
      }

      // remove last elements
      m_intGradVarGradsPntSet.pop_back();
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

vector< vector< RealVector > >& NSBR2BndFaceTermComputer::reconstructGivenPntsGrads(const vector< State* >& cellIntStates)
{
  // number of output points
  const CFuint nbrPnts = m_solPolyDerivCoefsPntSet.size();
  cf_assert(nbrPnts == m_gradsPntSet.size());

  // number of solution points/basis polynomials
  const CFuint nbrPolys = cellIntStates.size();
  for (CFuint iPnt = 0; iPnt < nbrPnts; ++iPnt)
  {
    for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
    {
      RealVector grad(m_dim);
      grad = 0.0;
      for (CFuint iDim = 0; iDim < m_dim; ++iDim)
      {
        cf_assert(nbrPolys == m_solPolyDerivCoefsPntSet[iPnt][iDim].size());
        for (CFuint iPoly = 0; iPoly < nbrPolys; ++iPoly)
        {
          grad[iDim] += m_solPolyDerivCoefsPntSet[iPnt][iDim][iPoly]*(*cellIntStates[iPoly])[iEq];
        }
      }

      // transform the gradient
      *m_intGradsPntSet[iPnt][iEq] = m_cellInvJacobMatrPntSet[iPnt]*grad;
    }
  }

  // compute internal gradient variable gradients
  m_navierStokesVarSet->setGradientVarGradients(m_intSolPntSetPtrs,
                                                m_intGradsPntSet,
                                                m_intGradVarGradsPntSet,
                                                m_intGradsPntSet.size());

  // compute ghost gradients
  m_bcStateComputer->computeGhostGradients(m_intGradVarGradsPntSet,
                                           m_ghostGradsPntSet,
                                           m_unitNormalPntSet,
                                           m_coordsPntSet);

  // average internal and ghost gradients
  for (CFuint iPnt = 0; iPnt < nbrPnts; ++iPnt)
  {
    for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
    {
      m_gradsPntSet[iPnt][iEq] = 0.5*(*m_intGradVarGradsPntSet[iPnt][iEq] + *m_ghostGradsPntSet[iPnt][iEq]);
    }
  }

  /// @todo: add lifting operators

  return m_gradsPntSet;
}

//////////////////////////////////////////////////////////////////////////////

void NSBR2BndFaceTermComputer::setup()
{
  CFAUTOTRACE;

  // call setup of the parent class
  BR2BndFaceTermComputer::setup();

  // dynamic_cast the diffusive varset to a navier-stokes varset
  m_navierStokesVarSet = m_diffusiveVarSet.d_castTo< NavierStokesVarSet >();

  // dynamic_cast the face diffusive flux computer
//   m_nsFaceDiffusiveFluxBR2 = m_faceDiffFluxComputer.d_castTo< NSFaceDiffusiveFluxBR2Approach >();

  // get the local spectral FD data
  vector< SpectralFDElementData* >& sdLocalData = getMethodData().getSDLocalData();
  cf_assert(sdLocalData.size() > 0);

  // number of flux points
  const CFuint nbrFlxPnts = sdLocalData[0]->getNbrOfFaceFlxPnts();

    // create gradients for flux points
  m_flxPntIntGradVarGrads  .resize(nbrFlxPnts);
  m_flxPntGhostGradVarGrads.resize(nbrFlxPnts);
  for (CFuint iFlx = 0; iFlx < nbrFlxPnts; ++iFlx)
  {
    for (CFuint iGrad = 0; iGrad < m_nbrEqs; ++iGrad)
    {
      m_flxPntIntGradVarGrads  [iFlx].push_back(new RealVector(m_dim));
      m_flxPntGhostGradVarGrads[iFlx].push_back(new RealVector(m_dim));
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void NSBR2BndFaceTermComputer::unsetup()
{
  CFAUTOTRACE;

  for (CFuint iFlx = 0; iFlx < m_flxPntIntGradVarGrads.size(); ++iFlx)
  {
    for (CFuint iGrad = 0; iGrad < m_flxPntIntGradVarGrads[iFlx].size(); ++iGrad)
    {
      deletePtr(m_flxPntIntGradVarGrads  [iFlx][iGrad]);
      deletePtr(m_flxPntGhostGradVarGrads[iFlx][iGrad]);
    }
    m_flxPntIntGradVarGrads  [iFlx].resize(0);
    m_flxPntGhostGradVarGrads[iFlx].resize(0);
  }
  m_flxPntIntGradVarGrads  .resize(0);
  m_flxPntGhostGradVarGrads.resize(0);

  // call unsetup of the parent class
  BR2BndFaceTermComputer::unsetup();
}

//////////////////////////////////////////////////////////////////////////////

  }  // namespace SpectralFD

}  // namespace COOLFluiD
