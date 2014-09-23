#include "Framework/MethodStrategyProvider.hh"

#include "SpectralFD/SpectralFDElementData.hh"

#include "SpectralFDNavierStokes/SpectralFDNavierStokes.hh"
#include "SpectralFDNavierStokes/NSIPBndFaceTermComputer.hh"

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
    NSIPBndFaceTermComputer,SpectralFDMethodData,BaseBndFaceTermComputer,SpectralFDNavierStokesModule >
NSIPBndFaceTermComputerProvider("NSIPBndFaceTermComputer");

//////////////////////////////////////////////////////////////////////////////

NSIPBndFaceTermComputer::NSIPBndFaceTermComputer(const std::string& name) :
  IPBndFaceTermComputer(name),
  m_navierStokesVarSet(CFNULL),
//   m_nsFaceDiffusiveFluxIP(CFNULL),
  m_flxPntIntGradVarGrads(),
  m_flxPntGhostGradVarGrads()
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

NSIPBndFaceTermComputer::~NSIPBndFaceTermComputer()
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

void NSIPBndFaceTermComputer::computeGhostGradients()
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
  m_navierStokesVarSet->setStateGradients(m_flxPntGhostRVSol,
                                          m_flxPntGhostGradVarGrads,
                                          m_flxPntGhostGrads,
                                          m_nbrFlxPnts[m_orient]);
}

//////////////////////////////////////////////////////////////////////////////

// void NSIPBndFaceTermComputer::computeDiffFaceTerm(RealVector& resUpdates)
// {
//   // compute the diffusive fluxes in the face flux points
//   m_nsFaceDiffusiveFluxIP->setBCStateComputer(m_bcStateComputer);
//   m_flxPntRiemannFlux = m_nsFaceDiffusiveFluxIP->computeBndDiffFlux(m_flxPntIntGrads,
//                                                                     m_flxPntIntRVSol,m_flxPntGhostRVSol,
//                                                                     m_faceInvCharLengths,m_unitNormalFlxPnts,
//                                                                     m_flxPntCoords,m_nbrFlxPnts[m_orient]);
//
//   // compute the actual face term
//   computeFaceTermFromFlxPntFluxes(resUpdates);
// }

//////////////////////////////////////////////////////////////////////////////

void NSIPBndFaceTermComputer::computeDiffFaceTermAndUpdateCoefContributions(RealVector& resUpdates,
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

void NSIPBndFaceTermComputer::setup()
{
  CFAUTOTRACE;

  // call setup of the parent class
  IPBndFaceTermComputer::setup();

  // dynamic_cast the diffusive varset to a navier-stokes varset
  m_navierStokesVarSet = m_diffusiveVarSet.d_castTo< NavierStokesVarSet >();

  // dynamic_cast the face diffusive flux computer
//   m_nsFaceDiffusiveFluxIP = m_faceDiffFluxComputer.d_castTo< NSFaceDiffusiveFluxIPApproach >();

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

void NSIPBndFaceTermComputer::unsetup()
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
  IPBndFaceTermComputer::unsetup();
}

//////////////////////////////////////////////////////////////////////////////

  }  // namespace SpectralFD

}  // namespace COOLFluiD
