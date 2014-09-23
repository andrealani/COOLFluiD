#include "Framework/MethodStrategyProvider.hh"

#include "SpectralFVNavierStokes/SpectralFVNavierStokes.hh"
#include "SpectralFVNavierStokes/NSQuadFreeBndFaceTermComputer.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Physics::NavierStokes;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace SpectralFV {

//////////////////////////////////////////////////////////////////////////////

Framework::MethodStrategyProvider<
    NSQuadFreeBndFaceTermComputer,SpectralFVMethodData,BaseBndFaceTermComputer,SpectralFVNavierStokesModule >
  NSQuadFreeBndFaceTermComputerProvider("NSQuadFreeBndFaceTermComputer");

//////////////////////////////////////////////////////////////////////////////

NSQuadFreeBndFaceTermComputer::NSQuadFreeBndFaceTermComputer(const std::string& name) :
  QuadFreeBndFaceTermComputer(name),
  m_navierStokesVarSet(CFNULL)
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

NSQuadFreeBndFaceTermComputer::~NSQuadFreeBndFaceTermComputer()
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

void NSQuadFreeBndFaceTermComputer::computeDiffFaceTermAndUpdateCoefContributions(RealVector& resUpdates,
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
  const CFreal muXSurfXSurfXRatio = m_navierStokesVarSet->getCurrDynViscosity()*m_surf*m_surf*m_cflConvDiffRatio;
  const CFreal rho = m_navierStokesVarSet->getDensity(*m_faceAvgSolInt);
  updateCoefContr = muXSurfXSurfXRatio/rho/m_cellVolume;
}

//////////////////////////////////////////////////////////////////////////////

void NSQuadFreeBndFaceTermComputer::setup()
{
  CFAUTOTRACE;

  QuadFreeBndFaceTermComputer::setup();

  // dynamic_cast the diffusive varset to a navier-stokes varset
  if (getMethodData().hasDiffTerm())
  {
    m_navierStokesVarSet = m_diffusiveVarSet.d_castTo< NavierStokesVarSet >();
  }
}

//////////////////////////////////////////////////////////////////////////////

void NSQuadFreeBndFaceTermComputer::unsetup()
{
  CFAUTOTRACE;


  QuadFreeBndFaceTermComputer::unsetup();
}

//////////////////////////////////////////////////////////////////////////////

  }  // namespace SpectralFV

}  // namespace COOLFluiD
