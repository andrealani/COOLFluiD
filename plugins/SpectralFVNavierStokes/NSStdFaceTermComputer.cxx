#include "Framework/MethodStrategyProvider.hh"

#include "SpectralFVNavierStokes/SpectralFVNavierStokes.hh"
#include "SpectralFVNavierStokes/NSStdFaceTermComputer.hh"

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
    NSStdFaceTermComputer,SpectralFVMethodData,BaseFaceTermComputer,SpectralFVNavierStokesModule >
  NSStdFaceTermComputerProvider("NSStdFaceTermComputer");

//////////////////////////////////////////////////////////////////////////////

NSStdFaceTermComputer::NSStdFaceTermComputer(const std::string& name) :
  StdFaceTermComputer(name),
  m_navierStokesVarSet(CFNULL)
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

NSStdFaceTermComputer::~NSStdFaceTermComputer()
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

void NSStdFaceTermComputer::computeDiffFaceTermAndUpdateCoefContributions(RealVector& resUpdates,
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
  const CFreal muXSurfXSurfXRatio = m_navierStokesVarSet->getCurrDynViscosity()*m_surf*m_surf*m_cflConvDiffRatio;

  const CFreal rhoL = m_navierStokesVarSet->getDensity(*m_faceAvgSol[LEFT ]);
  updateCoefContrL = muXSurfXSurfXRatio/rhoL/m_cellVolumes[LEFT ];

  const CFreal rhoR = m_navierStokesVarSet->getDensity(*m_faceAvgSol[RIGHT]);
  updateCoefContrR = muXSurfXSurfXRatio/rhoR/m_cellVolumes[RIGHT];
}

//////////////////////////////////////////////////////////////////////////////

void NSStdFaceTermComputer::setup()
{
  CFAUTOTRACE;

  StdFaceTermComputer::setup();

  // dynamic_cast the diffusive varset to a navier-stokes varset
  if (getMethodData().hasDiffTerm())
  {
    m_navierStokesVarSet = m_diffusiveVarSet.d_castTo< NavierStokesVarSet >();
  }
}

//////////////////////////////////////////////////////////////////////////////

void NSStdFaceTermComputer::unsetup()
{
  CFAUTOTRACE;

  StdFaceTermComputer::unsetup();
}

//////////////////////////////////////////////////////////////////////////////

  }  // namespace SpectralFV

}  // namespace COOLFluiD
