#include "Framework/MethodStrategyProvider.hh"

#include "SpectralFDNavierStokes/SpectralFDNavierStokes.hh"
#include "SpectralFDNavierStokes/NavierStokesBndFaceTermComputer.hh"

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
    NavierStokesBndFaceTermComputer,SpectralFDMethodData,BaseBndFaceTermComputer,SpectralFDNavierStokesModule >
NavierStokesBndFaceTermComputerProvider("NavierStokesBndFaceTermComputer");

//////////////////////////////////////////////////////////////////////////////

NavierStokesBndFaceTermComputer::NavierStokesBndFaceTermComputer(const std::string& name) :
  BaseBndFaceTermComputer(name),
  m_navierStokesVarSet(CFNULL)
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

NavierStokesBndFaceTermComputer::~NavierStokesBndFaceTermComputer()
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

void NavierStokesBndFaceTermComputer::computeDiffFaceTermAndUpdateCoefContributions(RealVector& resUpdates,
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
    /// @todo replace cell volume with Jacobian determinant at the flux point
  }
}

//////////////////////////////////////////////////////////////////////////////

void NavierStokesBndFaceTermComputer::setup()
{
  CFAUTOTRACE;

  // call setup of the parent class
  BaseBndFaceTermComputer::setup();

  // dynamic_cast the diffusive varset to a navier-stokes varset
  m_navierStokesVarSet = m_diffusiveVarSet.d_castTo< NavierStokesVarSet >();
}

//////////////////////////////////////////////////////////////////////////////

void NavierStokesBndFaceTermComputer::unsetup()
{
  CFAUTOTRACE;

  // call unsetup of the parent class
  BaseBndFaceTermComputer::unsetup();
}

//////////////////////////////////////////////////////////////////////////////

  }  // namespace SpectralFD

}  // namespace COOLFluiD
