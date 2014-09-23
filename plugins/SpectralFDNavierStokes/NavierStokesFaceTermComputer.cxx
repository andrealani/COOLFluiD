#include "Framework/MethodStrategyProvider.hh"

#include "SpectralFDNavierStokes/SpectralFDNavierStokes.hh"
#include "SpectralFDNavierStokes/NavierStokesFaceTermComputer.hh"

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
    NavierStokesFaceTermComputer,SpectralFDMethodData,BaseFaceTermComputer,SpectralFDNavierStokesModule >
NavierStokesFaceTermComputerProvider("NavierStokesFaceTermComputer");

//////////////////////////////////////////////////////////////////////////////

NavierStokesFaceTermComputer::NavierStokesFaceTermComputer(const std::string& name) :
  BaseFaceTermComputer(name),
  m_navierStokesVarSet(CFNULL)
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

NavierStokesFaceTermComputer::~NavierStokesFaceTermComputer()
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

void NavierStokesFaceTermComputer::computeDiffFaceTermAndUpdateCoefContributions(vector< RealVector>& resUpdates,
                                                                         vector< CFreal >& updateCoefContr)
{
  // compute the face term
  computeDiffFaceTerm(resUpdates);

  // compute the update coefficient contributions for the neighbouring cells
  cf_assert(updateCoefContr.size() == 2);
  const CFreal dynVisc = m_navierStokesVarSet->getCurrDynViscosity();
  for (CFuint iSide = 0; iSide < 2; ++iSide)
  {
    updateCoefContr[iSide] = 0.0;
    for (CFuint iFlx = 0; iFlx < m_nbrFlxPnts[m_orient]; ++iFlx)
    {
      const CFreal jacobXJacobXIntCoef = m_faceJacobVecAbsSizeFlxPnts[iFlx]*
                                         m_faceJacobVecAbsSizeFlxPnts[iFlx]*
                                         (*m_faceIntegrationCoefs)[iFlx]*m_cflConvDiffRatio;
      const CFreal rho = m_navierStokesVarSet->getDensity(*m_flxPntSol[iSide][iFlx]);
      updateCoefContr[iSide] += dynVisc*jacobXJacobXIntCoef/rho/m_cellVolumes[iSide];
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void NavierStokesFaceTermComputer::setup()
{
  CFAUTOTRACE;

  // call setup of the parent class
  BaseFaceTermComputer::setup();

  // dynamic_cast the diffusive varset to a navier-stokes varset
  m_navierStokesVarSet = m_diffusiveVarSet.d_castTo< NavierStokesVarSet >();
}

//////////////////////////////////////////////////////////////////////////////

void NavierStokesFaceTermComputer::unsetup()
{
  CFAUTOTRACE;

  // call unsetup of the parent class
  BaseFaceTermComputer::unsetup();
}

//////////////////////////////////////////////////////////////////////////////

  }  // namespace SpectralFD

}  // namespace COOLFluiD
