#include "Framework/MethodStrategyProvider.hh"

#include "SpectralFDNavierStokes/SpectralFDNavierStokes.hh"
#include "SpectralFDNavierStokes/NavierStokesVolTermComputer.hh"

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
    NavierStokesVolTermComputer,SpectralFDMethodData,BaseVolTermComputer,SpectralFDNavierStokesModule >
NavierStokesVolTermComputerProvider("NavierStokesVolTermComputer");

//////////////////////////////////////////////////////////////////////////////

NavierStokesVolTermComputer::NavierStokesVolTermComputer(const std::string& name) :
  BaseVolTermComputer(name),
  m_navierStokesVarSet(CFNULL)
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

NavierStokesVolTermComputer::~NavierStokesVolTermComputer()
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

void NavierStokesVolTermComputer::computeGradientVolumeTerm(vector< vector< RealVector > >& gradUpdates)
{
  // compute the gradient variables in the flux points
  m_navierStokesVarSet->setGradientVars(m_solRVInFlxPnts,m_gradVarsInFlxPnts,m_nbrFlxPnts);

  // compute the actual volume term contribution to the gradient
  computeGradVolTermFromFlxPntSol(gradUpdates);
}

//////////////////////////////////////////////////////////////////////////////

void NavierStokesVolTermComputer::setup()
{
  CFAUTOTRACE;

  // call setup of the parent class
  BaseVolTermComputer::setup();

  // dynamic_cast the diffusive varset to a navier-stokes varset
  m_navierStokesVarSet = m_diffusiveVarSet.d_castTo< NavierStokesVarSet >();
}

//////////////////////////////////////////////////////////////////////////////

void NavierStokesVolTermComputer::unsetup()
{
  CFAUTOTRACE;

  // call unsetup of the parent class
  BaseVolTermComputer::unsetup();
}

//////////////////////////////////////////////////////////////////////////////

  }  // namespace SpectralFD

}  // namespace COOLFluiD
