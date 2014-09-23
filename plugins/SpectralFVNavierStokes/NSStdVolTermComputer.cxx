#include "Framework/MethodStrategyProvider.hh"

#include "SpectralFVNavierStokes/SpectralFVNavierStokes.hh"
#include "SpectralFVNavierStokes/NSStdVolTermComputer.hh"

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
    NSStdVolTermComputer,SpectralFVMethodData,BaseVolTermComputer,SpectralFVNavierStokesModule >
  NSStdVolTermComputerProvider("NSStdVolTermComputer");

//////////////////////////////////////////////////////////////////////////////

NSStdVolTermComputer::NSStdVolTermComputer(const std::string& name) :
  StdVolTermComputer(name),
  m_navierStokesVarSet(CFNULL)
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

NSStdVolTermComputer::~NSStdVolTermComputer()
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

void NSStdVolTermComputer::computeGradientVolumeTerm(vector< vector< RealVector > >& gradUpdates)
{
  // compute the gradient variables in the flux points
  m_navierStokesVarSet->setGradientVars(m_solRVInFlxPnts,m_gradVarsInFlxPnts,m_nbrFlxPnts);

  // compute the actual volume term contribution to the gradient
  computeGradVolTermFromFlxPntSol(gradUpdates);
}

//////////////////////////////////////////////////////////////////////////////

void NSStdVolTermComputer::setup()
{
  CFAUTOTRACE;

  StdVolTermComputer::setup();

  // dynamic_cast the diffusive varset to a navier-stokes varset
  if (getMethodData().hasDiffTerm())
  {
    m_navierStokesVarSet = m_diffusiveVarSet.d_castTo< NavierStokesVarSet >();
  }
}

//////////////////////////////////////////////////////////////////////////////

void NSStdVolTermComputer::unsetup()
{
  CFAUTOTRACE;

  StdVolTermComputer::unsetup();
}

//////////////////////////////////////////////////////////////////////////////

  }  // namespace SpectralFV

}  // namespace COOLFluiD

