#include "Framework/MethodStrategyProvider.hh"

#include "SpectralFVNavierStokes/SpectralFVNavierStokes.hh"
#include "SpectralFVNavierStokes/NSQuadFreeVolTermComputer.hh"

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
    NSQuadFreeVolTermComputer,SpectralFVMethodData,BaseVolTermComputer,SpectralFVNavierStokesModule >
  NSQuadFreeVolTermComputerProvider("NSQuadFreeVolTermComputer");

//////////////////////////////////////////////////////////////////////////////

NSQuadFreeVolTermComputer::NSQuadFreeVolTermComputer(const std::string& name) :
  QuadFreeVolTermComputer(name),
  m_navierStokesVarSet(CFNULL)
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

NSQuadFreeVolTermComputer::~NSQuadFreeVolTermComputer()
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

void NSQuadFreeVolTermComputer::computeGradientVolumeTerm(vector< vector< RealVector > >& gradUpdates)
{
  // compute the gradient variables in the flux points
  m_navierStokesVarSet->setGradientVars(m_solRVInFlxPnts,m_gradVarsInFlxPnts,m_nbrFlxPnts);

  // compute the actual volume term contribution to the gradient
  computeGradVolTermFromFlxPntSol(gradUpdates);
}

//////////////////////////////////////////////////////////////////////////////

void NSQuadFreeVolTermComputer::setup()
{
  CFAUTOTRACE;

  QuadFreeVolTermComputer::setup();

  // dynamic_cast the diffusive varset to a navier-stokes varset
  if (getMethodData().hasDiffTerm())
  {
    m_navierStokesVarSet = m_diffusiveVarSet.d_castTo< NavierStokesVarSet >();
  }
}

//////////////////////////////////////////////////////////////////////////////

void NSQuadFreeVolTermComputer::unsetup()
{
  CFAUTOTRACE;

  QuadFreeVolTermComputer::unsetup();
}

//////////////////////////////////////////////////////////////////////////////

  }  // namespace SpectralFV

}  // namespace COOLFluiD

