#include "FiniteVolumeNavierStokes/FiniteVolumeNavierStokes.hh"
#include "NavierStokes/Euler2DVarSet.hh"
#include "FiniteVolumeNavierStokes/SubInletEuler2DTtPtAlpha.hh"
#include "Framework/MethodCommandProvider.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Physics::NavierStokes;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<SubInletEuler2DTtPtAlpha, CellCenterFVMData, FiniteVolumeNavierStokesModule> 
subInletEuler2DTtPtAlphaFVMCCProvider("SubInletEuler2DTtPtAlphaFVMCC");

//////////////////////////////////////////////////////////////////////////////

void SubInletEuler2DTtPtAlpha::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOption< CFreal >("Ttot","total temperature");
   options.addConfigOption< CFreal >("Ptot","total pressure");
   options.addConfigOption< CFreal >("alpha","alpha");
}

//////////////////////////////////////////////////////////////////////////////

SubInletEuler2DTtPtAlpha::SubInletEuler2DTtPtAlpha(const std::string& name) :
  FVMCC_BC(name),
  _varSet(CFNULL),
  _dataInnerState(),
  _dataGhostState()
{
   addConfigOptionsTo(this);
  _tTotal = 0.0;
   setParameter("Ttot",&_tTotal);

   _pTotal = 0.0;
   setParameter("Ptot",&_pTotal);

   _alpha = 0.0;
   setParameter("alpha",&_alpha);
}

//////////////////////////////////////////////////////////////////////////////

SubInletEuler2DTtPtAlpha::~SubInletEuler2DTtPtAlpha()
{
}

//////////////////////////////////////////////////////////////////////////////

void SubInletEuler2DTtPtAlpha::setGhostState(GeometricEntity *const face)
{
  State *const innerState = face->getState(0);
  State *const ghostState = face->getState(1);

  // set the physical data starting from the inner state
  _varSet->computePhysicalData(*innerState, _dataInnerState);

  const CFreal gamma = _varSet->getModel()->getGamma();
  const CFreal gammaMinus1 = gamma - 1.;
  const CFreal gammaDivGammaMinus1 = gamma/gammaMinus1;
  const CFreal R = _varSet->getModel()->getR();
  const CFreal u = _dataInnerState[EulerTerm::VX];
  const CFreal v = _dataInnerState[EulerTerm::VY];
  const CFreal uSqvSq = u*u + v*v;
  const CFreal pInnerState = _varSet->getModel()->getPressureFromState(_dataInnerState[EulerTerm::P]);
  const CFreal tInnerState = pInnerState / (R*_dataInnerState[EulerTerm::RHO]);
  const CFreal machInner = sqrt(uSqvSq/(gamma*R*tInnerState));
  const CFreal coefficient = 1.0 + 0.5*gammaMinus1*machInner*machInner;
  const CFreal tTotalInner = tInnerState*coefficient;
  const CFreal coeffPow = pow(coefficient, gamma/gammaMinus1);
  const CFreal pTotalInner = pInnerState*coeffPow;
  const CFreal tgAlphaInner = v/u;
  
  // ghost state quantities
  const CFreal tTotalGhost = 2.0*_tTotal - tTotalInner;
  const CFreal pTotalGhost = 2.0*_pTotal - pTotalInner;
  const CFreal tgAlphaGhost = 2.0*tan(_alpha) - tgAlphaInner;
  const CFreal machGhost = machInner;
  const CFreal tGhost = tTotalGhost / coefficient;
  
  // set the physical data for the ghost state
  const CFreal pressure = pTotalGhost/coeffPow;
  _dataGhostState[EulerTerm::P] = pressure - _varSet->getModel()->getPressInf();
  _dataGhostState[EulerTerm::RHO] = pressure/(R*tGhost); 
  _dataGhostState[EulerTerm::VX] = machGhost*sqrt(gamma*R*tGhost/(1.0 + tgAlphaGhost*tgAlphaGhost));
  _dataGhostState[EulerTerm::VY] = tgAlphaGhost*_dataGhostState[EulerTerm::VX];
  _dataGhostState[EulerTerm::V] = sqrt(_dataGhostState[EulerTerm::VX]*
				       _dataGhostState[EulerTerm::VX] +
				       _dataGhostState[EulerTerm::VY]*
				       _dataGhostState[EulerTerm::VY]);
  _dataGhostState[EulerTerm::H] = (gammaDivGammaMinus1*pressure + 0.5*_dataGhostState[EulerTerm::RHO]*
				   _dataGhostState[EulerTerm::V]*
				   _dataGhostState[EulerTerm::V])/_dataGhostState[EulerTerm::RHO];
  _dataGhostState[EulerTerm::A] = sqrt(gamma*pressure/_dataGhostState[EulerTerm::RHO]);
  _dataGhostState[EulerTerm::T] = tGhost;
  
  // set the ghost state starting from the physical data
  _varSet->computeStateFromPhysicalData(_dataGhostState, *ghostState);
}

//////////////////////////////////////////////////////////////////////////////

void SubInletEuler2DTtPtAlpha::setup()
{
  FVMCC_BC::setup();
	
  _varSet = getMethodData().getUpdateVar().d_castTo<Euler2DVarSet>();
  _varSet->getModel()->resizePhysicalData(_dataInnerState);
  _varSet->getModel()->resizePhysicalData(_dataGhostState);

  _pTotal/=_varSet->getModel()->getPressRef();
  _tTotal/=_varSet->getModel()->getTempRef();

}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD
