#include "FiniteVolumeNavierStokes/FiniteVolumeNavierStokes.hh"
#include "NavierStokes/Euler1DVarSet.hh"
#include "FiniteVolumeNavierStokes/SubInletEuler1DTtPt.hh"
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

MethodCommandProvider<SubInletEuler1DTtPt, CellCenterFVMData, FiniteVolumeNavierStokesModule> 
subInletEuler1DTtPtFVMCCProvider("SubInletEuler1DTtPtFVMCC");

//////////////////////////////////////////////////////////////////////////////

void SubInletEuler1DTtPt::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOption< CFreal >("Ttot","total temperature");
   options.addConfigOption< CFreal >("Ptot","total pressure");
}

//////////////////////////////////////////////////////////////////////////////

SubInletEuler1DTtPt::SubInletEuler1DTtPt(const std::string& name) :
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

}

//////////////////////////////////////////////////////////////////////////////

SubInletEuler1DTtPt::~SubInletEuler1DTtPt()
{
}

//////////////////////////////////////////////////////////////////////////////

void SubInletEuler1DTtPt::setGhostState(GeometricEntity *const face)
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
  const CFreal uSqvSq = u*u;
  const CFreal pInnerState = _dataInnerState[EulerTerm::P];
  const CFreal tInnerState = pInnerState / (R*_dataInnerState[EulerTerm::RHO]);
  const CFreal machInner = sqrt(uSqvSq/(gamma*R*tInnerState));
  const CFreal coefficient = 1.0 + 0.5*gammaMinus1*machInner*machInner;
  const CFreal tTotalInner = tInnerState*coefficient;
  const CFreal coeffPow = pow(coefficient, gamma/gammaMinus1);
  const CFreal pTotalInner = pInnerState*coeffPow;

  // ghost state quantities
  const CFreal tTotalGhost = 2.0*_tTotal - tTotalInner;
  const CFreal pTotalGhost = 2.0*_pTotal - pTotalInner;
  const CFreal machGhost = machInner;
  const CFreal tGhost = tTotalGhost / coefficient;

  // set the physical data for the ghost state
  _dataGhostState[EulerTerm::P] = pTotalGhost/coeffPow;
  _dataGhostState[EulerTerm::RHO] = _dataGhostState[EulerTerm::P]/(R*tGhost);
  _dataGhostState[EulerTerm::VX] = machGhost*sqrt(gamma*R*tGhost);						 
  _dataGhostState[EulerTerm::V] = sqrt(_dataGhostState[EulerTerm::VX]*
				 _dataGhostState[EulerTerm::VX]); 
  _dataGhostState[EulerTerm::H] = (gammaDivGammaMinus1*_dataGhostState[EulerTerm::P]
				   + 0.5*_dataGhostState[EulerTerm::RHO]*
				   _dataGhostState[EulerTerm::V]*
				   _dataGhostState[EulerTerm::V])/
  _dataGhostState[EulerTerm::RHO];
  _dataGhostState[EulerTerm::A] = sqrt(gamma*_dataGhostState[EulerTerm::P]/
				       _dataGhostState[EulerTerm::RHO]);
  
  _dataGhostState[EulerTerm::T] = tGhost;
  
  // set the ghost state starting from the physical data
  _varSet->computeStateFromPhysicalData(_dataGhostState, *ghostState);
 }

//////////////////////////////////////////////////////////////////////////////

void SubInletEuler1DTtPt::setup()
{
  FVMCC_BC::setup();
	
  _varSet = getMethodData().getUpdateVar().d_castTo<Euler1DVarSet>();
  _varSet->getModel()->resizePhysicalData(_dataInnerState);
  _varSet->getModel()->resizePhysicalData(_dataGhostState);

  _pTotal/=_varSet->getModel()->getPressRef();
  _tTotal/=_varSet->getModel()->getTempRef();

}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD
