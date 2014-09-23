#include "FiniteVolumeNavierStokes/FiniteVolumeNavierStokes.hh"
#include "NavierStokes/Euler2DVarSet.hh"
#include "FiniteVolumeNavierStokes/SubInletNSTurb2DTtPtAlpha.hh"
#include "Framework/MethodCommandProvider.hh"
#include "NavierStokes/NavierStokes2DVarSet.hh"

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

MethodCommandProvider<SubInletNSTurb2DTtPtAlpha, CellCenterFVMData, FiniteVolumeNavierStokesModule> 
subInletNSTurb2DTtPtAlphaFVMCCProvider("SubInletNSTurb2DTtPtAlphaFVMCC");

//////////////////////////////////////////////////////////////////////////////

void SubInletNSTurb2DTtPtAlpha::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOption< CFreal >("Ttot","total temperature");
   options.addConfigOption< CFreal >("Ptot","total pressure");
   options.addConfigOption< CFreal >("alpha","alpha");
}

//////////////////////////////////////////////////////////////////////////////

SubInletNSTurb2DTtPtAlpha::SubInletNSTurb2DTtPtAlpha(const std::string& name) :
  FVMCC_BC(name),
  _varSetTurb(CFNULL),
  _diffVarTurb(CFNULL),
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

SubInletNSTurb2DTtPtAlpha::~SubInletNSTurb2DTtPtAlpha()
{
}

//////////////////////////////////////////////////////////////////////////////

void SubInletNSTurb2DTtPtAlpha::setGhostState(GeometricEntity *const face)
{
  State *const innerState = face->getState(0);
  State *const ghostState = face->getState(1);

  // set the physical data starting from the inner state
  _varSetTurb->computePhysicalData(*innerState, _dataInnerState);

  const CFreal gamma = _varSetTurb->getModel()->getGamma();
  const CFreal gammaMinus1 = gamma - 1.;
  const CFreal gammaDivGammaMinus1 = gamma/gammaMinus1;
  const CFreal R = _varSetTurb->getModel()->getR();
  const CFreal u = _dataInnerState[EulerTerm::VX];
  const CFreal v = _dataInnerState[EulerTerm::VY];
  const CFreal uSqvSq = u*u + v*v;
  const CFreal pInnerState = _dataInnerState[EulerTerm::P];
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
  _dataGhostState[EulerTerm::P] = pTotalGhost/coeffPow;
  _dataGhostState[EulerTerm::RHO] = _dataGhostState[EulerTerm::P]/(R*tGhost);
  _dataGhostState[EulerTerm::VX] = machGhost*sqrt(gamma*R*tGhost/
						  (1.0 + tgAlphaGhost*tgAlphaGhost));
  _dataGhostState[EulerTerm::VY] = tgAlphaGhost*_dataGhostState[EulerTerm::VX];
  _dataGhostState[EulerTerm::V] = sqrt(_dataGhostState[EulerTerm::VX]*
				 _dataGhostState[EulerTerm::VX] +
				 _dataGhostState[EulerTerm::VY]*
				 _dataGhostState[EulerTerm::VY]);
  _dataGhostState[EulerTerm::H] = (gammaDivGammaMinus1*_dataGhostState[EulerTerm::P]
				   + 0.5*_dataGhostState[EulerTerm::RHO]*
				   _dataGhostState[EulerTerm::V]*
				   _dataGhostState[EulerTerm::V])/
    _dataGhostState[EulerTerm::RHO];
  _dataGhostState[EulerTerm::A] = sqrt(gamma*_dataGhostState[EulerTerm::P]/
				       _dataGhostState[EulerTerm::RHO]);
  
  _dataGhostState[EulerTerm::T] = tGhost;
    
//unused//  const CFuint iK = _varSetTurb->getModel()->getFirstScalarVar(0);
//unused//  const CFuint nbTurbVars = _varSetTurb->getModel()->getNbScalarVars(0);





    ///@todo this is only valid for k-Omega model
    //check the adimensionalisation - velocities are in (m/s), temperature in (Kelvin)

    //Values taken from: F.R. Menter: Two-Eq Eddy-Viscosity Turbulence Models for Engineering Applications (Aug 1994)
    const CFreal L = PhysicalModelStack::getActive()->getImplementor()->getRefLength() ;
    const CFreal Uinf = _dataGhostState[EulerTerm::V] ;
    //Wilcox BC
    _turbVars[1] = Uinf/L;
    //Menter BC
    //_turbVars[1] = 10.*Uinf/L;

    const CFreal pdim = _dataGhostState[EulerTerm::P] * _varSetTurb->getModel()->getPressRef();
    const CFreal Tdim = _dataGhostState[EulerTerm::P] / ( _dataGhostState[EulerTerm::RHO] * R ) * _varSetTurb->getModel()->getTempRef();
    
    const CFreal muInf = _diffVarTurb->getModel().getDynViscosityDim(pdim, Tdim)/
      (_diffVarTurb->getModel().getReferencePhysicalData())[NSTurbTerm::MU];
        
    //upper bound
    const CFreal muTurbInf = muInf/100.;
    //lower bound
    //const CFreal muTurbInf = muInf/100000.;

    _turbVars[0] = _turbVars[1] * muTurbInf ;


    _dataGhostState[4] = _turbVars[0];
    _dataGhostState[5] = _turbVars[1];

  // set the ghost state starting from the physical data
  _varSetTurb->computeStateFromPhysicalData(_dataGhostState, *ghostState);
 }

//////////////////////////////////////////////////////////////////////////////

void SubInletNSTurb2DTtPtAlpha::setup()
{
  FVMCC_BC::setup();
 
  _varSetTurb = getMethodData().getUpdateVar().d_castTo<ConvTurb2DVarSet>();
  cf_assert(_varSetTurb.isNotNull());

  _varSetTurb->getModel()->resizePhysicalData(_dataInnerState);
  _varSetTurb->getModel()->resizePhysicalData(_dataGhostState);

  if(_turbVars.size() == 0){
    _turbVars.resize(_varSetTurb->getModel()->getNbScalarVars(0));
  }

  _diffVarTurb = getMethodData().getDiffusiveVar().d_castTo<DiffTurb2DVarSet>();
  cf_assert(_diffVarTurb.isNotNull());

  _pTotal/=_varSetTurb->getModel()->getPressRef();
  _tTotal/=_varSetTurb->getModel()->getTempRef();

}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD
