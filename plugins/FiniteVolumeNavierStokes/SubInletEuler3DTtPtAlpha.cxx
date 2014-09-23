#include "FiniteVolumeNavierStokes/FiniteVolumeNavierStokes.hh"
#include "NavierStokes/Euler3DVarSet.hh"
#include "FiniteVolumeNavierStokes/SubInletEuler3DTtPtAlpha.hh"
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

MethodCommandProvider<SubInletEuler3DTtPtAlpha, CellCenterFVMData, FiniteVolumeNavierStokesModule> 
subInletEuler3DTtPtAlphaFVMCCProvider("SubInletEuler3DTtPtAlphaFVMCC");

//////////////////////////////////////////////////////////////////////////////

void SubInletEuler3DTtPtAlpha::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< CFreal >("Ttot","total temperature");
  options.addConfigOption< CFreal >("alphaxy","alphaxy");
  options.addConfigOption< CFreal >("Ptot","total pressure");
  options.addConfigOption< CFreal >("alphaxz","alphaxz");
  options.addConfigOption< vector<CFuint> >("xyzIDs","IDs for the x,y,z coordinates");
}

//////////////////////////////////////////////////////////////////////////////

SubInletEuler3DTtPtAlpha::SubInletEuler3DTtPtAlpha(const std::string& name) :
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

  _alpha1 = 0.0;
   setParameter("alphaxy",&_alpha1);

  _alpha2 = 0.0;
  setParameter("alphaxz",&_alpha2);
  
  _xyzIDs = vector<CFuint>();
  setParameter("xyzIDs",&_xyzIDs);
}

//////////////////////////////////////////////////////////////////////////////

SubInletEuler3DTtPtAlpha::~SubInletEuler3DTtPtAlpha()
{
}

//////////////////////////////////////////////////////////////////////////////

void SubInletEuler3DTtPtAlpha::setGhostState(GeometricEntity *const face)
{
  State *const innerState = face->getState(0);
  State *const ghostState = face->getState(1);

  // set the physical data starting from the inner state
  _varSet->computePhysicalData(*innerState, _dataInnerState);

  const CFreal gamma = _varSet->getModel()->getGamma();
  const CFreal gammaMinus1 = gamma - 1.;
  const CFreal gammaDivGammaMinus1 = gamma/gammaMinus1;
  const CFreal R = _varSet->getModel()->getR();
  const CFreal u = _dataInnerState[EulerTerm::V+_xyzIDs[XX]+1];
  const CFreal v = _dataInnerState[EulerTerm::V+_xyzIDs[YY]+1];
  const CFreal w = _dataInnerState[EulerTerm::V+_xyzIDs[ZZ]+1];
  const CFreal V2 = u*u + v*v + w*w;
  const CFreal pInnerState = _varSet->getModel()->
     getPressureFromState(_dataInnerState[EulerTerm::P]);
  const CFreal tInnerState = _dataInnerState[EulerTerm::T];
  const CFreal machInner = sqrt(V2/(gamma*R*tInnerState));
  const CFreal coefficient = 1.0 + 0.5*gammaMinus1*machInner*machInner;
  const CFreal tTotalInner = tInnerState*coefficient;
  const CFreal coeffPow = pow(coefficient, gamma/gammaMinus1);
  const CFreal pTotalInner = pInnerState*coeffPow;
  const CFreal tgAlpha1Inner = v/u;
  const CFreal tgAlpha2Inner = w/u;
    
  // ghost state quantities
  const CFreal tTotalGhost = 2.0*_tTotal - tTotalInner;
  const CFreal pTotalGhost = 2.0*_pTotal - pTotalInner;
  const CFreal tgAlpha1Ghost = 2.0*tan(_alpha1) - tgAlpha1Inner;
  const CFreal tgAlpha2Ghost = 2.0*tan(_alpha2) - tgAlpha2Inner;
  const CFreal machGhost = machInner;
  const CFreal tGhost = tTotalGhost / coefficient;
  
  const CFreal vxGhost = machGhost*sqrt
    (gamma*R*tGhost/(1.0 + tgAlpha1Ghost*tgAlpha1Ghost + tgAlpha2Ghost*tgAlpha2Ghost));
  const CFreal vyGhost = tgAlpha1Ghost*vxGhost;
  const CFreal vzGhost = tgAlpha2Ghost*vxGhost;
  const CFreal vGhost = sqrt(vxGhost*vxGhost + vyGhost*vyGhost + vzGhost*vzGhost);
  
  // set the physical data for the ghost state
  const CFreal pGhost = pTotalGhost/coeffPow;
  _dataGhostState[EulerTerm::P] = pGhost - _varSet->getModel()->getPressInf();
  _dataGhostState[EulerTerm::RHO] = _varSet->getModel()->getDensity(pGhost, tGhost);
  _dataGhostState[EulerTerm::A] = sqrt(gamma*R*tGhost);
  _dataGhostState[EulerTerm::VX] = machGhost*_dataGhostState[EulerTerm::A]/
    sqrt(1.0 + tgAlpha1Ghost*tgAlpha1Ghost + tgAlpha2Ghost*tgAlpha2Ghost);
  _dataGhostState[EulerTerm::V+_xyzIDs[XX]+1] = vxGhost;
  _dataGhostState[EulerTerm::V+_xyzIDs[YY]+1] = vyGhost;
  _dataGhostState[EulerTerm::V+_xyzIDs[ZZ]+1] = vzGhost;
  _dataGhostState[EulerTerm::V] = vGhost;
  _dataGhostState[EulerTerm::H] = (gammaDivGammaMinus1*pGhost
				   + 0.5*_dataGhostState[EulerTerm::RHO]*_dataGhostState[EulerTerm::V]*
				   _dataGhostState[EulerTerm::V])/_dataGhostState[EulerTerm::RHO];
  _dataGhostState[EulerTerm::E] = _dataGhostState[EulerTerm::H] - pGhost/_dataGhostState[EulerTerm::RHO];
  _dataGhostState[EulerTerm::T] = tGhost;
    
  // set the ghost state starting from the physical data
  _varSet->computeStateFromPhysicalData(_dataGhostState, *ghostState);
}

//////////////////////////////////////////////////////////////////////////////

void SubInletEuler3DTtPtAlpha::setup()
{
  FVMCC_BC::setup();
  
  _varSet = getMethodData().getUpdateVar().d_castTo<EulerVarSet>();

  _varSet->getModel()->resizePhysicalData(_dataInnerState);
  _varSet->getModel()->resizePhysicalData(_dataGhostState);
  
  cf_assert(_varSet->getModel()->getPressRef() > 0.);
  cf_assert(_varSet->getModel()->getTempRef() > 0.);

  _pTotal /= _varSet->getModel()->getPressRef();
  _tTotal /= _varSet->getModel()->getTempRef();

  if (_xyzIDs.size() == 0) {
    _xyzIDs.resize(DIM_3D);
    for (CFuint i = 0; i < DIM_3D; ++i) {
      _xyzIDs[i] = i;
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD
