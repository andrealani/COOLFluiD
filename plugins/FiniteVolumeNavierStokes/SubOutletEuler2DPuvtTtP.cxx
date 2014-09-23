#include "Framework/MethodCommandProvider.hh"

#include "NavierStokes/Euler2DVarSet.hh"

#include "FiniteVolumeNavierStokes/FiniteVolumeNavierStokes.hh"
#include "FiniteVolumeNavierStokes/SubOutletEuler2DPuvtTtP.hh"

#include <cmath>

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

MethodCommandProvider<SubOutletEuler2DPuvtTtP, CellCenterFVMData, FiniteVolumeNavierStokesModule>
subOutletEuler2DPuvtTtPFVMCCProvider("SubOutletEuler2DPuvtTtPFVMCC");

//////////////////////////////////////////////////////////////////////////////

void SubOutletEuler2DPuvtTtP::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< CFreal > ("Tt","total temperature");
}

//////////////////////////////////////////////////////////////////////////////

SubOutletEuler2DPuvtTtP::SubOutletEuler2DPuvtTtP(const std::string& name) :
  FVMCC_BC(name)
{
  addConfigOptionsTo(this);
  Tt_imposed = 0.0;
  setParameter("Tt",&Tt_imposed);
}

//////////////////////////////////////////////////////////////////////////////

SubOutletEuler2DPuvtTtP::~SubOutletEuler2DPuvtTtP()
{
}

//////////////////////////////////////////////////////////////////////////////

void SubOutletEuler2DPuvtTtP::setGhostState(GeometricEntity *const face)
{
  State *const innerState = face->getState(0);
  State *const ghostState = face->getState(1);

  const CFreal gamma = (getMethodData().getUpdateVar().d_castTo<Euler2DVarSet>())->getModel()->getGamma();
  const CFreal R = (getMethodData().getUpdateVar().d_castTo<Euler2DVarSet>())->getModel()->getR();
  const CFreal M_inner_state =  sqrt(((*innerState)[1]*(*innerState)[1]+(*innerState)[2]*(*innerState)[2])/((*innerState)[3]*gamma*R));
  const CFreal M_ghost = M_inner_state;
  const CFreal Tt_inner_state = (*innerState)[3]*(1+((gamma-1)/2)*M_inner_state*M_inner_state);
  const CFreal Tt_ghost = 2*Tt_imposed - Tt_inner_state;
  const CFreal T_ghost = Tt_ghost/(1+((gamma-1)/2)*M_ghost*M_ghost);
  const CFreal P_ghost = (*innerState)[0]*pow(T_ghost/(*innerState)[3],(gamma/(gamma-1)));
  (*ghostState)[0] = P_ghost;
  (*ghostState)[1] = (*innerState)[1];
  (*ghostState)[2] = (*innerState)[2];
  (*ghostState)[3] = T_ghost;
}

//////////////////////////////////////////////////////////////////////////////

void SubOutletEuler2DPuvtTtP::setup()
{

  FVMCC_BC::setup();

  Tt_imposed /= (getMethodData().getUpdateVar().d_castTo<Euler2DVarSet>())->getModel()->getTempRef();

  m_dim = PhysicalModelStack::getActive()->getDim();
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD
