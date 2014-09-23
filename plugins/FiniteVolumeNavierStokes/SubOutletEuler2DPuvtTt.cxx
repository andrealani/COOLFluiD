#include "Framework/MethodCommandProvider.hh"

#include "NavierStokes/Euler2DVarSet.hh"

#include "FiniteVolumeNavierStokes/FiniteVolumeNavierStokes.hh"
#include "FiniteVolumeNavierStokes/SubOutletEuler2DPuvtTt.hh"

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

MethodCommandProvider<SubOutletEuler2DPuvtTt, CellCenterFVMData, FiniteVolumeNavierStokesModule>
subOutletEuler2DPuvtTtFVMCCProvider("SubOutletEuler2DPuvtTtFVMCC");

//////////////////////////////////////////////////////////////////////////////

void SubOutletEuler2DPuvtTt::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< CFreal > ("Tt","total temperature");
}

//////////////////////////////////////////////////////////////////////////////

SubOutletEuler2DPuvtTt::SubOutletEuler2DPuvtTt(const std::string& name) :
  FVMCC_BC(name)
{
  addConfigOptionsTo(this);
  Tt_imposed = 0.0;
  setParameter("Tt",&Tt_imposed);
}

//////////////////////////////////////////////////////////////////////////////

SubOutletEuler2DPuvtTt::~SubOutletEuler2DPuvtTt()
{
}

//////////////////////////////////////////////////////////////////////////////

void SubOutletEuler2DPuvtTt::setGhostState(GeometricEntity *const face)
{
  State *const innerState = face->getState(0);
  State *const ghostState = face->getState(1);

  (*ghostState)[0] = (*innerState)[0];
  (*ghostState)[1] = (*innerState)[1];
  (*ghostState)[2] = (*innerState)[2];
  const CFreal gamma = (getMethodData().getUpdateVar().d_castTo<Euler2DVarSet>())->getModel()->getGamma();
  const CFreal R = (getMethodData().getUpdateVar().d_castTo<Euler2DVarSet>())->getModel()->getR();
  const CFreal M_inner_state =  sqrt(((*innerState)[1]*(*innerState)[1]+(*innerState)[2]*(*innerState)[2])/((*innerState)[3]*gamma*R));
  const CFreal Tt_inner_state = (*innerState)[3]*(1+((gamma-1)/2)*M_inner_state*M_inner_state);
  const CFreal Tt_ghost = 2*Tt_imposed - Tt_inner_state;
  const CFreal T_ghost = Tt_ghost/(1+((gamma-1)/2)*M_inner_state*M_inner_state);
  (*ghostState)[3] = T_ghost;
}

//////////////////////////////////////////////////////////////////////////////

void SubOutletEuler2DPuvtTt::setup()
{

  FVMCC_BC::setup();
  
  Tt_imposed /= (getMethodData().getUpdateVar().d_castTo<Euler2DVarSet>())->getModel()->getTempRef();

  m_dim = PhysicalModelStack::getActive()->getDim();
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD
