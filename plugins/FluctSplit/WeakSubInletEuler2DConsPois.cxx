#include "FluctSplit/FluctSplitNavierStokes.hh"
#include "WeakSubInletEuler2DConsPois.hh"
#include "InwardNormalsData.hh"
#include "Framework/PhysicalModel.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/NamespaceSwitcher.hh"
//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Physics::NavierStokes;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<WeakSubInletEuler2DConsPois, FluctuationSplitData, FluctSplitNavierStokesModule> WeakSubInletEuler2DConsPoisProvider("WeakSubInletEuler2DConsPois");

void WeakSubInletEuler2DConsPois::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOption< CFreal >("Ptot","total pressure");

}

//////////////////////////////////////////////////////////////////////////////

WeakSubInletEuler2DConsPois::WeakSubInletEuler2DConsPois(const std::string& name) :
  WeakBC2D(name),
  _varSet()
{
   addConfigOptionsTo(this);

  _pTotal = 1.0;
   setParameter("Ptot",&_pTotal);
}

//////////////////////////////////////////////////////////////////////////////

WeakSubInletEuler2DConsPois::~WeakSubInletEuler2DConsPois()
{
}

//////////////////////////////////////////////////////////////////////////////

void WeakSubInletEuler2DConsPois::setup()
{
  CFAUTOTRACE;

  WeakBC2D::setup();
_varSet->setup();

  _varSet->setup();

  m_tTotal /= _varSet->getModel()->getTempRef();
  _pTotal /= _varSet->getModel()->getPressRef();
}

//////////////////////////////////////////////////////////////////////////////

void WeakSubInletEuler2DConsPois::setGhostState(const State& state,
					    State& gstate)
{
 const CFreal gamma = _varSet->getModel()->getGamma();
 const CFreal gammaMinus1 = gamma - 1.;
// const CFreal R = _varSet->getModel()->getR();
  const Node& node = state.getCoordinates();

  CFreal y = node[YY];
  const CFreal u = (y+0.5)-0.22*((y+0.5)*(y+0.5) - (y+0.5));
  const CFreal v = 0;

// unused //   const CFreal vel2 = u*u + v*v;

// const Node& node0 = state.getCoordinates();
//    const CFreal p = gammaMinus1*(state[3] - 0.5*state[0]*vel2);
  // mach number is extrapolated from inside the domain
//   const CFreal mach = sqrt(vel2/(gamma*p/state[0]));
//  const CFreal coeffM = 1. + 0.5*gammaMinus1*mach*mach;
  //const CFreal ghostT = m_tTotal/coeffM;
  const CFreal ghostP = _pTotal;
CFreal rho = state[0];
// const CFreal ghostP = rho*R*ghostT;

  gstate[0] = rho;
  gstate[1] = rho*u;
  gstate[2] = rho*v;
  gstate[3] = ghostP/gammaMinus1 + 0.5*(gstate[1]*gstate[1] +
					gstate[2]*gstate[2])/rho;

//   CF_DEBUG_OBJ(gstate[0]);
//   CF_DEBUG_OBJ(gstate[1]);
//   CF_DEBUG_OBJ(gstate[2]);
//   CF_DEBUG_OBJ(gstate[3]);
//
// CF_DEBUG_OBJ(state[0]);
//   CF_DEBUG_OBJ(state[1]);
//   CF_DEBUG_OBJ(state[2]);
//   CF_DEBUG_OBJ(state[3]);

}

//////////////////////////////////////////////////////////////////////////////

void WeakSubInletEuler2DConsPois::configure ( Config::ConfigArgs& args )
{
  WeakBC2D::configure(args);


  std::string name = getMethodData().getNamespace();
  Common::SafePtr<Namespace> nsp = NamespaceSwitcher::getInstance().getNamespace(name);
  Common::SafePtr<PhysicalModel> physModel = PhysicalModelStack::getInstance().getEntryByNamespace(nsp);

  std::string varSetName = "Euler2DCons";
  _varSet.reset((Environment::Factory<ConvectiveVarSet>::getInstance().getProvider(varSetName)->
    create(physModel->getImplementor()->getConvectiveTerm())).d_castTo<Physics::NavierStokes::Euler2DCons>());

  cf_assert(_varSet.isNotNull());

}
//////////////////////////////////////////////////////////////////////////////

CFreal WeakSubInletEuler2DConsPois::getAngle(const State& state)
{
  return 0.;
}

//////////////////////////////////////////////////////////////////////////////

} // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
