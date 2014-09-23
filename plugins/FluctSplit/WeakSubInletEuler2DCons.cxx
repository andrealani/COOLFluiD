#include "FluctSplit/FluctSplitNavierStokes.hh"
#include "WeakSubInletEuler2DCons.hh"
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

MethodCommandProvider<WeakSubInletEuler2DCons, FluctuationSplitData, FluctSplitNavierStokesModule> weakSubInletEuler2DConsProvider("WeakSubInletEuler2DCons");

//////////////////////////////////////////////////////////////////////////////

void WeakSubInletEuler2DCons::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOption< CFreal >("Ttot","total temperature");
   options.addConfigOption< CFreal >("angle","angle");
   options.addConfigOption< CFreal >("Ptot","total pressure");
}

//////////////////////////////////////////////////////////////////////////////

WeakSubInletEuler2DCons::WeakSubInletEuler2DCons(const std::string& name) :
  WeakBC2D(name),
  _varSet()
{
   addConfigOptionsTo(this);
  _tTotal = 1.0;
   setParameter("Ttot",&_tTotal);

  _pTotal = 1.0;
   setParameter("Ptot",&_pTotal);

  _angle = 0.0;
   setParameter("angle",&_angle);

}

//////////////////////////////////////////////////////////////////////////////

WeakSubInletEuler2DCons::~WeakSubInletEuler2DCons()
{
}

//////////////////////////////////////////////////////////////////////////////

void WeakSubInletEuler2DCons::setup()
{
  CFAUTOTRACE;

  WeakBC2D::setup();

  _varSet->setup();

  _angle *= MathTools::MathConsts::CFrealPi()/180.;
  _tTotal /= _varSet->getModel()->getTempRef();
  _pTotal /= _varSet->getModel()->getPressRef();
}


//////////////////////////////////////////////////////////////////////////////

void WeakSubInletEuler2DCons::setGhostState(const State& state,
					    State& gstate)
{
 const CFreal gamma = _varSet->getModel()->getGamma();
 const CFreal gammaMinus1 = gamma - 1.;
 const CFreal R = _varSet->getModel()->getR();
  const CFreal u = state[1]/state[0];
  const CFreal v = state[2]/state[0];
  const CFreal vel2 = u*u + v*v;

//   const Node& node0 = state.getCoordinates();
  const CFreal p = gammaMinus1*(state[3] - 0.5*state[0]*vel2);
  // mach number is extrapolated from inside the domain
  const CFreal mach = sqrt(vel2/(gamma*p/state[0]));
  //CF_DEBUG_OBJ(mach);
  const CFreal coeffM = 1. + 0.5*gammaMinus1*mach*mach;
  const CFreal ghostT = _tTotal/coeffM;
  const CFreal ghostP = _pTotal/pow(coeffM, gamma/gammaMinus1);
  CFreal angle = getAngle(state);

  const CFreal tgAngle = tan(angle);
  const CFreal coeffAngle = 1. + tgAngle*tgAngle;

  

if ((tgAngle*angle) < 0.0)
{
  // density
  gstate[0] = ghostP/(R*ghostT);
  gstate[1] = -mach*sqrt(gamma*R*ghostT/coeffAngle)*gstate[0];
  gstate[2] = gstate[1]*tgAngle;
  gstate[3] = ghostP/gammaMinus1 + 0.5*(gstate[1]*gstate[1] +
					gstate[2]*gstate[2])/gstate[0];
}
else
{
  // density
  gstate[0] = ghostP/(R*ghostT);
  gstate[1] = mach*sqrt(gamma*R*ghostT/coeffAngle)*gstate[0];
  gstate[2] = gstate[1]*tgAngle;
  gstate[3] = ghostP/gammaMinus1 + 0.5*(gstate[1]*gstate[1] +
					gstate[2]*gstate[2])/gstate[0];
}
}

//////////////////////////////////////////////////////////////////////////////

void WeakSubInletEuler2DCons::configure ( Config::ConfigArgs& args )
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
CFreal WeakSubInletEuler2DCons::getAngle(const State& state)
{
  return _angle;
}
//////////////////////////////////////////////////////////////////////////////

} // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
