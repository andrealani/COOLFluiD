#include "FluctSplit/FluctSplitNavierStokes.hh"
#include "WeakSubOutletEuler2DCons.hh"
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

MethodCommandProvider < WeakSubOutletEuler2DCons,
                        FluctuationSplitData,
                        FluctSplitNavierStokesModule >
theWeakSubOutletEuler2DConsProvider("WeakSubOutletEuler2DCons");

//////////////////////////////////////////////////////////////////////////////

void WeakSubOutletEuler2DCons::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOption< CFreal >("P","static pressure");
}

//////////////////////////////////////////////////////////////////////////////

WeakSubOutletEuler2DCons::WeakSubOutletEuler2DCons(const std::string& name) :
  WeakBC2D(name),
  _varSet()
{
   addConfigOptionsTo(this);
  _pressure = 1.0;
   setParameter("P",&_pressure);
}

//////////////////////////////////////////////////////////////////////////////

WeakSubOutletEuler2DCons::~WeakSubOutletEuler2DCons()
{
}

//////////////////////////////////////////////////////////////////////////////

void WeakSubOutletEuler2DCons::setup()
{
  CFAUTOTRACE;

  WeakBC2D::setup();

  _varSet->setup();

  _pressure /= _varSet->getModel()->getPressRef();
}

//////////////////////////////////////////////////////////////////////////////

void WeakSubOutletEuler2DCons::setGhostState(const State& state,
					     State& gstate)
{
  gstate[0] = state[0];
  gstate[1] = state[1];
  gstate[2] = state[2];
  gstate[3] = _pressure/(_varSet->getModel()->getGamma() - 1.0) +
    0.5*(gstate[1]*gstate[1] + gstate[2]*gstate[2])/gstate[0];
}

//////////////////////////////////////////////////////////////////////////////

void WeakSubOutletEuler2DCons::configure ( Config::ConfigArgs& args )
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

} // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
