#include "FluctSplit/FluctSplitNavierStokes.hh"
#include "WeakSubOutletEuler3DConsHO.hh"
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

MethodCommandProvider<WeakSubOutletEuler3DConsHO, FluctuationSplitData, FluctSplitNavierStokesModule> weakSubOutletEuler3DConsHOProvider("WeakSubOutletEuler3DConsHO");

//////////////////////////////////////////////////////////////////////////////

void WeakSubOutletEuler3DConsHO::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOption< CFreal >("P","static pressure");
}

//////////////////////////////////////////////////////////////////////////////

WeakSubOutletEuler3DConsHO::WeakSubOutletEuler3DConsHO(const std::string& name) :
  WeakBC3DHO(name),
  _varSet()
{
   addConfigOptionsTo(this);
  _pressure = 1.0;
   setParameter("P",&_pressure);
}

//////////////////////////////////////////////////////////////////////////////

WeakSubOutletEuler3DConsHO::~WeakSubOutletEuler3DConsHO()
{
}

//////////////////////////////////////////////////////////////////////////////

void WeakSubOutletEuler3DConsHO::setup()
{
  CFAUTOTRACE;

  WeakBC3DHO::setup();
  _varSet->setup();

  _pressure /= _varSet->getModel()->getPressRef();
}

//////////////////////////////////////////////////////////////////////////////

void WeakSubOutletEuler3DConsHO::configure ( Config::ConfigArgs& args )
{
  WeakBC3DHO::configure(args);

  std::string name = getMethodData().getNamespace();
  Common::SafePtr<Namespace> nsp = NamespaceSwitcher::getInstance().getNamespace(name);
  Common::SafePtr<PhysicalModel> physModel = PhysicalModelStack::getInstance().getEntryByNamespace(nsp);

  std::string varSetName = "Euler3DCons";
  _varSet.reset((Environment::Factory<ConvectiveVarSet>::getInstance().getProvider(varSetName)->
    create(physModel->getImplementor()->getConvectiveTerm())).d_castTo<Physics::NavierStokes::Euler3DCons>());

  cf_assert(_varSet.isNotNull());

}

//////////////////////////////////////////////////////////////////////////////

void WeakSubOutletEuler3DConsHO::setGhostState(const State& state,
					     State& gstate)
{
  gstate[0] = state[0];
  gstate[1] = state[1];
  gstate[2] = state[2];
  gstate[3] = state[3];
  gstate[4] = _pressure/(_varSet->getModel()->getGamma() - 1.) +
    0.5*(gstate[1]*gstate[1] + gstate[2]*gstate[2] +
	 gstate[3]*gstate[3])/gstate[0];
}

//////////////////////////////////////////////////////////////////////////////

} // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
