#include "FluctSplit/FluctSplitSpaceTimeNavierStokes.hh"
#include "UnsteadyWeakSubOutletEuler2DCons.hh"
#include "InwardNormalsData.hh"
#include "Framework/PhysicalModel.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/SubSystemStatus.hh"
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

MethodCommandProvider<UnsteadyWeakSubOutletEuler2DCons, FluctuationSplitData, FluctSplitSpaceTimeNavierStokesModule> unsteadyWeakSubOutletEuler2DConsProvider("UnsteadyWeakSubOutletEuler2DCons");

//////////////////////////////////////////////////////////////////////////////

void UnsteadyWeakSubOutletEuler2DCons::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOption< std::vector<std::string> >("Vars","Definition of the Variables.");
   options.addConfigOption< std::vector<std::string> >("Def","Definition of the Function defining the static pressure.");
}

//////////////////////////////////////////////////////////////////////////////

UnsteadyWeakSubOutletEuler2DCons::UnsteadyWeakSubOutletEuler2DCons(const std::string& name) :
  WeakBC2D(name),
  _varSet(),
  _vFunction()
{
   addConfigOptionsTo(this);
  _functions = std::vector<std::string>();
   setParameter("Def",&_functions);

  _vars = std::vector<std::string>();
   setParameter("Vars",&_vars);

  _pressure = 1.0;
  //addConfigOption("P","static pressure",&_pressure);
}

//////////////////////////////////////////////////////////////////////////////

UnsteadyWeakSubOutletEuler2DCons::~UnsteadyWeakSubOutletEuler2DCons()
{
}

//////////////////////////////////////////////////////////////////////////////

void UnsteadyWeakSubOutletEuler2DCons::setup()
{
   WeakBC2D::setup();

  _varSet->setup();
}

//////////////////////////////////////////////////////////////////////////////

void UnsteadyWeakSubOutletEuler2DCons::configure ( Config::ConfigArgs& args )
{
  WeakBC2D::configure(args);

  _vFunction.setFunctions(_functions);
  _vFunction.setVariables(_vars);
  try {
    _vFunction.parse();
  }
  catch (Common::ParserException& e) {
    CFout << e.what() << "\n";
    throw; // retrow the exception to signal the error to the user
  }

  const std::string name = getMethodData().getNamespace();
  Common::SafePtr<Namespace> nsp = NamespaceSwitcher::getInstance
    (SubSystemStatusStack::getCurrentName()).getNamespace(name);
  Common::SafePtr<PhysicalModel> physModel = PhysicalModelStack::getInstance().getEntryByNamespace(nsp);

  const std::string varSetName = "Euler2DCons";
  _varSet.reset((Environment::Factory<ConvectiveVarSet>::getInstance().getProvider(varSetName)->
    create(physModel->getImplementor()->getConvectiveTerm())).d_castTo<Physics::NavierStokes::Euler2DCons>());

  cf_assert(_varSet.isNotNull());

}

//////////////////////////////////////////////////////////////////////////////

void UnsteadyWeakSubOutletEuler2DCons::setGhostState(const State& state,
                                             State& gstate)
{

  // evaluate the value of the pressure
  RealVector variables(PhysicalModelStack::getActive()->getDim()+1);
  RealVector temp(PhysicalModelStack::getActive()->getDim());

  temp = (state.getCoordinates());

  for (CFuint i = 0; i < temp.size();++i){
    variables[i] = temp[i];
    }
  variables[temp.size()] = SubSystemStatusStack::getActive()->getCurrentTimeDim();

  RealVector result(1);
  _vFunction.evaluate(variables,result);
  cf_assert (result.size() == 1);

  _pressure = result[0]/_varSet->getModel()->getPressRef();

  gstate[0] = state[0];
  gstate[1] = state[1];
  gstate[2] = state[2];
  gstate[3] = _pressure/(_varSet->getModel()->getGamma() - 1.) +
    0.5*(gstate[1]*gstate[1] + gstate[2]*gstate[2])/gstate[0];
}

//////////////////////////////////////////////////////////////////////////////

} // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
