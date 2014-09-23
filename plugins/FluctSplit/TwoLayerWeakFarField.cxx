#include "FluctSplit/FluctSplitSpaceTime.hh"
#include "TwoLayerWeakFarField.hh"
#include "InwardNormalsData.hh"
#include "Framework/PhysicalModel.hh"
#include "Framework/SubSystemStatus.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/NamespaceSwitcher.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<TwoLayerWeakFarField, FluctuationSplitData, FluctSplitSpaceTimeModule> TwoLayerWeakFarFieldProvider("TwoLayerWeakFarField");

//////////////////////////////////////////////////////////////////////////////

void TwoLayerWeakFarField::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOption< std::string >("InputVar","Input variables.");
   options.addConfigOption< std::vector<std::string> >("Vars","Definition of the Variables.");
   options.addConfigOption< std::vector<std::string> >("Def","Definition of the Functions.");
}

//////////////////////////////////////////////////////////////////////////////

TwoLayerWeakFarField::TwoLayerWeakFarField(const std::string& name) :
  TwoLayerWeakBC2D(name),
   _inputToUpdateVar(),
  _vFunction(),
  _variables(),
  _input(CFNULL),
  _dimState(CFNULL)
{
  addConfigOptionsTo(this);

  _functions = std::vector<std::string>();
  setParameter("Def",&_functions);

  _vars = std::vector<std::string>();
  setParameter("Vars",&_vars);

  _inputVarStr = "Null";
  setParameter("InputVar",&_inputVarStr);
}

//////////////////////////////////////////////////////////////////////////////

TwoLayerWeakFarField::~TwoLayerWeakFarField()
{
}

//////////////////////////////////////////////////////////////////////////////

void TwoLayerWeakFarField::setup()
{
  TwoLayerWeakBC2D::setup();

  const CFuint maxNbStatesInCell = MeshDataStack::getActive()->Statistics().getMaxNbStatesInCell();

  _variables.resize(PhysicalModelStack::getActive()->getDim() + 1);
  _inputToUpdateVar->setup(maxNbStatesInCell);

  _dimState = new State();
  _input = new State();
}

//////////////////////////////////////////////////////////////////////////////

void TwoLayerWeakFarField::unsetup()
{
  deletePtr(_dimState);
  deletePtr(_input);

  TwoLayerWeakBC2D::unsetup();
}


//////////////////////////////////////////////////////////////////////////////

void TwoLayerWeakFarField::setGhostState(const State& state,
                                                    State& gstate)
{
  // set the values in the ghost state (whose position coincides
  // with the boundary states)
  const CFuint dim = PhysicalModelStack::getActive()->getDim();
  const RealVector& temp = state.getCoordinates();

  for (CFuint i = 0; i < dim; ++i){
    _variables[i] = temp[i];
  }
  _variables[temp.size()] = SubSystemStatusStack::getActive()->getCurrentTimeDim();

  // evaluate first the input states in input variables
  _vFunction.evaluate(_variables, *_input);

  // transform to conservative variables
  *_dimState = *_inputToUpdateVar->transform(_input);

  // Adimensionalize the value
  Common::SafePtr<Framework::ConvectiveVarSet> updateVarSet = getMethodData().getUpdateVar();
  updateVarSet->setAdimensionalValues(*_dimState, gstate);

 }

//////////////////////////////////////////////////////////////////////////////

void TwoLayerWeakFarField::configure ( Config::ConfigArgs& args )
{
  TwoLayerWeakBC2D::configure(args);

  std::string name = getMethodData().getNamespace();
  Common::SafePtr<Namespace> nsp = NamespaceSwitcher::getInstance().getNamespace(name);
  Common::SafePtr<PhysicalModel> physModel = PhysicalModelStack::getInstance().getEntryByNamespace(nsp);

  _updateVarStr = getMethodData().getUpdateVarStr();

  // create the transformer from input to update variables
  if (_inputVarStr == "Null") {
    _inputVarStr = _updateVarStr;
  }
  std::string provider = VarSetTransformer::getProviderName
    (physModel->getNameImplementor(),
     _inputVarStr, _updateVarStr);

  _inputToUpdateVar = Environment::Factory<VarSetTransformer>::getInstance().getProvider(provider)->
    create(physModel->getImplementor());

  cf_assert(_inputToUpdateVar.isNotNull());

  _vFunction.setFunctions(_functions);
  _vFunction.setVariables(_vars);
  try {
    _vFunction.parse();
  }
  catch (Common::ParserException& e) {
    CFout << e.what() << "\n";
    throw; // retrow the exception to signal the error to the user
  }
}

//////////////////////////////////////////////////////////////////////////////

} // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
