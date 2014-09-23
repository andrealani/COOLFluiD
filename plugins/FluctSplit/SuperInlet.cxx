#include "FluctSplit/FluctSplit.hh"
#include "FluctSplit/SuperInlet.hh"

#include "Framework/MethodCommandProvider.hh"
#include "Framework/MeshData.hh"
#include "Framework/SubSystemStatus.hh"
#include "Framework/NamespaceSwitcher.hh"
#include "Config/PositiveLessThanOne.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Config;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<SuperInlet,
                      FluctuationSplitData,
                      FluctSplitModule>
superInletProvider("SuperInlet");

//////////////////////////////////////////////////////////////////////////////

void SuperInlet::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< std::string >("InputVar","Input variables.");
  options.addConfigOption< std::vector<std::string> >("Vars","Definition of the Variables.");
  options.addConfigOption< std::vector<std::string> >("Def","Definition of the Functions.");
  options.addConfigOption< std::string >("Condition","Definition of the condition to apply the BC.");
  options.addConfigOption< bool >("AdimensionalValues","Input adimensional values.");
  options.addConfigOption< std::vector<CFuint> >
    ("InteractiveVarIDs", "IDs of the variables that will be changed interactively.");
  options.addConfigOption< CFreal, DynamicOption< ValidateOption < PositiveLessThanOne > > >
    ("InteractiveFactor", "Factor to multiply the selected InteractiveVarIDs (should be < 1).");
}

//////////////////////////////////////////////////////////////////////////////

SuperInlet::SuperInlet(const std::string& name) :
  FluctuationSplitCom(name),
  socket_rhs("rhs"),
  socket_isUpdated("isUpdated"),
  socket_states("states"),
  m_checkCondition(false),
  _nbEquationsToSkip(0),
  _varSet(CFNULL),
  _inputToUpdateVar(),
  _input(CFNULL)
{
  addConfigOptionsTo(this);

  m_functions = std::vector<std::string>();
  setParameter("Def",&m_functions);

  m_vars = std::vector<std::string>();
  setParameter("Vars",&m_vars);

  m_conditionStr = "";
  setParameter("Condition",&m_conditionStr);

  _inputVarStr = "Null";
  setParameter("InputVar",&_inputVarStr);

  _inputAdimensionalValues = false;
  setParameter("AdimensionalValues",&_inputAdimensionalValues);

  _interVarIDs = std::vector<CFuint>();
  setParameter("InteractiveVarIDs",&_interVarIDs);

  _interFactor = 1.0;
  setParameter("InteractiveFactor",&_interFactor);
}

//////////////////////////////////////////////////////////////////////////////

SuperInlet::~SuperInlet()
{
}

//////////////////////////////////////////////////////////////////////////////

void SuperInlet::setup()
{
  CFAUTOTRACE;

  FluctuationSplitCom::setup();

  _input = new State();
  const CFuint maxNbStatesInCell = MeshDataStack::getActive()->Statistics().getMaxNbStatesInCell();

  _inputToUpdateVar->setup(maxNbStatesInCell);
  _varSet = getMethodData().getUpdateVar();
}

//////////////////////////////////////////////////////////////////////////////

void SuperInlet::unsetup()
{
  CFAUTOTRACE;

  deletePtr(_input);

  FluctuationSplitCom::unsetup();
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
SuperInlet::needsSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSink> > result;

  result.push_back(&socket_rhs);
  result.push_back(&socket_isUpdated);
  result.push_back(&socket_states);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

void SuperInlet::executeOnTrs()
{
  CFAUTOTRACE;

  SafePtr<TopologicalRegionSet> trs = getCurrentTRS();
  CFLogDebugMax( "SuperInlet::execute() called for TRS: "
  << trs->getName() << "\n");

  DataHandle<bool> isUpdated = socket_isUpdated.getDataHandle();
  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();
  DataHandle< CFreal> rhs = socket_rhs.getDataHandle();

  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  const CFuint dim = PhysicalModelStack::getActive()->getDim();

  bool isUnsteady(false);
  if(SubSystemStatusStack::getActive()->getDT() > 0.) isUnsteady = true;
  State inletState;
  RealVector variables(dim);
  if(isUnsteady) variables.resize(dim+1);

  Common::SafePtr< vector<CFuint> > const statesIdx = trs->getStatesInTrs();

  bool applyBC = true;
  State dimState;

  for (CFuint iState = 0; iState < statesIdx->size(); ++iState) {
    const CFuint stateID = (*statesIdx)[iState];
    State *const state = states[stateID];

    // find if we should apply the BC to this node
    if (m_checkCondition) {
      const CFreal applyBCvalue = m_condition.Eval(state->getCoordinates());
      applyBC = ((!isUpdated[stateID]) && (applyBCvalue > 0.0));
    }
    else {
      applyBC = (!isUpdated[stateID]);
    }

    // apply the BC to this node
    if (applyBC){

      //Set the values of the variables xyz + time
      for (CFuint i = 0; i < state->getCoordinates().size();++i){
        variables[i] = state->getCoordinates()[i];
      }
      if(isUnsteady) variables[state->getCoordinates().size()] = SubSystemStatusStack::getActive()->getCurrentTimeDim();

      //Evaluate the function
      m_vFunction.evaluate(variables,*_input);

      // if some interactive variable IDs are specified, multiply
      // those variables by the given factor
      if (_interVarIDs.size() > 0) {
	for (CFuint i = 0; i < _interVarIDs.size(); ++i) {
	  (*_input)[_interVarIDs[i]] *= _interFactor;
	}
      }

      //Set the state value
      if (_inputAdimensionalValues){
        (*state) = *_inputToUpdateVar->transform(_input);
      }
      else{
        dimState = *_inputToUpdateVar->transform(_input);
        _varSet->setAdimensionalValues(dimState, *state);
      }


      //Reset the rhs
      for (CFuint iEq = 0; iEq < (nbEqs - _nbEquationsToSkip); ++iEq) {
        rhs(stateID, iEq, nbEqs) = 0.;
      }

      isUpdated[stateID] = true; // flagging is important!!!!!
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void SuperInlet::configure ( Config::ConfigArgs& args )
{
  CFAUTOTRACE;

  FluctuationSplitCom::configure(args);

// create the transformer from input to conservative variables
  _updateVarStr = getMethodData().getUpdateVarStr();
  if (_inputVarStr == "Null") {
    _inputVarStr = _updateVarStr;
  }

  std::string name = getMethodData().getNamespace();
  Common::SafePtr<Namespace> nsp = NamespaceSwitcher::getInstance().getNamespace(name);
  Common::SafePtr<PhysicalModel> physModel = PhysicalModelStack::getInstance().getEntryByNamespace(nsp);
  std::string provider = VarSetTransformer::getProviderName
    (physModel->getNameImplementor(),
    _inputVarStr, _updateVarStr);

  _inputToUpdateVar = Environment::Factory<VarSetTransformer>::getInstance().getProvider(provider)
    ->create(physModel->getImplementor());

  cf_assert(_inputToUpdateVar.isNotNull());


  // configure the expression for the boundary values
  m_vFunction.setFunctions(m_functions);
  m_vFunction.setVariables(m_vars);
  try {
    m_vFunction.parse();
  }
  catch (Common::ParserException& e) {
    CFout << e.what() << "\n";
    throw; // retrow the exception to signal the error to the user
  }

  // configure the expression that controls if to apply or not the BC
  try {
    if (!m_conditionStr.empty()) { setCondition(); }
  }
  catch (Common::ParserException& e) {
    CFout << e.what() << "\n";
    throw; // retrow the exception to signal the error to the user
  }
}

//////////////////////////////////////////////////////////////////////////////

void SuperInlet::setCondition()
{
  // some sanity checks
  if (Common::StringOps::getWords(m_conditionStr).size() != 1)
    throwConditionException (FromHere(),"included spaces in definition");

  if (m_vars.empty())
    throwConditionException (FromHere(),"didnt define the variables");

  // join the variables into a single comma separated string
  std::string vars (m_vars[0]);
  for(CFuint i = 1; i < m_vars.size(); i++) {
      vars += ",";
      vars += m_vars[i];
  }

  m_condition.Parse(m_conditionStr,vars);
  if (m_condition.ErrorMsg() != 0) {
    throwConditionException (FromHere(), std::string(m_condition.ErrorMsg()));
  }
  m_checkCondition = true;
}

//////////////////////////////////////////////////////////////////////////////

void SuperInlet::throwConditionException (const Common::CodeLocation& where, const std::string& add)
{
    std::string msg("ParseError in SuperInlet::setCondition(): ");
    msg += add;
    msg += " Conditon function: " + m_conditionStr;
    msg += " Vars: ";
    for(CFuint i = 0; i < m_vars.size(); i++) {
         msg += m_vars[i];
         msg += " ";
    }
    throw Common::ParserException (where,msg);
}

//////////////////////////////////////////////////////////////////////////////


    } // namespace FluctSplit



} // namespace COOLFluiD
