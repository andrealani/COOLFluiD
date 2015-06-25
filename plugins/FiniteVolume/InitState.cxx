#include "FiniteVolume/FiniteVolume.hh"


#include "InitState.hh"
#include "Common/CFLog.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Common/BadValueException.hh"
#include "Framework/NamespaceSwitcher.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<InitState, CellCenterFVMData, FiniteVolumeModule>
initStateProvider("InitState");

//////////////////////////////////////////////////////////////////////////////

void InitState::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< std::vector<std::string> >("Vars","Definition of the Variables.");
  options.addConfigOption< std::vector<std::string> >("Def","Definition of the Functions.");
  options.addConfigOption< std::string >("InputVar","Input variables.");
  options.addConfigOption< bool >("AdimensionalValues","Flag to input adimensional values.");
}

//////////////////////////////////////////////////////////////////////////////

InitState::InitState(const std::string& name) :
  CellCenterFVMCom(name),
  //socket_states("states"),
  socket_normals("normals"),
  socket_states("states"),
  socket_gstates("gstates"),
  socket_nodes("nodes"),
  _varSet(CFNULL),
  _inputToUpdateVar(),
  _input()
{
   addConfigOptionsTo(this);
  _functions = std::vector<std::string>();
   setParameter("Def",&_functions);

  _vars = std::vector<std::string>();
   setParameter("Vars",&_vars);

   _inputVarStr = "Null";
   setParameter("InputVar",&_inputVarStr);

  _inputAdimensionalValues = false;
   setParameter("AdimensionalValues",&_inputAdimensionalValues);

}

//////////////////////////////////////////////////////////////////////////////

InitState::~InitState()
{
}

//////////////////////////////////////////////////////////////////////////////

void InitState::executeOnTrs()
{
  SafePtr<TopologicalRegionSet> trs = getCurrentTRS();
  CFLogDebugMax( "InitState::executeOnTrs() called for TRS: "
  << trs->getName() << "\n");

  if (trs->getName() != "InnerFaces") {
    throw BadValueException (FromHere(),"InitState not applied to InnerFaces!!!");
  }

  // this cannot be used for FV boundary faces because
  // ghost state and inner state could have the same local ID !!!
  SafePtr<vector<CFuint> > trsStates = trs->getStatesInTrs();
  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();
  
  std::vector<CFuint>::iterator itd;
  if(_inputAdimensionalValues)
  {
    for (itd = trsStates->begin(); itd != trsStates->end(); ++itd) {
      State* const currState = states[(*itd)];
      _vFunction.evaluate(currState->getCoordinates(), *_input);
      *currState = *_inputToUpdateVar->transform(_input);
    }
  }
  else
  {
    State dimState;
    for (itd = trsStates->begin(); itd != trsStates->end(); ++itd) {
      State* const currState = states[(*itd)];
      _vFunction.evaluate(currState->getCoordinates(), *_input);
      dimState = *_inputToUpdateVar->transform(_input);
      _varSet->setAdimensionalValues(dimState, *currState);
    }
  }

}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
InitState::needsSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSink> > result;

  //result.push_back(&socket_states);
  result.push_back(&socket_normals);
  result.push_back(&socket_states);
  result.push_back(&socket_gstates);
  result.push_back(&socket_nodes);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

void InitState::setup()
{
  CFAUTOTRACE;

  CellCenterFVMCom::setup();

  _input = new State();

  const CFuint maxNbStatesInCell = MeshDataStack::getActive()->Statistics().getMaxNbStatesInCell();
  _inputToUpdateVar->setup(maxNbStatesInCell);

  _varSet = getMethodData().getUpdateVar();
}

//////////////////////////////////////////////////////////////////////////////

void InitState::unsetup()
{
  CFAUTOTRACE;

  deletePtr(_input);

  CellCenterFVMCom::unsetup();
}

//////////////////////////////////////////////////////////////////////////////

void InitState::configure ( Config::ConfigArgs& args )
{
  CFAUTOTRACE;

  CellCenterFVMCom::configure(args);

  const std::string name = getMethodData().getNamespace();
  Common::SafePtr<Namespace> nsp = NamespaceSwitcher::getInstance
    (SubSystemStatusStack::getCurrentName()).getNamespace(name);
  Common::SafePtr<PhysicalModel> physModel =
    PhysicalModelStack::getInstance().getEntryByNamespace(nsp);
  
  // create the transformer from input to update variables
  if (_inputVarStr == "Null") {
    _inputVarStr = getMethodData().getUpdateVarStr();
  }
  
  const std::string provider = VarSetTransformer::getProviderName
    (physModel->getConvectiveName(), _inputVarStr, getMethodData().getUpdateVarStr());
  
  _inputToUpdateVar =
    Environment::Factory<VarSetTransformer>::getInstance().getProvider(provider)->
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

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD
