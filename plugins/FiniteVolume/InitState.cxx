#include "FiniteVolume/FiniteVolume.hh"
#include "FiniteVolume/InitState.hh"
#include "Common/CFLog.hh"
#include "Common/BadValueException.hh"
#include "Common/PEFunctions.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/NamespaceSwitcher.hh"
#include "Environment/FileHandlerInput.hh"
#include "Environment/SingleBehaviorFactory.hh"

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
  options.addConfigOption< std::string >("DefFileName","Name of file where Functions are defined.");
  options.addConfigOption< std::string >("InputVar","Input variables.");
  options.addConfigOption< bool >("AdimensionalValues","Flag to input adimensional values.");
}

//////////////////////////////////////////////////////////////////////////////

InitState::InitState(const std::string& name) :
  CellCenterFVMCom(name),
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

  _functionsFileName = "";
  setParameter("DefFileName",&_functionsFileName);
}
      
//////////////////////////////////////////////////////////////////////////////

InitState::~InitState()
{
}

//////////////////////////////////////////////////////////////////////////////

void InitState::executeOnTrs()
{  
  CFLog(VERBOSE, "InitState::executeOnTrs() => START\n");
  
  SafePtr<TopologicalRegionSet> trs = getCurrentTRS();
  CFLogDebugMax( "InitState::executeOnTrs() called for TRS: " << trs->getName() << "\n");
  
  if (trs->getName() != "InnerFaces") {
    throw BadValueException (FromHere(),"InitState not applied to InnerFaces!!!");
  }
  
  // this cannot be used for FV boundary faces because
  // ghost state and inner state could have the same local ID !!!
  SafePtr<vector<CFuint> > trsStates = trs->getStatesInTrs();
  cf_assert(trsStates->size() > 0);
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
      cf_assert((*itd) < states.size());
      State* const currState = states[(*itd)];
      _vFunction.evaluate(currState->getCoordinates(), *_input);
      dimState = *_inputToUpdateVar->transform(_input);
      _varSet->setAdimensionalValues(dimState, *currState);
    }
  }
  
  CFLog(VERBOSE, "InitState::executeOnTrs() => END\n");
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

  // number of functions detected must match the total number of equations
  cf_assert(_functions.size() == PhysicalModelStack::getActive()->getNbEq());
  
  _input = new State();

  const CFuint maxNbStatesInCell =
    MeshDataStack::getActive()->Statistics().getMaxNbStatesInCell();
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
  
  CFLog(VERBOSE, "InitState::configure() => provider for _inputToUpdateVar is "<< provider << "\n");
  
  _inputToUpdateVar =
    FACTORY_GET_PROVIDER(getFactoryRegistry(), VarSetTransformer, provider)->
    create(physModel->getImplementor());
  
  if (_functionsFileName != "" && _functions.size() == 0) {
    runSerial<void, InitState, &InitState::readFunctionsFile>(this, name);
  }
  
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

void InitState::readFunctionsFile()
{
  boost::filesystem::path fpath(_functionsFileName);

  CFLog(INFO, "InitState::readFunctionsFile() " << _functionsFileName << "\n");
  
  SelfRegistPtr<Environment::FileHandlerInput> fhandle = 
    Environment::SingleBehaviorFactory<Environment::FileHandlerInput>::getInstance().create();
  /* ifstream& inputFile = fhandle->open(fpath);
     CFuint nbVars  = 0; 
   inputFile >> nbVars;
   cf_assert(nbVars == nbEqs);
   _functions.resize(nbVars);
   
   // this assumes that the functions are written each in one line, one per variable
   for (CFuint iVar = 0; iVar < nbVars; ++iVar) {
   inputFile >> _functions[iVar];
   }
   fhandle->close();
*/
  
  string line = "";
  CFuint nbLines = 0;
  ifstream& inputFile1 = fhandle->open(fpath);

  if (inputFile1.is_open()) {
    while (!inputFile1.eof()) {
      getline (inputFile1,line);
      ++nbLines;
    }
  }
  fhandle->close();
    
  ifstream& inputFile2 = fhandle->open(fpath);
  string functionName = "";
  for (CFuint iLine = 0; iLine < nbLines-1; ++iLine) {
    getline(inputFile2, line);
    // if you find "\", remove it and update the string
    if (line.find('\\') != std::string::npos) {
      CFLog(VERBOSE, "BEFORE InitState::readFunctionsFile() => " << line << "\n");
      line.erase(remove(line.begin(), line.end(), '\\'), line.end());
      CFLog(VERBOSE, "AFTER  InitState::readFunctionsFile() => " << line << "\n");
      StringOps::trim(line);
      functionName += line;
    }
    else {
      StringOps::trim(line);
      functionName += line;
      _functions.push_back(functionName);
      CFLog(VERBOSE, "InitState::readFunctionsFile() => [" << functionName << "]\n");
      functionName = "";
    }
  }
  
  fhandle->close();
 }

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD
