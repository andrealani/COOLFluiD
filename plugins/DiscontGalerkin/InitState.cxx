

#include "Common/CFLog.hh"
#include "Framework/MethodCommandProvider.hh"
#include "DiscontGalerkin/DiscontGalerkin.hh"
#include "DiscontGalerkin/InitState.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace DiscontGalerkin {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<InitState, DiscontGalerkinSolverData, DiscontGalerkinModule>
initStateProvider("InitState");

//////////////////////////////////////////////////////////////////////////////

void InitState::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< std::vector<std::string> >("Vars","Definition of the Variables.");
  options.addConfigOption< std::vector<std::string> >("Def","Definition of the Functions.");
}

//////////////////////////////////////////////////////////////////////////////

InitState::InitState(const std::string& name) :
  DiscontGalerkinSolverCom(name),
  socket_states("states")
{
   addConfigOptionsTo(this);
  _functions = std::vector<std::string>();
   setParameter("Def",&_functions);

  _vars = std::vector<std::string>();
   setParameter("Vars",&_vars);
}

//////////////////////////////////////////////////////////////////////////////

InitState::~InitState()
{
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
InitState::needsSockets()
{

  std::vector<Common::SafePtr<BaseDataSocketSink> > result;

  result.push_back(&socket_states);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

void InitState::setup()
{

  CFAUTOTRACE;

  DiscontGalerkinSolverCom::setup();
}

//////////////////////////////////////////////////////////////////////////////

void InitState::unsetup()
{
  CFAUTOTRACE;

  DiscontGalerkinSolverCom::unsetup();
}

//////////////////////////////////////////////////////////////////////////////

void InitState::configure ( Config::ConfigArgs& args )
{
  CFAUTOTRACE;

  DiscontGalerkinSolverCom::configure(args);

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

void InitState::executeOnTrs()
{
  CFAUTOTRACE;

  SafePtr<TopologicalRegionSet> trs = getCurrentTRS();
  CFLogDebugMin( "InitState::executeOnTrs() called for TRS: " << trs->getName() << "\n");
  Common::SafePtr<std::vector<CFuint> > statesIdx = trs->getStatesInTrs();
  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();
    for (CFuint iState = 0; iState < statesIdx->size(); ++iState) 
    {
      State *const state = states[(*statesIdx)[iState]];
      _vFunction.evaluate(state->getCoordinates(),*state);
    }
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace DiscontGalerkin

} // namespace COOLFluiD
