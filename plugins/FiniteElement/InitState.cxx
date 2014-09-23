

#include "Common/CFLog.hh"
#include "Framework/MeshData.hh"
#include "Framework/MethodCommandProvider.hh"
#include "FiniteElement/FiniteElement.hh"
#include "FiniteElement/InitState.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteElement {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<InitState, FiniteElementMethodData, FiniteElementModule> initStateProvider("InitState");

//////////////////////////////////////////////////////////////////////////////

void InitState::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOption< std::vector<std::string> >("Vars","Definition of the Variables.");
   options.addConfigOption< std::vector<std::string> >("Def","Definition of the Functions.");
}

//////////////////////////////////////////////////////////////////////////////

InitState::InitState(const std::string& name) :
  FiniteElementMethodCom(name),
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

void InitState::configure ( Config::ConfigArgs& args )
{
  CFAUTOTRACE;

  FiniteElementMethodCom::configure(args);

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
  CFLogDebugMin( "FiniteElement::InitState::executeOnTrs() called for TRS: " << trs->getName() << "\n");

  DataHandle < Framework::State*, Framework::GLOBAL > states  = socket_states.getDataHandle();

  SafePtr<vector<CFuint> > trsStates = trs->getStatesInTrs();
  std::vector<CFuint>::iterator itd;
  for (itd = trsStates->begin(); itd != trsStates->end(); ++itd) {
    const RealVector& node = states[*itd]->getCoordinates();
      _vFunction.evaluate(node,*states[*itd]);
//     CFout << "Assigning state: " << states[*itd]->getLocalID() << " at coord: " << states[*itd]->getCoordinates() << " the value: " << *(states[*itd]) <<"\n";
  }
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

    } // namespace FiniteElement

  } // namespace Numerics

} // namespace COOLFluiD
