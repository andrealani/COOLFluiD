#include "Common/CFLog.hh"

#include "Framework/MeshData.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/SubSystemStatus.hh"
#include "Framework/NamespaceSwitcher.hh"

#include "UFEM/UFEM.hh"
#include "UFEM/RestartState.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace UFEM {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<RestartState, UFEMSolverData, UFEMPlugin> restartStateProvider("RestartState");

//////////////////////////////////////////////////////////////////////////////

void RestartState::defineConfigOptions(Config::OptionList& options)
{
}

//////////////////////////////////////////////////////////////////////////////

RestartState::RestartState(const std::string& name) :
  UFEMSolverCom(name),
  socket_states("states"),
  socket_pastStates("pastStates"),
  socket_pastpastStates("pastpastStates")
{
  CFAUTOTRACE;
  addConfigOptionsTo(this);
}

//////////////////////////////////////////////////////////////////////////////

RestartState::~RestartState()
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

void RestartState::configure ( Config::ConfigArgs& args )
{
  CFAUTOTRACE;
  UFEMSolverCom::configure(args);
}

//////////////////////////////////////////////////////////////////////////////

void RestartState::executeOnTrs()
{
  CFAUTOTRACE;

  SafePtr<TopologicalRegionSet> trs = getCurrentTRS();

  CFLogDebugMin( "UFEM::RestartState::executeOnTrs() called for TRS: " << trs->getName() << "\n");

cout << "UFEM::RestartState::executeOnTrs" << endl << flush;

  DataHandle<State*,GLOBAL > states  = socket_states.getDataHandle();
  DataHandle<State*> pastStates = socket_pastStates.getDataHandle();
  DataHandle<State*> pastpastStates = socket_pastpastStates.getDataHandle();

  SafePtr<vector<CFuint> > trsStates = trs->getStatesInTrs();

  for (std::vector<CFuint>::iterator itd = trsStates->begin(); itd != trsStates->end(); ++itd)
  {
    *pastStates[*itd]=*states[*itd];
    *pastpastStates[*itd]=*states[*itd];
/*
    cout << *states[*itd] << "\n" << flush;
    cout << *pastStates[*itd] << "\n" << flush;
    cout << *pastpastStates[*itd] << "\n" << flush;
*/
  }
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
RestartState::needsSockets()
{
  CFAUTOTRACE;
  std::vector<Common::SafePtr<BaseDataSocketSink> > result;

  result.push_back(&socket_states);
  result.push_back(&socket_pastStates);
  result.push_back(&socket_pastpastStates);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace UFEM

} // namespace COOLFluiD
