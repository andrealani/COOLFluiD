#include "FluctSplit/FluctSplitSpaceTime.hh"


#include "TwoLayerInitState.hh"
#include "Common/CFLog.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/MeshData.hh"
#include "Framework/NamespaceSwitcher.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<TwoLayerInitState, FluctuationSplitData, FluctSplitSpaceTimeModule> TwoLayerInitStateProvider("TwoLayerInitState");

//////////////////////////////////////////////////////////////////////////////

TwoLayerInitState::TwoLayerInitState(const std::string& name) :
  InitState(name),
  socket_interStates("interStates")
{
}

//////////////////////////////////////////////////////////////////////////////

TwoLayerInitState::~TwoLayerInitState()
{
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
TwoLayerInitState::needsSockets()
{

  std::vector<Common::SafePtr<BaseDataSocketSink> > result = InitState::needsSockets();

  result.push_back(&socket_interStates);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

void TwoLayerInitState::executeOnTrs()
{
  CFAUTOTRACE;

  CFLogDebugMax( "TwoLayerInitState::executeOnTrs() called for TRS: "
  << getCurrentTRS()->getName() << "\n");

  DataHandle<State*> interStates = socket_interStates.getDataHandle();
  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();

  Common::SafePtr<std::vector<CFuint> > statesIdx = getCurrentTRS()->getStatesInTrs();


  if(_inputAdimensionalValues)
  {
    for (CFuint iState = 0; iState < statesIdx->size(); ++iState) {
      const CFuint stateID = (*statesIdx)[iState];
      State *const state = states[stateID];
      _vFunction.evaluate(state->getCoordinates(),(*_input));
      *state = *_inputToUpdateVar->transform(_input);
      *interStates[stateID] = *state;
    }
  }
  else
  {
    State dimState;
    for (CFuint iState = 0; iState < statesIdx->size(); ++iState) {
      const CFuint stateID = (*statesIdx)[iState];
      State *const state = states[stateID];
      _vFunction.evaluate(state->getCoordinates(),(*_input));
      dimState = *_inputToUpdateVar->transform(_input);
      _varSet->setAdimensionalValues(dimState, *state);
      *interStates[stateID] = *state;
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD
