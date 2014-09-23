#include "FluctSplit/FluctSplitSpaceTime.hh"
#include "TwoLayerSuperInlet.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/SubSystemStatus.hh"
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

MethodCommandProvider<TwoLayerSuperInlet, FluctuationSplitData, FluctSplitSpaceTimeModule> TwoLayerSuperInletProvider("TwoLayerSuperInlet");

//////////////////////////////////////////////////////////////////////////////

TwoLayerSuperInlet::TwoLayerSuperInlet(const std::string& name) :
  SuperInlet(name),
  socket_interRhs("interRhs"),
  socket_interStates("interStates")
{
}

//////////////////////////////////////////////////////////////////////////////

TwoLayerSuperInlet::~TwoLayerSuperInlet()
{
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
TwoLayerSuperInlet::needsSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSink> > result = SuperInlet::needsSockets();

  result.push_back(&socket_interRhs);
  result.push_back(&socket_interStates);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

void TwoLayerSuperInlet::executeOnTrs()
{
  CFAUTOTRACE;

  SafePtr<TopologicalRegionSet> trs = getCurrentTRS();
  CFLogDebugMax( "SuperInlet::execute() called for TRS: "
  << trs->getName() << "\n");

  DataHandle<bool> isUpdated = socket_isUpdated.getDataHandle();
  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();
  DataHandle< CFreal> rhs = socket_rhs.getDataHandle();
  DataHandle<CFreal> interRhs = socket_interRhs.getDataHandle();
  DataHandle<State*> interStates = socket_interStates.getDataHandle();

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
    State *const intState = interStates[stateID];
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

      //Set the state value
      if (_inputAdimensionalValues){
        (*state) = *_inputToUpdateVar->transform(_input);
      }
      else{
        dimState = *_inputToUpdateVar->transform(_input);
        _varSet->setAdimensionalValues(dimState, *state);
      }

      (*intState) = (*state);
      //Reset the rhs
      for (CFuint iEq = 0; iEq < (nbEqs - _nbEquationsToSkip); ++iEq) {
        rhs(stateID, iEq, nbEqs) = 0.;
        interRhs(stateID, iEq, nbEqs) = 0.;
      }

      isUpdated[stateID] = true; // flagging is important!!!!!
    }
  }
}

//////////////////////////////////////////////////////////////////////////////


    } // namespace FluctSplit



} // namespace COOLFluiD
