#include "LUSGSMethod/LUSGSMethod.hh"
#include "LUSGSMethod/StdPrepare.hh"
#include "Framework/State.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace LUSGSMethod {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<StdPrepare, LUSGSIteratorData, LUSGSMethodModule> stdPrepareProvider("StdPrepare");


//////////////////////////////////////////////////////////////////////////////

StdPrepare::StdPrepare(std::string name) :
    LUSGSIteratorCom(name),
    socket_states("states"),
    socket_pastStates("pastStates"),
    socket_statesSetStateIDs("statesSetStateIDs"),
    socket_rhsCurrStatesSet("rhsCurrStatesSet"),
    socket_statesSetIdx("statesSetIdx"),
    m_resAux(CFNULL)
{
}

//////////////////////////////////////////////////////////////////////////////

void StdPrepare::execute()
{
  CFAUTOTRACE;

  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();
  DataHandle<State*> pastStates  = socket_pastStates.getDataHandle();

  // Set Initial States to current states
  //  CFuint nodeID;
  if(SubSystemStatusStack::getActive()->isSubIterationFirstStep()){
    for(CFuint i = 0; i < states.size(); ++i) {
      *(pastStates[i]) = *(states[i]);
    }
  }

  // Get the number of states sets
  DataHandle< vector< CFuint > > statesSetStateIDs = socket_statesSetStateIDs.getDataHandle();
  getMethodData().setNbrStatesSets(statesSetStateIDs.size());

  // Gets the rhs vectors
  DataHandle< CFreal > rhsCurrStatesSet = socket_rhsCurrStatesSet.getDataHandle();

  // Derefence and resize auxiliary rhs variable
  RealVector& resAux = *m_resAux;
  resAux.resize(rhsCurrStatesSet.size());

  // Set states set index to zero
  DataHandle< CFint > statesSetIdx = socket_statesSetIdx.getDataHandle();
  statesSetIdx[0] = -1;
  statesSetIdx[1] = 1; // --> compute update coefficients
}

//////////////////////////////////////////////////////////////////////////////

void StdPrepare::setup()
{
  CFAUTOTRACE;

  // get the auxiliary rhs variable
  m_resAux = getMethodData().getResAux();

}

//////////////////////////////////////////////////////////////////////////////

vector<SafePtr<BaseDataSocketSink> > StdPrepare::needsSockets()
{
  vector<SafePtr<BaseDataSocketSink> > result;

  result.push_back(&socket_states);
  result.push_back(&socket_pastStates);
  result.push_back(&socket_statesSetStateIDs);
  result.push_back(&socket_rhsCurrStatesSet);
  result.push_back(&socket_statesSetIdx);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace LUSGSMethod

  } // namespace Numerics

} // namespace COOLFluiD
