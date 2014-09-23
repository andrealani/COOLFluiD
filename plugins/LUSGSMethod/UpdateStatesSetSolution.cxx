#include "Framework/PhysicalModel.hh"
#include "Framework/State.hh"
#include "LUSGSMethod/LUSGSMethod.hh"
#include "LUSGSMethod/UpdateStatesSetSolution.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace LUSGSMethod {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<UpdateStatesSetSolution, LUSGSIteratorData, LUSGSMethodModule>
    updateStatesSetSolutionProvider("UpdateStatesSetSolution");

//////////////////////////////////////////////////////////////////////////////

UpdateStatesSetSolution::UpdateStatesSetSolution(std::string name) :
  LUSGSIteratorCom(name),
  socket_rhsCurrStatesSet("rhsCurrStatesSet"),
  socket_statesSetIdx("statesSetIdx"),
  socket_states("states"),
  socket_statesSetStateIDs("statesSetStateIDs"),
  socket_isStatesSetParUpdatable("isStatesSetParUpdatable"),
  m_nbrEqs()
{
}

//////////////////////////////////////////////////////////////////////////////

void UpdateStatesSetSolution::execute()
{
  // Gets current states set index
  DataHandle< CFint > statesSetIdx = socket_statesSetIdx.getDataHandle();
  const CFuint currStatesSetIdx = statesSetIdx[0];

  // Gets isStatesSetParUpdatable datahandle
  DataHandle< bool > isStatesSetParUpdatable = socket_isStatesSetParUpdatable.getDataHandle();

  if (isStatesSetParUpdatable[currStatesSetIdx])
  {
    // Gets the rhs vectors datahandle
    DataHandle< CFreal > rhsCurrStatesSet = socket_rhsCurrStatesSet.getDataHandle();

    // rhsCurrStatesSet is the temporary placeholder for the dU
    DataHandle<CFreal>& dU = rhsCurrStatesSet;

    // Gets state datahandle
    DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();

    // Gets the current states ID
    DataHandle< vector< CFuint > > statesSetStateIDs = socket_statesSetStateIDs.getDataHandle();
    const vector< CFuint >& currStatesIDs = statesSetStateIDs[currStatesSetIdx];
    const CFuint currNbrStates = currStatesIDs.size();

    // Updates the current states set
    CFuint resIdx = 0;
    for (CFuint iState = 0; iState < currNbrStates; ++iState)
    {
      const CFuint stateID = currStatesIDs[iState];
      State& currState = *states[stateID];
      for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq, ++resIdx)
      {
        currState[iEq] += dU[resIdx];
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void UpdateStatesSetSolution::setup()
{
  CFAUTOTRACE;

  // call setup of parent class
  LUSGSIteratorCom::setup();

  // get the number of equations
  m_nbrEqs = PhysicalModelStack::getActive()->getNbEq();
}

//////////////////////////////////////////////////////////////////////////////

vector<SafePtr<BaseDataSocketSink> > UpdateStatesSetSolution::needsSockets()
{
  vector<SafePtr<BaseDataSocketSink> > result;

  result.push_back(&socket_statesSetIdx);
  result.push_back(&socket_rhsCurrStatesSet);
  result.push_back(&socket_states);
  result.push_back(&socket_statesSetStateIDs);
  result.push_back(&socket_isStatesSetParUpdatable);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace LUSGSMethod

  } // namespace Numerics

} // namespace COOLFluiD
