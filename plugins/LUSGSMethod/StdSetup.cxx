#include "LUSGSMethod/LUSGSMethod.hh"
#include "StdSetup.hh"
#include "Framework/PhysicalModel.hh"
#include "Framework/MeshData.hh"
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

MethodCommandProvider<StdSetup, LUSGSIteratorData, LUSGSMethodModule> stdSetupProvider("StdSetup");

//////////////////////////////////////////////////////////////////////////////

StdSetup::StdSetup(const std::string& name) :
  LUSGSIteratorCom(name),
  socket_updateCoeff("updateCoeff"),
  socket_pastStates("pastStates"),
  socket_statesSetIdx("statesSetIdx"),
//   socket_bStatesNeighbors("bStatesNeighbors"),
  socket_states("states"),
  socket_pivotLUFactorization("pivotLUFactorization"),
  socket_rhs("rhs")
{
}

//////////////////////////////////////////////////////////////////////////////

StdSetup::~StdSetup()
{
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSource> >
StdSetup::providesSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSource> > result;

  result.push_back(&socket_updateCoeff);
  result.push_back(&socket_pastStates  );
  result.push_back(&socket_statesSetIdx);
  result.push_back(&socket_pivotLUFactorization);
  result.push_back(&socket_rhs);

//   result.push_back(&socket_bStatesNeighbors);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

vector<SafePtr<BaseDataSocketSink> > StdSetup::needsSockets()
{
  vector<SafePtr<BaseDataSocketSink> > result;

  result.push_back(&socket_states);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

void StdSetup::execute()
{
  CFAUTOTRACE;

  // get the states datahandle
  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();

  // get number of states and number of equations
  const CFuint nbStates = states.size();

  // resize updateCoeff
  DataHandle<CFreal> updateCoeff = socket_updateCoeff.getDataHandle();
  updateCoeff.resize(nbStates);

  // States at time step u_n
  DataHandle<State*> pastStates = socket_pastStates.getDataHandle();

  // Resize the past states datahandle and allocate the states
  if(pastStates.size() != nbStates)
  {
    pastStates.resize(nbStates);
    for (CFuint i = 0; i < nbStates; ++i) {
      pastStates[i] = new State();
    }
  }

  // resize statesSetIdx (size = 1)
  DataHandle< CFint > statesSetIdx = socket_statesSetIdx.getDataHandle();
  statesSetIdx.resize(2);
  statesSetIdx[1] = 1; // used to tell to whether to recompute the update coefficient
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace LUSGSMethod

  } // namespace Numerics

} // namespace COOLFluiD
