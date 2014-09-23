#include "MatrixStabilityMethodWriter/MatrixStabilityMethodWriter.hh"
#include "MatrixStabilityMethodWriter/SetStates.hh"
#include "Framework/State.hh"
#include "Framework/MeshData.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace MatrixStabilityMethodWriter {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<SetStates,
                      MatrixStabilityMethodData,
                      MatrixStabilityMethodWriterModule
                     >
setStatesProvider("SetStates");

//////////////////////////////////////////////////////////////////////////////

SetStates::SetStates(const std::string& name) :
    MatrixStabilityMethodCom(name),
    socket_states("states")
{
}

//////////////////////////////////////////////////////////////////////////////

void SetStates::execute()
{
  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();

  if (getMethodData().getSetAllStatesToZero())
  {
    const CFuint nbrStates = states.size();
    cf_assert(nbrStates == getMethodData().getNbrStates());
    for (CFuint i = 0; i < nbrStates; ++i)
    {
      cf_assert(states[i]->size() == 1);
      *states[i] = 0.0;
    }
  }
  else
  {
    *(states[getMethodData().getStateIdx()]) =
        getMethodData().getSetStateToZero() ? 0.0 : 1.0;
  }
}

//////////////////////////////////////////////////////////////////////////////

vector<SafePtr<BaseDataSocketSink> > SetStates::needsSockets()
{
  vector<SafePtr<BaseDataSocketSink> > result;

  result.push_back(&socket_states);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace MatrixStabilityMethodWriter

  } // namespace Numerics

} // namespace COOLFluiD
