#include "LUSGSMethod/LUSGSMethod.hh"
#include "LUSGSMethod/UpdateStatesSetIndex.hh"


//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace LUSGSMethod {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<UpdateStatesSetIndex, LUSGSIteratorData, LUSGSMethodModule> updateStatesSetIndexProvider("UpdateStatesSetIndex");


//////////////////////////////////////////////////////////////////////////////

UpdateStatesSetIndex::UpdateStatesSetIndex(std::string name) :
    LUSGSIteratorCom(name),
    socket_statesSetIdx("statesSetIdx")
{
}

//////////////////////////////////////////////////////////////////////////////

void UpdateStatesSetIndex::execute()
{
  // Get state index datahandle
  DataHandle< CFint > statesSetIdx = socket_statesSetIdx.getDataHandle();

  if (getMethodData().isForwardSweep())
  {
    ++statesSetIdx[0];
    if (getMethodData().getNbrStatesSets() <= statesSetIdx[0])
    {
      getMethodData().setStopSweep(true);
      statesSetIdx[0] = getMethodData().getNbrStatesSets();
    }
  }
  else
  {
    --statesSetIdx[0];
    if (-1 >= statesSetIdx[0])
    {
      getMethodData().setStopSweep(true);
      statesSetIdx[0] = -1;
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

vector<SafePtr<BaseDataSocketSink> > UpdateStatesSetIndex::needsSockets()
{
  vector<SafePtr<BaseDataSocketSink> > result;

  result.push_back(&socket_statesSetIdx);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace LUSGSMethod

  } // namespace Numerics

} // namespace COOLFluiD
