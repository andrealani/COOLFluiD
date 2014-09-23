#include "FiniteVolume/FiniteVolume.hh"
#include "StdALEPrepare.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/MeshData.hh"
#include "Framework/SubSystemStatus.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<StdALEPrepare, CellCenterFVMData, FiniteVolumeModule> StdALEPrepareProvider("StdALEPrepare");

//////////////////////////////////////////////////////////////////////////////

void StdALEPrepare::execute()
{
  CFLogDebugMin( "StdALEPrepare::execute()" << "\n");

  DataHandle < Framework::Node*, Framework::GLOBAL > nodes  = socket_nodes.getDataHandle();
  DataHandle<Node*> pastNodes  = socket_pastNodes.getDataHandle();

  // Set the pastNodes = Nodes
  cf_assert (nodes.size() == pastNodes.size());

  if(SubSystemStatusStack::getActive()->isSubIterationFirstStep()){
    for (CFuint i=0; i < nodes.size();++i){
      *pastNodes[i] = *nodes[i];
    }
  }
  else{
    //Go back to the previous situation before continuing
    for (CFuint i=0; i < nodes.size();++i){
      *nodes[i] = *pastNodes[i];
    }
  }


}

//////////////////////////////////////////////////////////////////////////////

vector<SafePtr<BaseDataSocketSink> >
StdALEPrepare::needsSockets()
{
  vector<SafePtr<BaseDataSocketSink> > result;

  result.push_back(&socket_nodes);
  result.push_back(&socket_pastNodes);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD
