#include "FiniteVolume/FiniteVolume.hh"
#include "BDF2ALEPrepare.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/MeshData.hh"
#include "Framework/SubSystemStatus.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<BDF2ALEPrepare, CellCenterFVMData, FiniteVolumeModule> BDF2ALEPrepareProvider("BDF2ALEPrepare");

//////////////////////////////////////////////////////////////////////////////

void BDF2ALEPrepare::execute()
{
  CFAUTOTRACE;

  CFLog(VERBOSE, "BDF2ALEPrepare::execute() => start\n");
  
  DataHandle < Framework::Node*, Framework::GLOBAL > nodes = socket_nodes.getDataHandle();

  DataHandle<Node*> pastNodes = socket_pastNodes.getDataHandle();

  DataHandle<Node*> pastPastNodes = socket_pastPastNodes.getDataHandle();

  // Set the pastNodes = Nodes
  // Set the pastPastNodes = pastNodes
  cf_assert (nodes.size() == pastNodes.size());
  cf_assert (pastNodes.size() == pastPastNodes.size());

  CFLog(VERBOSE, "BDF2ALEPrepare::execute() => isSubIterationFirstStep [" << 
	SubSystemStatusStack::getActive()->isSubIterationFirstStep() << "]\n");
  
  if(SubSystemStatusStack::getActive()->isSubIterationFirstStep()){
    //Back up the nodes positions before continuing`
    for (CFuint i=0; i < nodes.size();++i){
      *pastPastNodes[i] = *pastNodes[i];
      *pastNodes[i] = *nodes[i];
    }
  }
  else{
    //Go back to the previous situation before continuing
    for (CFuint i=0; i < nodes.size();++i){
      *nodes[i] = *pastNodes[i];
    }
  }
  
  CFLog(VERBOSE, "BDF2ALEPrepare::execute() => end\n");
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
BDF2ALEPrepare::needsSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSink> > result;

  result.push_back(&socket_nodes);
  result.push_back(&socket_pastNodes);
  result.push_back(&socket_pastPastNodes);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD
