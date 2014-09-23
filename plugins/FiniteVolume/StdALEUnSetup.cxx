#include "FiniteVolume/FiniteVolume.hh"
#include "StdALEUnSetup.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/MeshData.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<StdALEUnSetup, CellCenterFVMData, FiniteVolumeModule> StdALEUnSetupProvider("StdALEUnSetup");

//////////////////////////////////////////////////////////////////////////////

StdALEUnSetup::StdALEUnSetup(const std::string& name) :
  CellCenterFVMCom(name),
  socket_pastNodes("pastNodes"),
  socket_futureNodes("futureNodes")
{
}

//////////////////////////////////////////////////////////////////////////////

StdALEUnSetup::~StdALEUnSetup()
{
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
StdALEUnSetup::needsSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSink> > result;

  result.push_back(&socket_pastNodes);
  result.push_back(&socket_futureNodes);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

void StdALEUnSetup::execute()
{
  CFAUTOTRACE;

  DataHandle<Node*> pastNodes = socket_pastNodes.getDataHandle();

  // delete the nodes
  for(CFuint i = 0; i < pastNodes.size(); i++) {
    deletePtr(pastNodes[i]);
  }

  DataHandle<Node*> futureNodes = socket_futureNodes.getDataHandle();

  // delete the nodes
  for(CFuint i = 0; i < futureNodes.size(); i++) {
    deletePtr(futureNodes[i]);
  }
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD
