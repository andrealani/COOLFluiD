#include "FiniteVolume/FiniteVolume.hh"
#include "BDF2ALEUnSetup.hh"
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

MethodCommandProvider<BDF2ALEUnSetup, CellCenterFVMData, FiniteVolumeModule> BDF2ALEUnSetupProvider("BDF2ALEUnSetup");

//////////////////////////////////////////////////////////////////////////////

BDF2ALEUnSetup::BDF2ALEUnSetup(const std::string& name) :
  StdALEUnSetup(name),
  socket_pastPastNodes("pastPastNodes")
{
}

//////////////////////////////////////////////////////////////////////////////

BDF2ALEUnSetup::~BDF2ALEUnSetup()
{
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
BDF2ALEUnSetup::needsSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSink> > result = StdALEUnSetup::needsSockets();

  result.push_back(&socket_pastPastNodes);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

void BDF2ALEUnSetup::execute()
{
  CFAUTOTRACE;

  DataHandle<Node*> pastPastNodes = socket_pastPastNodes.getDataHandle();

  // delete the nodes
  for(CFuint i = 0; i < pastPastNodes.size(); i++) {
    deletePtr(pastPastNodes[i]);
  }

}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD
