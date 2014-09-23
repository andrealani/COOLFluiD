#include "FiniteVolume/FiniteVolume.hh"
#include "BDF2ALESetup.hh"
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

MethodCommandProvider<BDF2ALESetup, CellCenterFVMData, FiniteVolumeModule> BDF2ALESetupProvider("BDF2ALESetup");

//////////////////////////////////////////////////////////////////////////////

BDF2ALESetup::BDF2ALESetup(const std::string& name) :
  StdALESetup(name),
  socket_pastPastNodes("pastPastNodes"),
  socket_avNormals("avNormals"),
  socket_normals("normals")
{
}

//////////////////////////////////////////////////////////////////////////////

BDF2ALESetup::~BDF2ALESetup()
{
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSource> >
BDF2ALESetup::providesSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSource> > result =
    StdALESetup::providesSockets();

  result.push_back(&socket_pastPastNodes);
  result.push_back(&socket_avNormals);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

void BDF2ALESetup::execute()
{
  CFAUTOTRACE;

  StdALESetup::execute();

  // Create the storage for the normals averaged between n and n+1 (avNormals)
  DataHandle<CFreal> normals = socket_normals.getDataHandle();

  DataHandle<CFreal> avNormals = socket_avNormals.getDataHandle();
  avNormals.resize(normals.size());

  //
  //   // Create the storage for the normals averaged between n-1 and n (pastAvNormals)
  //   DataHandle<CFreal> pastAvNormals = socket_pastAvNormals.getDataHandle();
  //   pastAvNormals.resize(normals.size());

  // Create the storage for the nodes coordinates at time n-1
  DataHandle < Framework::Node*, Framework::GLOBAL > nodes = socket_nodes.getDataHandle();

  const CFuint nbNodes = nodes.size();

  DataHandle<Node*> pastNodes = socket_pastNodes.getDataHandle();
  DataHandle<Node*> pastPastNodes = socket_pastPastNodes.getDataHandle();
  pastPastNodes.resize(nbNodes);

  for (CFuint i = 0; i < nbNodes; ++i) {
    pastPastNodes[i] = new Node();
    pastPastNodes[i]->setIsOnMesh(false);
    *(pastPastNodes[i]) = *(pastNodes[i]);
  }

}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
BDF2ALESetup::needsSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSink> > result  =
    StdALESetup::needsSockets();

  result.push_back(&socket_normals);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD
