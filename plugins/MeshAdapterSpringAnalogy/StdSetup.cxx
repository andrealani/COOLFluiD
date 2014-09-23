#include "MeshAdapterSpringAnalogy/MeshAdapterSpringAnalogy.hh"
#include "StdSetup.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/SubSystemStatus.hh"
#include "Framework/MeshData.hh"
#include "Framework/PhysicalModel.hh"
#include "Framework/Node.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace MeshAdapterSpringAnalogy {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<StdSetup, SpringAnalogyData, MeshAdapterSpringAnalogyModule> StdSetupProvider("StdSetup");

//////////////////////////////////////////////////////////////////////////////

StdSetup::StdSetup(std::string name) : SpringAnalogyCom(name),
  socket_nodes("nodes"),
  socket_averageVector("averageVector"),
  socket_sumWeight("sumWeight"),
  socket_isMovable("isMovable"),
  socket_nodalDisplacements("nodalDisplacements")
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

  result.push_back(&socket_averageVector);
  result.push_back(&socket_sumWeight);
  result.push_back(&socket_isMovable);
  result.push_back(&socket_nodalDisplacements);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
StdSetup::needsSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSink> > result;

  result.push_back(&socket_nodes);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

void StdSetup::execute()
{
  SubSystemStatusStack::getActive()->setMovingMesh(true);

  DataHandle < Framework::Node*, Framework::GLOBAL > nodes = socket_nodes.getDataHandle();

  const CFuint nbNodes = nodes.size();

   // This is for the SpringAnalogy only

  DataHandle<RealVector> averageVector = socket_averageVector.getDataHandle();
  averageVector.resize(nbNodes);

   for (CFuint i = 0; i < nbNodes; ++i) {
     averageVector[i].resize(PhysicalModelStack::getActive()->getDim());
     averageVector[i] = 0.;
   }

  DataHandle<CFreal> sumWeight = socket_sumWeight.getDataHandle();
  sumWeight.resize(nbNodes);
  sumWeight = 0.;

  DataHandle<bool> isMovable = socket_isMovable.getDataHandle();
  isMovable.resize(nbNodes);
  isMovable = true;

  DataHandle<RealVector> nodalDisplacements = socket_nodalDisplacements.getDataHandle();
  nodalDisplacements.resize(nbNodes);
   for (CFuint i = 0; i < nbNodes; ++i) {
     nodalDisplacements[i].resize(PhysicalModelStack::getActive()->getDim());
     nodalDisplacements[i] = 0.;
   }

  createNodesListInTRS();

  flagBoundaryNodes();
}

//////////////////////////////////////////////////////////////////////////////

void StdSetup::createNodesListInTRS()
{
  std::vector< Common::SafePtr<TopologicalRegionSet> > trs = MeshDataStack::getActive()->getTrsList();
  std::vector< Common::SafePtr<TopologicalRegionSet> >::iterator itrs;
  for (itrs = trs.begin(); itrs != trs.end(); ++itrs) {
    (*itrs)->createNodesList();
  }
}

//////////////////////////////////////////////////////////////////////////////

void StdSetup::flagBoundaryNodes()
{
  DataHandle< bool> isMovable = socket_isMovable.getDataHandle();

  // set the list of faces
  std::vector<Common::SafePtr<TopologicalRegionSet> > trs =
    MeshDataStack::getActive()->getTrsList();

  std::vector< Common::SafePtr<TopologicalRegionSet> >::iterator itrs;
  for (itrs = trs.begin(); itrs != trs.end(); ++itrs) {
    if (((*itrs)->getName() != "InnerCells") && ((*itrs)->getName() != "InnerFaces")) {
      Common::SafePtr<std::vector<CFuint> > nodesInTRS = (*itrs)->getNodesInTrs();
      std::vector<CFuint>::iterator it;
      for (it = nodesInTRS->begin(); it != nodesInTRS->end(); ++it) {
        isMovable[(*it)] = false;
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace MeshAdapterSpringAnalogy

  } // namespace Numerics

} // namespace COOLFluiD
