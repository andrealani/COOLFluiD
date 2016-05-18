#include "FiniteVolume/FiniteVolume.hh"
#include "StdALESetup.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/MeshData.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<StdALESetup, CellCenterFVMData, FiniteVolumeModule> StdALESetupProvider("StdALESetup");

//////////////////////////////////////////////////////////////////////////////

StdALESetup::StdALESetup(const std::string& name) :
  CellCenterFVMCom(name),
  socket_pastNodes("pastNodes"),
  socket_futureNodes("futureNodes"),
  socket_pastVolumes("pastVolumes"),
  socket_nodes("nodes"),
  socket_states("states"),
  socket_volumes("volumes")

{
}

//////////////////////////////////////////////////////////////////////////////

StdALESetup::~StdALESetup()
{
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSource> >
StdALESetup::providesSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSource> > result;

  result.push_back(&socket_pastNodes);
  result.push_back(&socket_futureNodes);
  result.push_back(&socket_pastVolumes);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

void StdALESetup::execute()
{
  CFAUTOTRACE;

  // Create the storage for the nodes coordinates at time n (pastNodes) and n+1(futureNodes)
  DataHandle < Framework::Node*, Framework::GLOBAL > nodes = socket_nodes.getDataHandle();
  DataHandle<Node*> pastNodes = socket_pastNodes.getDataHandle();
  DataHandle<Node*> futureNodes = socket_futureNodes.getDataHandle();
  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();
  DataHandle<CFreal> pastVolumes = socket_pastVolumes.getDataHandle();
  DataHandle<CFreal> volumes = socket_volumes.getDataHandle();
  const CFuint nbNodes = nodes.size();

  bool foundPastNodes = true;
  if(pastNodes.size() != nbNodes)
  {
    foundPastNodes = false;
    pastNodes.resize(nbNodes);
    for (CFuint i = 0; i < nbNodes; ++i) {
      pastNodes[i] = new Node();
      pastNodes[i]->setIsOnMesh(false);
      // initialization of future nodes
      *pastNodes[i] = *nodes[i];
    }
  }
  
  futureNodes.resize(nbNodes);
  for (CFuint i = 0; i < nbNodes; ++i) {
    futureNodes[i] = new Node();
    futureNodes[i]->setIsOnMesh(false);
    // initialization of future nodes
    *futureNodes[i] = *nodes[i];
  }
  
  // Create the storage for the volume at time n (pastVolume)
  pastVolumes.resize(states.size());
  if(!foundPastNodes)
  {
    for (CFuint i=0; i < volumes.size();++i){
      pastVolumes[i] = volumes[i];
    }
  }
  else{

    vector<RealVector> backNodes(nodes.size());
    for(CFuint iNode=0; iNode<nodes.size(); ++iNode)
    {
      backNodes[iNode].resize(nodes[iNode]->size());
      backNodes[iNode] = *(nodes[iNode]);
    }

    //Recompute past Volumes
    Common::SafePtr<TopologicalRegionSet> cells = MeshDataStack::getActive()->
      getTrs("InnerCells");

    Common::SafePtr<GeometricEntityPool<TrsGeoWithNodesBuilder> >
      geoBuilder = getMethodData().getGeoWithNodesBuilder();

    TrsGeoWithNodesBuilder::GeoData& geoData = geoBuilder->getDataGE();
    geoData.trs = cells;

    const CFuint nbElems = cells->getLocalNbGeoEnts();
    for (CFuint iElem = 0; iElem < nbElems; ++iElem) {
      // build the GeometricEntity
      geoData.idx = iElem;
      GeometricEntity *const cell = geoBuilder->buildGE();

      //Compute the pastVolume
      pastVolumes[iElem] = cell->computeVolume();

      cf_assert(pastVolumes[iElem] > 0.);

      //release the GeometricEntity
      geoBuilder->releaseGE();
    }

    for(CFuint iNode=0; iNode<nodes.size(); ++iNode)
    {
      *(nodes[iNode]) = backNodes[iNode];
    }

  }

}

//////////////////////////////////////////////////////////////////////////////

vector<SafePtr<BaseDataSocketSink> >
StdALESetup::needsSockets()
{
  vector<SafePtr<BaseDataSocketSink> > result;

  result.push_back(&socket_nodes);
  result.push_back(&socket_states);
  result.push_back(&socket_volumes);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD
