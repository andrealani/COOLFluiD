#include "FluctSplit/FluctSplitSpaceTime.hh"
#include "StdALESetup.hh"
#include "Framework/MeshData.hh"
#include "Framework/MethodCommandProvider.hh"
#include "InwardNormalsData.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<StdALESetup, FluctuationSplitData, FluctSplitSpaceTimeModule> StdALESetupProvider("StdALESetup");

//////////////////////////////////////////////////////////////////////////////

StdALESetup::StdALESetup(const std::string& name) :
  FluctuationSplitCom(name),
  socket_pastNormals("pastNormals"),
  socket_pastNormalsData("pastNormalsData"),
  socket_pastCellVolume("pastCellVolume"),
  socket_cellSpeed("cellSpeed"),
  socket_pastNodes("pastNodes"),
  socket_normals("normals"),
  socket_nodes("nodes")
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

  result.push_back(&socket_pastNormals);
  result.push_back(&socket_pastNormalsData);
  result.push_back(&socket_pastCellVolume);
  result.push_back(&socket_cellSpeed);
  result.push_back(&socket_pastNodes);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
StdALESetup::needsSockets()
{

  std::vector<Common::SafePtr<BaseDataSocketSink> > result;

  result.push_back(&socket_normals);
  result.push_back(&socket_nodes);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

void StdALESetup::execute()
{
  CFAUTOTRACE;

  Common::SafePtr<TopologicalRegionSet> cells = MeshDataStack::getActive()->getTrs("InnerCells");
  const CFuint nbCells = cells->getLocalNbGeoEnts();

  DataHandle<InwardNormalsData*> normals = socket_normals.getDataHandle();

  /// @todo try to find a way to avoid this because it takes lots of memory
  DataHandle<InwardNormalsData*> pastNormals = socket_pastNormals.getDataHandle();

  pastNormals.resize(normals.size());

  for (CFuint i=0; i<nbCells;++i){
    pastNormals[i] = new InwardNormalsData(*normals[i]);
  }

  DataHandle<CFreal> pastNormalsData = socket_pastNormalsData.getDataHandle();
  pastNormalsData.resize(normals.size());

  DataHandle<CFreal> pastCellVolume = socket_pastCellVolume.getDataHandle();
  pastCellVolume.resize(nbCells);

  DataHandle<RealVector> cellSpeed = socket_cellSpeed.getDataHandle();
  cellSpeed.resize(nbCells);

  for (CFuint i = 0; i < nbCells; ++i) {
    cellSpeed[i].resize(PhysicalModelStack::getActive()->getDim());
  }

  DataHandle < Framework::Node*, Framework::GLOBAL > nodes = socket_nodes.getDataHandle();

  const CFuint nbNodes = nodes.size();

  DataHandle<Node*> pastNodes = socket_pastNodes.getDataHandle();
  pastNodes.resize(nbNodes);

  for (CFuint i = 0; i < nbNodes; ++i) {
    pastNodes[i] = new Node();
  }
}

//////////////////////////////////////////////////////////////////////////////

 } // namespace FluctSplit



} // namespace COOLFluiD
