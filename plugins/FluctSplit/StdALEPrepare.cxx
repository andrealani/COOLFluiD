#include "FluctSplit/FluctSplitSpaceTime.hh"
#include "StdALEPrepare.hh"
#include "Framework/MeshData.hh"
#include "Framework/MethodCommandProvider.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<StdALEPrepare, FluctuationSplitData, FluctSplitSpaceTimeModule> StdALEPrepareProvider("StdALEPrepare");

//////////////////////////////////////////////////////////////////////////////

StdALEPrepare::StdALEPrepare(const std::string& name) :
FluctuationSplitCom(name),
socket_nodes("nodes"),
socket_pastNodes("pastNodes"),
socket_pastCellVolume("pastCellVolume"),
socket_cellSpeed("cellSpeed")
{
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
StdALEPrepare::needsSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSink> > result;

  result.push_back(&socket_nodes);
  result.push_back(&socket_pastNodes);
  result.push_back(&socket_pastCellVolume);
  result.push_back(&socket_cellSpeed);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

void StdALEPrepare::execute()
{
  Common::SafePtr<TopologicalRegionSet> cells = MeshDataStack::getActive()->getTrs("InnerCells");

  DataHandle < Framework::Node*, Framework::GLOBAL > nodes  = socket_nodes.getDataHandle();
  DataHandle<Node*> pastNodes  = socket_pastNodes.getDataHandle();
  DataHandle<CFreal> pastCellVolume = socket_pastCellVolume.getDataHandle();
  DataHandle<RealVector> cellSpeed = socket_cellSpeed.getDataHandle();

  // Set the pastNodes = Nodes
  cf_assert (nodes.size() == pastNodes.size());
  for (CFuint i=0; i < nodes.size();++i){
    *pastNodes[i] = *nodes[i];
  }


  Common::SafePtr<GeometricEntityPool<StdTrsGeoBuilder> >
  geoBuilder = getMethodData().getStdTrsGeoBuilder();

  StdTrsGeoBuilder::GeoData& geoData = geoBuilder->getDataGE();
  geoData.trs = cells;
  const CFuint nbGeos = cells->getLocalNbGeoEnts();

  // Store the past cell volume into pastCellInfo
  for (CFuint iCell = 0; iCell < nbGeos; ++iCell) {
    CFLogDebugMed( "Computing iCell = " << iCell << "\n");

    // build the GeometricEntity
    geoData.idx = iCell;
    GeometricEntity& cell = *geoBuilder->buildGE();

    pastCellVolume[iCell] = cell.computeVolume();
    cellSpeed[iCell] = cell.computeCentroid();
    geoBuilder->releaseGE();
  }

}

//////////////////////////////////////////////////////////////////////////////

 } // namespace FluctSplit



} // namespace COOLFluiD
