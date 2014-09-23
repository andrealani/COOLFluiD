#include "MeshLaplacianSmoothing/MeshLaplacianSmoothing.hh"
#include "StdPrepare.hh"
#include "Framework/MeshData.hh"
#include "Framework/PhysicalModel.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/SubSystemStatus.hh"
#include "Common/SelfRegistPtr.hh"
#include "MathTools/MathConsts.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace MeshLaplacianSmoothing {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<StdPrepare, LaplacianSmoothingData, MeshLaplacianSmoothingModule> StdPrepareProvider("StdPrepare");

//////////////////////////////////////////////////////////////////////////////

StdPrepare::StdPrepare(const std::string& name) :
LaplacianSmoothingCom(name),
  socket_nodes("nodes"),
  socket_pastNodes("pastNodes"),
  socket_nodalDisplacements("nodalDisplacements"),
  socket_nodeToCellConnectivity("nodeToCellConnectivity"),
  socket_qualityNode("qualityNode")
{
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
StdPrepare::needsSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSink> > result;

  result.push_back(&socket_nodes);
  result.push_back(&socket_pastNodes);
  result.push_back(&socket_nodalDisplacements);
  result.push_back(&socket_nodeToCellConnectivity);
  result.push_back(&socket_qualityNode);

  return result;
}


//////////////////////////////////////////////////////////////////////////////

void StdPrepare::execute()
{
  CFAUTOTRACE;

  computeWorstQualityCells();
}

//////////////////////////////////////////////////////////////////////////////

void StdPrepare::computeWorstQualityCells()
{
  CFAUTOTRACE;

  Common::SelfRegistPtr<MeshTools::QualityCalculator> qualityComputer =
    Environment::Factory<MeshTools::QualityCalculator>::getInstance().getProvider("Concrete")->create("Concrete");

  cf_assert(qualityComputer.isNotNull());

  /// Get the datahandle of pointers to GeomEntity
  DataHandle< vector<CFuint> > nodeToCellConn = socket_nodeToCellConnectivity.getDataHandle();
  DataHandle< CFreal> qualityNode = socket_qualityNode.getDataHandle();

  //Loop over Nodes
  // set the list of faces
  SafePtr<TopologicalRegionSet> trs =
    MeshDataStack::getActive()->getTrs("InnerCells");


  SafePtr<vector<CFuint> > nodesInTRS = trs->getNodesInTrs();
  vector<CFuint>::iterator it;
  for (it = nodesInTRS->begin(); it != nodesInTRS->end(); ++it) {
    const CFuint localID = *it;
    vector<CFuint>& cellsID = nodeToCellConn[localID];

    Common::SafePtr<GeometricEntityPool<StdTrsGeoBuilder> >
    geoBuilder = getMethodData().getStdTrsGeoBuilder();

    StdTrsGeoBuilder::GeoData& geoData = geoBuilder->getDataGE();
    geoData.trs = trs;
//unused//    const CFuint nbCells = trs->getLocalNbGeoEnts();

    CFreal worseQuality = 0.;
    for (CFuint iCell=0; iCell < cellsID.size(); ++iCell) {
      // build the GeometricEntity
      geoData.idx = cellsID[iCell];
      GeometricEntity* currCell = geoBuilder->buildGE();

      CFreal quality = qualityComputer->computeQuality(currCell);
      if (quality > worseQuality) worseQuality = quality;
      if (quality < 0.) worseQuality = MathTools::MathConsts::CFrealMax();

      geoBuilder->releaseGE();
    }
    if(qualityNode[localID] !=0.) qualityNode[localID] = worseQuality;
//    cout << qualityNode.size() << " " << localID << " "<< qualityNode[localID]<< endl;
  }

  cout << "I am out"<< endl;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace MeshLaplacianSmoothing

  } // namespace Numerics

} // namespace COOLFluiD
