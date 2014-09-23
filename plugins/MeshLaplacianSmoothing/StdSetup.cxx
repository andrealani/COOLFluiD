#include "Framework/MethodCommandProvider.hh"
#include "Framework/SubSystemStatus.hh"
#include "Framework/PhysicalModel.hh"
#include "Framework/MeshData.hh"
#include "Framework/Node.hh"

#include "MeshLaplacianSmoothing/MeshLaplacianSmoothing.hh"
#include "MeshLaplacianSmoothing/StdSetup.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace MeshLaplacianSmoothing {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<StdSetup, LaplacianSmoothingData, MeshLaplacianSmoothingModule> StdSetupProvider("StdSetup");

//////////////////////////////////////////////////////////////////////////////

StdSetup::StdSetup(std::string name) : LaplacianSmoothingCom(name),
  socket_averageVector("averageVector"),
  socket_sumWeight("sumWeight"),
  socket_qualityNode("qualityNode"),
  socket_nodalDisplacements("nodalDisplacements"),
  socket_nodeToCellConnectivity("nodeToCellConnectivity"),
  socket_nodes("nodes")
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
  result.push_back(&socket_qualityNode);
  result.push_back(&socket_nodeToCellConnectivity);
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
  CFAUTOTRACE;

  SubSystemStatusStack::getActive()->setMovingMesh(true);

  DataHandle < Framework::Node*, Framework::GLOBAL > nodes = socket_nodes.getDataHandle();
  const CFuint nbNodes = nodes.size();

  DataHandle<RealVector> averageVector = socket_averageVector.getDataHandle();
  averageVector.resize(nbNodes);

  for (CFuint i = 0; i < nbNodes; ++i) {
    averageVector[i].resize(PhysicalModelStack::getActive()->getDim());
    averageVector[i] = 0.;
  }

  DataHandle<CFreal> sumWeight = socket_sumWeight.getDataHandle();
  sumWeight.resize(nbNodes);

  for (CFuint i = 0; i < nbNodes; ++i) {
    sumWeight[i] = 0.;
  }

  DataHandle<CFreal> qualityNode = socket_qualityNode.getDataHandle();
  qualityNode.resize(nbNodes);
  qualityNode = 1.0;

  DataHandle<RealVector> nodalDisplacements = socket_nodalDisplacements.getDataHandle();
  nodalDisplacements.resize(nbNodes);
   for (CFuint i = 0; i < nbNodes; ++i) {
     nodalDisplacements[i].resize(PhysicalModelStack::getActive()->getDim());
     nodalDisplacements[i] = 0.;
   }

  DataHandle<vector<CFuint> > nodeToCellConnectivity = socket_nodeToCellConnectivity.getDataHandle();
  nodeToCellConnectivity.resize(nbNodes);

  createNodesListInTRS();

  flagBoundaryNodes();

  createNodeToCellConnectivity();
}

//////////////////////////////////////////////////////////////////////////////

void StdSetup::createNodesListInTRS()
{
  std::vector<Common::SafePtr<TopologicalRegionSet> > trs =
    MeshDataStack::getActive()->getTrsList();
  std::vector< Common::SafePtr<TopologicalRegionSet> >::iterator itrs;
  for (itrs = trs.begin(); itrs != trs.end(); ++itrs) {
    (*itrs)->createNodesList();
  }
}

//////////////////////////////////////////////////////////////////////////////

void StdSetup::flagBoundaryNodes()
{
  DataHandle< CFreal> qualityNode = socket_qualityNode.getDataHandle();

  // set the list of faces
  std::vector< Common::SafePtr<TopologicalRegionSet> > trs = MeshDataStack::getActive()->getTrsList();

  std::vector< Common::SafePtr<TopologicalRegionSet> >::iterator itrs;
  for (itrs = trs.begin(); itrs != trs.end(); ++itrs) {
    if ((*itrs)->getName() != "InnerCells") {
      SafePtr<vector<CFuint> > nodesInTRS = (*itrs)->getNodesInTrs();
      vector<CFuint>::iterator it;
      for (it = nodesInTRS->begin(); it != nodesInTRS->end(); ++it) {
        //Set the quality to the ideal quality -> they are never moved
        qualityNode[*it] = 0.;
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void StdSetup::createNodeToCellConnectivity()
{
  DataHandle<vector<CFuint> > nodeToCellConnectivity = socket_nodeToCellConnectivity.getDataHandle();

  vector<Common::SafePtr<TopologicalRegionSet> > trs = MeshDataStack::getActive()->getTrsList();
  vector< Common::SafePtr<TopologicalRegionSet> >::iterator itrs;
  for (itrs = trs.begin(); itrs != trs.end(); ++itrs) {
    if ((*itrs)->getName() == "InnerCells") {
      CFLogDebugMin("Building Node to Cell Connectivity...\n");

      Common::SafePtr<std::vector<CFuint> > nodesInTrs = (*itrs)->getNodesInTrs();

      Common::SafePtr<GeometricEntityPool<StdTrsGeoBuilder> >
      geoBuilder = getMethodData().getStdTrsGeoBuilder();

      StdTrsGeoBuilder::GeoData& geoData = geoBuilder->getDataGE();
      geoData.trs = (*itrs);
      const CFuint nbCells = (*itrs)->getLocalNbGeoEnts();

      for (CFuint iCell = 0; iCell < nbCells; ++iCell) {
        CFLogDebugMed( "Computing iCell = " << iCell << "\n");

        // build the GeometricEntity
        geoData.idx = iCell;
        GeometricEntity& currCell = *geoBuilder->buildGE();
        const CFuint nbCellNodes = currCell.getNodes()->size();

        for(CFuint iNode=0; iNode < nbCellNodes; iNode++)
        {
          CFuint localID = (currCell.getNode(iNode))->getLocalID();
          nodeToCellConnectivity[localID].push_back(iCell);
        }

        geoBuilder->releaseGE();
      }
      CFLogDebugMin("Finished building Node to Cell Connectivity...\n");
    } // end of if is InnerCells
  } // end of loop over TRSs
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace MeshLaplacianSmoothing

  } // namespace Numerics

} // namespace COOLFluiD
