#include "Framework/MeshData.hh"
#include "Framework/DomainModel.hh"
// #include "Framework/SubSystemStatus.hh"
// 
#include "Framework/MethodCommandProvider.hh"
// #include "Framework/MapGeoToTrsAndIdx.hh"
// #include "Framework/NamespaceSwitcher.hh"

#include "FluctSplit/FluctSplit.hh"
#include "FluctSplit/MoveBoundary.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<MoveBoundary, FluctuationSplitData, FluctSplitModule> MoveBoundaryProvider("MoveBoundary");

//////////////////////////////////////////////////////////////////////////////

MoveBoundary::MoveBoundary(const std::string& name) :
  FluctuationSplitCom(name),
  socket_states("states"),
  socket_nodes("nodes"),
  socket_isBState("isBState")
{
}

//////////////////////////////////////////////////////////////////////////////

MoveBoundary::~MoveBoundary()
{
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSource> >
MoveBoundary::providesSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSource> > result;
  return result;
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
MoveBoundary::needsSockets()
{

  std::vector<Common::SafePtr<BaseDataSocketSink> > result;

  result.push_back(&socket_isBState);
  result.push_back(&socket_states);
  result.push_back(&socket_nodes);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

void MoveBoundary::configure ( Config::ConfigArgs& args )
{
  CFAUTOTRACE;

  FluctuationSplitCom::configure(args);

//   //Loop over the TRS's and add the "TRSName" + "-boundaryNormals" datasocketsink to the m_dynamicSockets
//   const std::string name = getMethodData().getNamespace();
//   Common::SafePtr<Namespace> nsp = NamespaceSwitcher::getInstance().getNamespace(name);
//   Common::SafePtr<MeshData> meshData = MeshDataStack::getInstance().getEntryByNamespace(nsp);
//
//   vector<std::string> trsList = meshData->getTRSNameList();
//   const CFint nbTRSs = trsList.size();
//
//   for (CFint iTRS = 0; iTRS < nbTRSs; ++iTRS) {
//     const std::string trsName = trsList[iTRS];
//
//     if (trsName != "PartitionFaces" &&
//         trsName != "InnerCells" &&
//         trsName != "InnerFaces")
//     {
//       const std::string socketName = trsName + "-boundaryNormals";
//       const bool isEssential = false;
//       m_dynamicSockets.createSocketSink<const CFreal*>(socketName,isEssential);
//     }
//   }
}

//////////////////////////////////////////////////////////////////////////////

void MoveBoundary::executeOnTrs()
{
  CFAUTOTRACE;

  RealVector mapCoordFaceNode(DIM_1D);

  Common::SafePtr<TopologicalRegionSet> trs = getCurrentTRS();
  Common::SafePtr< vector<CFuint> > const statesIdx = trs->getStatesInTrs();

  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();
  DataHandle < Framework::Node*, Framework::GLOBAL > nodes  = socket_nodes.getDataHandle();

  // get the domain model
  SafePtr< DomainModel > domModel = MeshDataStack::getActive()->getDomainModel();

  // builder of faces
  SafePtr< GeometricEntityPool<StdTrsGeoBuilder> > geoBuilder =
    getMethodData().getStdTrsGeoBuilder();

  StdTrsGeoBuilder::GeoData& geoData = geoBuilder->getDataGE();
  geoData.trs = trs;

  // get number of TRs in this TRS
  const CFuint nbTRs = trs->getNbTRs();

  // loop over the TRs
  for (DomainModel::TRidx iTR = 0; iTR < nbTRs; ++iTR)
  {
    // get the TR
    SafePtr< TopologicalRegion > bcTR = trs->getTopologicalRegion(iTR);

    // get TR global index
    const std::string trKey = trs->getName() + StringOps::to_str(iTR);
    DomainModel::TRidx trGlobalIdx = domModel->getTRGlobalIdx(trKey);

    // loop over faces in this TR
    const CFuint nbrFaces = bcTR->getLocalNbGeoEnts();
    for (CFuint iFace = 0; iFace < nbrFaces; ++iFace)
    {
        // get face index in the TRS
        const CFuint faceIdx = bcTR->getGeoIDInTrs(iFace);

        // build the face GeometricEntity
        geoData.idx = faceIdx;
        GeometricEntity* face = geoBuilder->buildGE();

        // get face nodes
        vector< Node* >& nodes = *(face->getNodes());
        // move all nodes in face
        for (vector<Node*>::iterator itr = nodes.begin(), stop = nodes.end(); itr != stop; ++itr)
        {
          Node& this_node = (**itr);
          // guess for the mapped coordinates
          mapCoordFaceNode = 0.5;
          // project into the boundary and get the true mapped coordinates
          domModel->computeParamCoord(trGlobalIdx, this_node, mapCoordFaceNode);
          // compute the moved coordinates from the mapped corrdinates
          domModel->computeCoord(trGlobalIdx, mapCoordFaceNode, this_node);
        }

        // get face states
        vector< State* >& states = *(face->getStates());
        // move all nodes in face
        for (vector<State*>::iterator itr = states.begin(), stop = states.end(); itr != stop; ++itr)
        {
          Node& this_node = (*itr)->getCoordinates();
          // guess for the mapped coordinates
          mapCoordFaceNode = 0.5;
          // project into the boundary and get the true mapped coordinates
          domModel->computeParamCoord(trGlobalIdx, this_node, mapCoordFaceNode);
          // compute the moved coordinates from the mapped corrdinates
          domModel->computeCoord(trGlobalIdx, mapCoordFaceNode, this_node);
        }

        // release the face
        geoBuilder->releaseGE();
     } // loop faces
  } // loop TR's
}

//////////////////////////////////////////////////////////////////////////////

 } // namespace FluctSplit



} // namespace COOLFluiD
