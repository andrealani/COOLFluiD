#include "FiniteVolume/FiniteVolume.hh"
#include "FiniteVolume/StdALEUpdate.hh"
#include "FiniteVolume/FVMCC_PolyRec.hh"

#include "Framework/ComputeDummyStates.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/MeshData.hh"
#include "Framework/VolumeCalculator.hh"
#include "Framework/Node.hh"
#include "Framework/SetElementStateCoord.hh"
#include "Framework/ComputeFaceNormalsFVMCC.hh"
#include "Framework/SubSystemStatus.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<StdALEUpdate, CellCenterFVMData, FiniteVolumeModule> StdALEUpdateProvider("StdALEUpdate");

//////////////////////////////////////////////////////////////////////////////

StdALEUpdate::StdALEUpdate(const std::string& name) :
  CellCenterFVMCom(name),
  socket_nodes("nodes"),
  socket_pastNodes("pastNodes"),
  socket_futureNodes("futureNodes"),
  socket_states("states"),
  socket_isOutward("isOutward"),
  socket_normals("normals"),
  socket_faceAreas("faceAreas"),
  socket_volumes("volumes"),
  socket_pastVolumes("pastVolumes"),
  socket_gstates("gstates")
{
}

//////////////////////////////////////////////////////////////////////////////

void StdALEUpdate::execute()
{

  backupFutureNodes();

  if(SubSystemStatusStack::getActive()->isSubIterationFirstStep()){
    backupCellVolume();
  }

  ///@todo doing this here is dangerous because for the source terms,
  ///you would need the volume computed at intermediate location
  ///while here you are with nodes=futureNodes
  updateCellVolume();

  computeIntermediateNodes();

  resetIsOutward();
  updateNormalsData();
  updateFaceAreas();
  updateReconstructor();

  //!->To modify the Ghost nodes, we need the normals -> after updateNormalsData
  modifyOffMeshNodes();

 }

//////////////////////////////////////////////////////////////////////////////

void StdALEUpdate::updateReconstructor()
{
  Common::SafePtr<FVMCC_PolyRec> polyRec = getMethodData().getPolyReconstructor().d_castTo<FVMCC_PolyRec>();

  polyRec->updateWeights();

}

//////////////////////////////////////////////////////////////////////////////

void StdALEUpdate::backupFutureNodes()
{
  DataHandle < Framework::Node*, Framework::GLOBAL > nodes  = socket_nodes.getDataHandle();

  DataHandle<Node*> futureNodes = socket_futureNodes.getDataHandle();

  // Set the futureNodes equal to the displaced Nodes
  cf_assert (nodes.size() == futureNodes.size());
  for (CFuint i=0; i < nodes.size();++i){
    *futureNodes[i] = *nodes[i];
  }
}

//////////////////////////////////////////////////////////////////////////////

void StdALEUpdate::computeIntermediateNodes()
{
  DataHandle < Framework::Node*, Framework::GLOBAL > nodes = socket_nodes.getDataHandle();
  DataHandle<Node*> pastNodes = socket_pastNodes.getDataHandle();
  DataHandle<Node*> futureNodes =socket_futureNodes.getDataHandle();

  // Set the intermediate nodes r_int = 0.5*(r_n + r_n+1)
  /// @todo This is valid only for Cranck-Nicholson!!!!!!

  cf_assert (nodes.size() == futureNodes.size());
  for (CFuint i=0; i < nodes.size();++i){
    *nodes[i] = 0.5 * (*pastNodes[i]);
    *nodes[i] += 0.5 * (*futureNodes[i]);
  }

}

//////////////////////////////////////////////////////////////////////////////

void StdALEUpdate::modifyOffMeshNodes()
{
  /// The cell-center nodes and ghost nodes have NOT been moved during the mesh adaptation!!!
  //
  /// Modify the element state coordinates

  SafePtr<vector<ElementTypeData> > elementType =
    MeshDataStack::getActive()->getElementTypeData();

  CFuint nbElemTypes = elementType->size();

  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();

  SafePtr<TopologicalRegionSet> cells = MeshDataStack::getActive()->
    getTrs("InnerCells");

  SafePtr<GeometricEntityPool<TrsGeoWithNodesBuilder> > geoBuilder =
    getMethodData().getGeoWithNodesBuilder();

  TrsGeoWithNodesBuilder::GeoData& geoData = geoBuilder->getDataGE();
  geoData.trs = cells;

  vector<State*> eState(1);

  CFuint elemID = 0;
  for (CFuint iType = 0; iType < nbElemTypes; ++iType) {

    ///@todo change this for other order!!
    const std::string elemName = CFGeoShape::Convert::to_str((*elementType)[iType].getGeoShape()) + "LagrangeP1LagrangeP0";

    SelfRegistPtr<SetElementStateCoord> setStateCoord = 
      FACTORY_GET_PROVIDER(getFactoryRegistry(), SetElementStateCoord, elemName)->create();
    
    const CFuint nbElemPerType = (*elementType)[iType].getNbElems();
    for (CFuint iElem = 0; iElem < nbElemPerType; ++iElem, ++elemID) {
      // build the cell
      geoData.idx = elemID;
      GeometricEntity *const currCell = geoBuilder->buildGE();

      const vector<Node*>& eNodes = *currCell->getNodes();
      eState[0] = states[elemID];

      setStateCoord->update(eNodes,eState);

      // release the cell
      geoBuilder->releaseGE();
    }
  }

  ///@todo Do the same for the ghost states
  // compute the dummy (ghost) states
  ComputeDummyStates computeDummyStates;
  computeDummyStates.setDataSockets(socket_normals, socket_gstates,socket_states, socket_nodes);
  computeDummyStates.setup();
  computeDummyStates.updateAllDummyStates();

}

//////////////////////////////////////////////////////////////////////////////

void StdALEUpdate::resetIsOutward()
{
  DataHandle<CFint> isOutward = socket_isOutward.getDataHandle();

  for(CFuint i=0; i<isOutward.size(); i++)
  {
    isOutward[i] = -1;
  }
}

//////////////////////////////////////////////////////////////////////////////

void StdALEUpdate::updateNormalsData()
{
  CFAUTOTRACE;

  Common::SafePtr<vector<ElementTypeData> > elemTypes =
    MeshDataStack::getActive()->getElementTypeData();

  SafePtr<DataSocketSink<CFreal> > sinkNormalsPtr = &socket_normals;
  SafePtr<DataSocketSink<CFint> > sinkIsOutwardPtr = &socket_isOutward;

  for (CFuint iType = 0; iType < elemTypes->size(); ++iType) {

    const CFuint geoOrder = (*elemTypes)[iType].getGeoOrder();
    const std::string elemName = (*elemTypes)[iType].getShape() + CFPolyOrder::Convert::to_str(geoOrder);
    const string cname = "Face" + elemName;
    Common::SelfRegistPtr<ComputeNormals> computeFaceNormals =
      FACTORY_GET_PROVIDER(getFactoryRegistry(), ComputeNormals, cname)->create();
    
    const CFuint firstElem = (*elemTypes)[iType].getStartIdx();
    const CFuint lastElem  = (*elemTypes)[iType].getEndIdx();
    
    SelfRegistPtr<ComputeFaceNormalsFVMCC> faceNormalsComputer =
      computeFaceNormals.d_castTo<ComputeFaceNormalsFVMCC>();
    
    faceNormalsComputer->setSockets(sinkNormalsPtr, sinkIsOutwardPtr);
    (*faceNormalsComputer)(firstElem, lastElem);
  }
  
}

//////////////////////////////////////////////////////////////////////////////

void StdALEUpdate::updateFaceAreas()
{
  const CFuint totalNbFaces = MeshDataStack::getActive()->Statistics().getNbFaces();

  DataHandle<CFreal> normals = socket_normals.getDataHandle();

  const CFuint nbDim = PhysicalModelStack::getActive()->getDim();
  RealVector faceNormal(nbDim);
  for (CFuint iFace = 0; iFace < totalNbFaces; ++iFace) {
    const CFuint startID = iFace*nbDim;
    //Update the normals
    for (CFuint iDim = 0; iDim < nbDim; ++iDim) {
      faceNormal[iDim] = normals[startID + iDim];
    }

    //Compute-update the face Area
    /// @todo This is valid only for Cranck-Nicholson!!!!!! 0.5* (norm2+ pastFaceNormlas)
    (socket_faceAreas.getDataHandle())[iFace] = faceNormal.norm2();
  }
}

//////////////////////////////////////////////////////////////////////////////

void StdALEUpdate::backupCellVolume()
{
  DataHandle<CFreal> volumes = socket_volumes.getDataHandle();
  DataHandle<CFreal> pastVolumes = socket_pastVolumes.getDataHandle();

  for (CFuint i=0; i < volumes.size();++i){
    pastVolumes[i] = volumes[i];
  }
}

//////////////////////////////////////////////////////////////////////////////

void StdALEUpdate::updateCellVolume()
{
  DataHandle<CFreal> volumes = socket_volumes.getDataHandle();

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
    volumes[iElem] = cell->computeVolume();

    //release the GeometricEntity
    geoBuilder->releaseGE();
  }
}

//////////////////////////////////////////////////////////////////////////////

vector<SafePtr<BaseDataSocketSink> >
StdALEUpdate::needsSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSink> > result;

  result.push_back(&socket_nodes);
  result.push_back(&socket_pastNodes);
  result.push_back(&socket_futureNodes);
  result.push_back(&socket_states);
  result.push_back(&socket_isOutward);
  result.push_back(&socket_normals);
  result.push_back(&socket_faceAreas);
  result.push_back(&socket_volumes);
  result.push_back(&socket_pastVolumes);
  result.push_back(&socket_gstates);
  return result;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD
