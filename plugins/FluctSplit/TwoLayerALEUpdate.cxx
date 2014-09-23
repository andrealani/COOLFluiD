#include "FluctSplit/FluctSplitSpaceTime.hh"
#include "TwoLayerALEUpdate.hh"
#include "Framework/MeshData.hh"
#include "Framework/ComputeNormals.hh"
#include "MathTools/RealMatrix.hh"
#include "InwardNormalsData.hh"
#include "ComputeInwardNormals.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/SubSystemStatus.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<TwoLayerALEUpdate, FluctuationSplitData, FluctSplitSpaceTimeModule> TwoLayerALEUpdateProvider("TwoLayerALEUpdate");

//////////////////////////////////////////////////////////////////////////////

TwoLayerALEUpdate::TwoLayerALEUpdate(const std::string& name) :
  FluctuationSplitCom(name),
  socket_normals("normals"),
  socket_pastNormals("pastNormals"),
  socket_interNormals("interNormals"),
  socket_cellSpeed("cellSpeed")
{
}

//////////////////////////////////////////////////////////////////////////////

TwoLayerALEUpdate::~TwoLayerALEUpdate()
{
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
TwoLayerALEUpdate::needsSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSink> > result;

  result.push_back(&socket_normals);
  result.push_back(&socket_pastNormals);
  result.push_back(&socket_interNormals);
  result.push_back(&socket_cellSpeed);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

void TwoLayerALEUpdate::execute()
{
  computeCellSpeed();
  savePastNormals();
  updateNormalsData();
  computeInterNormalsData();
 }

//////////////////////////////////////////////////////////////////////////////

void TwoLayerALEUpdate::computeInterNormalsData()
{
  CFLogDebugMin( "TwoLayerALEUpdate::computeInterNormalsData() => begin" << "\n");

  SafePtr<vector<ElementTypeData> > elemTypes =
    MeshDataStack::getActive()->getElementTypeData();

  const CFuint nbElemTypes = elemTypes->size();

  for (CFuint iType = 0; iType < nbElemTypes; ++iType) {

    const CFuint geoOrder = (*elemTypes)[iType].getGeoOrder();
    const std::string elemName = (*elemTypes)[iType].getShape() + CFPolyOrder::Convert::to_str( geoOrder );

    SelfRegistPtr<ComputeNormals> computeFaceNormals = Environment::Factory<ComputeNormals>::getInstance().
      getProvider("Inward" + elemName)->create();

    SelfRegistPtr<ComputeInwardNormals> dPtr = computeFaceNormals.d_castTo<ComputeInwardNormals>();
    // set the data sockets
    dPtr->setNormalsSockets(&socket_normals);
    dPtr->setPastNormalsSockets(&socket_pastNormals);
    dPtr->setInterNormalsSockets(&socket_interNormals);

    const CFuint firstElem = (*elemTypes)[iType].getStartIdx();
    const CFuint lastElem  = (*elemTypes)[iType].getEndIdx();

    // Compute the weighted average between the normals at time n and n+1
    (*dPtr).average(firstElem, lastElem, iType);
  }

  CFLogDebugMin( "TwoLayerALEUpdate::computeInterNormalsData() => end" << "\n");
}

//////////////////////////////////////////////////////////////////////////////

void TwoLayerALEUpdate::computeCellSpeed()
{

  const CFreal dt = SubSystemStatusStack::getActive()->getDT();

  // Compute the speed of the centroid of the cell
  Common::SafePtr<TopologicalRegionSet> cells = MeshDataStack::getActive()->getTrs("InnerCells");

  DataHandle<RealVector> cellSpeed = socket_cellSpeed.getDataHandle();

  // Store the cell speed into cellSpeed datahandle
  Common::SafePtr<GeometricEntityPool<StdTrsGeoBuilder> >
  geoBuilder = getMethodData().getStdTrsGeoBuilder();

  StdTrsGeoBuilder::GeoData& geoData = geoBuilder->getDataGE();
  geoData.trs = cells;
  const CFuint nbGeos = cells->getLocalNbGeoEnts();

  // Store the cell speed into cellSpeed datahandle
  for (CFuint iCell = 0; iCell < nbGeos; ++iCell) {
    CFLogDebugMed( "Computing cell speed for iCell = " << iCell << "\n");
    // build the GeometricEntity
    geoData.idx = iCell;
    GeometricEntity& cell = *geoBuilder->buildGE();

    cellSpeed[iCell] -= cell.computeCentroid();
    cellSpeed[iCell] /= dt;
    cellSpeed[iCell] *= -1.;
    geoBuilder->releaseGE();
  }
}

//////////////////////////////////////////////////////////////////////////////

void TwoLayerALEUpdate::savePastNormals()
{
  CFLogDebugMin( "TwoLayerALEUpdate::savePastNormals() => begin" << "\n");

  DataHandle<InwardNormalsData*> normals = socket_normals.getDataHandle();

  DataHandle<InwardNormalsData*> pastNormals = socket_pastNormals.getDataHandle();

  cf_assert(normals.size() == pastNormals.size());
  for (CFuint i=0; i<normals.size(); ++i){
    *pastNormals[i] = *normals[i];
  }

  CFLogDebugMin( "TwoLayerALEUpdate::savePastNormals() => end" << "\n");
}

//////////////////////////////////////////////////////////////////////////////

void TwoLayerALEUpdate::updateNormalsData()
{
  CFLogDebugMin( "TwoLayerALEUpdate::updateNormalsData() => begin" << "\n");

  SafePtr<vector<ElementTypeData> > elemTypes =
    MeshDataStack::getActive()->getElementTypeData();

  const CFuint nbElemTypes = elemTypes->size();

  for (CFuint iType = 0; iType < nbElemTypes; ++iType) {

    const CFuint geoOrder = (*elemTypes)[iType].getGeoOrder();
    const std::string elemName = (*elemTypes)[iType].getShape() +  CFPolyOrder::Convert::to_str( geoOrder );

    SelfRegistPtr<ComputeNormals> computeFaceNormals = Environment::Factory<ComputeNormals>::getInstance().
      getProvider("Inward" + elemName)->create();

    SelfRegistPtr<ComputeInwardNormals> dPtr = computeFaceNormals.d_castTo<ComputeInwardNormals>();
    // set the data sockets
    dPtr->setNormalsSockets(&socket_normals);

    const CFuint firstElem = (*elemTypes)[iType].getStartIdx();
    const CFuint lastElem  = (*elemTypes)[iType].getEndIdx();

    (*dPtr).update(firstElem, lastElem, iType);
  }

  CFLogDebugMin( "TwoLayerALEUpdate::updateNormalsData() => end" << "\n");
}

//////////////////////////////////////////////////////////////////////////////

 } // namespace FluctSplit



} // namespace COOLFluiD
