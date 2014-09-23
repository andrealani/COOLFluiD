#include "MathTools/RealMatrix.hh"

#include "Framework/MeshData.hh"
#include "Framework/SubSystemStatus.hh"
#include "Framework/ComputeNormals.hh"
#include "Framework/MethodCommandProvider.hh"

#include "FluctSplit/InwardNormalsData.hh"
#include "FluctSplit/ComputeInwardNormals.hh"
#include "FluctSplit/FluctSplitSpaceTime.hh"
#include "FluctSplit/StdALEUpdate.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<StdALEUpdate, FluctuationSplitData, FluctSplitSpaceTimeModule> StdALEUpdateProvider("StdALEUpdate");

//////////////////////////////////////////////////////////////////////////////

StdALEUpdate::StdALEUpdate(const std::string& name) :
  FluctuationSplitCom(name),
  socket_normals("normals"),
  socket_pastNormals("pastNormals"),
  socket_cellSpeed("cellSpeed")
{
}

//////////////////////////////////////////////////////////////////////////////

StdALEUpdate::~StdALEUpdate()
{
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
StdALEUpdate::needsSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSink> > result;

  result.push_back(&socket_normals);
  result.push_back(&socket_pastNormals);
  result.push_back(&socket_cellSpeed);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

void StdALEUpdate::execute()
{
  computeCellSpeed();
  savePastNormals();
  updateNormalsData();
 }

//////////////////////////////////////////////////////////////////////////////

void StdALEUpdate::computeCellSpeed()
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

void StdALEUpdate::savePastNormals()
{
  CFLogDebugMin( "StdALEUpdate::savePastNormals() => begin" << "\n");

  DataHandle<InwardNormalsData*> normals = socket_normals.getDataHandle();

  DataHandle<InwardNormalsData*> pastNormals = socket_pastNormals.getDataHandle();

  cf_assert(normals.size() == pastNormals.size());
  for (CFuint i=0; i<normals.size(); ++i){
    *pastNormals[i] = *normals[i];
  }

  CFLogDebugMin( "StdALEUpdate::savePastNormals() => end" << "\n");
}

//////////////////////////////////////////////////////////////////////////////

void StdALEUpdate::updateNormalsData()
{
  CFLogDebugMin( "StdALEUpdate::updateNormalsData() => begin" << "\n");

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

    const CFuint firstCell = (*elemTypes)[iType].getStartIdx();
    const CFuint lastCell  = (*elemTypes)[iType].getEndIdx();

    (*dPtr).update(firstCell, lastCell, iType);
  }

  CFLogDebugMin( "StdALEUpdate::updateNormalsData() => end" << "\n");
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD
