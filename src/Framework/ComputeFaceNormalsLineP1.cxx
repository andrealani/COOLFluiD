#include "Framework/ComputeFaceNormalsLineP1.hh"
#include "Framework/MeshData.hh"
#include "Framework/PhysicalModel.hh"
#include "Environment/ObjectProvider.hh"
#include "Framework/TrsGeoWithNodesBuilder.hh"
#include "Framework/GeometricEntityPool.hh"
#include "Framework/PhysicalModel.hh"
#include "Framework/Framework.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<ComputeFaceNormalsLineP1,
			    ComputeNormals,
			    FrameworkLib>
computeFaceNormalsLineP1Provider("FaceLineP1");

//////////////////////////////////////////////////////////////////////////////

void ComputeFaceNormalsLineP1::operator() (const CFuint& iFirstElem,
					   const CFuint& iLastElem,
					   const CFuint& iType)
{
  SafePtr<TopologicalRegionSet> cells =
    MeshDataStack::getActive()->getTrs("InnerCells");

  SafePtr<ConnectivityTable<CFuint> > cellFaces =
    MeshDataStack::getActive()->getConnectivity("cellFaces");

  DataHandle<CFreal> normals = socket_normals->getDataHandle();
  DataHandle<CFint> isOutward = socket_isOutward->getDataHandle();

  GeometricEntityPool<TrsGeoWithNodesBuilder> geoBuilder;
  geoBuilder.setup();

  TrsGeoWithNodesBuilder::GeoData& geoData = geoBuilder.getDataGE();
  geoData.trs = cells;

  for (CFuint iElem = iFirstElem; iElem < iLastElem; ++iElem) {
    geoData.idx = iElem;
    GeometricEntity& cell = *geoBuilder.buildGE();
    cf_assert(cell.getNodes()->size() == 2);

    // we assume the normals like this: <-- 0------1 -->  if x goes -->
    CFuint faceID = (*cellFaces)(iElem, 0);
    if (isOutward[faceID] == -1) {
      normals[faceID] = -1.0;
      isOutward[faceID] = iElem;
    }

    faceID = (*cellFaces)(iElem, 1);
    if (isOutward[faceID] == -1) {
      normals[faceID] = 1.0;
      isOutward[faceID] = iElem;
    }
    geoBuilder.releaseGE();
  }
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
