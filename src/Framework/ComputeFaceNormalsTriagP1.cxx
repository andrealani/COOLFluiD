#include "Framework/ComputeFaceNormalsTriagP1.hh"
#include "Framework/MeshData.hh"
#include "Environment/ObjectProvider.hh"
#include "Framework/TrsGeoWithNodesBuilder.hh"
#include "Framework/GeometricEntityPool.hh"
#include "Framework/Framework.hh"
#include "Framework/PhysicalModel.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  
  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<ComputeFaceNormalsTriagP1,
			    ComputeNormals,
			    FrameworkLib>
computeFaceNormalsTriagP1Provider("FaceTriagP1");

//////////////////////////////////////////////////////////////////////////////

void ComputeFaceNormalsTriagP1::operator() (const CFuint& iFirstElem,
                                            const CFuint& iLastElem,
                                            const CFuint& iType)
{
  SafePtr<TopologicalRegionSet> cells =
    MeshDataStack::getActive()->getTrs("InnerCells");

  SafePtr<ConnectivityTable<CFuint> > cellFaces =
    MeshDataStack::getActive()->getConnectivity("cellFaces");

  DataHandle<CFreal> normals = socket_normals->getDataHandle();
  DataHandle<CFint> isOutward = socket_isOutward->getDataHandle();

  const CFuint dim = PhysicalModelStack::getActive()->getDim();

  GeometricEntityPool<TrsGeoWithNodesBuilder> geoBuilder;
  geoBuilder.setup();
  
  TrsGeoWithNodesBuilder::GeoData& geoData = geoBuilder.getDataGE();
  geoData.trs = cells;
  
  for (CFuint iElem = iFirstElem; iElem < iLastElem; ++iElem) {
    geoData.idx = iElem;
    GeometricEntity& cell = *geoBuilder.buildGE();
    const vector<Node*> *const nodes = cell.getNodes();

    cf_assert(nodes->size() == 3);

    CFuint faceID = (*cellFaces)(iElem, 0);
    CFuint faceCompID = faceID*dim;

    if (isOutward[faceID] == -1) {
      normals[faceCompID] = - (*(*nodes)[0])[1] + (*(*nodes)[1])[1];
      normals[faceCompID + 1] = - (*(*nodes)[1])[0] + (*(*nodes)[0])[0];

      isOutward[faceID] = iElem;
    }

    faceID = (*cellFaces)(iElem, 1);
    faceCompID = faceID*dim;

    if (isOutward[faceID] == -1) {
      normals[faceCompID] =  - (*(*nodes)[1])[1] + (*(*nodes)[2])[1];
      normals[faceCompID + 1] = - (*(*nodes)[2])[0] + (*(*nodes)[1])[0];

      isOutward[faceID] = iElem;
    }

    faceID = (*cellFaces)(iElem, 2);
    faceCompID = faceID*dim;

    if (isOutward[faceID] == -1) {
      normals[faceCompID] = - (*(*nodes)[2])[1] + (*(*nodes)[0])[1];
      normals[faceCompID + 1] = - (*(*nodes)[0])[0] + (*(*nodes)[2])[0];

      isOutward[faceID] = iElem;
    }
    geoBuilder.releaseGE();
  }
 }

//////////////////////////////////////////////////////////////////////////////

    } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
