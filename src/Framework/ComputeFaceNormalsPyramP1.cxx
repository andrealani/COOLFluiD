#include "Framework/ComputeFaceNormalsPyramP1.hh"
#include "Framework/MeshData.hh"
#include "Framework/GeometricEntity.hh"
#include "Environment/ObjectProvider.hh"
#include "Framework/TrsGeoWithNodesBuilder.hh"
#include "Framework/GeometricEntityPool.hh"
#include "Framework/NormalsCalculator.hh"
#include "Framework/Framework.hh"
#include "Framework/PhysicalModel.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<ComputeFaceNormalsPyramP1,
			    ComputeNormals,
			    FrameworkLib>
computeFaceNormalsPyramP1Provider("FacePyramP1");

//////////////////////////////////////////////////////////////////////////////

void ComputeFaceNormalsPyramP1::operator() (const CFuint& iFirstElem,
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

  NormalsCalculator calc;

  GeometricEntityPool<TrsGeoWithNodesBuilder> geoBuilder;
  geoBuilder.setup();

  TrsGeoWithNodesBuilder::GeoData& geoData = geoBuilder.getDataGE();
  geoData.trs = cells;

  const CFuint nbNodesInCell = cells->getNbNodesInGeo(iFirstElem);
  const CFuint nbDim = PhysicalModelStack::getActive()->getDim();
  cf_assert(nbNodesInCell == 5);
  cf_assert(nbDim == DIM_3D);
  //   const CFuint nbFacesInCell = 5;

  for (CFuint iElem = iFirstElem; iElem < iLastElem; ++iElem) {
    geoData.idx = iElem;
    GeometricEntity& cell = *geoBuilder.buildGE();
    const vector<Node*> *const nodes = cell.getNodes();

    // put coordinates in the coordinate matrix
    _coordinates.setRow(*(*nodes)[0],0);
    _coordinates.setRow(*(*nodes)[1],1);
    _coordinates.setRow(*(*nodes)[2],2);
    _coordinates.setRow(*(*nodes)[3],3);
    _coordinates.setRow(*(*nodes)[4],4);

    // computation of Face normals
    calc.computePyramNormals(_coordinates, _faceNormals);

    // Face 0321
    CFuint faceID = (*cellFaces)(iElem, 0);
    CFuint faceCompID = faceID*dim;

    if (isOutward[faceID] == -1) {
      for (CFuint i = 0; i < dim; ++i) {
        normals[faceCompID + i] = _faceNormals(0,i);
      }

      isOutward[faceID] = iElem;
    }

    // Face 014
    faceID = (*cellFaces)(iElem, 1);
    faceCompID = faceID*dim;

    if (isOutward[faceID] == -1) {
      for (CFuint i = 0; i < dim; ++i) {
        normals[faceCompID + i] = _faceNormals(1,i);
      }

      isOutward[faceID] = iElem;
    }

   // Face 124
    faceID = (*cellFaces)(iElem, 2);
    faceCompID = faceID*dim;

    if (isOutward[faceID] == -1) {
      for (CFuint i = 0; i < dim; ++i) {
        normals[faceCompID + i] = _faceNormals(2,i);
      }

      isOutward[faceID] = iElem;
    }

    // Face 234
    faceID = (*cellFaces)(iElem, 3);
    faceCompID = faceID*dim;

    if (isOutward[faceID] == -1) {
      for (CFuint i = 0; i < dim; ++i) {
        normals[faceCompID + i] = _faceNormals(3,i);
      }

      isOutward[faceID] = iElem;
    }

    // Face 043
    faceID = (*cellFaces)(iElem, 4);
    faceCompID = faceID*dim;

    if (isOutward[faceID] == -1) {
      for (CFuint i = 0; i < dim; ++i) {
        normals[faceCompID + i] = _faceNormals(4,i);
      }

      isOutward[faceID] = iElem;
    }
    geoBuilder.releaseGE();
  }
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
