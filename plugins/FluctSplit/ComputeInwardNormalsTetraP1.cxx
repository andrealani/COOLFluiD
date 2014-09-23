#include "Environment/ObjectProvider.hh"
#include "Framework/NormalsCalculator.hh"
#include "Framework/SubSystemStatus.hh"
#include "Framework/MeshData.hh"

#include "FluctSplit/FluctSplit.hh"
#include "FluctSplit/ComputeInwardNormalsTetraP1.hh"
#include "FluctSplit/InwardNormalsData.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<ComputeInwardNormalsTetraP1,
               ComputeNormals,
               FluctSplitModule>
computeInwardNormalsTetraP1Provider("InwardTetraP1");

//////////////////////////////////////////////////////////////////////////////

ComputeInwardNormalsTetraP1::ComputeInwardNormalsTetraP1() :
ComputeInwardNormals()
{
}

//////////////////////////////////////////////////////////////////////////////

void ComputeInwardNormalsTetraP1::operator() (const CFuint& iFirstCell,
                                              const CFuint& iLastCell,
                                              const CFuint& iType)
{
  NormalsCalculator calc;

  //set the IDS mapping a given state with the opposite faceID
  // in the local (within the element) connectivity
  vector<CFuint> ids(4);
  ids[0] = 2; // node 0 opposite to face 2
  ids[1] = 3; // node 1 opposite to face 3
  ids[2] = 1; // node 2 opposite to face 1
  ids[3] = 0; // node 3 opposite to face 0

  InwardNormalsData::setStateToFaceID(iType, ids);

  Common::SafePtr<TopologicalRegionSet> cells =
    MeshDataStack::getActive()->getTrs("InnerCells");

  DataHandle<InwardNormalsData*> normals = socket_normals->getDataHandle();

  DataHandle<CFreal> normalsData = socket_normalsData->getDataHandle();

  DataHandle<CFuint> tempSize = socket_tempSize->getDataHandle();

  const CFuint nbNodesInCell = cells->getNbNodesInGeo(iFirstCell);
  const CFuint nbDim = PhysicalModelStack::getActive()->getDim();
  cf_assert(nbNodesInCell == 4);
  cf_assert(nbDim == DIM_3D);

  const CFuint nbFacesInCell = 4;

  RealMatrix faceNormals(nbFacesInCell, DIM_3D);
  RealVector faceAreas(nbFacesInCell);
  RealMatrix coordinates(nbNodesInCell,DIM_3D);
  RealVector tmp(DIM_3D);
  CFreal* nodalNormals = CFNULL;
  CFreal* nodalAreas = CFNULL;

  Common::SafePtr<GeometricEntityPool<StdTrsGeoBuilder> >
  geoBuilder = getStdTrsGeoBuilder();

  StdTrsGeoBuilder::GeoData& geoData = geoBuilder->getDataGE();
  geoData.trs = cells;

  for (CFuint iCell = iFirstCell; iCell < iLastCell; ++iCell) {

    // build the GeometricEntity
    geoData.idx = iCell;

    GeometricEntity& currCell = *geoBuilder->buildGE();

    const vector<Node*>* const nodes = currCell.getNodes();
    cf_assert(nodes != CFNULL);
    cf_assert(nodes->size() == nbNodesInCell);

    // put coordinates in the coordinate matrix
    coordinates.setRow(*(*nodes)[0],0);
    coordinates.setRow(*(*nodes)[1],1);
    coordinates.setRow(*(*nodes)[2],2);
    coordinates.setRow(*(*nodes)[3],3);

    // computation of Face normals

    calc.computeTetraNormals(coordinates,faceNormals);

    // computation of face areas
    /**
     * @todo this can be substituted for
     *       \code areas[iNode] = faceNormals.getRow(iNode).norm2() \endcode
     *       but it is expensive due to temporary hidden allocation.
     */
    for (CFuint iFace = 0; iFace < nbFacesInCell; ++iFace) {
      CFreal sum = 0.;
      for (CFuint iDim = 0; iDim < nbDim; ++iDim) {
        sum += faceNormals(iFace, iDim) * faceNormals(iFace, iDim);
      }
      faceAreas[iFace] = sqrt(sum);
    }

    CFreal *const startFaceData = &normalsData[tempSize[0]];
    for (CFuint i = 0; i < faceNormals.size(); ++i) {
      normalsData[tempSize[0]++] = faceNormals[i];
    }

    CFreal *const startFaceArea = &normalsData[tempSize[0]];
    for (CFuint i = 0; i < faceAreas.size(); ++i) {
      normalsData[tempSize[0]++] = faceAreas[i];
    }

    normals[iCell] = new InwardNormalsData(startFaceData,
                                           startFaceArea,
                                           nodalNormals,
                                           nodalAreas,
                                           nbFacesInCell,
                                           iType);
    ///release the GeometricEntity
    geoBuilder->releaseGE();

  }
}

//////////////////////////////////////////////////////////////////////////////


void ComputeInwardNormalsTetraP1::update(const CFuint& iFirstCell,
                                         const CFuint& iLastCell,
                                         const CFuint& iType)
{
  NormalsCalculator calc;

  //set the IDS mapping a given state with the opposite faceID
  // in the local (within the element) connectivity
  vector<CFuint> ids(4);
  ids[0] = 2; // node 0 opposite to face 2
  ids[1] = 3; // node 1 opposite to face 3
  ids[2] = 1; // node 2 opposite to face 1
  ids[3] = 0; // node 3 opposite to face 0

  InwardNormalsData::setStateToFaceID(iType, ids);

  Common::SafePtr<TopologicalRegionSet> cells =
    MeshDataStack::getActive()->getTrs("InnerCells");

  DataHandle<InwardNormalsData*> normals = socket_normals->getDataHandle();

  const CFuint nbNodesInCell = cells->getNbNodesInGeo(iFirstCell);
  const CFuint nbDim = PhysicalModelStack::getActive()->getDim();
  cf_assert(nbNodesInCell == 4);
  cf_assert(nbDim == DIM_3D);

  const CFuint nbFacesInCell = 4;

  RealMatrix faceNormals(nbFacesInCell, DIM_3D);
  RealVector faceAreas(nbFacesInCell);
  RealMatrix coordinates(nbNodesInCell,DIM_3D);
  RealVector tmp(DIM_3D);

  Common::SafePtr<GeometricEntityPool<StdTrsGeoBuilder> >
  geoBuilder = getStdTrsGeoBuilder();

  StdTrsGeoBuilder::GeoData& geoData = geoBuilder->getDataGE();
  geoData.trs = cells;

  for (CFuint iCell = iFirstCell; iCell < iLastCell; ++iCell) {

    // build the GeometricEntity
    geoData.idx = iCell;

    GeometricEntity& currCell = *geoBuilder->buildGE();

    const vector<Node*>* const nodes = currCell.getNodes();
    cf_assert(nodes != CFNULL);
    cf_assert(nodes->size() == nbNodesInCell);

    // put coordinates in the coordinate matrix
    coordinates.setRow(*(*nodes)[0],0);
    coordinates.setRow(*(*nodes)[1],1);
    coordinates.setRow(*(*nodes)[2],2);
    coordinates.setRow(*(*nodes)[3],3);

    // computation of Face normals

    calc.computeTetraNormals(coordinates,faceNormals);

    // computation of face areas
    /**
     * @todo this can be substituted for
     *       \code areas[iNode] = faceNormals.getRow(iNode).norm2() \endcode
     *       but it is expensive due to temporary hidden allocation.
     */
    for (CFuint iFace = 0; iFace < nbFacesInCell; ++iFace) {
      CFreal sum = 0.;
      for (CFuint iDim = 0; iDim < nbDim; ++iDim) {
        sum += faceNormals(iFace, iDim) * faceNormals(iFace, iDim);
      }
      faceAreas[iFace] = sqrt(sum);
    }

    normals[iCell]->setFaceNormals(faceNormals);
    normals[iCell]->setFaceAreas(faceAreas);
    normals[iCell]->setTypeID(iType);

    ///release the GeometricEntity
    geoBuilder->releaseGE();
  }
}

//////////////////////////////////////////////////////////////////////////////

void ComputeInwardNormalsTetraP1::average(const CFuint& iFirstCell,
                                          const CFuint& iLastCell,
                                          const CFuint& iType)
{
///@todo here we should average the current and pastNormals instead of recomputing the normals at the intermediate
  NormalsCalculator calc;

   //set the IDS mapping a given state with the opposite faceID
   // in the local (within the element) connectivity
  vector<CFuint> ids(4);
  ids[0] = 2; // node 0 opposite to face 2
  ids[1] = 3; // node 1 opposite to face 3
  ids[2] = 1; // node 2 opposite to face 1
  ids[3] = 0; // node 3 opposite to face 0

  InwardNormalsData::setStateToFaceID(iType, ids);

  DataHandle<InwardNormalsData*> normals = socket_normals->getDataHandle();

  DataHandle<InwardNormalsData*> pastNormals = socket_pastNormals->getDataHandle();

  DataHandle<InwardNormalsData*> interNormals = socket_interNormals->getDataHandle();

  Common::SafePtr<TopologicalRegionSet> cells =
    MeshDataStack::getActive()->getTrs("InnerCells");

  const CFuint nbNodesInCell = cells->getNbNodesInGeo(iFirstCell);

  cf_assert(nbNodesInCell == 4);
  cf_assert(PhysicalModelStack::getActive()->getDim() == DIM_3D);

  const CFuint nbFacesInCell = 4;

   RealMatrix faceNormals(nbFacesInCell, DIM_3D);
   RealVector faceAreas(nbFacesInCell);
   RealMatrix coordinates(nbNodesInCell,DIM_3D);
   RealVector tmp(DIM_3D);

  const CFreal weight1 = SubSystemStatusStack::getActive()->getInnerDTRatio(0);
  const CFreal weight2 = SubSystemStatusStack::getActive()->getInnerDTRatio(1);

  for (CFuint iCell = iFirstCell; iCell < iLastCell; ++iCell) {

    faceNormals(0,0) = pastNormals[iCell]->getFaceNormComp(0,0)*weight2 + normals[iCell]->getFaceNormComp(0,0)*weight1;
    faceNormals(0,1) = pastNormals[iCell]->getFaceNormComp(0,1)*weight2 + normals[iCell]->getFaceNormComp(0,1)*weight1;
    faceNormals(0,2) = pastNormals[iCell]->getFaceNormComp(0,2)*weight2 + normals[iCell]->getFaceNormComp(0,2)*weight1;

    faceNormals(1,0) = pastNormals[iCell]->getFaceNormComp(1,0)*weight2 + normals[iCell]->getFaceNormComp(1,0)*weight1;
    faceNormals(1,1) = pastNormals[iCell]->getFaceNormComp(1,1)*weight2 + normals[iCell]->getFaceNormComp(1,1)*weight1;
    faceNormals(1,2) = pastNormals[iCell]->getFaceNormComp(1,2)*weight2 + normals[iCell]->getFaceNormComp(1,2)*weight1;

    faceNormals(2,0) = pastNormals[iCell]->getFaceNormComp(2,0)*weight2 + normals[iCell]->getFaceNormComp(2,0)*weight1;
    faceNormals(2,1) = pastNormals[iCell]->getFaceNormComp(2,1)*weight2 + normals[iCell]->getFaceNormComp(2,1)*weight1;
    faceNormals(2,2) = pastNormals[iCell]->getFaceNormComp(2,2)*weight2 + normals[iCell]->getFaceNormComp(2,2)*weight1;

    faceNormals(3,0) = pastNormals[iCell]->getFaceNormComp(3,0)*weight2 + normals[iCell]->getFaceNormComp(3,0)*weight1;
    faceNormals(3,1) = pastNormals[iCell]->getFaceNormComp(3,1)*weight2 + normals[iCell]->getFaceNormComp(3,1)*weight1;
    faceNormals(3,2) = pastNormals[iCell]->getFaceNormComp(3,2)*weight2 + normals[iCell]->getFaceNormComp(3,2)*weight1;

    for (CFuint iFace = 0; iFace < nbFacesInCell; ++iFace)
    {
      CFreal area = 0.;
      for (CFuint iDim = 0; iDim < DIM_3D; ++iDim) {
        area += faceNormals(iFace,iDim)*faceNormals(iFace, iDim);
      }
      faceAreas[iFace] = sqrt(area);
    }

    interNormals[iCell]->setFaceNormals(faceNormals);
    interNormals[iCell]->setFaceAreas(faceAreas);
    interNormals[iCell]->setTypeID(iType);

  }

}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
