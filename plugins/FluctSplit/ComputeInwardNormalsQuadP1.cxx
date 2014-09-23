

#include "ComputeInwardNormalsQuadP1.hh"
#include "Framework/MeshData.hh"
#include "InwardNormalsData.hh"
#include "Environment/ObjectProvider.hh"

#include "FluctSplit/FluctSplit.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<ComputeInwardNormalsQuadP1,
           ComputeNormals,
               FluctSplitModule>
computeInwardNormalsQuadP1Provider("InwardQuadP1");

//////////////////////////////////////////////////////////////////////////////

void ComputeInwardNormalsQuadP1::operator() (const CFuint& iFirstCell,
                                             const CFuint& iLastCell,
                                             const CFuint& iType)
{
  DataHandle<InwardNormalsData*> normals = socket_normals->getDataHandle();

  DataHandle<CFreal> normalsData = socket_normalsData->getDataHandle();

  DataHandle<CFuint> tempSize = socket_tempSize->getDataHandle();

  Common::SafePtr<TopologicalRegionSet> cells =
    MeshDataStack::getActive()->getTrs("InnerCells");

  const CFuint nbNodesInCell = cells->getNbNodesInGeo(iFirstCell);
  const CFuint nbDim = PhysicalModelStack::getActive()->getDim();
  cf_assert(nbNodesInCell == 4);
  cf_assert(nbDim == DIM_2D);

  const CFuint nbFacesInCell = 4;

  RealMatrix faceNormals(nbFacesInCell, DIM_2D);
  RealVector faceAreas(0.0, nbFacesInCell);
  RealMatrix nodalNormals(nbNodesInCell, DIM_2D);
  RealVector nodalAreas(nbNodesInCell);

#ifndef NDEBUG
  std::string filename("QuadP1Normals.dat");
  ofstream nfile(filename.c_str());
#endif

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

        // compute face normals
    faceNormals(0,0) = (*(*nodes)[1])[1] - (*(*nodes)[0])[1];
    faceNormals(0,1) = (*(*nodes)[0])[0] - (*(*nodes)[1])[0];

    faceNormals(1,0) = (*(*nodes)[2])[1] - (*(*nodes)[1])[1];
    faceNormals(1,1) = (*(*nodes)[1])[0] - (*(*nodes)[2])[0];

    faceNormals(2,0) = (*(*nodes)[3])[1] - (*(*nodes)[2])[1];
    faceNormals(2,1) = (*(*nodes)[2])[0] - (*(*nodes)[3])[0];

    faceNormals(3,0) = (*(*nodes)[0])[1] - (*(*nodes)[3])[1];
    faceNormals(3,1) = (*(*nodes)[3])[0] - (*(*nodes)[0])[0];

    for (CFuint iFace = 0; iFace < nbFacesInCell; ++iFace) {
      CFreal area = 0.;
      for (CFuint iDim = 0; iDim < DIM_2D; ++iDim) {
        area += faceNormals(iFace,iDim)*faceNormals(iFace, iDim);
      }
      faceAreas[iFace] = sqrt(area);
    }

    // computation of nodal normals
    nodalNormals(0,0) = (*(*nodes)[1])[1] - (*(*nodes)[3])[1];
    nodalNormals(0,1) = (*(*nodes)[3])[0] - (*(*nodes)[1])[0];

    nodalNormals(1,0) = (*(*nodes)[2])[1] - (*(*nodes)[0])[1];
    nodalNormals(1,1) = (*(*nodes)[0])[0] - (*(*nodes)[2])[0];

    nodalNormals(2,0) = - nodalNormals(0,0);
    nodalNormals(2,1) = - nodalNormals(0,1);

    nodalNormals(3,0) = - nodalNormals(1,0);
    nodalNormals(3,1) = - nodalNormals(1,1);

    for (CFuint iNode = 0; iNode < nbNodesInCell; ++iNode) {

      /**
       * @todo this can be substituted for
       *       \code areas[iNode] = nodalNormals.getRow(iNode).norm2() \endcode
       *       but it is expensive due to temporary hidden allocation.
       */
      CFreal sum = 0.;
      for (CFuint iDim = 0; iDim < nbDim; ++iDim) {
        sum += nodalNormals(iNode, iDim) * nodalNormals(iNode, iDim);
      }
      nodalAreas[iNode] = sqrt(sum);
    }

#ifndef NDEBUG
    nfile << "Cell: "     << "\n" << iCell        << "\n";
    nfile << "Nodes: "    << "\n";
    for(CFuint iNode = 0; iNode < nbNodesInCell; ++iNode){
       nfile << (*(*nodes)[iNode]) << "\n";
    }
    nfile << "NNormals: " << "\n" << nodalNormals ;
    nfile << "FNormals: " << "\n" << faceNormals  ;
    nfile << "NAreas: "   << "\n" << nodalAreas   << "\n";
    nfile << "FAreas: "   << "\n" << faceAreas    << "\n" << "\n";
#endif

    CFreal *const startFaceData = &normalsData[tempSize[0]];
    for (CFuint i = 0; i < faceNormals.size(); ++i) {
      normalsData[tempSize[0]++] = faceNormals[i];
    }

    CFreal *const startFaceArea = &normalsData[tempSize[0]];
    for (CFuint i = 0; i < faceAreas.size(); ++i) {
      normalsData[tempSize[0]++] = faceAreas[i];
    }

    CFreal *const startNodeData = &normalsData[tempSize[0]];
    for (CFuint i = 0; i < nodalNormals.size(); ++i) {
      normalsData[tempSize[0]++] = nodalNormals[i];
    }

    CFreal *const startNodeArea = &normalsData[tempSize[0]];
    for (CFuint i = 0; i < nodalAreas.size(); ++i) {
      normalsData[tempSize[0]++] = nodalAreas[i];
    }

    normals[iCell] = new InwardNormalsData(startFaceData,
                                           startFaceArea,
                                           startNodeData,
                                           startNodeArea,
                                           nbFacesInCell,
                                           iType);

    ///release the GeometricEntity
    geoBuilder->releaseGE();
  }
}

//////////////////////////////////////////////////////////////////////////////

void ComputeInwardNormalsQuadP1::update(const CFuint& iFirstCell,
                                        const CFuint& iLastCell,
                                        const CFuint& iType)
{
  DataHandle<InwardNormalsData*> normals = socket_normals->getDataHandle();

  Common::SafePtr<TopologicalRegionSet> cells =
    MeshDataStack::getActive()->getTrs("InnerCells");

  const CFuint nbNodesInCell = cells->getNbNodesInGeo(iFirstCell);
  const CFuint nbDim = PhysicalModelStack::getActive()->getDim();
  cf_assert(nbNodesInCell == 4);
  cf_assert(nbDim == DIM_2D);

  const CFuint nbFacesInCell = 4;

  RealMatrix faceNormals(nbFacesInCell, DIM_2D);
  RealVector faceAreas(0.0, nbFacesInCell);
  RealMatrix nodalNormals(nbNodesInCell, DIM_2D);
  RealVector nodalAreas(nbNodesInCell);

#ifndef NDEBUG
  std::string filename("QuadP1Normals.dat");
  ofstream nfile(filename.c_str());
#endif

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

    // compute face normals
    faceNormals(0,0) = (*(*nodes)[1])[1] - (*(*nodes)[0])[1];
    faceNormals(0,1) = (*(*nodes)[0])[0] - (*(*nodes)[1])[0];

    faceNormals(1,0) = (*(*nodes)[2])[1] - (*(*nodes)[1])[1];
    faceNormals(1,1) = (*(*nodes)[1])[0] - (*(*nodes)[2])[0];

    faceNormals(2,0) = (*(*nodes)[3])[1] - (*(*nodes)[2])[1];
    faceNormals(2,1) = (*(*nodes)[2])[0] - (*(*nodes)[3])[0];

    faceNormals(3,0) = (*(*nodes)[0])[1] - (*(*nodes)[3])[1];
    faceNormals(3,1) = (*(*nodes)[3])[0] - (*(*nodes)[0])[0];

    for (CFuint iFace = 0; iFace < nbFacesInCell; ++iFace) {
      CFreal area = 0.;
      for (CFuint iDim = 0; iDim < DIM_2D; ++iDim) {
        area += faceNormals(iFace,iDim)*faceNormals(iFace, iDim);
      }
      faceAreas[iFace] = sqrt(area);
    }

    // computation of nodal normals
    nodalNormals(0,0) = (*(*nodes)[1])[1] - (*(*nodes)[3])[1];
    nodalNormals(0,1) = (*(*nodes)[3])[0] - (*(*nodes)[1])[0];

    nodalNormals(1,0) = (*(*nodes)[2])[1] - (*(*nodes)[0])[1];
    nodalNormals(1,1) = (*(*nodes)[0])[0] - (*(*nodes)[2])[0];

    nodalNormals(2,0) = - nodalNormals(0,0);
    nodalNormals(2,1) = - nodalNormals(0,1);

    nodalNormals(3,0) = - nodalNormals(1,0);
    nodalNormals(3,1) = - nodalNormals(1,1);

    for (CFuint iNode = 0; iNode < nbNodesInCell; ++iNode) {

      /**
       * @todo this can be substituted for
       *       \code areas[iNode] = nodalNormals.getRow(iNode).norm2() \endcode
       *       but it is expensive due to temporary hidden allocation.
       */
      CFreal sum = 0.;
      for (CFuint iDim = 0; iDim < nbDim; ++iDim) {
        sum += nodalNormals(iNode, iDim) * nodalNormals(iNode, iDim);
      }
      nodalAreas[iNode] = sqrt(sum);
    }

#ifndef NDEBUG
    nfile << "Cell: "     << "\n" << iCell        << "\n";
    nfile << "Nodes: "    << "\n";
    for(CFuint iNode = 0; iNode < nbNodesInCell; ++iNode){
       nfile << (*(*nodes)[iNode]) << "\n";
    }
    nfile << "NNormals: " << "\n" << nodalNormals ;
    nfile << "FNormals: " << "\n" << faceNormals  ;
    nfile << "NAreas: "   << "\n" << nodalAreas   << "\n";
    nfile << "FAreas: "   << "\n" << faceAreas    << "\n" << "\n";
#endif


    normals[iCell]->setFaceNormals(faceNormals);
    normals[iCell]->setFaceAreas(faceAreas);
    normals[iCell]->setNodalNormals(nodalNormals);
    normals[iCell]->setNodalAreas(nodalAreas);
    normals[iCell]->setTypeID(iType);

    ///release the GeometricEntity
    geoBuilder->releaseGE();
  }
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
