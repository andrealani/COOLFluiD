

#include "ComputeInwardNormalsPyramP1.hh"
#include "Framework/ComputeNormals.hh"
#include "Framework/MeshData.hh"
#include "InwardNormalsData.hh"
#include "Environment/ObjectProvider.hh"
#include "Framework/NormalsCalculator.hh"

#include "FluctSplit/FluctSplit.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<ComputeInwardNormalsPyramP1,
               ComputeNormals,
               FluctSplitModule>
computeInwardNormalsPyramP1Provider("InwardPyramP1");

//////////////////////////////////////////////////////////////////////////////

void ComputeInwardNormalsPyramP1::operator() (const CFuint& iFirstCell,
                                              const CFuint& iLastCell,
                                              const CFuint& iType)
{
  NormalsCalculator calc;

  Common::SafePtr<TopologicalRegionSet> cells =
    MeshDataStack::getActive()->getTrs("InnerCells");

  DataHandle<InwardNormalsData*> normals = socket_normals->getDataHandle();

  DataHandle<CFreal> normalsData = socket_normalsData->getDataHandle();

  DataHandle<CFuint> tempSize = socket_tempSize->getDataHandle();

  const CFuint nbNodesInCell = cells->getNbNodesInGeo(iFirstCell);
  const CFuint nbDim = PhysicalModelStack::getActive()->getDim();
  cf_assert(nbNodesInCell == 5);
  cf_assert(nbDim == DIM_3D);

  const CFuint nbFacesInCell = 5;

  RealMatrix nodalNormals(nbNodesInCell, DIM_3D);
  RealVector nodalAreas  (nbNodesInCell);

  RealMatrix faceNormals(nbFacesInCell, DIM_3D);
  RealVector faceAreas  (nbFacesInCell);

  RealMatrix coordinates(nbNodesInCell,DIM_3D);

  RealVector tmp(DIM_3D);

#ifndef NDEBUG
  std::string filename("PyramP1normals.dat");
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

    // put coordinates in the coordinate matrix
    coordinates.setRow(*(*nodes)[0],0);
    coordinates.setRow(*(*nodes)[1],1);
    coordinates.setRow(*(*nodes)[2],2);
    coordinates.setRow(*(*nodes)[3],3);
    coordinates.setRow(*(*nodes)[4],4);

    // computation of Face normals

    calc.computePyramNormals(coordinates,faceNormals);

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

    // computation of nodal faces

    // NN0 = - F2 - F3
    tmp = faceNormals.getRow<RealVector>(2) + faceNormals.getRow<RealVector>(3);
    nodalNormals.setRow(tmp,0);

    // NN1 = - F4 - F3
    tmp = faceNormals.getRow<RealVector>(4) + faceNormals.getRow<RealVector>(3);
    nodalNormals.setRow(tmp,1);

    // NN2 = - F4 - F1
    tmp = faceNormals.getRow<RealVector>(4) + faceNormals.getRow<RealVector>(1);
    nodalNormals.setRow(tmp,2);

    // NN3 = - F2 - F1
    tmp = faceNormals.getRow<RealVector>(2) + faceNormals.getRow<RealVector>(1);
    nodalNormals.setRow(tmp,3);

    // NN4 = - F0
    nodalNormals.setRow(faceNormals.getRow<RealVector>(0),4);

    nodalNormals *= -1.;

    // computation of nodal areas
    /**
     * @todo this can be substituted for
     *       \code areas[iNode] = nodalNormals.getRow<RealVector>(iNode).norm2() \endcode
     *       but it is expensive due to temporary hidden allocation.
     */
    for (CFuint iNode = 0; iNode < nbNodesInCell; ++iNode) {
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

void ComputeInwardNormalsPyramP1::update(const CFuint& iFirstCell,
                                         const CFuint& iLastCell,
                                         const CFuint& iType)
{
  NormalsCalculator calc;

  Common::SafePtr<TopologicalRegionSet> cells =
    MeshDataStack::getActive()->getTrs("InnerCells");


  DataHandle<InwardNormalsData*> normals = socket_normals->getDataHandle();
  const CFuint nbDim = PhysicalModelStack::getActive()->getDim();
  const CFuint nbNodesInCell = cells->getNbNodesInGeo(iFirstCell);

  cf_assert(nbNodesInCell == 5);
  cf_assert(nbDim == DIM_3D);

  const CFuint nbFacesInCell = 5;

  RealMatrix nodalNormals(nbNodesInCell, DIM_3D);
  RealVector nodalAreas  (nbNodesInCell);

  RealMatrix faceNormals(nbFacesInCell, DIM_3D);
  RealVector faceAreas  (nbFacesInCell);

  RealMatrix coordinates(nbNodesInCell,DIM_3D);

  RealVector tmp(DIM_3D);

#ifndef NDEBUG
  std::string filename("PyramP1normals.dat");
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

    // put coordinates in the coordinate matrix
    coordinates.setRow(*(*nodes)[0],0);
    coordinates.setRow(*(*nodes)[1],1);
    coordinates.setRow(*(*nodes)[2],2);
    coordinates.setRow(*(*nodes)[3],3);
    coordinates.setRow(*(*nodes)[4],4);

    // computation of Face normals

    calc.computePyramNormals(coordinates,faceNormals);

    // computation of face areas
    /**
     * @todo this can be substituted for
     *       \code areas[iNode] = faceNormals.getRow<RealVector>(iNode).norm2() \endcode
     *       but it is expensive due to temporary hidden allocation.
     */
    for (CFuint iFace = 0; iFace < nbFacesInCell; ++iFace) {
      CFreal sum = 0.;
      for (CFuint iDim = 0; iDim < nbDim; ++iDim) {
        sum += faceNormals(iFace, iDim) * faceNormals(iFace, iDim);
      }
      faceAreas[iFace] = sqrt(sum);
    }

    // computation of nodal faces

    // NN0 = - F2 - F3
    tmp = faceNormals.getRow<RealVector>(2) + faceNormals.getRow<RealVector>(3);
    nodalNormals.setRow(tmp,0);

    // NN1 = - F4 - F3
    tmp = faceNormals.getRow<RealVector>(4) + faceNormals.getRow<RealVector>(3);
    nodalNormals.setRow(tmp,1);

    // NN2 = - F4 - F1
    tmp = faceNormals.getRow<RealVector>(4) + faceNormals.getRow<RealVector>(1);
    nodalNormals.setRow(tmp,2);

    // NN3 = - F2 - F1
    tmp = faceNormals.getRow<RealVector>(2) + faceNormals.getRow<RealVector>(1);
    nodalNormals.setRow(tmp,3);

    // NN4 = - F0
    nodalNormals.setRow(faceNormals.getRow<RealVector>(0),4);

    nodalNormals *= -1.;

    // computation of nodal areas
    /**
     * @todo this can be substituted for
     *       \code areas[iNode] = nodalNormals.getRow<RealVector>(iNode).norm2() \endcode
     *       but it is expensive due to temporary hidden allocation.
     */
    for (CFuint iNode = 0; iNode < nbNodesInCell; ++iNode) {
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
