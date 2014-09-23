

#ifndef NDEBUG
#include <fstream>
#endif

#include "ComputeInwardNormalsTriagP1.hh"
#include "InwardNormalsData.hh"
#include "Framework/MeshData.hh"
#include "Environment/ObjectProvider.hh"

#include "Framework/SubSystemStatus.hh"
#include "FluctSplit/FluctSplit.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<ComputeInwardNormalsTriagP1,
           ComputeNormals,
               FluctSplitModule>
computeInwardNormalsTriagP1Provider("InwardTriagP1");

//////////////////////////////////////////////////////////////////////////////

ComputeInwardNormalsTriagP1::ComputeInwardNormalsTriagP1() :
ComputeInwardNormals()
{
}

//////////////////////////////////////////////////////////////////////////////

void ComputeInwardNormalsTriagP1::operator() (const CFuint& iFirstCell,
                                              const CFuint& iLastCell,
                                              const CFuint& iType)
 {
   //set the IDS mapping a given state with the opposite faceID
   // in the local (within the element) connectivity
   vector<CFuint> ids(3);
   ids[0] = 1; // node 0 opposite to face 1
   ids[1] = 2; // node 1 opposite to face 2
   ids[2] = 0; // node 2 opposite to face 0

   InwardNormalsData::setStateToFaceID(iType, ids);

   Common::SafePtr<TopologicalRegionSet> cells =
    MeshDataStack::getActive()->getTrs("InnerCells");

   DataHandle<InwardNormalsData*> normals = socket_normals->getDataHandle();

   DataHandle<CFreal> normalsData = socket_normalsData->getDataHandle();

   DataHandle<CFuint> tempSize = socket_tempSize->getDataHandle();

   cf_assert(cells->getNbNodesInGeo(iFirstCell) == 3);
   cf_assert(PhysicalModelStack::getActive()->getDim() == DIM_2D);

   const CFuint nbFacesInCell = 3;
   const CFuint dim = PhysicalModelStack::getActive()->getDim();

   RealMatrix faceNormals(nbFacesInCell, dim);
   RealVector faceAreas(0.0, nbFacesInCell);
   CFreal* nodalNormals = CFNULL;
   CFreal* nodalAreas = CFNULL;

//#ifndef NDEBUG
//  std::string filename("TriagP1normals.dat");
//  ofstream nfile(filename.c_str());
//#endif

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
    cf_assert(nodes->size() == cells->getNbNodesInGeo(iFirstCell));

     faceNormals(0,0) = (*(*nodes)[0])[1] - (*(*nodes)[1])[1];
     faceNormals(0,1) = (*(*nodes)[1])[0] - (*(*nodes)[0])[0];

     faceNormals(1,0) = (*(*nodes)[1])[1] - (*(*nodes)[2])[1];
     faceNormals(1,1) = (*(*nodes)[2])[0] - (*(*nodes)[1])[0];

     faceNormals(2,0) = (*(*nodes)[2])[1] - (*(*nodes)[0])[1];
     faceNormals(2,1) = (*(*nodes)[0])[0] - (*(*nodes)[2])[0];

     faceNormals *= -1.;

     for (CFuint iFace = 0; iFace < nbFacesInCell; ++iFace) {
       CFreal area = 0.;
       for (CFuint iDim = 0; iDim < dim; ++iDim) {
   area += faceNormals(iFace,iDim)*faceNormals(iFace, iDim);
       }
       faceAreas[iFace] = sqrt(area);
     }

/*#ifndef NDEBUG
    nfile << "Cell: "     << "\n" << iCell        << "\n";
    nfile << "Nodes: "    << "\n";
    for(CFuint iNode = 0; iNode < cells->getNbNodesInGeo(iFirstCell); ++iNode){
      nfile << (*(*nodes)[iNode]) << "\n";
    }
    nfile << "NNormals: " << "\n" << nodalNormals ;
    nfile << "FNormals: " << "\n" << faceNormals  ;
    nfile << "NAreas: "   << "\n" << nodalAreas   << "\n";
    nfile << "FAreas: "   << "\n" << faceAreas    << "\n" << "\n";
#endif*/

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
    // release the GeometricEntity
    geoBuilder->releaseGE();
   }
 }

//////////////////////////////////////////////////////////////////////////////

void ComputeInwardNormalsTriagP1::update(const CFuint& iFirstCell,
                                         const CFuint& iLastCell,
                                         const CFuint& iType)
{
   //set the IDS mapping a given state with the opposite faceID
   // in the local (within the element) connectivity
   vector<CFuint> ids(3);
   ids[0] = 1; // node 0 opposite to face 1
   ids[1] = 2; // node 1 opposite to face 2
   ids[2] = 0; // node 2 opposite to face 0

   InwardNormalsData::setStateToFaceID(iType, ids);

   Common::SafePtr<TopologicalRegionSet> cells =
    MeshDataStack::getActive()->getTrs("InnerCells");

   DataHandle<InwardNormalsData*> normals = socket_normals->getDataHandle();

   cf_assert(cells->getNbNodesInGeo(iFirstCell) == 3);
   cf_assert(PhysicalModelStack::getActive()->getDim() == DIM_2D);

   const CFuint nbFacesInCell = 3;
   const CFuint dim = PhysicalModelStack::getActive()->getDim();

   RealMatrix faceNormals(nbFacesInCell, dim);
   RealVector faceAreas(0.0, nbFacesInCell);

#ifndef NDEBUG
  std::string filename("TriagP1normals.dat");
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
    cf_assert(nodes->size() == cells->getNbNodesInGeo(iFirstCell));

    faceNormals(0,0) = (*(*nodes)[0])[1] - (*(*nodes)[1])[1];
    faceNormals(0,1) = (*(*nodes)[1])[0] - (*(*nodes)[0])[0];

    faceNormals(1,0) = (*(*nodes)[1])[1] - (*(*nodes)[2])[1];
    faceNormals(1,1) = (*(*nodes)[2])[0] - (*(*nodes)[1])[0];

    faceNormals(2,0) = (*(*nodes)[2])[1] - (*(*nodes)[0])[1];
    faceNormals(2,1) = (*(*nodes)[0])[0] - (*(*nodes)[2])[0];

    faceNormals *= -1.;

    for (CFuint iFace = 0; iFace < nbFacesInCell; ++iFace)
    {
      CFreal area = 0.;
      for (CFuint iDim = 0; iDim < dim; ++iDim) {
        area += faceNormals(iFace,iDim)*faceNormals(iFace, iDim);
      }
      faceAreas[iFace] = sqrt(area);
    }

#ifndef NDEBUG
    nfile << "Cell: "     << "\n" << iCell        << "\n";
    nfile << "Nodes: "    << "\n";
    for(CFuint iNode = 0; iNode < cells->getNbNodesInGeo(iFirstCell); ++iNode){
      nfile << (*(*nodes)[iNode]) << "\n";
    }
    nfile << "FNormals: " << "\n" << faceNormals  ;
    nfile << "FAreas: "   << "\n" << faceAreas    << "\n" << "\n";
#endif

    normals[iCell]->setFaceNormals(faceNormals);
    normals[iCell]->setFaceAreas(faceAreas);
    normals[iCell]->setTypeID(iType);

    // release the GeometricEntity
    geoBuilder->releaseGE();
   }
 }

//////////////////////////////////////////////////////////////////////////////

void ComputeInwardNormalsTriagP1::average(const CFuint& iFirstCell,
                                          const CFuint& iLastCell,
                                          const CFuint& iType)
{
   //set the IDS mapping a given state with the opposite faceID
   // in the local (within the element) connectivity
   vector<CFuint> ids(3);
   ids[0] = 1; // node 0 opposite to face 1
   ids[1] = 2; // node 1 opposite to face 2
   ids[2] = 0; // node 2 opposite to face 0

   InwardNormalsData::setStateToFaceID(iType, ids);

   Common::SafePtr<TopologicalRegionSet> cells =
    MeshDataStack::getActive()->getTrs("InnerCells");

   DataHandle<InwardNormalsData*> normals = socket_normals->getDataHandle();

   DataHandle<InwardNormalsData*> pastNormals = socket_pastNormals->getDataHandle();

   DataHandle<InwardNormalsData*> interNormals = socket_interNormals->getDataHandle();

   cf_assert(cells->getNbNodesInGeo(iFirstCell) == 3);
   cf_assert(PhysicalModelStack::getActive()->getDim() == DIM_2D);

   const CFuint nbFacesInCell = 3;
   const CFuint dim = PhysicalModelStack::getActive()->getDim();
   const CFreal weight1 = SubSystemStatusStack::getActive()->getInnerDTRatio(0);
   const CFreal weight2 = SubSystemStatusStack::getActive()->getInnerDTRatio(1);

   RealMatrix faceNormals(nbFacesInCell, dim);
   RealVector faceAreas(0.0, nbFacesInCell);

#ifndef NDEBUG
  std::string filename("TriagP1normals.dat");
  ofstream nfile(filename.c_str());
#endif

//   Common::SafePtr<GeometricEntityPool<StdTrsGeoBuilder> >
//   geoBuilder = getStdTrsGeoBuilder();
//
//   StdTrsGeoBuilder::GeoData& geoData = geoBuilder->getDataGE();
//   geoData.trs = cells;

  for (CFuint iCell = iFirstCell; iCell < iLastCell; ++iCell) {

    // build the GeometricEntity
//    geoData.idx = iCell;

//    GeometricEntity& currCell = *geoBuilder->buildGE();

    faceNormals(0,0) = pastNormals[iCell]->getFaceNormComp(0,0)*weight2 + normals[iCell]->getFaceNormComp(0,0)*weight1;
    faceNormals(0,1) = pastNormals[iCell]->getFaceNormComp(0,1)*weight2 + normals[iCell]->getFaceNormComp(0,1)*weight1;

    faceNormals(1,0) = pastNormals[iCell]->getFaceNormComp(1,0)*weight2 + normals[iCell]->getFaceNormComp(1,0)*weight1;
    faceNormals(1,1) = pastNormals[iCell]->getFaceNormComp(1,1)*weight2 + normals[iCell]->getFaceNormComp(1,1)*weight1;

    faceNormals(2,0) = pastNormals[iCell]->getFaceNormComp(2,0)*weight2 + normals[iCell]->getFaceNormComp(2,0)*weight1;
    faceNormals(2,1) = pastNormals[iCell]->getFaceNormComp(2,1)*weight2 + normals[iCell]->getFaceNormComp(2,1)*weight1;

    for (CFuint iFace = 0; iFace < nbFacesInCell; ++iFace)
    {
      CFreal area = 0.;
      for (CFuint iDim = 0; iDim < dim; ++iDim) {
        area += faceNormals(iFace,iDim)*faceNormals(iFace, iDim);
      }
      faceAreas[iFace] = sqrt(area);
    }

#ifndef NDEBUG
    nfile << "Elem: "     << "\n" << iCell        << "\n";
    nfile << "FNormals: " << "\n" << faceNormals  ;
    nfile << "FAreas: "   << "\n" << faceAreas    << "\n" << "\n";
#endif

    interNormals[iCell]->setFaceNormals(faceNormals);
    interNormals[iCell]->setFaceAreas(faceAreas);
    interNormals[iCell]->setTypeID(iType);

    // release the GeometricEntity
//     geoBuilder->releaseGE();
   }
 }

//////////////////////////////////////////////////////////////////////////////


    } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
