#include "ComputeInwardNormalsHexaP1.hh"
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

Environment::ObjectProvider<ComputeInwardNormalsHexaP1,
               ComputeNormals,
               FluctSplitModule>
computeInwardNormalsHexaP1Provider("InwardHexaP1");

//////////////////////////////////////////////////////////////////////////////

void ComputeInwardNormalsHexaP1::operator() (const CFuint& iFirstCell,
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
  cf_assert(nbNodesInCell == 8);
  cf_assert(nbDim == DIM_3D);

  const CFuint nbFacesInCell = 6;

#ifndef NDEBUG
  //  std::string filename("HexaP1normals.dat");
  //  ofstream nfile(filename.c_str());
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

    // put _coordinates in the coordinate matrix
    _coordinates.setRow(*(*nodes)[0],0);
    _coordinates.setRow(*(*nodes)[1],1);
    _coordinates.setRow(*(*nodes)[2],2);
    _coordinates.setRow(*(*nodes)[3],3);
    _coordinates.setRow(*(*nodes)[4],4);
    _coordinates.setRow(*(*nodes)[5],5);
    _coordinates.setRow(*(*nodes)[6],6);
    _coordinates.setRow(*(*nodes)[7],7);

    // computation of Face normals

    calc.computeHexaNormals(_coordinates,_faceNormals);

    // computation of face areas
    /**
     * @todo this can be substituted for
     *       \code areas[iNode] = _faceNormals.getRow(iNode).norm2() \endcode
     *       but it is expensive due to temporary hidden allocation.
     */
    for (CFuint iFace = 0; iFace < nbFacesInCell; ++iFace) {
      CFreal sum = 0.;
      for (CFuint iDim = 0; iDim < nbDim; ++iDim) {
  sum += _faceNormals(iFace, iDim) * _faceNormals(iFace, iDim);
      }
      _faceAreas[iFace] = sqrt(sum);
    }

    // computation of nodal normals

    // NN0 = F0 + F2 + F5
    _tmp = _faceNormals.getRow<RealVector>(0) + _faceNormals.getRow<RealVector>(2);
    _tmp += _faceNormals.getRow<RealVector>(5);
    _nodalNormals.setRow(_tmp,0);

    // NN1 = F0 + F2 + F3
    _tmp = _faceNormals.getRow<RealVector>(0) + _faceNormals.getRow<RealVector>(2);
    _tmp += _faceNormals.getRow<RealVector>(3);
    _nodalNormals.setRow(_tmp,1);

    // NN2 = F0 + F3 + F4
    _tmp = _faceNormals.getRow<RealVector>(0) + _faceNormals.getRow<RealVector>(3);
    _tmp += _faceNormals.getRow<RealVector>(4);
    _nodalNormals.setRow(_tmp,2);

    // NN3 = F0 + F4 + F5
    _tmp = _faceNormals.getRow<RealVector>(0) + _faceNormals.getRow<RealVector>(4);
    _tmp += _faceNormals.getRow<RealVector>(5);
    _nodalNormals.setRow(_tmp,3);

    // NN4 = F1 + F2 + F5
    _tmp = _faceNormals.getRow<RealVector>(1) + _faceNormals.getRow<RealVector>(2);
    _tmp += _faceNormals.getRow<RealVector>(5);
    _nodalNormals.setRow(_tmp,4);

    // NN5 = F1 + F2 + F3
    _tmp = _faceNormals.getRow<RealVector>(1) + _faceNormals.getRow<RealVector>(2);
    _tmp += _faceNormals.getRow<RealVector>(3);
    _nodalNormals.setRow(_tmp,5);

    // NN6 = F1 + F3 + F4
    _tmp = _faceNormals.getRow<RealVector>(1) + _faceNormals.getRow<RealVector>(3);
    _tmp += _faceNormals.getRow<RealVector>(4);
    _nodalNormals.setRow(_tmp,6);

    // NN7 = F1 + F4 + F5
    _tmp = _faceNormals.getRow<RealVector>(1) + _faceNormals.getRow<RealVector>(4);
    _tmp += _faceNormals.getRow<RealVector>(5);
    _nodalNormals.setRow(_tmp,7);

   // computation of nodal areas
    /**
     * @todo this can be substituted for
     *       \code areas[iNode] = _nodalNormals.getRow<RealVector>(iNode).norm2() \endcode
     */
    for (CFuint iNode = 0; iNode < nbNodesInCell; ++iNode) {
      CFreal sum = 0.;
      for (CFuint iDim = 0; iDim < nbDim; ++iDim) {
  sum += _nodalNormals(iNode, iDim) * _nodalNormals(iNode, iDim);
      }
      _nodalAreas[iNode] = sqrt(sum);
    }

#ifndef NDEBUG
    /*nfile << "Cell: "     << "\n" << iCell        << "\n";
      nfile << "Nodes: "    << "\n";
      for(CFuint iNode = 0; iNode < nbNodesInCell; ++iNode){
      nfile << (*(*nodes)[iNode]) << "\n";
      }
      nfile << "NNormals: " << "\n" << _nodalNormals ;
      nfile << "FNormals: " << "\n" << _faceNormals  ;
      nfile << "NAreas: "   << "\n" << _nodalAreas   << "\n";
      nfile << "FAreas: "   << "\n" << _faceAreas    << "\n" << "\n";*/
#endif

    CFreal *const startFaceData = &normalsData[tempSize[0]];
    for (CFuint i = 0; i < _faceNormals.size(); ++i) {
      normalsData[tempSize[0]++] = _faceNormals[i];
    }

    CFreal *const startFaceArea = &normalsData[tempSize[0]];
    for (CFuint i = 0; i < _faceAreas.size(); ++i) {
      normalsData[tempSize[0]++] = _faceAreas[i];
    }

    CFreal *const startNodeData = &normalsData[tempSize[0]];
    for (CFuint i = 0; i < _nodalNormals.size(); ++i) {
      normalsData[tempSize[0]++] = _nodalNormals[i];
    }

    CFreal *const startNodeArea = &normalsData[tempSize[0]];
    for (CFuint i = 0; i < _nodalAreas.size(); ++i) {
      normalsData[tempSize[0]++] = _nodalAreas[i];
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

void ComputeInwardNormalsHexaP1::update(const CFuint& iFirstCell,
                                        const CFuint& iLastCell,
                                        const CFuint& iType)
{
  NormalsCalculator calc;

  Common::SafePtr<TopologicalRegionSet> cells =
    MeshDataStack::getActive()->getTrs("InnerCells");

  DataHandle<InwardNormalsData*> normals = socket_normals->getDataHandle();
  const CFuint nbDim = PhysicalModelStack::getActive()->getDim();
  const CFuint nbNodesInCell = cells->getNbNodesInGeo(iFirstCell);

  cf_assert(nbNodesInCell == 8);
  cf_assert(nbDim == DIM_3D);

  const CFuint nbFacesInCell = 6;

#ifndef NDEBUG
  //  std::string filename("HexaP1normals.dat");
  // ofstream nfile(filename.c_str());
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

    // put _coordinates in the coordinate matrix
    _coordinates.setRow(*(*nodes)[0],0);
    _coordinates.setRow(*(*nodes)[1],1);
    _coordinates.setRow(*(*nodes)[2],2);
    _coordinates.setRow(*(*nodes)[3],3);
    _coordinates.setRow(*(*nodes)[4],4);
    _coordinates.setRow(*(*nodes)[5],5);
    _coordinates.setRow(*(*nodes)[6],6);
    _coordinates.setRow(*(*nodes)[7],7);

    // computation of Face normals

    calc.computeHexaNormals(_coordinates,_faceNormals);

    // computation of face areas
    /**
     * @todo this can be substituted for
     *       \code areas[iNode] = _faceNormals.getRow<RealVector>(iNode).norm2() \endcode
     *       but it is expensive due to temporary hidden allocation.
     */
    for (CFuint iFace = 0; iFace < nbFacesInCell; ++iFace) {
      CFreal sum = 0.;
      for (CFuint iDim = 0; iDim < nbDim; ++iDim) {
  sum += _faceNormals(iFace, iDim) * _faceNormals(iFace, iDim);
      }
      _faceAreas[iFace] = sqrt(sum);
    }

    // computation of nodal normals

    // NN0 = F0 + F2 + F5
    _tmp = _faceNormals.getRow<RealVector>(0) + _faceNormals.getRow<RealVector>(2);
    _tmp += _faceNormals.getRow<RealVector>(5);
    _nodalNormals.setRow(_tmp,0);

    // NN1 = F0 + F2 + F3
    _tmp = _faceNormals.getRow<RealVector>(0) + _faceNormals.getRow<RealVector>(2);
    _tmp += _faceNormals.getRow<RealVector>(3);
    _nodalNormals.setRow(_tmp,1);

    // NN2 = F0 + F3 + F4
    _tmp = _faceNormals.getRow<RealVector>(0) + _faceNormals.getRow<RealVector>(3);
    _tmp += _faceNormals.getRow<RealVector>(4);
    _nodalNormals.setRow(_tmp,2);

    // NN3 = F0 + F4 + F5
    _tmp = _faceNormals.getRow<RealVector>(0) + _faceNormals.getRow<RealVector>(4);
    _tmp += _faceNormals.getRow<RealVector>(5);
    _nodalNormals.setRow(_tmp,3);

    // NN4 = F1 + F2 + F5
    _tmp = _faceNormals.getRow<RealVector>(1) + _faceNormals.getRow<RealVector>(2);
    _tmp += _faceNormals.getRow<RealVector>(5);
    _nodalNormals.setRow(_tmp,4);

    // NN5 = F1 + F2 + F3
    _tmp = _faceNormals.getRow<RealVector>(1) + _faceNormals.getRow<RealVector>(2);
    _tmp += _faceNormals.getRow<RealVector>(3);
    _nodalNormals.setRow(_tmp,5);

    // NN6 = F1 + F3 + F4
    _tmp = _faceNormals.getRow<RealVector>(1) + _faceNormals.getRow<RealVector>(3);
    _tmp += _faceNormals.getRow<RealVector>(4);
    _nodalNormals.setRow(_tmp,6);

    // NN7 = F1 + F4 + F5
    _tmp = _faceNormals.getRow<RealVector>(1) + _faceNormals.getRow<RealVector>(4);
    _tmp += _faceNormals.getRow<RealVector>(5);
    _nodalNormals.setRow(_tmp,7);

   // computation of nodal areas
    /**
     * @todo this can be substituted for
     *       \code areas[iNode] = _nodalNormals.getRow<RealVector>(iNode).norm2() \endcode
     */
    for (CFuint iNode = 0; iNode < nbNodesInCell; ++iNode) {
      CFreal sum = 0.;
      for (CFuint iDim = 0; iDim < nbDim; ++iDim) {
  sum += _nodalNormals(iNode, iDim) * _nodalNormals(iNode, iDim);
      }
      _nodalAreas[iNode] = sqrt(sum);
    }

#ifndef NDEBUG
    /*    nfile << "Cell: "     << "\n" << iCell        << "\n";
	  nfile << "Nodes: "    << "\n";
	  for(CFuint iNode = 0; iNode < nbNodesInCell; ++iNode){
	  nfile << (*(*nodes)[iNode]) << "\n";
	  }
	  nfile << "NNormals: " << "\n" << _nodalNormals ;
	  nfile << "FNormals: " << "\n" << _faceNormals  ;
	  nfile << "NAreas: "   << "\n" << _nodalAreas   << "\n";
	  nfile << "FAreas: "   << "\n" << _faceAreas    << "\n" << "\n";*/
#endif

    normals[iCell]->setFaceNormals(_faceNormals);
    normals[iCell]->setFaceAreas(_faceAreas);
    normals[iCell]->setNodalNormals(_nodalNormals);
    normals[iCell]->setNodalAreas(_nodalAreas);
    normals[iCell]->setTypeID(iType);

    ///release the GeometricEntity
    geoBuilder->releaseGE();
  }
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
