#include <ctime>
#include <set>


#include "Environment/ObjectProvider.hh"

#include "Common/BadValueException.hh"
#include "Framework/MeshData.hh"
#include "Framework/LocalConnectionData.hh"
#include "Framework/GeometricEntityProvider.hh"
#include "Framework/Face.hh"
#include "Framework/Cell.hh"
#include "Framework/PhysicalModel.hh"
#include "Framework/BadFormatException.hh"
#include "Framework/MapGeoToTrsAndIdx.hh"
#include "FluctSplit/FluctSplit.hh"
#include "FluctSplit/RDS_MeshDataBuilder.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<RDS_MeshDataBuilder,
               MeshDataBuilder,
               FluctSplitModule,
               1>
RDS_MeshDataBuilderProvider("RDS");

//////////////////////////////////////////////////////////////////////////////

RDS_MeshDataBuilder::RDS_MeshDataBuilder(const std::string& name) :
  MeshDataBuilder(name)
{
}

//////////////////////////////////////////////////////////////////////////////

RDS_MeshDataBuilder::~RDS_MeshDataBuilder()
{
}

//////////////////////////////////////////////////////////////////////////////

void RDS_MeshDataBuilder::releaseMemory()
{
  MeshDataBuilder::releaseMemory();
}

//////////////////////////////////////////////////////////////////////////////

void RDS_MeshDataBuilder::setMaxNbStatesInCell()
{
  Common::SafePtr<TopologicalRegionSet> cells =
    MeshDataStack::getActive()->getTrs("InnerCells");

  CFuint maxNbStates = 0;
  const CFuint nbGeos = cells->getLocalNbGeoEnts();
  for (CFuint iGeo = 0; iGeo < nbGeos; ++iGeo) {
    if (cells->getNbStatesInGeo(iGeo) > maxNbStates) {
      maxNbStates = cells->getNbStatesInGeo(iGeo);
    }
  }

  MeshDataStack::getActive()->Statistics().setMaxNbStatesInCell(maxNbStates);
}

//////////////////////////////////////////////////////////////////////////////

void RDS_MeshDataBuilder::setMaxNbNodesInCell()
{
  Common::SafePtr<TopologicalRegionSet> cells =
    MeshDataStack::getActive()->getTrs("InnerCells");

  CFuint maxNbNodes = 0;
  const CFuint nbGeos = cells->getLocalNbGeoEnts();
  for (CFuint iGeo = 0; iGeo < nbGeos; ++iGeo) {
    if (cells->getNbNodesInGeo(iGeo) > maxNbNodes) {
      maxNbNodes = cells->getNbNodesInGeo(iGeo);
    }
  }

  MeshDataStack::getActive()->Statistics().setMaxNbNodesInCell(maxNbNodes);
}

//////////////////////////////////////////////////////////////////////////////

void RDS_MeshDataBuilder::setMaxNbFacesInCell()
{
  SafePtr<vector<ElementTypeData> > elementType =
    getCFmeshData().getElementTypeData();

  CFuint maxNbFaces = 0;
  for (CFuint iType = 0; iType < getNbElementTypes(); ++iType) {
    CFuint nbFaces = LocalConnectionData::getInstance().getNbFacesInShape
      ((*elementType)[iType].getGeoShape());
    if (nbFaces > maxNbFaces) {
      maxNbFaces = nbFaces;
    }
  }

  MeshDataStack::getActive()->Statistics().setMaxNbFacesInCell(maxNbFaces);
}

//////////////////////////////////////////////////////////////////////////////

void RDS_MeshDataBuilder::setMapGeoToTrs()
{
  CFuint nbBFaces = 0;
  vector<Common::SafePtr<TopologicalRegionSet> > trs = MeshDataStack::getActive()->getTrsList();
  for (CFuint i = 0; i < trs.size(); ++i) {
    SafePtr<TopologicalRegionSet> currTrs = trs[i];
    if (currTrs->getName() != "InnerCells") {
      nbBFaces += currTrs->getLocalNbGeoEnts();
    }
  }
  // set the number of faces which is the number of boundary faces in MeshData
  MeshDataStack::getActive()->Statistics().setNbFaces(nbBFaces);

  MapGeoToTrsAndIdx* mapGeoToTrs = new MapGeoToTrsAndIdx();
  const CFuint nbFaces = MeshDataStack::getActive()->Statistics().getNbFaces();
  mapGeoToTrs->resize(nbFaces);

  const bool isOnBoundary = true;
  for (CFuint i = 0; i < trs.size(); ++i) {
    SafePtr<TopologicalRegionSet> currTrs = trs[i];
    if (currTrs->getName() != "InnerCells") {
      const CFuint nbTrsFaces = currTrs->getLocalNbGeoEnts();
      for (CFuint iFace = 0; iFace < nbTrsFaces; ++iFace) {
        const CFuint faceID = currTrs->getLocalGeoID(iFace);
        mapGeoToTrs->setMappingData(faceID, currTrs, iFace, isOnBoundary);
      }
    }
  }

  MeshDataStack::getActive()->storeMapGeoToTrs("MapFacesToTrs", mapGeoToTrs);

}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
