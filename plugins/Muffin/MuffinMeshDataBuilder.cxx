#include "Environment/ObjectProvider.hh"
#include "Framework/MapGeoToTrsAndIdx.hh"
#include "Framework/BadFormatException.hh"
#include "Framework/MeshData.hh"
#include "Framework/PhysicalModel.hh"
#include "Muffin/MuffinModule.hh"
#include "Muffin/MuffinMeshDataBuilder.hh"
#include "Muffin/Muffin.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace COOLFluiD::Framework;

namespace COOLFluiD {
  namespace Muffin {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider< MuffinMeshDataBuilder,MeshDataBuilder,MuffinModule,1 >
  camusMeshDataBuilderProvider("Muffin");

//////////////////////////////////////////////////////////////////////////////

MuffinMeshDataBuilder::MuffinMeshDataBuilder(const std::string& name) :
  MeshDataBuilder(name)
{
}

//////////////////////////////////////////////////////////////////////////////

MuffinMeshDataBuilder::~MuffinMeshDataBuilder()
{
}

//////////////////////////////////////////////////////////////////////////////

void MuffinMeshDataBuilder::releaseMemory()
{
  MeshDataBuilder::releaseMemory();
}

//////////////////////////////////////////////////////////////////////////////

void MuffinMeshDataBuilder::setMaxNbStatesInCell()
{
  MeshDataStack::getActive()->Statistics().setMaxNbStatesInCell(
    PhysicalModelStack::getActive()->getNbEq()+1 );
}

//////////////////////////////////////////////////////////////////////////////

void MuffinMeshDataBuilder::setMaxNbNodesInCell()
{
  MeshDataStack::getActive()->Statistics().setMaxNbNodesInCell(
    PhysicalModelStack::getActive()->getDim()+1 );
}

//////////////////////////////////////////////////////////////////////////////

void MuffinMeshDataBuilder::setMaxNbFacesInCell()
{
  MeshDataStack::getActive()->Statistics().setMaxNbFacesInCell(
    PhysicalModelStack::getActive()->getDim()+1 );
}

//////////////////////////////////////////////////////////////////////////////

void MuffinMeshDataBuilder::setMapGeoToTrs()
{
  CFuint nbBFaces = 0;
  std::vector< Common::SafePtr< TopologicalRegionSet > > trs =
    MeshDataStack::getActive()->getTrsList();
  for (CFuint i = 0; i < trs.size(); ++i) {
    Common::SafePtr< TopologicalRegionSet > currTrs = trs[i];
    if (currTrs->hasTag("face"))
    {
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
    Common::SafePtr< TopologicalRegionSet > currTrs = trs[i];
    if (currTrs->hasTag("face"))
    {
      const CFuint nbTrsFaces = currTrs->getLocalNbGeoEnts();
      for (CFuint iFace = 0; iFace < nbTrsFaces; ++iFace) {
        const CFuint faceID = currTrs->getLocalGeoID(iFace);
        mapGeoToTrs->setMappingData
          (faceID, currTrs, iFace, isOnBoundary);
      }
    }
  }

  MeshDataStack::getActive()->storeMapGeoToTrs("MapFacesToTrs", mapGeoToTrs);

}

//////////////////////////////////////////////////////////////////////////////

  }  // end of namespace Muffin
}  // end of namespace COOLFluiD

