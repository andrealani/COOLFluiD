#include <ctime>
#include <set>

#include "Common/BadValueException.hh"

#include "Environment/ObjectProvider.hh"

#include "Framework/MeshData.hh"
#include "Framework/LocalConnectionData.hh"
#include "Framework/GeometricEntityProvider.hh"
#include "Framework/Face.hh"
#include "Framework/Cell.hh"
#include "Framework/PhysicalModel.hh"
#include "Framework/BadFormatException.hh"
#include "Framework/MapGeoToTrsAndIdx.hh"

#include "UFEM/UFEM.hh"
#include "UFEM/UFEM_MeshDataBuilder.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace UFEM {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<UFEM_MeshDataBuilder, MeshDataBuilder, UFEMPlugin, 1>
aUFEM_MeshDataBuilderProvider("UFEM");

//////////////////////////////////////////////////////////////////////////////

UFEM_MeshDataBuilder::UFEM_MeshDataBuilder(const std::string& name) :
  MeshDataBuilder(name)
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

UFEM_MeshDataBuilder::~UFEM_MeshDataBuilder()
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

void UFEM_MeshDataBuilder::releaseMemory()
{
  CFAUTOTRACE;
  MeshDataBuilder::releaseMemory();
}

//////////////////////////////////////////////////////////////////////////////

void UFEM_MeshDataBuilder::setMaxNbStatesInCell()
{
  CFAUTOTRACE;
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

void UFEM_MeshDataBuilder::setMaxNbNodesInCell()
{
  CFAUTOTRACE;
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

void UFEM_MeshDataBuilder::setMaxNbFacesInCell()
{
  CFAUTOTRACE;
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

void UFEM_MeshDataBuilder::setMapGeoToTrs()
{
  CFAUTOTRACE;
  CFuint nbBFaces = 0;
  vector<Common::SafePtr<TopologicalRegionSet> > trs = MeshDataStack::getActive()->getTrsList();
  for (CFuint i = 0; i < trs.size(); ++i) {
    SafePtr< TopologicalRegionSet > currTrs = trs[i];
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
    SafePtr<TopologicalRegionSet> currTrs = trs[i];
    if (currTrs->hasTag("face"))
    {
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

    } // namespace UFEM

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
