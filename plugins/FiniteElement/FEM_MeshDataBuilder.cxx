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
#include "FiniteElement/FiniteElement.hh"
#include "FiniteElement/FEM_MeshDataBuilder.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteElement {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<FEM_MeshDataBuilder,
               MeshDataBuilder,
         FiniteElementModule,
               1>
fem_MeshDataBuilderProvider("FiniteElement");

//////////////////////////////////////////////////////////////////////////////

FEM_MeshDataBuilder::FEM_MeshDataBuilder(const std::string& name) :
  MeshDataBuilder(name)
{
}

//////////////////////////////////////////////////////////////////////////////

FEM_MeshDataBuilder::~FEM_MeshDataBuilder()
{
}

//////////////////////////////////////////////////////////////////////////////

void FEM_MeshDataBuilder::releaseMemory()
{
  MeshDataBuilder::releaseMemory();
}

//////////////////////////////////////////////////////////////////////////////

void FEM_MeshDataBuilder::setMaxNbStatesInCell()
{

  CFuint maxNbStates = 0;

  // get the TRS of Cells
  std::vector< Common::SafePtr<TopologicalRegionSet> > cells = MeshDataStack::getActive()->getTrsList();
  for(CFuint iTrs=0; iTrs < cells.size(); ++iTrs)
  {
    if(cells[iTrs]->hasTag("inner")){
      const CFuint nbGeos = cells[iTrs]->getLocalNbGeoEnts();
      for (CFuint iGeo = 0; iGeo < nbGeos; ++iGeo) {
        if (cells[iTrs]->getNbStatesInGeo(iGeo) > maxNbStates) {
          maxNbStates = cells[iTrs]->getNbStatesInGeo(iGeo);
        }
      }
    }
  }

  MeshDataStack::getActive()->Statistics().setMaxNbStatesInCell(maxNbStates);
}

//////////////////////////////////////////////////////////////////////////////

void FEM_MeshDataBuilder::setMaxNbNodesInCell()
{
  CFuint maxNbNodes = 0;

  // get the TRS of Cells
  std::vector< Common::SafePtr<TopologicalRegionSet> > cells = MeshDataStack::getActive()->getTrsList();
  for(CFuint iTrs=0; iTrs < cells.size(); ++iTrs)
  {
    if(cells[iTrs]->hasTag("inner")){

      const CFuint nbGeos = cells[iTrs]->getLocalNbGeoEnts();
      for (CFuint iGeo = 0; iGeo < nbGeos; ++iGeo) {
        if (cells[iTrs]->getNbNodesInGeo(iGeo) > maxNbNodes) {
          maxNbNodes = cells[iTrs]->getNbNodesInGeo(iGeo);
        }
      }
    }
  }

  MeshDataStack::getActive()->Statistics().setMaxNbNodesInCell(maxNbNodes);
}

//////////////////////////////////////////////////////////////////////////////

void FEM_MeshDataBuilder::setMaxNbFacesInCell()
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

void FEM_MeshDataBuilder::setMapGeoToTrs()
{
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

    } // namespace FiniteElement

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
