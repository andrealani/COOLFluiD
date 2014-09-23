#include "Framework/Storage.hh"
#include "CreateBoundaryNodalNormals.hh"
#include "InwardNormalsData.hh"
#include "Framework/MeshData.hh"
#include "Framework/PhysicalModel.hh"
#include "Common/CFMap.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

CreateBoundaryNodalNormals::CreateBoundaryNodalNormals
(SafePtr<GeometricEntityPool<StdTrsGeoBuilder> > geoBuilder) :
  _geoBuilder(geoBuilder),
  socket_normals("Null"),
  socket_faceNeighCell("Null")
{
}

//////////////////////////////////////////////////////////////////////////////

CreateBoundaryNodalNormals::~CreateBoundaryNodalNormals()
{
}

//////////////////////////////////////////////////////////////////////////////

void CreateBoundaryNodalNormals::create
(const std::vector<Common::SafePtr<Framework::TopologicalRegionSet> >& trsList,
 std::vector< std::vector<RealVector> >& bcNormals)
{
  const CFuint nbTrs = trsList.size();
  for (CFuint iTrs = 0; iTrs < nbTrs; ++iTrs) {
    createOnTrs(trsList[iTrs], bcNormals[iTrs]);
  }
}

//////////////////////////////////////////////////////////////////////////////

void CreateBoundaryNodalNormals::setDataSockets(
                                  DataSocketSink< InwardNormalsData* > normalsSocket,
                                  DataSocketSink< Common::Trio<CFuint, CFuint, Common::SafePtr<Framework::TopologicalRegionSet> > > faceNeighCellSocket )
{

  socket_normals = normalsSocket;
  socket_faceNeighCell = faceNeighCellSocket;
}

//////////////////////////////////////////////////////////////////////////////

void CreateBoundaryNodalNormals::createOnTrs
(SafePtr<TopologicalRegionSet> trs,
 vector<RealVector>& bcNormals)
{
  CFAUTOTRACE;

  CFLog(VERBOSE, "Creating boundary nodal normals for TRS : " << trs->getName() << "\n");

  SafePtr<vector<CFuint> > statesSet = trs->getStatesInTrs();
  const CFuint nbStates = statesSet->size();
  CFMap<CFuint, CFuint> mapStateID(nbStates);
  CFuint maxID = 0;
  CFuint countState = 0;

  //create mapping State* <-> ID of the state
  vector<CFuint>::const_iterator it;
  for (it = statesSet->begin(); it != statesSet->end(); ++it) {
    mapStateID.insert(*it, countState);
    ++countState;
    if ((*it) > maxID) {
      maxID = (*it);
    }
  }
  mapStateID.sortKeys();

  //allocate memory for the boundary normals
  //if maxID == N, nb of possible IDs is N+1 (0 is included)
  bcNormals.resize(maxID + 1);

  const CFuint dim = PhysicalModelStack::getActive()->getDim();

  vector<CFuint>::const_iterator itd;
  for (itd = statesSet->begin(); itd != statesSet->end(); ++itd) {
    const CFuint stateID = (*itd);
    bcNormals[stateID].resize(dim);
    bcNormals[stateID] = 0.;
  }

  // temporary vector
  RealVector nxyz(0., dim);

  // temporary storage for the sum of the normals referencing each node
  vector<RealVector> normals(nbStates);
  for (CFuint i = 0; i < normals.size(); ++i) {
    normals[i].resize(dim);
    normals[i] = 0.;
  }

  // temporary storage for the abs of the normals
  RealVector absN(0., nbStates);

  // handle to the normals
  DataHandle< InwardNormalsData*> normalsData = socket_normals.getDataHandle();

  // handle to the neighbor cell
  DataHandle<Common::Trio<CFuint, CFuint, Common::SafePtr<Framework::TopologicalRegionSet> > >
    faceNeighCell = socket_faceNeighCell.getDataHandle();

  StdTrsGeoBuilder::GeoData& geoData = _geoBuilder->getDataGE();
  geoData.trs = trs;

  const CFuint nbFaces = trs->getLocalNbGeoEnts();
  for (CFuint iFace = 0; iFace < nbFaces; ++iFace) {
    geoData.idx = iFace;
    GeometricEntity *const currFace = _geoBuilder->buildGE();

    const CFuint faceID = currFace->getID();
    const CFuint cellTrsID = faceNeighCell[faceID].first;
    const CFuint iFaceLocal = faceNeighCell[faceID].second;
    SafePtr<TopologicalRegionSet> cellTrs = faceNeighCell[faceID].third;
    const CFuint cellLocalID = cellTrs->getLocalGeoID(cellTrsID);
    for (CFuint iDim = 0; iDim < dim; ++iDim) {
      nxyz[iDim] = -normalsData[cellLocalID]->
	getFaceNormComp(iFaceLocal, iDim);
    }

    // loop over all the states in the face to distribute contributions
    // to the corresponding nodal normal from this face
    const vector<State*> *const states = currFace->getStates();
    const CFuint nbStatesInFace = states->size();
    for (CFuint iState = 0; iState < nbStatesInFace; ++iState) {
      const CFuint stateID = mapStateID.find
	((*states)[iState]->getLocalID());
      normals[stateID] += nxyz;
      absN[stateID]    += normalsData[cellLocalID]->getAreaFace(iFaceLocal);
    }

    // release the face
    _geoBuilder->releaseGE();
  }

  for (CFuint i = 0; i < nbStates; ++i) {
    normals[i] /= absN[i];
    bcNormals[(*statesSet)[i]] = normals[i];
  }
}

//////////////////////////////////////////////////////////////////////////////

 } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
