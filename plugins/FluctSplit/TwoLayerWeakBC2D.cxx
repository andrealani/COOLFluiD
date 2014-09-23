#include "TwoLayerWeakBC2D.hh"
#include "InwardNormalsData.hh"
#include "Framework/SubSystemStatus.hh"
#include "Framework/MeshData.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

TwoLayerWeakBC2D::TwoLayerWeakBC2D(const std::string& name) :
  WeakBC(name),
  socket_interRhs("interRhs"),
  socket_pastNormals("pastNormals",false),    // only needed if moving
  socket_interNormals("interNormals",false),  // only needed if moving
  socket_nodes("nodes",false),                // only needed if moving
  socket_pastNodes("pastNodes",false),        // only needed if moving
  socket_pastStates("pastStates"),
  socket_interStates("interStates"),
  socket_isUpdated("isUpdated"),
  m_dt(),
  m_layer(0),
  m_wallSpeed()
{
}

//////////////////////////////////////////////////////////////////////////////

TwoLayerWeakBC2D::~TwoLayerWeakBC2D()
{
}

//////////////////////////////////////////////////////////////////////////////

void TwoLayerWeakBC2D::setup()
{
  WeakBC::setup();

  m_wallSpeed.resize(PhysicalModelStack::getActive()->getDim());
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
TwoLayerWeakBC2D::needsSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSink> > result = WeakBC::needsSockets();

  result.push_back(&socket_interRhs);
  result.push_back(&socket_pastNormals);
  result.push_back(&socket_interNormals);
  result.push_back(&socket_nodes);
  result.push_back(&socket_pastNodes);
  result.push_back(&socket_pastStates);
  result.push_back(&socket_interStates);
  result.push_back(&socket_isUpdated);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

void TwoLayerWeakBC2D::executeOnTrs()
{
  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  const CFreal oEminusAlpha = 1. - m_alpha;
  // states forming the ghost cell
  vector<State*> states(2);
  RealVector flux0(nbEqs);
  RealVector flux1(nbEqs);

  // get the data handle for the states
  DataHandle<State*,GLOBAL> allStates = socket_states.getDataHandle();
  DataHandle< CFreal> interRhs = socket_interRhs.getDataHandle();
  DataHandle< CFreal> rhs = socket_rhs.getDataHandle();

  // get the data handle for the past states
  DataHandle<State*> pastStates = socket_pastStates.getDataHandle();

  // get the data handle for the intermidiate states
  DataHandle<State*> interStates = socket_interStates.getDataHandle();

  DataHandle<bool> isUpdated = socket_isUpdated.getDataHandle();

  // create two ghost states (no local ID is needed for these states)
  State* gstate0 = new State();
  State* gstate1 = new State();

  State* state0 = new State();
  State* state1 = new State();

  bool isOnMesh = false;
  Node* coord0 = new Node(isOnMesh);
  Node* coord1 = new Node(isOnMesh);

  Common::SafePtr<Framework::TopologicalRegionSet> const faces = getCurrentTRS();

  Common::SafePtr<GeometricEntityPool<StdTrsGeoBuilder> >
  geoBuilder = getMethodData().getStdTrsGeoBuilder();

  StdTrsGeoBuilder::GeoData& geoData = geoBuilder->getDataGE();
  geoData.trs = faces;
  const CFuint nbFaces = faces->getLocalNbGeoEnts();

  for (CFuint iFace = 0; iFace < nbFaces; ++iFace) {
    CFLogDebugMed( "Computing iFace = " << iFace << "\n");

    // build the GeometricEntity
    geoData.idx = iFace;
    GeometricEntity& currFace = *geoBuilder->buildGE();

    const CFuint stateID0 = currFace.getState(0)->getLocalID();
    const CFuint stateID1 = currFace.getState(1)->getLocalID();

    /// Cell K1
    m_layer = 0;

    // Set dt1
    m_dt = SubSystemStatusStack::getActive()->getInnerDT(m_layer);

    // set the face normal
    if (SubSystemStatusStack::getActive()->isMovingMesh()){
      setMovingFaceNormal(currFace.getID());
    }
    else{
      setFaceNormal(currFace.getID());
    }

    // Unsteady, state=0.5*(pastState+interState)
    for (CFuint j=0; j<nbEqs ; ++j)
    {
      (*state0)[j] = 0.5*((*interStates[stateID0])[j] + (*pastStates[stateID0])[j]);
      (*state1)[j] = 0.5*((*interStates[stateID1])[j] + (*pastStates[stateID1])[j]);
    }

    RealVector stateCoord0 = ((*interStates[stateID0]).getCoordinates()) ;
    stateCoord0 += ((*pastStates[stateID0]).getCoordinates());
    stateCoord0 *= 0.5;
    RealVector stateCoord1 = ((*interStates[stateID1]).getCoordinates()) ;
    stateCoord1 += ((*pastStates[stateID1]).getCoordinates());
    stateCoord1 *= 0.5;

    *coord0 = stateCoord0;
    *coord1 = stateCoord1;

    (*state0).setSpaceCoordinates(coord0);
    (*state1).setSpaceCoordinates(coord1);

    // set the appropriate values in the ghost states
    setGhostState(*state0, *gstate0);
    setGhostState(*state1, *gstate1);

    // set ghost state == first state
    // corresponding B-state == second state
    states[0] = gstate0;
    states[1] = state0;

    // compute the corresponding fluctuation and the update coeff
    computeFlux(states, flux0);

    // set ghost state == first state
    // corresponding B-state == second state
    states[0] = gstate1;
    states[1] = state1;

    // compute the corresponding fluctuation and the update coeff
    computeFlux(states, flux1);

    // Multiply by dt1
    flux0 *= m_dt;
    flux1 *= m_dt;

    for (CFuint iEq = 0; iEq < nbEqs; ++iEq)
    {
      // distribute corrections to the boundary nodes
      if (!isUpdated[stateID0])
      {
        interRhs(stateID0, iEq, nbEqs) +=
          (m_alpha*flux0[iEq] + oEminusAlpha*flux1[iEq]);

        m_flagState[stateID0] = true;
      }

      if (!isUpdated[stateID1])
      {
        interRhs(stateID1, iEq, nbEqs) +=
            (m_alpha*flux1[iEq] + oEminusAlpha*flux0[iEq]);

        m_flagState[stateID1] = true;
      }
    }

    /// Cell K2
    m_layer = 1;
    // Set dt1
    m_dt = SubSystemStatusStack::getActive()->getInnerDT(m_layer);

    // set the face normal
    if (SubSystemStatusStack::getActive()->isMovingMesh())
    {
      setMovingFaceNormal(currFace.getID());
    }
    else
    {
      setFaceNormal(currFace.getID());
    }

    // Unsteady, state=0.5*(pastState+interState)
    for (CFuint j=0; j<nbEqs ; ++j)
    {
      (*state0)[j] = 0.5*((*interStates[stateID0])[j] + (*allStates[stateID0])[j]);
      (*state1)[j] = 0.5*((*interStates[stateID1])[j] + (*allStates[stateID1])[j]);
    }

    stateCoord0 = ((*interStates[stateID0]).getCoordinates()) ;
    stateCoord0 += ((*allStates[stateID0]).getCoordinates());
    stateCoord0 *= 0.5;

    stateCoord1 = ((*interStates[stateID1]).getCoordinates()) ;
    stateCoord1 += ((*allStates[stateID1]).getCoordinates());
    stateCoord1 *= 0.5;

    *coord0 = stateCoord0;
    *coord1 = stateCoord1;

    (*state0).setSpaceCoordinates(coord0);
    (*state1).setSpaceCoordinates(coord1);

    // set the appropriate values in the ghost states
    setGhostState(*state0, *gstate0);
    setGhostState(*state1, *gstate1);

    // set ghost state == first state
    // corresponding B-state == second state
    states[0] = gstate0;
    states[1] = state0;
    // compute the corresponding fluctuation and the update coeff
    computeFlux(states, flux0);

    // set ghost state == first state
    // corresponding B-state == second state
    states[0] = gstate1;
    states[1] = state1;

    // compute the corresponding fluctuation and the update coeff
    computeFlux(states, flux1);

    // Multiply by dt2
    flux0 *= m_dt;
    flux1 *= m_dt;

    for (CFuint iEq = 0; iEq < nbEqs; ++iEq)
    {
      // distribute corrections to the boundary nodes
      if (!isUpdated[stateID0])
      {
        rhs(stateID0, iEq, nbEqs) +=
          (m_alpha*flux0[iEq] + oEminusAlpha*flux1[iEq]);

        m_flagState[stateID0] = true;
      }

      if (!isUpdated[stateID1])
      {
        rhs(stateID1, iEq, nbEqs) +=
            (m_alpha*flux1[iEq] + oEminusAlpha*flux0[iEq]);

        m_flagState[stateID1] = true;
      }
    }
  //Release the geometric entity
  geoBuilder->releaseGE();
  }

  // this could be done in a more efficient way ...
  const CFuint nbStates = m_flagState.size();
  for (CFuint i = 0; i < nbStates; ++i) {
    if (m_flagState[i] == true) {
      isUpdated[i] = true;
    }
  }

  deletePtr(coord0);
  deletePtr(coord1);
  deletePtr(state0);
  deletePtr(state1);
  deletePtr(gstate0);
  deletePtr(gstate1);
 }

//////////////////////////////////////////////////////////////////////////////

void TwoLayerWeakBC2D::setMovingFaceNormal(const CFuint faceID)
{
  DataHandle< InwardNormalsData*> interNormals = socket_interNormals.getDataHandle();
  DataHandle< InwardNormalsData*> pastNormals = socket_pastNormals.getDataHandle();
  DataHandle< InwardNormalsData*> normals = socket_normals.getDataHandle();
  DataHandle<Common::Trio<CFuint, CFuint, Common::SafePtr<Framework::TopologicalRegionSet> > >
    faceNeighCell = socket_faceNeighCell.getDataHandle();

  const CFuint cellTrsID = faceNeighCell[faceID].first;
  const CFuint iFaceLocal = faceNeighCell[faceID].second;
  Common::SafePtr<TopologicalRegionSet> cellTrs = faceNeighCell[faceID].third;
  const CFuint cellLocalID = cellTrs->getLocalGeoID(cellTrsID);

  const CFuint dim = PhysicalModelStack::getActive()->getDim();
  for (CFuint iDim = 0; iDim < dim; ++iDim) {
  // external normal is already inward for ghost state
  // The normal is the average of the past and current
  // normals...
  if (m_layer == 0){
    m_faceNormal[iDim] = -0.5*(interNormals[cellLocalID]->
			      getFaceNormComp(iFaceLocal, iDim)+
			      pastNormals[cellLocalID]->getFaceNormComp
			      (iFaceLocal, iDim));
  }

  if (m_layer == 1) {
    m_faceNormal[iDim] = -0.5*(normals[cellLocalID]->getFaceNormComp
			      (iFaceLocal, iDim)+interNormals[cellLocalID]->
			      getFaceNormComp(iFaceLocal, iDim));
  }
  }

}

//////////////////////////////////////////////////////////////////////////////

} // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
