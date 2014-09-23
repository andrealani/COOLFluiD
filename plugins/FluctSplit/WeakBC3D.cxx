#include "WeakBC3D.hh"
#include "InwardNormalsData.hh"
#include "Framework/MeshData.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

WeakBC3D::WeakBC3D(const std::string& name) :
  WeakBC(name),
  socket_isUpdated("isUpdated"),
  socket_pastStates("pastStates", false),
  m_fluxes()
{
}

//////////////////////////////////////////////////////////////////////////////

WeakBC3D::~WeakBC3D()
{
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
WeakBC3D::needsSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSink> > result = WeakBC::needsSockets();

  result.push_back(&socket_isUpdated);
  result.push_back(&socket_pastStates);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

void WeakBC3D::setup()
{
  WeakBC::setup();

  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  // size the fluxes
  // triangular faces ONLY for  now
  m_fluxes.resize(3);
  for (CFuint i = 0; i < m_fluxes.size(); ++i) {
    m_fluxes[i].resize(nbEqs);
  }
  
  
  state0new = new State();
  state1new = new State();
  state2new = new State();
}

//////////////////////////////////////////////////////////////////////////////

void WeakBC3D::executeOnTrs()
{
  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  const CFreal third = 1./3.;
  const CFreal halfOEminAlphaThird = 0.5*(1. - m_alpha)*third;
  const CFreal alphaThird = m_alpha*third;
  bool isUnsteady = false;	
  CFreal dt = SubSystemStatusStack::getActive()->getDT();
  if(dt < 0.) dt = 1.;
  else isUnsteady = true;
  
  // states forming the ghost cell
  // @todo this works only for triangular boundary faces
  vector<State*> states(2);

  // get the data handle for the states
  DataHandle<State*,GLOBAL> allStates = socket_states.getDataHandle();
  // get the data handle for the is Updated flags
  DataHandle<bool> isUpdated = socket_isUpdated.getDataHandle();
  // get the data handle for the rhs
  DataHandle<CFreal> rhs = socket_rhs.getDataHandle();

  // create two ghost states (no local ID is needed for these states)
  State* gstate0 = new State();
  State* gstate1 = new State();
  State* gstate2 = new State();

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

    ///@todo extend this for quads...
    cf_assert(currFace.nbNodes() == 3);

    // set the face normal
    setFaceNormal(currFace.getID());

    State *const state0 = currFace.getState(0);
    State *const state1 = currFace.getState(1);
    State *const state2 = currFace.getState(2);

    const CFuint stateID0 = state0->getLocalID();
    const CFuint stateID1 = state1->getLocalID();
    const CFuint stateID2 = state2->getLocalID();



        if (isUnsteady){
           DataHandle<State*> pastStates = socket_pastStates.getDataHandle();

             // Unsteady, state=0.5*(pastState+interState)

               for (CFuint j=0; j<nbEqs ; ++j){

               (*state0new)[j] = 0.5*((*state0)[j] + (*pastStates[stateID0])[j]);
		(*state1new)[j] = 0.5*((*state1)[j] + (*pastStates[stateID1])[j]);
		(*state2new)[j] = 0.5*((*state2)[j] + (*pastStates[stateID2])[j]);
           }
        }
        else {
           // clone the states in order to use the indexes

          state0new->clone(*state0);
	  state1new->clone(*state1);
	  state2new->clone(*state2);
          }

    // set the appropriate values in the ghost states
    setGhostState(*state0new, *gstate0);
    setGhostState(*state1new, *gstate1);
    setGhostState(*state2new, *gstate2);

    // set ghost state == first state
    // corresponding B-state == second state
    states[0] = gstate0;
    states[1] = state0;
    // compute the corresponding fluctuation and the update coeff
    computeFlux(states, m_fluxes[0]);

    // set ghost state == first state
    // corresponding B-state == second state
    states[0] = gstate1;
    states[1] = state1;
    // compute the corresponding fluctuation and the update coeff
    computeFlux(states, m_fluxes[1]);

    // set ghost state == first state
    // corresponding B-state == second state
    states[0] = gstate2;
    states[1] = state2;
    // compute the corresponding fluctuation and the update coeff
    computeFlux(states, m_fluxes[2]);

    for (CFuint iEq = 0; iEq < nbEqs; ++iEq) {
      // distribute corrections to the boundary nodes
      if (!isUpdated[stateID0]) {
        rhs(stateID0, iEq, nbEqs) +=
          alphaThird*m_fluxes[0][iEq] + halfOEminAlphaThird*
          (m_fluxes[1][iEq] + m_fluxes[2][iEq]);

        m_flagState[stateID0] = true;
      }

      if (!isUpdated[stateID1]) {
        rhs(stateID1, iEq, nbEqs) +=
          alphaThird*m_fluxes[1][iEq] + halfOEminAlphaThird*
          (m_fluxes[0][iEq] + m_fluxes[2][iEq]);

        m_flagState[stateID1] = true;
      }

      if (!isUpdated[stateID2]) {
        rhs(stateID2, iEq, nbEqs) +=
          alphaThird*m_fluxes[2][iEq] + halfOEminAlphaThird*
          (m_fluxes[0][iEq] + m_fluxes[1][iEq]);

        m_flagState[stateID2] = true;
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

  deletePtr(gstate0);
  deletePtr(gstate1);
  deletePtr(gstate2);
}

//////////////////////////////////////////////////////////////////////////////

} // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
