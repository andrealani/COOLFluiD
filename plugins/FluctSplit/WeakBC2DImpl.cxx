#include "WeakBC2DImpl.hh"
#include "InwardNormalsData.hh"
#include "Framework/LSSMatrix.hh"
#include "Framework/BlockAccumulator.hh"
#include "Framework/MeshData.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

WeakBC2DImpl::WeakBC2DImpl(const std::string& name) :
  WeakBC(name),
  socket_isUpdated("isUpdated"),
  socket_pastStates("pastStates", false)
{
}

//////////////////////////////////////////////////////////////////////////////

WeakBC2DImpl::~WeakBC2DImpl()
{
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
WeakBC2DImpl::needsSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSink> > result = WeakBC::needsSockets();

  result.push_back(&socket_isUpdated);
  result.push_back(&socket_pastStates);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

void WeakBC2DImpl::unsetup()
{

  for ( CFuint iState = 0; iState < m_gstate.size(); ++iState){
    deletePtr(m_gstate[iState]);
  }

  WeakBC::unsetup();
}

//////////////////////////////////////////////////////////////////////////////

void WeakBC2DImpl::setup()
{
  WeakBC::setup();

  Common::SafePtr<GeometricEntityPool<StdTrsGeoBuilder> >
  geoBuilder = getMethodData().getStdTrsGeoBuilder();

  StdTrsGeoBuilder::GeoData& geoData = geoBuilder->getDataGE();

  std::vector< Common::SafePtr<TopologicalRegionSet> >& trss = getTrsList();

  vector< SafePtr<TopologicalRegionSet> >::iterator itrs = trss.begin();
  CFuint max_nb_states = 0 ;
  for (; itrs != trss.end(); ++itrs)
  {
    SafePtr<TopologicalRegionSet> faces = *itrs;

    geoData.trs = faces;

    const CFuint nbFaces = faces->getLocalNbGeoEnts();

    for (CFuint iFace = 0; iFace < nbFaces; ++iFace)
    {
      // build the GeometricEntity
      geoData.idx = iFace;
      GeometricEntity& currFace = *geoBuilder->buildGE();
      const CFuint nb_state = currFace.nbStates();
      max_nb_states = std::max(max_nb_states, nb_state);
  geoBuilder->releaseGE();
    }
  }
  // allocation of the flux and Jacobian flux
  m_flux.resize(max_nb_states);
  m_fluxJacob.resize(max_nb_states);

  // create nbstate ghost states (no local ID is needed for these states)
  m_gstate.resize(max_nb_states);

  m_state.resize(max_nb_states);
  const bool isGhost = true;
  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
   for ( CFuint iState = 0; iState < max_nb_states; ++iState){
      m_flux[iState].resize(nbEqs);
      m_fluxJacob[iState].resize(nbEqs, nbEqs);
      m_gstate[iState] = new State(m_flux[iState], isGhost);
    }
}

//////////////////////////////////////////////////////////////////////////////

void WeakBC2DImpl::executeOnTrs()
{
  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  const CFreal oEminusAlpha = 1. - m_alpha;

  SafePtr<LinearSystemSolver> lss =
    getMethodData().getLinearSystemSolver()[0];

  SafePtr<LSSMatrix> jacobMatrix = lss->getMatrix();

  vector<State*> states(2);

  // get the data handle for the states
  DataHandle<State*,GLOBAL> allStates = socket_states.getDataHandle();
  // get the data handle for the is Updated flags
  DataHandle<bool> isUpdated = socket_isUpdated.getDataHandle();
  // get the data handle for the rhs
  DataHandle<CFreal> rhs = socket_rhs.getDataHandle();
  bool isUnsteady = false;

  CFreal dt = SubSystemStatusStack::getActive()->getDT();
  if(dt < 0.) dt = 1.;
  else isUnsteady = true;

// unused // const bool isGhost = true;

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
    const CFuint nb_state = currFace.nbStates();

    // block accumulator nb_state*nb_state
    auto_ptr<BlockAccumulator> acc(lss->createBlockAccumulator(nb_state, nb_state, nbEqs));

    // set the face normal
    setFaceNormal(currFace.getID());


    for ( CFuint iState = 0; iState < nb_state; ++iState){
      m_state[iState] = currFace.getState(iState);
      acc->setRowColIndex(iState, m_state[iState]->getLocalID());
    CFuint stateID = currFace.getState(iState)->getLocalID();
     if (isUnsteady){
      DataHandle<State*> pastStates = socket_pastStates.getDataHandle();

       // Unsteady, state=0.5*(pastState+interState)
       for (CFuint j=0; j<nbEqs ; ++j)
	 (*m_state[iState])[j] = 0.5*((*allStates[stateID])[j] + (*pastStates[stateID])[j]);
     }
       else {
	 // clone the states in order to use the indexes
	 
	 m_state[iState]->clone(*allStates[stateID]);
	 m_gstate[iState]->clone(*allStates[stateID]);	
       }
     }
    for ( CFuint iState = 0; iState < nb_state; ++iState){
      // set ghost state == first state
      // corresponding B-state == second state
      states[0] = m_gstate[iState];
      states[1] = m_state[iState];

      // compute the corresponding fluctuation and the update coeff
      computeFluxAndJacob(states, m_flux[iState], m_fluxJacob[iState]);

      m_fluxJacob[iState] *= -1.;
    }

    for (CFuint iEq = 0; iEq < nbEqs; ++iEq) {
      // distribute corrections to the boundary nodes
      for (CFuint iState = 0 ; iState < nb_state; ++iState) {

        CFuint stateID = currFace.getState(iState)->getLocalID();
        CFreal sum_flux_not_iState = 0.;

        for (CFuint jState = 0; jState < nb_state ; ++jState){
          if (jState != iState)
            sum_flux_not_iState += m_flux[jState][iEq];
        }
        if (!isUpdated[stateID]) {
          rhs(stateID, iEq, nbEqs) +=
                     (m_alpha*m_flux[iState][iEq] + oEminusAlpha*sum_flux_not_iState);

          m_flagState[stateID] = true;
        }
      }
    }

if (getMethodData().doComputeJacobian()) {

    // compute the jacobian contributions for the involved states
    for (CFuint iEq = 0; iEq < nbEqs; ++iEq) {
      for (CFuint jEq = 0; jEq < nbEqs; ++jEq) {
        // distribute corrections to the boundary nodes
        for (CFuint iState = 0; iState < nb_state; ++iState){
          CFuint stateID = currFace.getState(iState)->getLocalID();
          if (!isUpdated[stateID] && m_state[iState]->isParUpdatable()) {
            acc->setValue(iState, iState, iEq, jEq,m_alpha*m_fluxJacob[iState](iEq, jEq));

              for (CFuint jState = 0; jState < nb_state; ++jState){
                CFuint stateID_j = currFace.getState(jState)->getLocalID();
                if ((iState != jState) && !isUpdated[stateID_j])
                  acc->setValue(iState, jState, iEq, jEq, oEminusAlpha*m_fluxJacob[jState](iEq, jEq));
              }
            }
          }
        }
      }

    // add the contributions to the jacobian
    jacobMatrix->addValues(*acc);

    acc->reset();
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

}

//////////////////////////////////////////////////////////////////////////////

} // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
