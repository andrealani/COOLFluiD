#include "WeakBC3DImpl.hh"
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

WeakBC3DImpl::WeakBC3DImpl(const std::string& name) :
  WeakBC(name),
  socket_isUpdated("isUpdated"),
  socket_pastStates("pastStates", false),
  m_fluxes(),
  m_fluxJacobs()
{
}

//////////////////////////////////////////////////////////////////////////////

WeakBC3DImpl::~WeakBC3DImpl()
{
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
WeakBC3DImpl::needsSockets()
{

  std::vector<Common::SafePtr<BaseDataSocketSink> > result = WeakBC::needsSockets();

  result.push_back(&socket_isUpdated);
  result.push_back(&socket_pastStates);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

void WeakBC3DImpl::setup()
{
  WeakBC::setup();

  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  // size the fluxes
  // only triangular faces for now
  m_fluxes.resize(3);
  for (CFuint i = 0; i < m_fluxes.size(); ++i) {
    m_fluxes[i].resize(nbEqs);
  }

  // size the flux jacobian3s
  m_fluxJacobs.resize(3);
  for (CFuint i = 0; i < m_fluxJacobs.size(); ++i) {
    m_fluxJacobs[i].resize(nbEqs, nbEqs);
  }

  state0new = new State();
  state1new = new State();
  state2new = new State();
}

//////////////////////////////////////////////////////////////////////////////

void WeakBC3DImpl::executeOnTrs()
{
  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  const CFreal third = 1./3.;
  const CFreal halfOEminAlphaThird = 0.5*(1. - m_alpha)*third;
  const CFreal alphaThird = m_alpha*third;
  bool isUnsteady = false;	
  CFreal dt = SubSystemStatusStack::getActive()->getDT();
  if(dt < 0.) dt = 1.;
  else isUnsteady = true;

  

  SafePtr<LinearSystemSolver> lss =
    getMethodData().getLinearSystemSolver()[0];

  SafePtr<LSSMatrix> jacobMatrix = lss->getMatrix();

  // block accumulator 3*3
  /// @todo it is assumed that only triangular faces are used
  auto_ptr<BlockAccumulator> acc(lss->createBlockAccumulator(3, 3, nbEqs));

  // get the data handle for the states
  DataHandle<State*,GLOBAL> allStates = socket_states.getDataHandle();

  // get the data handle for the isUpdated flags
  DataHandle<bool> isUpdated = socket_isUpdated.getDataHandle();

  // get the data handle for the isUpdated flags
  DataHandle<CFreal> rhs = socket_rhs.getDataHandle();

  RealVector temp(0.0, nbEqs);
  const bool isGhost = true;
  // create two ghost states (no local ID is needed for these states)
  State* gstate0 = new State(temp, isGhost);
  State* gstate1 = new State(temp, isGhost);
  State* gstate2 = new State(temp, isGhost);

  vector<State*> states(2);
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

    // set the face normal
    setFaceNormal(currFace.getID());

    ///@todo extend this for quads...
    cf_assert(currFace.nbStates() == 3);

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


    // set the row - column
    acc->setRowColIndex(0, stateID0);
    acc->setRowColIndex(1, stateID1);
    acc->setRowColIndex(2, stateID2);

    // set ghost state == first state
    // corresponding B-state == second state
    states[0] = gstate0;
    states[1] = state0new;
    // compute the corresponding fluctuation and the update coeff
    computeFluxAndJacob(states, m_fluxes[0], m_fluxJacobs[0]);

    // set ghost state == first state
    // corresponding B-state == second state
    states[0] = gstate1;
    states[1] = state1new;
    // compute the corresponding fluctuation and the update coeff
    computeFluxAndJacob(states, m_fluxes[1], m_fluxJacobs[1]);

    // set ghost state == first state
    // corresponding B-state == second state
    states[0] = gstate2;
    states[1] = state2new;
    // compute the corresponding fluctuation and the update coeff
    computeFluxAndJacob(states, m_fluxes[2], m_fluxJacobs[2]);

    m_fluxJacobs[0] *= -1.;
    m_fluxJacobs[1] *= -1.;
    m_fluxJacobs[2] *= -1.;

    // distribute contributions to the nodes
    for (CFuint iEq = 0; iEq < nbEqs; ++iEq) {
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

    // compute the jacobian contributions for the involved states
///@todo Need to check if the node contributing to a row should be
///      checked to be updated (Nadege)
if (getMethodData().doComputeJacobian()) {

    for (CFuint iEq = 0; iEq < nbEqs; ++iEq) {
      for (CFuint jEq = 0; jEq < nbEqs; ++jEq) {

        if (!isUpdated[stateID0] && state0->isParUpdatable()) {
          acc->setValue(0, 0, iEq, jEq,
                        alphaThird*m_fluxJacobs[0](iEq, jEq));

          acc->setValue(0, 1, iEq, jEq,
                        halfOEminAlphaThird*m_fluxJacobs[1](iEq, jEq));

          acc->setValue(0, 2, iEq, jEq,
                        halfOEminAlphaThird*m_fluxJacobs[2](iEq, jEq));
        }


        if (!isUpdated[stateID1] && state1->isParUpdatable()) {
          acc->setValue(1, 0, iEq, jEq,
                        halfOEminAlphaThird*m_fluxJacobs[0](iEq, jEq));

          acc->setValue(1, 1, iEq, jEq,
                        alphaThird*m_fluxJacobs[1](iEq, jEq));

          acc->setValue(1, 2, iEq, jEq,
                        halfOEminAlphaThird*m_fluxJacobs[2](iEq, jEq));
        }


        if (!isUpdated[stateID2] && state2->isParUpdatable()) {
          acc->setValue(2, 0, iEq, jEq,
                        halfOEminAlphaThird*m_fluxJacobs[0](iEq, jEq));

          acc->setValue(2, 1, iEq, jEq,
                        halfOEminAlphaThird*m_fluxJacobs[1](iEq, jEq));

          acc->setValue(2, 2, iEq, jEq,
                        alphaThird*m_fluxJacobs[2](iEq, jEq));
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

  deletePtr(gstate0);
  deletePtr(gstate1);
  deletePtr(gstate2);
}

//////////////////////////////////////////////////////////////////////////////

} // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
