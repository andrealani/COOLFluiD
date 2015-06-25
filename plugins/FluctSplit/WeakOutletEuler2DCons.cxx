#include "FluctSplit/FluctSplitNavierStokes.hh"
#include "WeakOutletEuler2DCons.hh"
#include "InwardNormalsData.hh"
#include "Framework/PhysicalModel.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/MeshData.hh"
#include "Framework/NamespaceSwitcher.hh"
//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Physics::NavierStokes;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<WeakOutletEuler2DCons, FluctuationSplitData, FluctSplitNavierStokesModule> weakOutletEuler2DConsProvider("WeakOutletEuler2DCons");

//////////////////////////////////////////////////////////////////////////////

void WeakOutletEuler2DCons::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOption< CFreal >("T","total temperature");
}

//////////////////////////////////////////////////////////////////////////////

WeakOutletEuler2DCons::WeakOutletEuler2DCons(const std::string& name) :
  WeakBC2D(name),
  socket_isUpdated("isUpdated"),
  m_varSet()
{
   addConfigOptionsTo(this);
  m_totTemperature = 0.0;
   setParameter("T",&m_totTemperature);
}

//////////////////////////////////////////////////////////////////////////////

WeakOutletEuler2DCons::~WeakOutletEuler2DCons()
{
}

//////////////////////////////////////////////////////////////////////////////

void WeakOutletEuler2DCons::setup()
{
  WeakBC2D::setup();

  m_varSet->setup();

  m_totTemperature /= m_varSet->getModel()->getTempRef();
}

//////////////////////////////////////////////////////////////////////////////

void WeakOutletEuler2DCons::executeOnTrs()
{
  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  const CFreal oEminusAlpha = 1. - m_alpha;
  // states forming the ghost cell
  vector<State*> states(2);
  RealVector flux0(nbEqs);
  RealVector flux1(nbEqs);

  // get the data handle for the states
  DataHandle<State*,GLOBAL> allStates = socket_states.getDataHandle();

  DataHandle<bool> isUpdated = socket_isUpdated.getDataHandle();
  DataHandle<CFreal> rhs = socket_rhs.getDataHandle();

  // create two ghost states (no local ID is needed for these states)
  State* gstate0 = new State();
  State* gstate1 = new State();

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

    State *const state0 = currFace.getState(0);
    State *const state1 = currFace.getState(1);
    const CFuint stateID0 = state0->getLocalID();
    const CFuint stateID1 = state1->getLocalID();

    const CFreal gamma = m_varSet->getModel()->getGamma();
    const CFreal gammaMinus1 = gamma - 1.;

    // compute the third eigenvalue for the first face state
    CFreal invRho = 1./(*state0)[0];
    CFreal u = (*state0)[1]*invRho;
    CFreal v = (*state0)[2]*invRho;
    CFreal Vn = m_faceNormal[0]*u + m_faceNormal[1]*v;
    CFreal p = gammaMinus1*((*state0)[3] - 0.5*(*state0)[0]*(u*u+v*v));
    CFreal a = sqrt(gamma*p*invRho);
    const CFreal lamda30 = Vn + a;

    // if the third eigenvalue is positive apply subsonic outlet BC
    // to calcutate the ghost state and the corresponding flux
    if (lamda30 > 0.0) {
      // set the appropriate values in the ghost states
      setGhostState(*state0, *gstate0);
      // set ghost state == first state
      // corresponding B-state == second state
      states[0] = gstate0;
      states[1] = state0;
      // compute the corresponding fluctuation and the update coeff
      computeFlux(states, flux0);
    }

    // compute the third eigenvalue for the second face state
    invRho = 1./(*state1)[0];
    u = (*state1)[1]*invRho;
    v = (*state1)[2]*invRho;
    Vn = m_faceNormal[0]*u + m_faceNormal[1]*v;
    p = gammaMinus1*((*state1)[3] - 0.5*(*state1)[0]*(u*u+v*v));
    a = sqrt(gamma*p*invRho);
    const CFreal lamda31 = Vn + a;

    // if the third eigenvalue is positive apply subsonic outlet BC
    // to calcutate the ghost state and the corresponding flux
    if (lamda31 > 0.0) {
      setGhostState(*state1, *gstate1);
      // set ghost state == first state
      // corresponding B-state == second state
      states[0] = gstate1;
      states[1] = state1;
      // compute the corresponding fluctuation and the update coeff
      computeFlux(states, flux1);
    }

    for (CFuint iEq = 0; iEq < nbEqs; ++iEq) {

      // distribute corrections to the boundary nodes
      // but only if they haven't been updated yet
      if (!isUpdated[stateID0] && lamda30 > 0.0) {
        rhs(stateID0, iEq, nbEqs) += m_alpha*flux0[iEq];

        if (lamda31 > 0.0) {
          rhs(stateID0, iEq, nbEqs) += oEminusAlpha*flux1[iEq];
        }
        m_flagState[stateID0] = true;
      }

      if (!isUpdated[stateID1] && lamda31 > 0.0) {
        rhs(stateID1, iEq, nbEqs) += m_alpha*flux1[iEq];

        if (lamda30 > 0.0) {
          rhs(stateID1, iEq, nbEqs) += oEminusAlpha*flux0[iEq];
        }

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

  deletePtr(gstate0);
  deletePtr(gstate1);
 }

//////////////////////////////////////////////////////////////////////////////

void WeakOutletEuler2DCons::setGhostState(const State& state,
                                          State& gstate)
{
  // Mach, static pressure, alpha are extrapolated from inside
 const CFreal gamma = m_varSet->getModel()->getGamma();
 const CFreal gammaMinus1 = gamma - 1.;
 const CFreal R = m_varSet->getModel()->getR();
 const CFreal u = state[1]/state[0];
 const CFreal v = state[2]/state[0];
 const CFreal vel2 = u*u + v*v;
  // pressure is extrapolated from inside the domain
  const CFreal p = gammaMinus1*(state[3] - 0.5*state[0]*vel2);
  // mach number is extrapolated from inside the domain
  const CFreal mach = sqrt(vel2/(gamma*p/state[0]));
  const CFreal coeffM = 1. + 0.5*gammaMinus1*mach*mach;
  const CFreal ghostT = m_totTemperature/coeffM;
  const CFreal ghostRho = p/(ghostT*R);
  // alpha is extrapolated from inside the domain
  const CFreal tgAngle = tan(v/u);
  const CFreal coeffAngle = 1. + tgAngle*tgAngle;
  const CFreal ghostU = mach*sqrt(gamma*R*ghostT/coeffAngle);
  const CFreal ghostV = ghostU*tgAngle;

  gstate[0] = ghostRho;
  gstate[1] = ghostRho*ghostU;
  gstate[2] = ghostRho*ghostV;
  gstate[3] = p/gammaMinus1 + 0.5*(gstate[1]*gstate[1] +
                                   gstate[2]*gstate[2])/gstate[0];
 }

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
WeakOutletEuler2DCons::needsSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSink> > result = WeakBC::needsSockets();

  result.push_back(&socket_isUpdated);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

void WeakOutletEuler2DCons::configure ( Config::ConfigArgs& args )
{
  WeakBC2D::configure(args);

  std::string name = getMethodData().getNamespace();
  Common::SafePtr<Namespace> nsp = NamespaceSwitcher::getInstance
    (SubSystemStatusStack::getCurrentName()).getNamespace(name);
  Common::SafePtr<PhysicalModel> physModel = PhysicalModelStack::getInstance().getEntryByNamespace(nsp);
  
  std::string varSetName = "Euler2DCons";
  m_varSet.reset((Environment::Factory<ConvectiveVarSet>::getInstance().getProvider(varSetName)->
    create(physModel->getImplementor()->getConvectiveTerm())).d_castTo<Physics::NavierStokes::Euler2DCons>());

  cf_assert(m_varSet.isNotNull());

}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
