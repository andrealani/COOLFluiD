#include "FluctSplit/FluctSplitSpaceTimeNavierStokes.hh"
#include "UnsteadyWeakSlipWallEulerHO2DCons.hh"
#include "InwardNormalsData.hh"
#include "Framework/CFL.hh"
#include "Framework/SubSystemStatus.hh"
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

MethodCommandProvider<UnsteadyWeakSlipWallEulerHO2DCons, FluctuationSplitData, FluctSplitSpaceTimeNavierStokesModule> UnsteadyWeakSlipWallEulerHO2DConsProvider("UnsteadyWeakSlipWallEulerHO2DCons");

//////////////////////////////////////////////////////////////////////////////

void UnsteadyWeakSlipWallEulerHO2DCons::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOption< CFreal >("alpha","distribution coefficient");
}

//////////////////////////////////////////////////////////////////////////////

UnsteadyWeakSlipWallEulerHO2DCons::UnsteadyWeakSlipWallEulerHO2DCons(const std::string& name) :
  FluctuationSplitCom(name),
  _varSet(),
  socket_rhs("rhs"),
  socket_normals("normals"),
  socket_nodes("nodes"),
  socket_faceNeighCell("faceNeighCell"),
  socket_states("states"),
  socket_pastStates("pastStates"),
  socket_interStates("interStates"),
  socket_isUpdated("isUpdated"),
  _dt(),
  _wallSpeed(),
  _flagState(),
  m_fluxn1(),
  m_fluxn12(),
  m_fluxn(),
  m_flux()
{
   addConfigOptionsTo(this);
  _alpha = 0.66;
   setParameter("alpha",&_alpha);

}

//////////////////////////////////////////////////////////////////////////////

UnsteadyWeakSlipWallEulerHO2DCons::~UnsteadyWeakSlipWallEulerHO2DCons()
{
}

//////////////////////////////////////////////////////////////////////////////

void UnsteadyWeakSlipWallEulerHO2DCons::setup()
{
  _wallSpeed.resize(PhysicalModelStack::getActive()->getDim());
  _varSet->setup();
  // get the data handle for the states
  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();

  _flagState.resize(states.size());
  _flagState = false;

  m_fluxn1.resize(3); //There are 3 nodes per face because we are considering a P2 element
  m_fluxn12.resize(3); //There are 3 nodes per face because we are considering a P2 element
  m_fluxn.resize(3); //There are 3 nodes per face because we are considering a P2 element
  m_flux.resize(3); //There are 3 nodes per face because we are considering a P2 element

  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  for (CFuint iState = 0; iState<3; ++iState ){
    m_fluxn1[iState].resize(nbEqs);
    m_fluxn12[iState].resize(nbEqs);
    m_fluxn[iState].resize(nbEqs);
    m_flux[iState].resize(nbEqs);
  }

}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
UnsteadyWeakSlipWallEulerHO2DCons::needsSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSink> > result;

  result.push_back(&socket_rhs);
  result.push_back(&socket_normals);
  result.push_back(&socket_nodes);
  result.push_back(&socket_faceNeighCell);
  result.push_back(&socket_states);
  result.push_back(&socket_pastStates);
  result.push_back(&socket_interStates);
  result.push_back(&socket_isUpdated);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

void UnsteadyWeakSlipWallEulerHO2DCons::configure ( Config::ConfigArgs& args )
{
  FluctuationSplitCom::configure(args);

  const std::string name = getMethodData().getNamespace();
  Common::SafePtr<Namespace> nsp = NamespaceSwitcher::getInstance
    (SubSystemStatusStack::getCurrentName()).getNamespace(name);
  Common::SafePtr<PhysicalModel> physModel = PhysicalModelStack::getInstance().getEntryByNamespace(nsp);
  
  const std::string varSetName = "Euler2DCons";
  _varSet.reset((Environment::Factory<ConvectiveVarSet>::getInstance().getProvider(varSetName)->
		 create(physModel->getImplementor()->getConvectiveTerm())).d_castTo<Physics::NavierStokes::Euler2DCons>());
  
  cf_assert(_varSet.isNotNull());

}

//////////////////////////////////////////////////////////////////////////////

void UnsteadyWeakSlipWallEulerHO2DCons::executeOnTrs()
{
// The first iteration needs a special function because at this point
  // only two layer are available althought we need three.

  CFuint nbIter = SubSystemStatusStack::getActive()->getNbIter();
  if(nbIter == 1) docomputeFirstBC();
  else docomputeBC();
}

//////////////////////////////////////////////////////////////////////////
void UnsteadyWeakSlipWallEulerHO2DCons::docomputeFirstBC()
{
  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  const CFuint dim = PhysicalModelStack::getActive()->getDim();
  const CFreal oEminusAlpha = 1. - _alpha;
  RealVector normal(dim);

  _dt = SubSystemStatusStack::getActive()->getDT();

  // get the data handle for the rhs
  DataHandle< CFreal> rhs = socket_rhs.getDataHandle();

  // get the data handle for the states
  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();

  // get the data handle for the past states
  DataHandle<State*> interStates = socket_interStates.getDataHandle();

  // get the data handle the update of the states
  DataHandle<bool> isUpdated = socket_isUpdated.getDataHandle();

  Common::SafePtr<GeometricEntityPool<StdTrsGeoBuilder> >
    geoBuilder = getMethodData().getStdTrsGeoBuilder();

  StdTrsGeoBuilder::GeoData& geoData = geoBuilder->getDataGE();
  geoData.trs = getCurrentTRS();

  const CFuint nbFaces = getCurrentTRS()->getLocalNbGeoEnts();
  for (CFuint iFace = 0; iFace < nbFaces; ++iFace) {
    geoData.idx = iFace;
    GeometricEntity *const currFace = geoBuilder->buildGE();

    // set the face normal
    setFaceNormal(currFace->getID(), normal);

    const CFuint stateID0 = currFace->getState(0)->getLocalID();
    const CFuint stateID1 = currFace->getState(1)->getLocalID();
    const CFuint stateID2 = currFace->getState(2)->getLocalID();

    // compute the normal fluxes corrections for both the states
    // of this cell
    computeNormalFlux((*states[stateID0]), normal, m_fluxn1[0]);
    computeNormalFlux((*states[stateID1]), normal, m_fluxn1[1]);
    computeNormalFlux((*states[stateID2]), normal, m_fluxn1[2]);
    computeNormalFlux((*interStates[stateID0]), normal, m_fluxn[0]);
    computeNormalFlux((*interStates[stateID1]), normal, m_fluxn[1]);
    computeNormalFlux((*interStates[stateID2]), normal, m_fluxn[2]);

    m_flux[0] = _dt*0.5*(m_fluxn1[0]+m_fluxn[0]);
    m_flux[1] = _dt*0.5*(m_fluxn1[1]+m_fluxn[1]);
    m_flux[2] = _dt*0.5*(m_fluxn1[2]+m_fluxn[2]);

    // distribute contributions to the two nodes
    for (CFuint iEq = 0; iEq < nbEqs; ++iEq) {
      if (!isUpdated[stateID0]){
        rhs(stateID0, iEq, nbEqs) -= _alpha*(1.0/6.0)*m_flux[0][iEq] + (oEminusAlpha/2.0)*(1.0/6.0)*m_flux[1][iEq] + (oEminusAlpha/2.0)*(2.0/3.0)*m_flux[2][iEq];

        _flagState[stateID0] = true;
      }

      if (!isUpdated[stateID1]){
        rhs(stateID1, iEq, nbEqs) -= _alpha*(1.0/6.0)*m_flux[1][iEq] + (oEminusAlpha/2.0)*(1.0/6.0)*m_flux[0][iEq] + (oEminusAlpha/2.0)*(2.0/3.0)*m_flux[2][iEq];

        _flagState[stateID1] = true;
      }

      if (!isUpdated[stateID2]){
        rhs(stateID2, iEq, nbEqs) -= _alpha*(2.0/3.0)*m_flux[2][iEq] + (oEminusAlpha/2.0)*(1.0/6.0)*m_flux[1][iEq] + (oEminusAlpha/2.0)*(1.0/6.0)*m_flux[0][iEq];

        _flagState[stateID2] = true;
      }
}

    // release the face
    geoBuilder->releaseGE();
  }
}

//////////////////////////////////////////////////////////////////////////
void UnsteadyWeakSlipWallEulerHO2DCons::docomputeBC()
{
  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  const CFuint dim = PhysicalModelStack::getActive()->getDim();
  const CFreal oEminusAlpha = 1. - _alpha;
  RealVector normal(dim);


  const CFreal dt1 = SubSystemStatusStack::getActive()->getDT();
  const CFreal dt0 =SubSystemStatusStack::getActive()->getPreviousDT();
  const CFreal q = dt1/dt0;
  const CFreal Jpast  = -dt1*(q*q)/(6.0*(1.0 + q));
  const CFreal Jinter  = dt1*(3.0 + q)/6.0;
  const CFreal Jpresent  = dt1*(3.0 + 2.0*q)/(6.0*(1.0 + q));

  // get the data handle for the rhs
  DataHandle< CFreal> rhs = socket_rhs.getDataHandle();

  // get the data handle for the states
  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();

  // get the data handle for the past states
  DataHandle<State*> pastStates = socket_pastStates.getDataHandle();

  // get the data handle for the past states
  DataHandle<State*> interStates = socket_interStates.getDataHandle();

  // get the data handle the update of the states
  DataHandle<bool> isUpdated = socket_isUpdated.getDataHandle();

  Common::SafePtr<GeometricEntityPool<StdTrsGeoBuilder> >
    geoBuilder = getMethodData().getStdTrsGeoBuilder();

  StdTrsGeoBuilder::GeoData& geoData = geoBuilder->getDataGE();
  geoData.trs = getCurrentTRS();

  const CFuint nbFaces = getCurrentTRS()->getLocalNbGeoEnts();
  for (CFuint iFace = 0; iFace < nbFaces; ++iFace) {
    geoData.idx = iFace;
    GeometricEntity *const currFace = geoBuilder->buildGE();

    // set the face normal
    setFaceNormal(currFace->getID(), normal);

    const CFuint stateID0 = currFace->getState(0)->getLocalID();
    const CFuint stateID1 = currFace->getState(1)->getLocalID();
    const CFuint stateID2 = currFace->getState(2)->getLocalID();

    // compute the normal fluxes corrections for both the states
    // of this cell
    computeNormalFlux((*states[stateID0]), normal, m_fluxn1[0]);
    computeNormalFlux((*states[stateID1]), normal, m_fluxn1[1]);
    computeNormalFlux((*states[stateID2]), normal, m_fluxn1[2]);
    computeNormalFlux((*interStates[stateID0]), normal, m_fluxn12[0]);
    computeNormalFlux((*interStates[stateID1]), normal, m_fluxn12[1]);
    computeNormalFlux((*interStates[stateID2]), normal, m_fluxn12[2]);
    computeNormalFlux((*pastStates[stateID0]), normal, m_fluxn[0]);
    computeNormalFlux((*pastStates[stateID1]), normal, m_fluxn[1]);
    computeNormalFlux((*pastStates[stateID2]), normal, m_fluxn[2]);

    m_flux[0] = Jpast*m_fluxn[0] + Jinter*m_fluxn12[0] + Jpresent*m_fluxn1[0];
    m_flux[1] = Jpast*m_fluxn[1] + Jinter*m_fluxn12[1] + Jpresent*m_fluxn1[1];
    m_flux[2] = Jpast*m_fluxn[2] + Jinter*m_fluxn12[2] + Jpresent*m_fluxn1[2];

    // distribute contributions to the two nodes
    for (CFuint iEq = 0; iEq < nbEqs; ++iEq) {
      if (!isUpdated[stateID0]){
        rhs(stateID0, iEq, nbEqs) -= _alpha*(1.0/6.0)*m_flux[0][iEq] + (oEminusAlpha/2.0)*(1.0/6.0)*m_flux[1][iEq] + (oEminusAlpha/2.0)*(2.0/3.0)*m_flux[2][iEq];

        _flagState[stateID0] = true;
      }

      if (!isUpdated[stateID1]){
        rhs(stateID1, iEq, nbEqs) -= _alpha*(1.0/6.0)*m_flux[1][iEq] + (oEminusAlpha/2.0)*(1.0/6.0)*m_flux[0][iEq] + (oEminusAlpha/2.0)*(2.0/3.0)*m_flux[2][iEq];

        _flagState[stateID1] = true;
      }

      if (!isUpdated[stateID2]){
        rhs(stateID2, iEq, nbEqs) -= _alpha*(2.0/3.0)*m_flux[2][iEq] + (oEminusAlpha/2.0)*(1.0/6.0)*m_flux[1][iEq] + (oEminusAlpha/2.0)*(1.0/6.0)*m_flux[0][iEq];

        _flagState[stateID2] = true;
      }
}

    // release the face
    geoBuilder->releaseGE();
  }
}


//////////////////////////////////////////////////////////////////////////////

void UnsteadyWeakSlipWallEulerHO2DCons::computeNormalFlux(
               const RealVector& state,
               const RealVector& normal,
               RealVector& flux) const
 {
   const CFreal nx = normal[0];
   const CFreal ny = normal[1];

   const CFreal rho = state[0];
   const CFreal rhoU = state[1];
   const CFreal rhoV = state[2];
   const CFreal rhoE = state[3];
   const CFreal u = rhoU/rho;
   const CFreal v = rhoV/rho;
   // unused // const CFreal E = rhoE/rho;
   const CFreal un = u*nx + v*ny;
   const CFreal gamma = _varSet->getModel()->getGamma();
   const CFreal halfGammaMinus1 = 0.5*(gamma - 1.);
   // unused //   const CFreal gammaMinus1 = _varSet.getModel()->getGammaMinus1();
   const CFreal vel2 = u*u + v*v;
   const CFreal rhoV2 = rho*vel2;
   const CFreal rhoH = gamma*rhoE - halfGammaMinus1*rhoV2;

   // Compute Flux for non-moving wall
   flux[0] = rho*un;
   flux[1] = un*rhoU;
   flux[2] = un*rhoV;
   flux[3] = un*rhoH;


 }

//////////////////////////////////////////////////////////////////////////////

void UnsteadyWeakSlipWallEulerHO2DCons::setFaceNormal(const CFuint faceID,
					RealVector& normal)
{

  DataHandle<Common::Trio<CFuint, CFuint, Common::SafePtr<Framework::TopologicalRegionSet> > >
    faceNeighCell = socket_faceNeighCell.getDataHandle();
  DataHandle< InwardNormalsData*> normals = socket_normals.getDataHandle();

  const CFuint cellTrsID = faceNeighCell[faceID].first;
  const CFuint iFaceLocal = faceNeighCell[faceID].second;
  Common::SafePtr<TopologicalRegionSet> cellTrs = faceNeighCell[faceID].third;
  const CFuint cellLocalID = cellTrs->getLocalGeoID(cellTrsID);

  // If the wall is moving, the normal is the average of
  // the normals at the two instants

    normal[0] = - normals[cellLocalID]->getFaceNormComp(iFaceLocal, 0);
    normal[1] = - normals[cellLocalID]->getFaceNormComp(iFaceLocal, 1);

}

//////////////////////////////////////////////////////////////////////////////

  } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
