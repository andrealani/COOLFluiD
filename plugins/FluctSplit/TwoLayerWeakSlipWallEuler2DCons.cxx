#include "FluctSplit/FluctSplitSpaceTimeNavierStokes.hh"
#include "TwoLayerWeakSlipWallEuler2DCons.hh"
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

MethodCommandProvider<TwoLayerWeakSlipWallEuler2DCons, FluctuationSplitData, FluctSplitSpaceTimeNavierStokesModule> twoLayerWeakSlipWallEuler2DConsProvider("TwoLayerWeakSlipWallEuler2DCons");

//////////////////////////////////////////////////////////////////////////////

void TwoLayerWeakSlipWallEuler2DCons::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOption< CFreal >("alpha","distribution coefficient");
}

//////////////////////////////////////////////////////////////////////////////

TwoLayerWeakSlipWallEuler2DCons::TwoLayerWeakSlipWallEuler2DCons(const std::string& name) :
  FluctuationSplitCom(name),
  _varSet(),
  socket_rhs("rhs"),
  socket_interRhs("interRhs"),
  socket_normals("normals"),
  socket_interNormals("interNormals", false), // only used if moving
  socket_pastNormals("pastNormals", false), // only used if moving
  socket_nodes("nodes",false), // only used if moving
  socket_pastNodes("pastNodes", false),     // only used if moving
  socket_faceNeighCell("faceNeighCell"),
  socket_states("states"),
  socket_interStates("interStates"),
  socket_pastStates("pastStates"),
  socket_isUpdated("isUpdated"),
  _dt(),
  _layer(0),
  _wallSpeed(),
  _flagState()
{
   addConfigOptionsTo(this);
  _alpha = 1.0;
   setParameter("alpha",&_alpha);

}

//////////////////////////////////////////////////////////////////////////////

TwoLayerWeakSlipWallEuler2DCons::~TwoLayerWeakSlipWallEuler2DCons()
{
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
TwoLayerWeakSlipWallEuler2DCons::needsSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSink> > result;

  result.push_back(&socket_rhs);
  result.push_back(&socket_interRhs);
  result.push_back(&socket_normals);
  result.push_back(&socket_interNormals);
  result.push_back(&socket_pastNormals);
  result.push_back(&socket_nodes);
  result.push_back(&socket_pastNodes);
  result.push_back(&socket_faceNeighCell);
  result.push_back(&socket_states);
  result.push_back(&socket_interStates);
  result.push_back(&socket_pastStates);
  result.push_back(&socket_isUpdated);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

void TwoLayerWeakSlipWallEuler2DCons::setup()
{

  // get the data handle for the states
  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();

  _flagState.resize(states.size());
  _flagState = false;

  _varSet->setup();
  _wallSpeed.resize(PhysicalModelStack::getActive()->getDim());
 }

//////////////////////////////////////////////////////////////////////////////

void TwoLayerWeakSlipWallEuler2DCons::configure ( Config::ConfigArgs& args )
{
  FluctuationSplitCom::configure(args);

  std::string name = getMethodData().getNamespace();
  Common::SafePtr<Namespace> nsp = NamespaceSwitcher::getInstance().getNamespace(name);
  Common::SafePtr<PhysicalModel> physModel = PhysicalModelStack::getInstance().getEntryByNamespace(nsp);

  std::string varSetName = "Euler2DCons";
  _varSet.reset((Environment::Factory<ConvectiveVarSet>::getInstance().getProvider(varSetName)->
    create(physModel->getImplementor()->getConvectiveTerm())).d_castTo<Physics::NavierStokes::Euler2DCons>());

  cf_assert(_varSet.isNotNull());

}


//////////////////////////////////////////////////////////////////////////////

void TwoLayerWeakSlipWallEuler2DCons::executeOnTrs()
{
  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  const CFuint dim = PhysicalModelStack::getActive()->getDim();
  const CFreal oEminusAlpha = 1. - _alpha;

  RealVector normal(dim);
  RealVector flux0(0.0, nbEqs);
  RealVector flux1(0.0, nbEqs);

  // get the data handle for the rhs
  DataHandle< CFreal> rhs = socket_rhs.getDataHandle();

  // get the data handle for the intermediate rhs
  DataHandle< CFreal> interRhs = socket_interRhs.getDataHandle();

  // get the data handle for the past states
  DataHandle<State*> pastStates = socket_pastStates.getDataHandle();

  // get the data handle for the past states
  DataHandle<State*> interStates = socket_interStates.getDataHandle();

  // get the data handle for the states
  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();

  // get the data handle the update of the states
  DataHandle<bool> isUpdated = socket_isUpdated.getDataHandle();

 Common::SafePtr<GeometricEntityPool<StdTrsGeoBuilder> >
    geoBuilder = getMethodData().getStdTrsGeoBuilder();

  StdTrsGeoBuilder::GeoData& geoData = geoBuilder->getDataGE();
  geoData.trs = getCurrentTRS();

  State state0;
  State state1;

  const CFuint nbFaces = getCurrentTRS()->getLocalNbGeoEnts();
  for (CFuint iFace = 0; iFace < nbFaces; ++iFace) {
    geoData.idx = iFace;
    GeometricEntity *const currFace = geoBuilder->buildGE();
    const CFuint stateID0 = currFace->getState(0)->getLocalID();
    const CFuint stateID1 = currFace->getState(1)->getLocalID();

    // if Moving Mesh, compute the wall velocity
    for (CFuint iDim=0; iDim < dim; ++iDim){
      _wallSpeed[iDim] = 0.;
    }
    if (SubSystemStatusStack::getActive()->isMovingMesh()){
      DataHandle < Framework::Node*, Framework::GLOBAL > nodes = socket_nodes.getDataHandle();
      DataHandle<Node*> pastNodes = socket_pastNodes.getDataHandle();

      _dt = SubSystemStatusStack::getActive()->getDT();
      for (CFuint iDim=0; iDim < dim; ++iDim){
	_wallSpeed[iDim]  = (*nodes[stateID0])[iDim]-(*pastNodes[stateID0])[iDim];
	_wallSpeed[iDim] += (*nodes[stateID1])[iDim]-(*pastNodes[stateID1])[iDim];
	_wallSpeed[iDim] *= 0.5;
	_wallSpeed[iDim] *= _dt;
      }
    }

    /// Cell K1
      _layer = 0;
      // Set dt1
      _dt = SubSystemStatusStack::getActive()->getInnerDT(_layer);
      // set the face normal
      setFaceNormal(currFace->getID(), normal);

      // Unsteady, state=0.5*(pastState+interState)
      for (CFuint j=0; j<nbEqs ; ++j){
	state0[j] = 0.5*((*interStates[stateID0])[j] + (*pastStates[stateID0])[j]);
	state1[j] = 0.5*((*interStates[stateID1])[j] + (*pastStates[stateID1])[j]);
      }

      // compute the normal fluxes corrections for both states
      // of this cell
      computeNormalFlux(state0, normal, flux0);
      computeNormalFlux(state1, normal, flux1);

      // distribute contributions to the two nodes
      for (CFuint iEq = 0; iEq < nbEqs; ++iEq) {
        if (!isUpdated[stateID0]) {
          interRhs(stateID0, iEq, nbEqs) -=
            0.5*(_alpha*flux0[iEq] + oEminusAlpha*flux1[iEq]);

          _flagState[stateID0] = true;
        }

        if (!isUpdated[stateID1]) {
          interRhs(stateID1, iEq, nbEqs) -=
            0.5*(_alpha*flux1[iEq] + oEminusAlpha*flux0[iEq]);

          _flagState[stateID1] = true;
        }

      }

      /// Cell K2
      _layer = 1;
      // Set dt2
      _dt = SubSystemStatusStack::getActive()->getInnerDT(_layer);

      // set the face normal
      setFaceNormal(currFace->getID(), normal);

      // Unsteady, state=0.5*(pastState+State)
      for (CFuint j=0; j<nbEqs ; ++j){
          state0[j] = 0.5*((*interStates[stateID0])[j] + (*states[stateID0])[j]);
          state1[j] = 0.5*((*interStates[stateID1])[j] + (*states[stateID1])[j]);
        }

      // compute the normal fluxes corrections for both states
      // of this cell
      computeNormalFlux(state0, normal, flux0);
      computeNormalFlux(state1, normal, flux1);

      for (CFuint iEq = 0; iEq < nbEqs; ++iEq) {
        if (!isUpdated[stateID0]) {
          rhs(stateID0, iEq, nbEqs) -=
            0.5*(_alpha*flux0[iEq] + oEminusAlpha*flux1[iEq]);

          _flagState[stateID0] = true;
        }

        if (!isUpdated[stateID1]) {
          rhs(stateID1, iEq, nbEqs) -=
            0.5*(_alpha*flux1[iEq] + oEminusAlpha*flux0[iEq]);

          _flagState[stateID1] = true;
        }
      }

    // release the face
    geoBuilder->releaseGE();
  }
}

//////////////////////////////////////////////////////////////////////////////

void TwoLayerWeakSlipWallEuler2DCons::computeNormalFlux(
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
   const CFreal un = u*nx + v*ny;
   const CFreal gamma = _varSet->getModel()->getGamma();
   const CFreal halfGammaMinus1 = 0.5*(gamma - 1.);
   const CFreal vel2 = u*u + v*v;
   const CFreal rhoV2 = rho*vel2;
   const CFreal rhoH = gamma*rhoE - halfGammaMinus1*rhoV2;
   // unused // const CFreal uuvv =  0.5*(u*u + v*v);
   // unused // const CFreal p = gammaMinus1*(rhoE - rho*uuvv);
   const CFreal velocity = nx*_wallSpeed[0] + ny*_wallSpeed[1];

   // Compute Flux for non-moving wall
   flux[0] = rho*un;
   flux[1] = un*rhoU;
   flux[2] = un*rhoV;
   flux[3] = un*rhoH;

   if (SubSystemStatusStack::getActive()->isMovingMesh()){
    // Modification of the flux due to the movement of the wall
    flux[0] += rho*velocity;
    flux[1] += rhoU*velocity;
    flux[2] += rhoV*velocity;
    flux[3] += rhoH*velocity;
    }

   // Scale with dt
   flux[0] *= _dt;
   flux[1] *= _dt;
   flux[2] *= _dt;
   flux[3] *= _dt;

 }

//////////////////////////////////////////////////////////////////////////////

void TwoLayerWeakSlipWallEuler2DCons::setFaceNormal
(const CFuint faceID, RealVector& normal)
{
  DataHandle< InwardNormalsData*> normals = socket_normals.getDataHandle();
  DataHandle<Common::Trio<CFuint, CFuint, Common::SafePtr<Framework::TopologicalRegionSet> > >
    faceNeighCell = socket_faceNeighCell.getDataHandle();

  const CFuint cellTrsID = faceNeighCell[faceID].first;
  const CFuint iFaceLocal = faceNeighCell[faceID].second;
  Common::SafePtr<TopologicalRegionSet> cellTrs = faceNeighCell[faceID].third;
  const CFuint cellLocalID = cellTrs->getLocalGeoID(cellTrsID);

  // If the wall is moving, the normal is the average of
  // the normals at the two instants
  if (SubSystemStatusStack::getActive()->isMovingMesh()){
    if (_layer == 0){
      DataHandle< InwardNormalsData*> pastNormals = socket_pastNormals.getDataHandle();
      DataHandle< InwardNormalsData*> interNormals = socket_interNormals.getDataHandle();
      normal[0] = -0.5*(interNormals[cellLocalID]->getFaceNormComp(iFaceLocal, 0)+pastNormals[cellLocalID]->getFaceNormComp(iFaceLocal, 0));
      normal[1] = -0.5*(interNormals[cellLocalID]->getFaceNormComp(iFaceLocal, 1)+pastNormals[cellLocalID]->getFaceNormComp(iFaceLocal, 1));
    }
    if (_layer == 1){
      DataHandle< InwardNormalsData*> interNormals = socket_interNormals.getDataHandle();
      normal[0] = -0.5*(normals[cellLocalID]->getFaceNormComp(iFaceLocal, 0)+interNormals[cellLocalID]->getFaceNormComp(iFaceLocal, 0));
      normal[1] = -0.5*(normals[cellLocalID]->getFaceNormComp(iFaceLocal, 1)+interNormals[cellLocalID]->getFaceNormComp(iFaceLocal, 1));
    }
  }
  else{
    normal[0] = -normals[cellLocalID]->getFaceNormComp(iFaceLocal, 0);
    normal[1] = -normals[cellLocalID]->getFaceNormComp(iFaceLocal, 1);
  }
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
