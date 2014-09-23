#include "FluctSplit/FluctSplitSpaceTimeNavierStokes.hh"
#include "UnsteadyWeakSlipWallEuler3DCons.hh"
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

MethodCommandProvider<UnsteadyWeakSlipWallEuler3DCons, FluctuationSplitData, FluctSplitSpaceTimeNavierStokesModule> unsteadyWeakSlipWallEuler3DConsProvider("UnsteadyWeakSlipWallEuler3DCons");

//////////////////////////////////////////////////////////////////////////////

void UnsteadyWeakSlipWallEuler3DCons::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOption< CFreal >("alpha","distribution coefficient");
}

//////////////////////////////////////////////////////////////////////////////

UnsteadyWeakSlipWallEuler3DCons::UnsteadyWeakSlipWallEuler3DCons(const std::string& name) :
  FluctuationSplitCom(name),
  socket_rhs("rhs"),
  socket_normals("normals"),
  socket_pastNormals("pastNormals", false), // only used if moving
  socket_nodes("nodes",false),
  socket_pastNodes("pastNodes", false),     // only used if moving
  socket_faceNeighCell("faceNeighCell"),
  socket_states("states"),
  socket_pastStates("pastStates"),
  socket_isUpdated("isUpdated"),
  _varSet(),
  _fluxes(),
  _dt(),
  _wallSpeed(),
  _flagState()
{
   addConfigOptionsTo(this);
  _alpha = 0.66;
   setParameter("alpha",&_alpha);
}

//////////////////////////////////////////////////////////////////////////////

UnsteadyWeakSlipWallEuler3DCons::~UnsteadyWeakSlipWallEuler3DCons()
{
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
UnsteadyWeakSlipWallEuler3DCons::needsSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSink> > result;

  result.push_back(&socket_rhs);
  result.push_back(&socket_normals);
  result.push_back(&socket_pastNormals);
  result.push_back(&socket_nodes);
  result.push_back(&socket_pastNodes);
  result.push_back(&socket_faceNeighCell);
  result.push_back(&socket_states);
  result.push_back(&socket_pastStates);
  result.push_back(&socket_isUpdated);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

void UnsteadyWeakSlipWallEuler3DCons::setup()
{
  _wallSpeed.resize(PhysicalModelStack::getActive()->getDim());
  _varSet->setup();
  // size the fluxes
  _fluxes.resize(MeshDataStack::getActive()->Statistics().getMaxNbStatesInCell()-1);
  for (CFuint i = 0; i < _fluxes.size(); ++i) {
    _fluxes[i].resize(PhysicalModelStack::getActive()->getNbEq());
  }

  // get the data handle for the states
  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();

  _flagState.resize(states.size());
  _flagState = false;
}

//////////////////////////////////////////////////////////////////////////////

void UnsteadyWeakSlipWallEuler3DCons::configure ( Config::ConfigArgs& args )
{
  FluctuationSplitCom::configure(args);

  std::string name = getMethodData().getNamespace();
  Common::SafePtr<Namespace> nsp = NamespaceSwitcher::getInstance().getNamespace(name);
  Common::SafePtr<PhysicalModel> physModel = PhysicalModelStack::getInstance().getEntryByNamespace(nsp);

  std::string varSetName = "Euler3DCons";
  _varSet.reset((Environment::Factory<ConvectiveVarSet>::getInstance().getProvider(varSetName)->
    create(physModel->getImplementor()->getConvectiveTerm())).d_castTo<Physics::NavierStokes::Euler3DCons>());

  cf_assert(_varSet.isNotNull());

}

//////////////////////////////////////////////////////////////////////////////

void UnsteadyWeakSlipWallEuler3DCons::executeOnTrs()
{
  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  const CFuint dim = PhysicalModelStack::getActive()->getDim();
  const CFreal oEminusAlpha = 1. - _alpha;
  const CFreal third = 1./3.;
  _dt = SubSystemStatusStack::getActive()->getDT();
  RealVector normal(0.0,dim);

  // get the data handle for the rhs
  DataHandle< CFreal> rhs = socket_rhs.getDataHandle();

  // get the data handle for the states
  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();

  // get the data handle for the past states
  DataHandle<State*> pastStates = socket_pastStates.getDataHandle();

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

    const vector<State*> *const statesInFace =
      currFace->getStates();
    const CFuint nbStatesInFace = statesInFace->size();

    ///@todo do this more general for any number of nodes
    cf_assert(nbStatesInFace == 3);

    State state0;
    State state1;
    State state2;
    const CFuint stateID0 = (*statesInFace)[0]->getLocalID();
    const CFuint stateID1 = (*statesInFace)[1]->getLocalID();
    const CFuint stateID2 = (*statesInFace)[2]->getLocalID();

    // Unsteady, state=0.5*(state+pastState)
    for (CFuint j=0; j<nbEqs ; ++j){
      state0[j] = 0.5*((*states[stateID0])[j] + (*pastStates[stateID0])[j]);
      state1[j] = 0.5*((*states[stateID1])[j] + (*pastStates[stateID1])[j]);
      state2[j] = 0.5*((*states[stateID2])[j] + (*pastStates[stateID2])[j]);
    }

    // if Moving Mesh, compute the wall velocity
    for (CFuint iDim=0; iDim < dim; ++iDim){
      _wallSpeed[iDim] = 0.;
    }
    if (SubSystemStatusStack::getActive()->isMovingMesh()){
      DataHandle < Framework::Node*, Framework::GLOBAL > nodes = socket_nodes.getDataHandle();
      DataHandle<Node*> pastNodes = socket_pastNodes.getDataHandle();

      for (CFuint iDim=0; iDim < dim; ++iDim){
	_wallSpeed[iDim]  = (*nodes[stateID0])[iDim]-(*pastNodes[stateID0])[iDim];
	_wallSpeed[iDim] += (*nodes[stateID1])[iDim]-(*pastNodes[stateID1])[iDim];
	_wallSpeed[iDim] += (*nodes[stateID2])[iDim]-(*pastNodes[stateID2])[iDim];
	_wallSpeed[iDim] *= third;
	_wallSpeed[iDim] /= _dt;
      }
    }

    for (CFuint is = 0; is < nbStatesInFace; ++is) {
      // compute the normal fluxes corrections for both the states
      // of this cell
      computeNormalFlux(*(*statesInFace)[is],
			normal,
			_fluxes[is]);
    }


    // distribute contributions to the nodes
    for (CFuint iEq = 0; iEq < nbEqs; ++iEq) {
      if (!isUpdated[stateID0]) {
	rhs(stateID0, iEq, nbEqs) -=
	  third*(_alpha*_fluxes[0][iEq] + oEminusAlpha*0.5*
		 (_fluxes[1][iEq] + _fluxes[2][iEq]));
	_flagState[stateID0] = true;
      }

      if (!isUpdated[stateID1]) {
	rhs(stateID1, iEq, nbEqs) -=
	  third*(_alpha*_fluxes[1][iEq] + oEminusAlpha*0.5*
		 (_fluxes[0][iEq] + _fluxes[2][iEq]));
	_flagState[stateID1] = true;
      }

      if (!isUpdated[stateID2]) {
	rhs(stateID2, iEq, nbEqs) -=
	  third*(_alpha*_fluxes[2][iEq] + oEminusAlpha*0.5*
		 (_fluxes[0][iEq] + _fluxes[1][iEq]));
	_flagState[stateID2] = true;
      }
    }

    // release the face
    geoBuilder->releaseGE();
  }
}

//////////////////////////////////////////////////////////////////////////////

void UnsteadyWeakSlipWallEuler3DCons::computeNormalFlux(
               const RealVector& state,
               const RealVector& normal,
               RealVector& flux) const
 {
   const CFreal nx = normal[0];
   const CFreal ny = normal[1];
   const CFreal nz = normal[2];

   const CFreal rho = state[0];
   const CFreal u = state[1]/rho;
   const CFreal v = state[2]/rho;
   const CFreal w = state[3]/rho;
   const CFreal un = u*nx + v*ny + w*nz;
   const CFreal gamma = _varSet->getModel()->getGamma();
   const CFreal halfGammaMinus1 = 0.5*(gamma - 1.);
   // unused // const CFreal gammaMinus1 = _varSet.getModel()->getGammaMinus1();
   const CFreal rhoV2 = rho*(u*u + v*v +w*w);
   const CFreal velocity = nx*_wallSpeed[0] + ny*_wallSpeed[1] + nz*_wallSpeed[2];
   const CFreal rhoH = gamma*state[4] - halfGammaMinus1*rhoV2;

   // Compute Flux for non-moving wall
   flux[0] = rho*un;
   flux[1] = un*state[1];
   flux[2] = un*state[2];
   flux[3] = un*state[3];
   flux[4] = un*rhoH;

   if (SubSystemStatusStack::getActive()->isMovingMesh()){
    // Modification of the flux due to the movement of the wall
    flux[0] -= state[0]*velocity;
    flux[1] -= state[1]*velocity;
    flux[2] -= state[2]*velocity;
    flux[3] -= state[3]*velocity;
    flux[4] -= rhoH*velocity;
    }

   // Scale with dt
   flux[0] *= _dt;
   flux[1] *= _dt;
   flux[2] *= _dt;
   flux[3] *= _dt;
   flux[4] *= _dt;

 }

//////////////////////////////////////////////////////////////////////////////

void UnsteadyWeakSlipWallEuler3DCons::setFaceNormal(const CFuint faceID,
						    RealVector& normal)
{
  DataHandle< InwardNormalsData*> normals = socket_normals.getDataHandle();
  DataHandle< InwardNormalsData*> pastNormals = socket_pastNormals.getDataHandle();
  DataHandle<Common::Trio<CFuint, CFuint, Common::SafePtr<Framework::TopologicalRegionSet> > >
    faceNeighCell = socket_faceNeighCell.getDataHandle();

  const CFuint cellTrsID = faceNeighCell[faceID].first;
  const CFuint iFaceLocal = faceNeighCell[faceID].second;
  Common::SafePtr<TopologicalRegionSet> cellTrs = faceNeighCell[faceID].third;
  const CFuint cellLocalID = cellTrs->getLocalGeoID(cellTrsID);

  // If the wall is moving, the normal is the average of
  // the normals at the two instants
  if (SubSystemStatusStack::getActive()->isMovingMesh()){
    DataHandle< InwardNormalsData*> pastNormals = socket_pastNormals.getDataHandle();
    normal[0] = -0.5*(normals[cellLocalID]->getFaceNormComp(iFaceLocal, 0)+pastNormals[cellLocalID]->getFaceNormComp(iFaceLocal, 0));
    normal[1] = -0.5*(normals[cellLocalID]->getFaceNormComp(iFaceLocal, 1)+pastNormals[cellLocalID]->getFaceNormComp(iFaceLocal, 1));
    normal[2] = -0.5*(normals[cellLocalID]->getFaceNormComp(iFaceLocal, 2)+pastNormals[cellLocalID]->getFaceNormComp(iFaceLocal, 2));

  }
  else{
    normal[0] = -normals[cellLocalID]->getFaceNormComp(iFaceLocal, 0);
    normal[1] = -normals[cellLocalID]->getFaceNormComp(iFaceLocal, 1);
    normal[2] = -normals[cellLocalID]->getFaceNormComp(iFaceLocal, 2);
  }
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
