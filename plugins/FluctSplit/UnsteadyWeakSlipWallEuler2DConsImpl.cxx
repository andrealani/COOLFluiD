#include "FluctSplit/FluctSplitNavierStokes.hh"
#include "UnsteadyWeakSlipWallEuler2DConsImpl.hh"
#include "FluctSplit/InwardNormalsData.hh"
#include "Framework/CFL.hh"
#include "Framework/SubSystemStatus.hh"
#include "Framework/PhysicalModel.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/MeshData.hh"
#include "Framework/LSSMatrix.hh"
#include "Framework/BlockAccumulator.hh"
#include "Framework/NamespaceSwitcher.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Physics::NavierStokes;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<UnsteadyWeakSlipWallEuler2DConsImpl, FluctuationSplitData, FluctSplitNavierStokesModule> unsteadyWeakSlipWallEuler2DConsImplProvider("UnsteadyWeakSlipWallEuler2DConsImpl");

//////////////////////////////////////////////////////////////////////////////

void UnsteadyWeakSlipWallEuler2DConsImpl::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOption< CFreal >("alpha","distribution coefficient");
}

//////////////////////////////////////////////////////////////////////////////

UnsteadyWeakSlipWallEuler2DConsImpl::UnsteadyWeakSlipWallEuler2DConsImpl(const std::string& name) :
  FluctuationSplitCom(name),
  socket_rhs("rhs"),
  socket_states("states"),
  socket_pastStates("pastStates"),
  socket_isUpdated("isUpdated"),
  socket_normals("normals"),
  socket_faceNeighCell("faceNeighCell"),
  socket_pastNormals("pastNormals",false),
  socket_nodes("nodes",false),
  socket_pastNodes("pastNodes",false),
  _varSet(),
  _dt(),
  _wallSpeed(),
  _im0(),
  _in0(),
  _im1(),
  _in1(),
  _flagState()
{
   addConfigOptionsTo(this);
  _alpha = 0.66;
   setParameter("alpha",&_alpha);

}

//////////////////////////////////////////////////////////////////////////////

UnsteadyWeakSlipWallEuler2DConsImpl::~UnsteadyWeakSlipWallEuler2DConsImpl()
{
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
UnsteadyWeakSlipWallEuler2DConsImpl::needsSockets()
{

  std::vector<Common::SafePtr<BaseDataSocketSink> > result;

  result.push_back(&socket_rhs);
  result.push_back(&socket_states);
  result.push_back(&socket_pastStates);
  result.push_back(&socket_isUpdated);
  result.push_back(&socket_normals);
  result.push_back(&socket_faceNeighCell);
  result.push_back(&socket_pastNormals);
  result.push_back(&socket_nodes);
  result.push_back(&socket_pastNodes);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

void UnsteadyWeakSlipWallEuler2DConsImpl::setup()
{

  _wallSpeed.resize(PhysicalModelStack::getActive()->getDim());
  _im0.resize(PhysicalModelStack::getActive()->getNbEq());
  _in0.resize(PhysicalModelStack::getActive()->getNbEq());
  _im1.resize(PhysicalModelStack::getActive()->getNbEq());
  _in1.resize(PhysicalModelStack::getActive()->getNbEq());

  _varSet->setup();
  // get the data handle for the states
  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();

  _flagState.resize(states.size());
  _flagState = false;

}

//////////////////////////////////////////////////////////////////////////////

void UnsteadyWeakSlipWallEuler2DConsImpl::configure ( Config::ConfigArgs& args )
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

void UnsteadyWeakSlipWallEuler2DConsImpl::executeOnTrs()
{

  DataHandle<bool> isUpdated = socket_isUpdated.getDataHandle();
  DataHandle<CFreal> rhs = socket_rhs.getDataHandle();
  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();
  DataHandle<State*> pastStates = socket_pastStates.getDataHandle();

  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  const CFuint dim = PhysicalModelStack::getActive()->getDim();
  const CFreal halfOEminusAlpha = 0.5*(1. - _alpha);
  const CFreal halfAlpha = 0.5*_alpha;
  RealVector normal(dim);
  RealVector flux0(0.0, nbEqs);
  RealVector flux1(0.0, nbEqs);
  RealMatrix fluxJacob0(nbEqs, nbEqs, 0.0);
  RealMatrix fluxJacob1(nbEqs, nbEqs, 0.0);
  _dt = SubSystemStatusStack::getActive()->getDT();

  SafePtr<LinearSystemSolver> lss =
    getMethodData().getLinearSystemSolver()[0];

  SafePtr<LSSMatrix> jacobMatrix = lss->getMatrix();

  // block accumulator 2*2
  BlockAccumulator* acc = getMethodData().getLinearSystemSolver()[0]->
    createBlockAccumulator(2, 2, nbEqs);

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

    State state0;
    State state1;
    const CFuint stateID0 = currFace->getState(0)->getLocalID();
    const CFuint stateID1 = currFace->getState(1)->getLocalID();

    // Unsteady, state=0.5*(state+pastState)
    for (CFuint j=0; j<nbEqs ; ++j){
      state0[j] = 0.5*((*states[stateID0])[j] + (*pastStates[stateID0])[j]);
      state1[j] = 0.5*((*states[stateID1])[j] + (*pastStates[stateID1])[j]);
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
	_wallSpeed[iDim] *= 0.5;
	_wallSpeed[iDim] /= _dt;
      }
    }

    // set the row - column
    acc->setRowColIndex(0, currFace->getState(0)->getLocalID());
    acc->setRowColIndex(1, currFace->getState(1)->getLocalID());

    // compute the normal fluxes corrections for both the states
    // of this cell
    computeNormalFluxAndJacob(state0, normal, flux0, fluxJacob0);
    computeNormalFluxAndJacob(state1, normal, flux1, fluxJacob1);


    // distribute contributions to the two nodes
    for (CFuint iEq = 0; iEq < nbEqs; ++iEq) {
      rhs(stateID0, iEq, nbEqs) -=
	halfAlpha*flux0[iEq] + halfOEminusAlpha*flux1[iEq];

      rhs(stateID1, iEq, nbEqs) -=
	halfAlpha*flux1[iEq] + halfOEminusAlpha*flux0[iEq];
    }

    // compute the jacobian contributions for the involved states
    for (CFuint iEq = 0; iEq < nbEqs; ++iEq) {
      for (CFuint jEq = 0; jEq < nbEqs; ++jEq) {

	if (state0.isParUpdatable()) {
	  acc->setValue(0, 0, iEq, jEq,
			halfAlpha*fluxJacob0(iEq, jEq));

	  acc->setValue(0, 1, iEq, jEq,
			halfOEminusAlpha*fluxJacob1(iEq, jEq));
	}

	if (state1.isParUpdatable()) {
	  acc->setValue(1, 0, iEq, jEq,
			halfOEminusAlpha*fluxJacob0(iEq, jEq));

	  acc->setValue(1, 1, iEq, jEq,
			halfAlpha*fluxJacob1(iEq, jEq));
	}
      }
    }

    // add the contributions to the jacobian
    jacobMatrix->addValues(*acc);

    acc->reset();

    // release the face
    geoBuilder->releaseGE();
  }
}

//////////////////////////////////////////////////////////////////////////////

void UnsteadyWeakSlipWallEuler2DConsImpl::computeNormalFluxAndJacob(
               const RealVector& state,
               const RealVector& normal,
               RealVector& flux,
               RealMatrix& fluxJacob) const
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
  const CFreal gammaMinus1 = gamma - 1.;
  const CFreal halfGammaMinus1 = 0.5*gammaMinus1;
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

  if (SubSystemStatusStack::getActive()->isMovingMesh())
  {
    // Modification of the flux due to the movement of the wall
    flux[0] -= rho*velocity;
    flux[1] -= rhoU*velocity;
    flux[2] -= rhoV*velocity;
    flux[3] -= rhoH*velocity;
  }

  // Scale with dt
  flux[0] *= _dt;
  flux[1] *= _dt;
  flux[2] *= _dt;
  flux[3] *= _dt;

  // the analytical jacobian of the normal fluxes
  fluxJacob(0,1) = _dt*nx;
  fluxJacob(0,2) = _dt*ny;
  fluxJacob(1,0) = -_dt*u*un;
  fluxJacob(1,1) = _dt*(un + u*nx);
  fluxJacob(1,2) = _dt*u*ny;
  fluxJacob(2,0) = -_dt*v*un;
  fluxJacob(2,1) = _dt*v*nx;
  fluxJacob(2,2) = _dt*(un + v*ny);
  fluxJacob(3,0) = _dt*un*(halfGammaMinus1*vel2 - rhoH/rho);
  fluxJacob(3,1) = _dt*(rhoH/rho*nx - un*u*gammaMinus1);
  fluxJacob(3,2) = _dt*(rhoH/rho*ny - un*v*gammaMinus1);
  fluxJacob(3,3) = _dt*(gamma*un);

  if (SubSystemStatusStack::getActive()->isMovingMesh())
  {
    // Modification of the jacobian of the flux
    // due to the movement of the wall
    fluxJacob(0,0) -= _dt*velocity;
    fluxJacob(1,1) -= _dt*velocity;
    fluxJacob(2,2) -= _dt*velocity;
    fluxJacob(3,3) -= _dt*gamma*velocity;
  }
 }

//////////////////////////////////////////////////////////////////////////////

void UnsteadyWeakSlipWallEuler2DConsImpl::setFaceNormal
(const CFuint faceID,
 RealVector& normal)
{
  // get the data handle for the states
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
    DataHandle< InwardNormalsData*> pastNormals = socket_pastNormals.getDataHandle();
    normal[0] = -0.5*(normals[cellLocalID]->getFaceNormComp(iFaceLocal, 0) + pastNormals[cellLocalID]->getFaceNormComp(iFaceLocal, 0));
    normal[1] = -0.5*(normals[cellLocalID]->getFaceNormComp(iFaceLocal, 1) + pastNormals[cellLocalID]->getFaceNormComp(iFaceLocal, 1));
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
