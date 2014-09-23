#include "FluctSplit/FluctSplitSpaceTimeNavierStokes.hh"
#include "UnsteadyWeakSlipWallEuler2D.hh"
#include "InwardNormalsData.hh"
#include "Framework/CFL.hh"
#include "Framework/SubSystemStatus.hh"
#include "Framework/PhysicalModel.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/MeshData.hh"
#include "Framework/NamespaceSwitcher.hh"
#include "NavierStokes/Euler2DVarSet.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Physics::NavierStokes;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<UnsteadyWeakSlipWallEuler2D, FluctuationSplitData, 
		      FluctSplitSpaceTimeNavierStokesModule> 
unsteadyWeakSlipWallEuler2DProvider("UnsteadyWeakSlipWallEuler2D");

MethodCommandProvider<UnsteadyWeakSlipWallEuler2D, FluctuationSplitData, 
		      FluctSplitSpaceTimeNavierStokesModule> 
unsteadyWeakSlipWallEuler2DconsProvider("UnsteadyWeakSlipWallEuler2DCons");

//////////////////////////////////////////////////////////////////////////////

UnsteadyWeakSlipWallEuler2D::UnsteadyWeakSlipWallEuler2D(const std::string& name) :
  WeakSlipWallEuler2D(name),
  socket_pastNormals("pastNormals", false), // only used if moving
  socket_nodes("nodes"),
  socket_pastNodes("pastNodes", false),     // only used if moving
  socket_pastStates("pastStates"),
  _dt(),
  _wallSpeed()
{
}
      
//////////////////////////////////////////////////////////////////////////////

UnsteadyWeakSlipWallEuler2D::~UnsteadyWeakSlipWallEuler2D()
{
}

//////////////////////////////////////////////////////////////////////////////

void UnsteadyWeakSlipWallEuler2D::setup()
{
  WeakSlipWallEuler2D::setup();
  
  _wallSpeed.resize(PhysicalModelStack::getActive()->getDim());
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
UnsteadyWeakSlipWallEuler2D::needsSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSink> > result = 
    WeakSlipWallEuler2D::needsSockets();
  
  result.push_back(&socket_pastNormals);
  result.push_back(&socket_nodes);
  result.push_back(&socket_pastNodes);
  result.push_back(&socket_pastStates);
  return result;
}
      
//////////////////////////////////////////////////////////////////////////////

void UnsteadyWeakSlipWallEuler2D::executeOnTrs()
{
  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  const CFuint dim = PhysicalModelStack::getActive()->getDim();
  const CFreal oEminusAlpha = 1. - _alpha;
  RealVector normal(dim);
  RealVector flux0(0.0, nbEqs);
  RealVector flux1(0.0, nbEqs);
  _dt = SubSystemStatusStack::getActive()->getDT();

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

  State state0;
  State state1;

  const CFuint nbFaces = getCurrentTRS()->getLocalNbGeoEnts();
  for (CFuint iFace = 0; iFace < nbFaces; ++iFace) {
    geoData.idx = iFace;
    GeometricEntity *const currFace = geoBuilder->buildGE();

    // set the face normal
    setFaceNormal(currFace->getID(), normal);

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

    // compute the normal fluxes corrections for both the states
    // of this cell
    computeNormalFlux(state0, normal, flux0);
    computeNormalFlux(state1, normal, flux1);

    // distribute contributions to the two nodes
    for (CFuint iEq = 0; iEq < nbEqs; ++iEq) {
      if (!isUpdated[stateID0]) {
	rhs(stateID0, iEq, nbEqs) -= 0.5*(_alpha*flux0[iEq] + oEminusAlpha*flux1[iEq]);
	_flagState[stateID0] = true;
      }
      
      if (!isUpdated[stateID1]) {
	rhs(stateID1, iEq, nbEqs) -= 0.5*(_alpha*flux1[iEq] + oEminusAlpha*flux0[iEq]);
	_flagState[stateID1] = true;
      }
    }

    // release the face
    geoBuilder->releaseGE();
  }
}

//////////////////////////////////////////////////////////////////////////////

void UnsteadyWeakSlipWallEuler2D::computeNormalFlux(const Framework::State& state,
						    const RealVector& normal,
						    RealVector& flux)
{
  State& ss = *(const_cast<State*>(&state));
  _varSet->computePhysicalData(ss, _physicalData);
  
  const CFreal un = _physicalData[EulerTerm::VX]*normal[XX] +
    _physicalData[EulerTerm::VY]*normal[YY];
  const CFreal rho = _physicalData[EulerTerm::RHO];
  const CFreal rhoUn = rho*un;
  
  flux[0] = rhoUn;
  flux[1] = rhoUn*_physicalData[EulerTerm::VX];
  flux[2] = rhoUn*_physicalData[EulerTerm::VY];
  flux[3] = rhoUn*_physicalData[EulerTerm::H];
  
  if (SubSystemStatusStack::getActive()->isMovingMesh()){
    const CFreal velocity = normal[XX]*_wallSpeed[0] + normal[YY]*_wallSpeed[1];
    // Modification of the flux due to the movement of the wall
    flux[0] -= rho*velocity;
    flux[1] -= rho*_physicalData[EulerTerm::VX]*velocity;
    flux[2] -= rho*_physicalData[EulerTerm::VY]*velocity;
    flux[3] -= rho*_physicalData[EulerTerm::H]*velocity;
  }
  
  // Scale with dt
  flux[0] *= _dt;
  flux[1] *= _dt;
  flux[2] *= _dt;
  flux[3] *= _dt;
}

//////////////////////////////////////////////////////////////////////////////

void UnsteadyWeakSlipWallEuler2D::setFaceNormal(const CFuint faceID,
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
  if (SubSystemStatusStack::getActive()->isMovingMesh()){
    DataHandle< InwardNormalsData*> pastNormals = socket_pastNormals.getDataHandle();
    normal[0] = -0.5*(normals[cellLocalID]->getFaceNormComp(iFaceLocal, 0) + pastNormals[cellLocalID]->getFaceNormComp(iFaceLocal, 0));
    normal[1] = -0.5*(normals[cellLocalID]->getFaceNormComp(iFaceLocal, 1) + pastNormals[cellLocalID]->getFaceNormComp(iFaceLocal, 1));
  }
  else{
    normal[0] = - normals[cellLocalID]->getFaceNormComp(iFaceLocal, 0);
    normal[1] = - normals[cellLocalID]->getFaceNormComp(iFaceLocal, 1);
  }
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace FluctSplit

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
