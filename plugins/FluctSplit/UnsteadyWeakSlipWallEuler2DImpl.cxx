#include "FluctSplit/FluctSplitSpaceTimeNavierStokes.hh"
#include "UnsteadyWeakSlipWallEuler2DImpl.hh"
#include "InwardNormalsData.hh"
#include "Framework/PhysicalModel.hh"
#include "Framework/MethodCommandProvider.hh"
#include "NavierStokes/Euler2DVarSet.hh"
#include "Framework/MeshData.hh"
#include "Framework/SubSystemStatus.hh"
#include "Framework/BlockAccumulator.hh"
#include "Framework/LSSMatrix.hh"

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

MethodCommandProvider<UnsteadyWeakSlipWallEuler2DImpl,
		      FluctuationSplitData,
		      FluctSplitSpaceTimeNavierStokesModule>
UnsteadyWeakSlipWallEuler2DImplProvider("UnsteadyWeakSlipWallEuler2DImpl");

//////////////////////////////////////////////////////////////////////////////

UnsteadyWeakSlipWallEuler2DImpl::UnsteadyWeakSlipWallEuler2DImpl(const std::string& name) :
  WeakSlipWall2DImpl(name),
  socket_pastStates("pastStates"),
  socket_pastNormals("pastNormals",false),
  socket_nodes("nodes",false),
  socket_pastNodes("pastNodes",false),
  _varSet(CFNULL)
{
}

//////////////////////////////////////////////////////////////////////////////

UnsteadyWeakSlipWallEuler2DImpl::~UnsteadyWeakSlipWallEuler2DImpl()
{
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
UnsteadyWeakSlipWallEuler2DImpl::needsSockets()
{

  std::vector<Common::SafePtr<BaseDataSocketSink> > result = WeakSlipWall2DImpl::needsSockets();

  result.push_back(&socket_pastStates);
  result.push_back(&socket_pastNormals);
  result.push_back(&socket_nodes);
  result.push_back(&socket_pastNodes);

  return result;
}


//////////////////////////////////////////////////////////////////////////////

void UnsteadyWeakSlipWallEuler2DImpl::setup()
{
  WeakSlipWall2DImpl::setup();

  _varSet = getMethodData().getUpdateVar().d_castTo<Euler2DVarSet>();

  // set up the physical data
  _varSet->getModel()->resizePhysicalData(_physicalData);

  _wallSpeed.resize(PhysicalModelStack::getActive()->getDim());
}

//////////////////////////////////////////////////////////////////////////////

void UnsteadyWeakSlipWallEuler2DImpl::executeOnTrs()
{

  DataHandle<bool> isUpdated = socket_isUpdated.getDataHandle();
  DataHandle<CFreal> rhs = socket_rhs.getDataHandle();
  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();
  DataHandle<State*> pastStates = socket_pastStates.getDataHandle();

  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  const CFuint dim = PhysicalModelStack::getActive()->getDim();
  const CFreal dt = SubSystemStatusStack::getActive()->getDT();

  RealVector normal(dim);


  SafePtr<LinearSystemSolver> lss =
    getMethodData().getLinearSystemSolver()[0];

  SafePtr<LSSMatrix> jacobMatrix = lss->getMatrix();

  // block accumulator 2*2


  Common::SafePtr<GeometricEntityPool<StdTrsGeoBuilder> >
    geoBuilder = getMethodData().getStdTrsGeoBuilder();

  StdTrsGeoBuilder::GeoData& geoData = geoBuilder->getDataGE();
  geoData.trs = getCurrentTRS();

  const CFuint nbFaces = getCurrentTRS()->getLocalNbGeoEnts();
  for (CFuint iFace = 0; iFace < nbFaces; ++iFace) {
    geoData.idx = iFace;
    GeometricEntity *const currFace = geoBuilder->buildGE();
    const CFuint nb_state = currFace->nbStates();
    const CFreal OEminusAlpha_over_nb_sate = (1.0/nb_state)*(1. - m_alpha);
    const CFreal Alpha_over_nb_sate = (1.0/nb_state)*m_alpha;
    // set the face normal
    setFaceNormal(currFace->getID(), normal);

    State state0;
    State state1;
    const CFuint stateID0 = currFace->getState(0)->getLocalID();
    const CFuint stateID1 = currFace->getState(1)->getLocalID();

    // Unsteady, state=0.5*(state+pastState)
    for (CFuint iEq=0; iEq<nbEqs ; ++iEq){
      state0[iEq] = 0.5*((*states[stateID0])[iEq] + (*pastStates[stateID0])[iEq]);
      state1[iEq] = 0.5*((*states[stateID1])[iEq] + (*pastStates[stateID1])[iEq]);
    }

    // if Moving Mesh, compute the wall velocity
    if (SubSystemStatusStack::getActive()->isMovingMesh()){
      DataHandle < Framework::Node*, Framework::GLOBAL > nodes = socket_nodes.getDataHandle();
      DataHandle<Node*> pastNodes = socket_pastNodes.getDataHandle();

      for (CFuint iDim=0; iDim < dim; ++iDim){
	_wallSpeed[iDim]  = (*nodes[stateID0])[iDim]-(*pastNodes[stateID0])[iDim];
	_wallSpeed[iDim] += (*nodes[stateID1])[iDim]-(*pastNodes[stateID1])[iDim];
	_wallSpeed[iDim] *= 0.5;
	_wallSpeed[iDim] /= dt;
      }
    }
    else{
      for (CFuint iDim=0; iDim < dim; ++iDim){
        _wallSpeed[iDim] = 0.;
      }
    }


     auto_ptr<BlockAccumulator> acc(lss->createBlockAccumulator(nb_state, nb_state, nbEqs));


    for ( CFuint iState = 0; iState < nb_state; ++iState){
      m_state[iState] = currFace->getState(iState);
      acc->setRowColIndex(iState, m_state[iState]->getLocalID());
    }


    // compute the normal fluxes corrections for both the states
    // of this cell
    computeNormalFluxAndJacob(state0, normal, m_flux[0], m_fluxJacob[0]);
    computeNormalFluxAndJacob(state1, normal, m_flux[1], m_fluxJacob[1]);

    // distribute contributions to the two nodes

    for (CFuint iEq = 0; iEq < nbEqs; ++iEq) {
      for (CFuint iState = 0; iState < nb_state; ++iState){
//         const CFuint stateID = currFace->getState(iState)->getLocalID();
        CFreal sum_flux_not_iState = 0.;
        for (CFuint jState = 0; jState < nb_state ; ++jState){
          if (jState != iState)
           sum_flux_not_iState += m_flux[jState][iEq];
        }
        rhs(m_state[iState]->getLocalID(), iEq, nbEqs) -=
                  Alpha_over_nb_sate*m_flux[iState][iEq] + OEminusAlpha_over_nb_sate*sum_flux_not_iState;
      }
    }

    // compute the jacobian contributions for the involved states

    for (CFuint iEq = 0; iEq < nbEqs; ++iEq) {
      for (CFuint jEq = 0; jEq < nbEqs; ++jEq) {
        for (CFuint iState = 0; iState < nb_state; ++iState){
//           const CFuint stateID = currFace->getState(iState)->getLocalID();
          if (m_state[iState]->isParUpdatable()) {
            acc->setValue(iState, iState, iEq, jEq,Alpha_over_nb_sate*m_fluxJacob[iState](iEq, jEq));
            for (CFuint jState = 0; jState < nb_state; ++jState){
              if (iState != jState)
                acc->setValue(iState, jState, iEq, jEq, OEminusAlpha_over_nb_sate*m_fluxJacob[jState](iEq, jEq));
              }
            }
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

void UnsteadyWeakSlipWallEuler2DImpl::setFaceNormal
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
    normal[XX] = -0.5*(normals[cellLocalID]->getFaceNormComp(iFaceLocal, 0) + pastNormals[cellLocalID]->getFaceNormComp(iFaceLocal, 0));
    normal[YY] = -0.5*(normals[cellLocalID]->getFaceNormComp(iFaceLocal, 1) + pastNormals[cellLocalID]->getFaceNormComp(iFaceLocal, 1));
  }
  else{
    normal[XX] = -normals[cellLocalID]->getFaceNormComp(iFaceLocal, 0);
    normal[YY] = -normals[cellLocalID]->getFaceNormComp(iFaceLocal, 1);
  }
 }

//////////////////////////////////////////////////////////////////////////////

void UnsteadyWeakSlipWallEuler2DImpl::computeNormalFluxAndJacob
(const State& state,
 const RealVector& normal,
 RealVector& flux,
 RealMatrix& fluxJacob)
{
  State& ss = *(const_cast<State*>(&state));
  _varSet->computePhysicalData(ss, _physicalData);

  const CFreal nx = normal[XX];
  const CFreal ny = normal[YY];
  const CFreal rho = _physicalData[EulerTerm::RHO];
  const CFreal u = _physicalData[EulerTerm::VX];
  const CFreal v = _physicalData[EulerTerm::VY];
  const CFreal un = u*nx + v*ny;
  const CFreal vel2 = u*u + v*v;
  const CFreal H = _physicalData[EulerTerm::H];
  const CFreal gamma = _varSet->getModel()->getGamma();
  const CFreal gammaMinus1 = gamma - 1.;
  const CFreal halfGammaMinus1 = 0.5*gammaMinus1;
  const CFreal dt = SubSystemStatusStack::getActive()->getDT();
  const CFreal velocity = nx*_wallSpeed[0] + ny*_wallSpeed[1];

  flux[0] = rho*un;
  flux[1] = un*rho*u;
  flux[2] = un*rho*v;
  flux[3] = un*rho*H;

  if (SubSystemStatusStack::getActive()->isMovingMesh())
  {
    // Modification of the flux due to the movement of the wall
    flux[0] -= rho*velocity;
    flux[1] -= rho*u*velocity;
    flux[2] -= rho*v*velocity;
    flux[3] -= rho*H*velocity;
  }

  // Scale with dt
  flux[0] *= dt;
  flux[1] *= dt;
  flux[2] *= dt;
  flux[3] *= dt;

   if (!getMethodData().isResidualTransformationNeeded()) {
     // the analytical jacobian of the normal fluxes
     fluxJacob(0,1) = nx;
     fluxJacob(0,2) = ny;
     fluxJacob(1,0) = -u*un;
     fluxJacob(1,1) = un + u*nx;
     fluxJacob(1,2) = u*ny;
     fluxJacob(2,0) = -v*un;
     fluxJacob(2,1) = v*nx;
     fluxJacob(2,2) = un + v*ny;
     fluxJacob(3,0) = un*(halfGammaMinus1*vel2 - H);
     fluxJacob(3,1) = H*nx - un*u*gammaMinus1;
     fluxJacob(3,2) = H*ny - un*v*gammaMinus1;
     fluxJacob(3,3) = _varSet->getModel()->getGamma()*un;

     if (SubSystemStatusStack::getActive()->isMovingMesh())
     {
       fluxJacob *= dt;

       // Modification of the jacobian of the flux
       // due to the movement of the wall
       fluxJacob(0,0) -= dt*velocity;
       fluxJacob(1,1) -= dt*velocity;
       fluxJacob(2,2) -= dt*velocity;
       fluxJacob(3,3) -= dt*gamma*velocity;
     }
   }
   else {
     // the analytical jacobian of the normal fluxes
     _tJacob(0,1) = nx;
     _tJacob(0,2) = ny;
     _tJacob(1,0) = -u*un;
     _tJacob(1,1) = un + u*nx;
     _tJacob(1,2) = u*ny;
     _tJacob(2,0) = -v*un;
     _tJacob(2,1) = v*nx;
     _tJacob(2,2) = un + v*ny;
     _tJacob(3,0) = un*(halfGammaMinus1*vel2 - H);
     _tJacob(3,1) = H*nx - un*u*gammaMinus1;
     _tJacob(3,2) = H*ny - un*v*gammaMinus1;
     _tJacob(3,3) = _varSet->getModel()->getGamma()*un;

     if (SubSystemStatusStack::getActive()->isMovingMesh())
     {
       _tJacob *= dt;
       // Modification of the jacobian of the flux
       // due to the movement of the wall
       _tJacob(0,0) -= dt*velocity;
       _tJacob(1,1) -= dt*velocity;
       _tJacob(2,2) -= dt*velocity;
       _tJacob(3,3) -= dt*gamma*velocity;
     }

     // set the transformation from update to solution in update
     SafePtr<VarSetMatrixTransformer> updateToSolInUpdate =
       getMethodData().getUpdateToSolutionInUpdateMatTrans();

     updateToSolInUpdate->setMatrix(state);
     const RealMatrix& tMatrix = *updateToSolInUpdate->getMatrix();

     fluxJacob = _tJacob*tMatrix;
   }
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
