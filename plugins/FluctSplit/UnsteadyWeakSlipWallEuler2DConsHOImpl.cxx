#include "FluctSplit/FluctSplitNavierStokes.hh"
#include "UnsteadyWeakSlipWallEuler2DConsHOImpl.hh"
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

MethodCommandProvider<UnsteadyWeakSlipWallEuler2DConsHOImpl, FluctuationSplitData, FluctSplitNavierStokesModule> UnsteadyWeakSlipWallEuler2DConsHOImplProvider("UnsteadyWeakSlipWallEuler2DConsHOImpl");

//////////////////////////////////////////////////////////////////////////////

void UnsteadyWeakSlipWallEuler2DConsHOImpl::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOption< CFreal >("alpha","distribution coefficient");
}

//////////////////////////////////////////////////////////////////////////////

UnsteadyWeakSlipWallEuler2DConsHOImpl::UnsteadyWeakSlipWallEuler2DConsHOImpl(const std::string& name) :
  FluctuationSplitCom(name),
  socket_rhs("rhs"),
  socket_states("states"),
  socket_pastStates("pastStates"),
  socket_interStates("interStates"),
  socket_isUpdated("isUpdated"),
  socket_normals("normals"),
  socket_faceNeighCell("faceNeighCell"),
  socket_nodes("nodes",false),
  _varSet(),
  _dt(),
  _wallSpeed(),
  _im0(),
  _in0(),
  _im1(),
  _in1(),
  _flagState(),
  m_fluxJacobn1(),
  m_fluxJacobn12(),
  m_fluxJacobn(),
  m_fluxJacob()
{              
   addConfigOptionsTo(this);
   m_alpha = 0.66;
   setParameter("alpha",&m_alpha);

}

//////////////////////////////////////////////////////////////////////////////

UnsteadyWeakSlipWallEuler2DConsHOImpl::~UnsteadyWeakSlipWallEuler2DConsHOImpl()
{
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
UnsteadyWeakSlipWallEuler2DConsHOImpl::needsSockets()
{

  std::vector<Common::SafePtr<BaseDataSocketSink> > result;

  result.push_back(&socket_rhs);
  result.push_back(&socket_states);
  result.push_back(&socket_pastStates);
  result.push_back(&socket_isUpdated);
  result.push_back(&socket_normals);
  result.push_back(&socket_faceNeighCell);
  result.push_back(&socket_interStates);
  result.push_back(&socket_nodes);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

void UnsteadyWeakSlipWallEuler2DConsHOImpl::setup()
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

  // Flux of n-1
  m_fluxn1.resize(3); //There are 3 nodes per face because we are considering a P2 element
  // Flux of n-2
  m_fluxn12.resize(3); //There are 3 nodes per face because we are considering a P2 element
  // Flux of n
  m_fluxn.resize(3); //There are 3 nodes per face because we are considering a P2 element
  // Total Flux
  m_flux.resize(3); //There are 3 nodes per face because we are considering a P2 element
  // Jacobian flux of n-1
  m_fluxJacobn1.resize(3); //There are 3 nodes per face because we are considering a P2 element
  // Jacobian Flux of n-2
  m_fluxJacobn12.resize(3); //There are 3 nodes per face because we are considering a P2 element
  // Jacobian Flux of n
  m_fluxJacobn.resize(3); //There are 3 nodes per face because we are considering a P2 element
  // Total Jacobian Flux
  m_fluxJacob.resize(3); //There are 3 nodes per face because we are considering a P2 element


  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  for (CFuint iState = 0; iState<3; ++iState ){
    m_fluxn1[iState].resize(nbEqs);
    m_fluxn12[iState].resize(nbEqs);
    m_fluxn[iState].resize(nbEqs);
    m_flux[iState].resize(nbEqs);
    m_fluxJacobn1[iState].resize(nbEqs,nbEqs);
    m_fluxJacobn12[iState].resize(nbEqs,nbEqs);
    m_fluxJacobn[iState].resize(nbEqs,nbEqs);
    m_fluxJacob[iState].resize(nbEqs,nbEqs);
  }
}

//////////////////////////////////////////////////////////////////////////////

void UnsteadyWeakSlipWallEuler2DConsHOImpl::configure ( Config::ConfigArgs& args )
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

void UnsteadyWeakSlipWallEuler2DConsHOImpl::executeOnTrs()
{


// The first iteration needs a special function because at this point
  // only two layer are available althought we need three.

  CFuint nbIter = SubSystemStatusStack::getActive()->getNbIter();
  if(nbIter == 1) docomputeFirstBC();
  else docomputeBC();
}

//////////////////////////////////////////////////////////////////////////////
void UnsteadyWeakSlipWallEuler2DConsHOImpl::docomputeFirstBC(){
  DataHandle<bool> isUpdated = socket_isUpdated.getDataHandle();
  DataHandle<CFreal> rhs = socket_rhs.getDataHandle();
  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();
  DataHandle<State*> interStates = socket_interStates.getDataHandle();
  
  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  const CFuint dim = PhysicalModelStack::getActive()->getDim();
  RealVector normal(dim);
  CFreal halfoEminusAlpha = (1.0-m_alpha)/2.0;
  CFreal onesixth = (1.0/6.0);
  CFreal twothird = (2.0/3.0);
  _dt = SubSystemStatusStack::getActive()->getDT();
  
  SafePtr<LinearSystemSolver> lss =
    getMethodData().getLinearSystemSolver()[0];

  SafePtr<LSSMatrix> jacobMatrix = lss->getMatrix();
  

  
  Common::SafePtr<GeometricEntityPool<StdTrsGeoBuilder> >
    geoBuilder = getMethodData().getStdTrsGeoBuilder();

  StdTrsGeoBuilder::GeoData& geoData = geoBuilder->getDataGE();
  geoData.trs = getCurrentTRS();
  
  const CFuint nbFaces = getCurrentTRS()->getLocalNbGeoEnts();
  for (CFuint iFace = 0; iFace < nbFaces; ++iFace) {
    geoData.idx = iFace;
    GeometricEntity *const currFace = geoBuilder->buildGE();

    const CFuint nb_state = currFace->nbStates();


     auto_ptr<BlockAccumulator> acc(lss->createBlockAccumulator(nb_state, nb_state, nbEqs));
for ( CFuint iState = 0; iState < nb_state; ++iState){
      
      acc->setRowColIndex(iState, currFace->getState(iState)->getLocalID());
    }
    // set the face normal
    setFaceNormal(currFace->getID(), normal);
  
    const CFuint stateID0 = currFace->getState(0)->getLocalID();
    const CFuint stateID1 = currFace->getState(1)->getLocalID();
    const CFuint stateID2 = currFace->getState(2)->getLocalID();
    

    // compute the normal fluxes corrections for both the states
    // of this cell
    computeNormalFluxAndJacob((*states[stateID0]), normal, m_fluxn[0], m_fluxJacobn[0]);
    computeNormalFluxAndJacob((*states[stateID1]), normal, m_fluxn[1], m_fluxJacobn[1]);
    computeNormalFluxAndJacob((*states[stateID2]), normal, m_fluxn[2], m_fluxJacobn[2]);
    computeNormalFluxAndJacob((*interStates[stateID0]), normal, m_fluxn1[0], m_fluxJacobn1[0]);
    computeNormalFluxAndJacob((*interStates[stateID1]), normal, m_fluxn1[1], m_fluxJacobn1[1]);
    computeNormalFluxAndJacob((*interStates[stateID2]), normal, m_fluxn1[2], m_fluxJacobn1[2]);
  
    m_flux[0] = _dt*0.5*(m_fluxn1[0]+m_fluxn[0]);
    m_flux[1] = _dt*0.5*(m_fluxn1[1]+m_fluxn[1]);
    m_flux[2] = _dt*0.5*(m_fluxn1[2]+m_fluxn[2]);
  
    m_fluxJacob[0] = _dt*0.5*(m_fluxJacobn1[0]+m_fluxJacobn[0]);
    m_fluxJacob[1] = _dt*0.5*(m_fluxJacobn1[1]+m_fluxJacobn[1]);
    m_fluxJacob[2] = _dt*0.5*(m_fluxJacobn1[2]+m_fluxJacobn[2]);
      

 // distribute contributions to the three nodes
    for (CFuint iEq = 0; iEq < nbEqs; ++iEq) {
      if (!isUpdated[stateID0]){
        rhs(stateID0, iEq, nbEqs) -= m_alpha*onesixth*m_flux[0][iEq] + halfoEminusAlpha*onesixth*m_flux[1][iEq] + halfoEminusAlpha*twothird*m_flux[2][iEq];

        _flagState[stateID0] = true;
      }

      if (!isUpdated[stateID1]){
        rhs(stateID1, iEq, nbEqs) -= m_alpha*onesixth*m_flux[1][iEq] + halfoEminusAlpha*onesixth*m_flux[0][iEq] + halfoEminusAlpha*twothird*m_flux[2][iEq];

        _flagState[stateID1] = true;
      }

      if (!isUpdated[stateID2]){
        rhs(stateID2, iEq, nbEqs) -= m_alpha*twothird*m_flux[2][iEq] + halfoEminusAlpha*onesixth*m_flux[1][iEq] + halfoEminusAlpha*onesixth*m_flux[0][iEq];

        _flagState[stateID2] = true;
      }

    }
  if (getMethodData().doComputeJacobian()) {
    // compute the jacobian contributions for the involved states

    for (CFuint iEq = 0; iEq < nbEqs; ++iEq) {
      for (CFuint jEq = 0; jEq < nbEqs; ++jEq) {
        if (!isUpdated[stateID0]){
          acc->setValue(0, 0, iEq, jEq,            m_alpha*onesixth*m_fluxJacob[0](iEq, jEq));
          acc->setValue(0, 1, iEq, jEq, halfoEminusAlpha*onesixth*m_fluxJacob[1](iEq, jEq));
          acc->setValue(0, 2, iEq, jEq, halfoEminusAlpha*twothird*m_fluxJacob[2](iEq, jEq));
        }

        if (!isUpdated[stateID1]){
          acc->setValue(1, 0, iEq, jEq, halfoEminusAlpha*onesixth*m_fluxJacob[0](iEq, jEq));
          acc->setValue(1, 1, iEq, jEq,            m_alpha*onesixth*m_fluxJacob[1](iEq, jEq));
          acc->setValue(1, 2, iEq, jEq, halfoEminusAlpha*twothird*m_fluxJacob[2](iEq, jEq));

        }

        if (!isUpdated[stateID2]){
          acc->setValue(2, 0, iEq, jEq, halfoEminusAlpha*onesixth*m_fluxJacob[0](iEq, jEq));
          acc->setValue(2, 1, iEq, jEq, halfoEminusAlpha*onesixth*m_fluxJacob[1](iEq, jEq));
          acc->setValue(2, 2, iEq, jEq,            m_alpha*twothird*m_fluxJacob[2](iEq, jEq));
        }


          }
        }
    // add the contributions to the jacobian
    jacobMatrix->addValues(*acc);
    acc->reset();
  }
    // release the face
    geoBuilder->releaseGE();
  }
}


//////////////////////////////////////////////////////////////////////////////
void UnsteadyWeakSlipWallEuler2DConsHOImpl::docomputeBC(){
  DataHandle<bool> isUpdated = socket_isUpdated.getDataHandle();
  DataHandle<CFreal> rhs = socket_rhs.getDataHandle();
  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();
  DataHandle<State*> interStates = socket_interStates.getDataHandle();
  DataHandle<State*> pastStates = socket_pastStates.getDataHandle();

  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  const CFuint dim = PhysicalModelStack::getActive()->getDim();
  RealVector normal(dim);
  CFreal halfoEminusAlpha = (1.0-m_alpha)/2.0;
  CFreal onesixth = (1.0/6.0);
  CFreal twothird = (2.0/3.0);
  const CFreal dt1 = SubSystemStatusStack::getActive()->getDT();
  const CFreal dt0 =SubSystemStatusStack::getActive()->getPreviousDT();
  const CFreal q = dt1/dt0;
  const CFreal Jpast  = -dt1*(q*q)/(6.0*(1.0 + q)); 
  const CFreal Jinter  = dt1*(3.0 + q)/6.0;
  const CFreal Jpresent  = dt1*(3.0 + 2.0*q)/(6.0*(1.0 + q)); 

  SafePtr<LinearSystemSolver> lss =
    getMethodData().getLinearSystemSolver()[0];

  SafePtr<LSSMatrix> jacobMatrix = lss->getMatrix();

  // block accumulator 2*2
  BlockAccumulator* acc = getMethodData().getLinearSystemSolver()[0]->
    createBlockAccumulator(3, 3, nbEqs);

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
  

    // set the row - column
    acc->setRowColIndex(0, currFace->getState(0)->getLocalID());
    acc->setRowColIndex(1, currFace->getState(1)->getLocalID());
    acc->setRowColIndex(2, currFace->getState(2)->getLocalID());

    // compute the normal fluxes corrections for both the states
    // of this cell
    computeNormalFluxAndJacob((*states[stateID0]), normal, m_fluxn[0], m_fluxJacobn[0]);
    computeNormalFluxAndJacob((*states[stateID1]), normal, m_fluxn[1], m_fluxJacobn[1]);
    computeNormalFluxAndJacob((*states[stateID2]), normal, m_fluxn[2], m_fluxJacobn[2]);
    computeNormalFluxAndJacob((*interStates[stateID0]), normal, m_fluxn1[0], m_fluxJacobn1[0]);
    computeNormalFluxAndJacob((*interStates[stateID1]), normal, m_fluxn1[1], m_fluxJacobn1[1]);
    computeNormalFluxAndJacob((*interStates[stateID2]), normal, m_fluxn1[2], m_fluxJacobn1[2]);
    computeNormalFluxAndJacob((*pastStates[stateID0]), normal, m_fluxn12[0], m_fluxJacobn12[0]);
    computeNormalFluxAndJacob((*pastStates[stateID1]), normal, m_fluxn12[1], m_fluxJacobn12[1]);
    computeNormalFluxAndJacob((*pastStates[stateID2]), normal, m_fluxn12[2], m_fluxJacobn12[2]);

    m_flux[0] = (Jpast*m_fluxn12[0]+Jinter*m_fluxn1[0]+Jpresent*m_fluxn[0]);
    m_flux[1] = (Jpast*m_fluxn12[1]+Jinter*m_fluxn1[1]+Jpresent*m_fluxn[1]);
    m_flux[2] = (Jpast*m_fluxn12[2]+Jinter*m_fluxn1[2]+Jpresent*m_fluxn[2]);

    m_fluxJacob[0] = (Jinter*m_fluxJacobn1[0]+Jpresent*m_fluxJacobn[0]+Jpast*m_fluxJacobn12[0]);
    m_fluxJacob[1] = (Jinter*m_fluxJacobn1[1]+Jpresent*m_fluxJacobn[1]+Jpast*m_fluxJacobn12[1]);
    m_fluxJacob[2] = (Jinter*m_fluxJacobn1[2]+Jpresent*m_fluxJacobn[2]+Jpast*m_fluxJacobn12[2]);
    

 // distribute contributions to the three nodes
    for (CFuint iEq = 0; iEq < nbEqs; ++iEq) {
      if (!isUpdated[stateID0]){
        rhs(stateID0, iEq, nbEqs) -= m_alpha*onesixth*m_flux[0][iEq] + halfoEminusAlpha*onesixth*m_flux[1][iEq] + halfoEminusAlpha*twothird*m_flux[2][iEq];

        _flagState[stateID0] = true;
      }

      if (!isUpdated[stateID1]){
        rhs(stateID1, iEq, nbEqs) -= m_alpha*onesixth*m_flux[1][iEq] + halfoEminusAlpha*onesixth*m_flux[0][iEq] + halfoEminusAlpha*twothird*m_flux[2][iEq];

        _flagState[stateID1] = true;
      }

      if (!isUpdated[stateID2]){
        rhs(stateID2, iEq, nbEqs) -= m_alpha*twothird*m_flux[2][iEq] + halfoEminusAlpha*onesixth*m_flux[1][iEq] + halfoEminusAlpha*onesixth*m_flux[0][iEq];

        _flagState[stateID2] = true;
      }

    }

    // compute the jacobian contributions for the involved states

    for (CFuint iEq = 0; iEq < nbEqs; ++iEq) {
      for (CFuint jEq = 0; jEq < nbEqs; ++jEq) {
        if (!isUpdated[stateID0]){
          acc->setValue(0, 0, iEq, jEq,            m_alpha*onesixth*m_fluxJacob[0](iEq, jEq));
          acc->setValue(0, 1, iEq, jEq, halfoEminusAlpha*onesixth*m_fluxJacob[1](iEq, jEq));
          acc->setValue(0, 2, iEq, jEq, halfoEminusAlpha*twothird*m_fluxJacob[2](iEq, jEq));
        }

        if (!isUpdated[stateID1]){
          acc->setValue(1, 0, iEq, jEq, halfoEminusAlpha*onesixth*m_fluxJacob[0](iEq, jEq));
          acc->setValue(1, 1, iEq, jEq,            m_alpha*onesixth*m_fluxJacob[1](iEq, jEq));
          acc->setValue(1, 2, iEq, jEq, halfoEminusAlpha*twothird*m_fluxJacob[2](iEq, jEq));

        }

        if (!isUpdated[stateID2]){
          acc->setValue(2, 0, iEq, jEq, halfoEminusAlpha*onesixth*m_fluxJacob[0](iEq, jEq));
          acc->setValue(2, 1, iEq, jEq, halfoEminusAlpha*onesixth*m_fluxJacob[1](iEq, jEq));
          acc->setValue(2, 2, iEq, jEq,            m_alpha*twothird*m_fluxJacob[2](iEq, jEq));
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

void UnsteadyWeakSlipWallEuler2DConsHOImpl::computeNormalFluxAndJacob(
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

  // Compute Flux for non-moving wall
  flux[0] = rho*un;
  flux[1] = un*rhoU;
  flux[2] = un*rhoV;
  flux[3] = un*rhoH;


  // the analytical jacobian of the normal fluxes
  fluxJacob(0,1) = nx;
  fluxJacob(0,2) = ny;
  fluxJacob(1,0) = -u*un;
  fluxJacob(1,1) = (un + u*nx);
  fluxJacob(1,2) = u*ny;
  fluxJacob(2,0) = -v*un;
  fluxJacob(2,1) = v*nx;
  fluxJacob(2,2) = (un + v*ny);
  fluxJacob(3,0) = un*(halfGammaMinus1*vel2 - rhoH/rho);
  fluxJacob(3,1) = (rhoH/rho*nx - un*u*gammaMinus1);
  fluxJacob(3,2) = (rhoH/rho*ny - un*v*gammaMinus1);
  fluxJacob(3,3) = (gamma*un);

 }

//////////////////////////////////////////////////////////////////////////////

void UnsteadyWeakSlipWallEuler2DConsHOImpl::setFaceNormal
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

  
    normal[0] = -normals[cellLocalID]->getFaceNormComp(iFaceLocal, 0);
    normal[1] = -normals[cellLocalID]->getFaceNormComp(iFaceLocal, 1);
 }

//////////////////////////////////////////////////////////////////////////////

  } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
