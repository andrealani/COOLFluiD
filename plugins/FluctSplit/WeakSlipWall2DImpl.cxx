#include "WeakSlipWall2DImpl.hh"
#include "InwardNormalsData.hh"
#include "Framework/PhysicalModel.hh"
#include "Framework/BlockAccumulator.hh"
#include "Framework/LSSMatrix.hh"
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

void WeakSlipWall2DImpl::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOption< CFreal >("alpha","distribution coefficient");
}

//////////////////////////////////////////////////////////////////////////////

WeakSlipWall2DImpl::WeakSlipWall2DImpl(const std::string& name) :
  FluctuationSplitCom(name),
  socket_rhs("rhs"),
  socket_states("states"),
  socket_nodes("nodes"),
  socket_isUpdated("isUpdated"),
  socket_normals("normals"),
  socket_faceNeighCell("faceNeighCell"),
  _physicalData(),
  _im0(),
  _in0(),
  _im1(),
  _in1(),
  _flagState(),
  _tJacob()
{
   addConfigOptionsTo(this);
  m_alpha = 1.0;
   setParameter("alpha",&m_alpha);
}

//////////////////////////////////////////////////////////////////////////////

WeakSlipWall2DImpl::~WeakSlipWall2DImpl()
{
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
WeakSlipWall2DImpl::needsSockets()
{

  std::vector<Common::SafePtr<BaseDataSocketSink> > result;

  result.push_back(&socket_rhs);
  result.push_back(&socket_states);
  result.push_back(&socket_isUpdated);
  result.push_back(&socket_normals);
  result.push_back(&socket_faceNeighCell);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

void WeakSlipWall2DImpl::setup()
{

  _im0.resize(PhysicalModelStack::getActive()->getNbEq());
  _in0.resize(PhysicalModelStack::getActive()->getNbEq());
  _im1.resize(PhysicalModelStack::getActive()->getNbEq());
  _in1.resize(PhysicalModelStack::getActive()->getNbEq());
  _tJacob.resize(PhysicalModelStack::getActive()->getNbEq(), PhysicalModelStack::getActive()->getNbEq());
  _tJacob = 0.;

  // get the data handle for the states
  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();
  _flagState.resize(states.size());
  _flagState = false;

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

  m_state.resize(max_nb_states);
//   const bool isGhost = true;
  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
   for ( CFuint iState = 0; iState < max_nb_states; ++iState){
      m_flux[iState].resize(nbEqs);
      m_fluxJacob[iState].resize(nbEqs, nbEqs);
    }
 }

//////////////////////////////////////////////////////////////////////////////

void WeakSlipWall2DImpl::executeOnTrs()
{

  DataHandle<bool> isUpdated = socket_isUpdated.getDataHandle();
  DataHandle<CFreal> rhs = socket_rhs.getDataHandle();

  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  const CFuint dim = PhysicalModelStack::getActive()->getDim();


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


     auto_ptr<BlockAccumulator> acc(lss->createBlockAccumulator(nb_state, nb_state, nbEqs));


 for ( CFuint iState = 0; iState < nb_state; ++iState){
      m_state[iState] = currFace->getState(iState);
      acc->setRowColIndex(iState, m_state[iState]->getLocalID());
    }


    // compute the normal fluxes corrections for both the states
    // of this cell
    for (CFuint iState = 0; iState < nb_state; ++iState)
        computeNormalFluxAndJacob((*m_state[iState]), normal, m_flux[iState], m_fluxJacob[iState]);

    // distribute contributions to the two nodes

    for (CFuint iEq = 0; iEq < nbEqs; ++iEq) {
      for (CFuint iState = 0; iState < nb_state; ++iState){
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

void WeakSlipWall2DImpl::setFaceNormal
(const CFuint faceID, RealVector& normal)
{
  // get the data handle for the states
  DataHandle< InwardNormalsData*> normals = socket_normals.getDataHandle();
  DataHandle<Common::Trio<CFuint, CFuint, Common::SafePtr<Framework::TopologicalRegionSet> > >
    faceNeighCell = socket_faceNeighCell.getDataHandle();

  const CFuint cellTrsID = faceNeighCell[faceID].first;
  const CFuint iFaceLocal = faceNeighCell[faceID].second;
  Common::SafePtr<TopologicalRegionSet> cellTrs = faceNeighCell[faceID].third;
  const CFuint cellLocalID = cellTrs->getLocalGeoID(cellTrsID);

  normal[XX] = -normals[cellLocalID]->getFaceNormComp(iFaceLocal, 0);
  normal[YY] = -normals[cellLocalID]->getFaceNormComp(iFaceLocal, 1);
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace FluctSplit

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
