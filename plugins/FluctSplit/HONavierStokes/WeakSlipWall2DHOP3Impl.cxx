#include "FluctSplit/HONavierStokes/WeakSlipWall2DHOP3Impl.hh"
#include "FluctSplit/InwardNormalsData.hh"
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

void WeakSlipWall2DHOP3Impl::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOption< CFreal >("alpha","distribution coefficient");
   options.addConfigOption< bool >("ExactNormNaca0012","Use the exact normals for a NACA0012");
}

//////////////////////////////////////////////////////////////////////////////

WeakSlipWall2DHOP3Impl::WeakSlipWall2DHOP3Impl(const std::string& name) :
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

  m_exact_norm = false;
  setParameter("ExactNormNaca0012",&m_exact_norm);

}

//////////////////////////////////////////////////////////////////////////////

WeakSlipWall2DHOP3Impl::~WeakSlipWall2DHOP3Impl()
{
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
WeakSlipWall2DHOP3Impl::needsSockets()
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

void WeakSlipWall2DHOP3Impl::setup()
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

// vector of nodes of the face 
  m_nodes.resize(max_nb_states);
  m_state.resize(max_nb_states);
//   const bool isGhost = true;
  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
   for ( CFuint iState = 0; iState < max_nb_states; ++iState){
      m_flux[iState].resize(nbEqs);
      m_fluxJacob[iState].resize(nbEqs, nbEqs);
    }
 }

//////////////////////////////////////////////////////////////////////////////

void WeakSlipWall2DHOP3Impl::executeOnTrs()
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

  /// vector of nodes of the face
  std::vector<Framework::Node*> m_nodes;
  CFreal thirdoEminusAlpha =  (1.0-m_alpha)/3.0;
  CFreal oneight = (1.0/8.0);
  CFreal threeeight= (3.0/8.0);

  const CFuint nbFaces = getCurrentTRS()->getLocalNbGeoEnts();
  for (CFuint iFace = 0; iFace < nbFaces; ++iFace) {
    geoData.idx = iFace;
    GeometricEntity *const currFace = geoBuilder->buildGE();
    const CFuint nb_state = currFace->nbStates();

    // set the face normal
    setFaceNormal(currFace->getID(), normal);


     auto_ptr<BlockAccumulator> acc(lss->createBlockAccumulator(nb_state, nb_state, nbEqs));


 for ( CFuint iState = 0; iState < nb_state; ++iState){
      m_state[iState] = currFace->getState(iState);
      acc->setRowColIndex(iState, m_state[iState]->getLocalID());
    }


    const CFuint dim = PhysicalModelStack::getActive()->getDim();
    m_qdNormal.resize(dim);
    // compute the normal fluxes corrections for both the states
    // of this cell
    for (CFuint iState = 0; iState < nb_state; ++iState){
      // If m_exact_norm is true then we use the true normal of the NACA0012 to compute
      // the correction flux that we want at the wall
      if (m_exact_norm){
      Node& node0 = (*m_state[iState]).getCoordinates();

      // Coordinate of the vertex
      const CFreal x0 = node0[XX];
      const CFreal y0 = node0[YY];

      if ((y0 >= 0.0) && (x0 != 0.0) && (x0 < 1.0)){
        m_qdNormal[XX] = -0.6*(( 0.2969/(2.0*sqrt(x0))) - 0.1260 - 0.3516*2.0*x0 + 3.0*0.2843*x0*x0 - 0.1015*4.0*x0*x0*x0 );
        m_qdNormal[YY] = 1.0;
      }
      else if ((y0 <= 0.0) && (x0 != 0.0) && (x0 < 1.0)){
        m_qdNormal[XX] = -0.6*(( 0.2969/(2.0*sqrt(x0))) - 0.1260 - 0.3516*2.0*x0 + 3.0*0.2843*x0*x0 - 0.1015*4.0*x0*x0*x0 );
        m_qdNormal[YY] = -1.0;
      }
      else if ((x0 == 0.0)){
        m_qdNormal[XX] = -1.0;
        m_qdNormal[YY] = 0.0;
      }
      else {m_qdNormal = normal;}
      m_qdNormal /= m_qdNormal.norm2();
      m_qdNormal *= normal.norm2();
      }
      
      else {m_qdNormal = normal; }
      computeNormalFluxAndJacob((*m_state[iState]), m_qdNormal, m_flux[iState], m_fluxJacob[iState]);
  }

    // distribute contributions to the two nodes
    for (CFuint iEq = 0; iEq < nbEqs; ++iEq) {
      if (!isUpdated[m_state[0]->getLocalID()]){
        rhs(m_state[0]->getLocalID(), iEq, nbEqs) -= m_alpha*oneight*m_flux[0][iEq] + thirdoEminusAlpha*(oneight*m_flux[1][iEq] + threeeight*m_flux[2][iEq] + threeeight*m_flux[3][iEq]);

      }

      if (!isUpdated[m_state[1]->getLocalID()]){
        rhs(m_state[1]->getLocalID(), iEq, nbEqs) -= m_alpha*oneight*m_flux[1][iEq] + thirdoEminusAlpha*(oneight*m_flux[0][iEq] + threeeight*m_flux[2][iEq] + threeeight*m_flux[3][iEq]);

      }

      if (!isUpdated[m_state[2]->getLocalID()]){
        rhs(m_state[2]->getLocalID(), iEq, nbEqs) -= m_alpha*threeeight*m_flux[2][iEq] + thirdoEminusAlpha*(oneight*m_flux[1][iEq] + oneight*m_flux[0][iEq] + threeeight*m_flux[3][iEq]);

      }

  if (!isUpdated[m_state[3]->getLocalID()]){
        rhs(m_state[3]->getLocalID(), iEq, nbEqs) -= m_alpha*threeeight*m_flux[3][iEq] + thirdoEminusAlpha*(oneight*m_flux[1][iEq] + oneight*m_flux[0][iEq] + threeeight*m_flux[2][iEq]);

      }
   }
    // compute the jacobian contributions for the involved states

    for (CFuint iEq = 0; iEq < nbEqs; ++iEq) {
      for (CFuint jEq = 0; jEq < nbEqs; ++jEq) {
        if (!isUpdated[m_state[0]->getLocalID()]){
          acc->setValue(0, 0, iEq, jEq,            m_alpha*oneight*m_fluxJacob[0](iEq, jEq));
          acc->setValue(0, 1, iEq, jEq, thirdoEminusAlpha*oneight*m_fluxJacob[1](iEq, jEq));
          acc->setValue(0, 2, iEq, jEq, thirdoEminusAlpha*threeeight*m_fluxJacob[2](iEq, jEq));
          acc->setValue(0, 3, iEq, jEq, thirdoEminusAlpha*threeeight*m_fluxJacob[3](iEq, jEq));
        }

        if (!isUpdated[m_state[1]->getLocalID()]){
          acc->setValue(1, 0, iEq, jEq, thirdoEminusAlpha*oneight*m_fluxJacob[0](iEq, jEq));
          acc->setValue(1, 1, iEq, jEq,            m_alpha*oneight*m_fluxJacob[1](iEq, jEq));
          acc->setValue(1, 2, iEq, jEq, thirdoEminusAlpha*threeeight*m_fluxJacob[2](iEq, jEq));
          acc->setValue(1, 3, iEq, jEq, thirdoEminusAlpha*threeeight*m_fluxJacob[3](iEq, jEq));
        }

        if (!isUpdated[m_state[2]->getLocalID()]){
          acc->setValue(2, 0, iEq, jEq, thirdoEminusAlpha*oneight*m_fluxJacob[0](iEq, jEq));
          acc->setValue(2, 1, iEq, jEq, thirdoEminusAlpha*oneight*m_fluxJacob[1](iEq, jEq));
          acc->setValue(2, 2, iEq, jEq,            m_alpha*threeeight*m_fluxJacob[2](iEq, jEq));
          acc->setValue(2, 3, iEq, jEq, thirdoEminusAlpha*threeeight*m_fluxJacob[3](iEq, jEq));
        }

      if (!isUpdated[m_state[3]->getLocalID()]){
          acc->setValue(3, 0, iEq, jEq, thirdoEminusAlpha*oneight*m_fluxJacob[0](iEq, jEq));
          acc->setValue(3, 1, iEq, jEq, thirdoEminusAlpha*oneight*m_fluxJacob[1](iEq, jEq));
          acc->setValue(3, 2, iEq, jEq, thirdoEminusAlpha*threeeight*m_fluxJacob[2](iEq, jEq));
          acc->setValue(3, 3, iEq, jEq,            m_alpha*threeeight*m_fluxJacob[3](iEq, jEq));
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

void WeakSlipWall2DHOP3Impl::setFaceNormal
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
