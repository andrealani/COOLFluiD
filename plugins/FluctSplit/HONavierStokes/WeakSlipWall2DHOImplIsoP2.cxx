#include "FluctSplit/HONavierStokes/WeakSlipWall2DHOImplIsoP2.hh"
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

void WeakSlipWall2DHOImplIsoP2::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOption< CFreal >("alpha","distribution coefficient");
}

//////////////////////////////////////////////////////////////////////////////

WeakSlipWall2DHOImplIsoP2::WeakSlipWall2DHOImplIsoP2(const std::string& name) :
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

WeakSlipWall2DHOImplIsoP2::~WeakSlipWall2DHOImplIsoP2()
{
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
WeakSlipWall2DHOImplIsoP2::needsSockets()
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

void WeakSlipWall2DHOImplIsoP2::setup()
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

  // vector of states of the face
    m_states.resize(max_nb_states);
// vector of nodes of the face 
  m_nodes.resize(max_nb_states);
  // allocation of the flux and Jacobian flux
  m_flux.resize(max_nb_states);
  m_fluxJacob.resize(max_nb_states);

//   const bool isGhost = true;
  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
   for ( CFuint iState = 0; iState < max_nb_states; ++iState){
      m_flux[iState].resize(nbEqs);
      m_fluxJacob[iState].resize(nbEqs, nbEqs);
    }


  const CFuint dim = PhysicalModelStack::getActive()->getDim();
  m_qdNormal.resize(dim);
 }

//////////////////////////////////////////////////////////////////////////////

void WeakSlipWall2DHOImplIsoP2::executeOnTrs()
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

  CFreal halfoEminusAlpha = (1.0-m_alpha)/2.0;
  CFreal onesixth = (1.0/6.0);
  CFreal twothird = (2.0/3.0);
  const CFuint nbFaces = getCurrentTRS()->getLocalNbGeoEnts();
  for (CFuint iFace = 0; iFace < nbFaces; ++iFace) {
    geoData.idx = iFace;
    GeometricEntity *const currFace = geoBuilder->buildGE();
    const CFuint nb_state = currFace->nbStates();

    // set the face normal
  //  setFaceNormal(currFace->getID(), normal);


     auto_ptr<BlockAccumulator> acc(lss->createBlockAccumulator(nb_state, nb_state, nbEqs));

 for ( CFuint iState = 0; iState < nb_state; ++iState){
      m_states[iState] = currFace->getState(iState);
      m_nodes[iState] = currFace->getNode(iState);
      acc->setRowColIndex(iState, m_states[iState]->getLocalID());
    }


      ///The normals are rescaled to unit length!

      ///Normal at local node 0 (first end node)
      m_CP2N.ComputeBNormal(m_nodes,0.0,m_qdNormal);
      //Jacobian at local node 0
      const CFreal J0 = m_qdNormal.norm2();
      m_qdNormal = -1.0/J0 * m_qdNormal;

      //CFreal x0 = (*m_nodes[0])[XX];
      //CFreal y0 = (*m_nodes[0])[YY];

      //if ((y0 >= 0.0) && (x0 != 0.0))
     // {
      //  m_qdNormal[XX] = -0.6*(( 0.2969/(2.0*sqrt(x0))) - 0.1260 - 0.3516*2.0*x0 + 3.0*0.2843*x0*x0 - 0.1015*4.0*x0*x0*x0 );
      //  m_qdNormal[YY] = 1.0;
     // }
    //else if ((y0 <= 0.0) && (x0 != 0.0))
    //  {
     //   m_qdNormal[XX] = -0.6*(( 0.2969/(2.0*sqrt(x0))) - 0.1260 - 0.3516*2.0*x0 + 3.0*0.2843*x0*x0 - 0.1015*4.0*x0*x0*x0 );
//        m_qdNormal[YY] = -1.0;
  //    }
   // else if ((x0 == 0.0))
     // {
       // m_qdNormal[XX] = -1.0;
       // m_qdNormal[YY] = 0.0;
     // }


     // m_qdNormal /= m_qdNormal.norm2();


      computeNormalFluxAndJacob(*m_states[0], m_qdNormal, m_flux[0], m_fluxJacob[0]);

      ///Normal at local node 1 (second end node)
      m_CP2N.ComputeBNormal(m_nodes,1.0,m_qdNormal);
      //Jacobian at local node 1
      const CFreal J1 = m_qdNormal.norm2();
      m_qdNormal = -1.0/J1 * m_qdNormal;

     // CFreal x1 = (*m_nodes[1])[XX];
//      CFreal y1 = (*m_nodes[1])[YY];

  //    if ((y1 >= 0.0) && (x1 != 0.0))
    //  {
      //  m_qdNormal[XX] = -0.6*(( 0.2969/(2.0*sqrt(x1))) - 0.1260 - 0.3516*2.0*x1 + 3.0*0.2843*x1*x1 - 0.1015*4.0*x1*x1*x1 );
        //m_qdNormal[YY] = 1.0;
     // }
    //else if ((y1 <= 0.0) && (x1 != 0.0))
      //{
        //m_qdNormal[XX] = -0.6*(( 0.2969/(2.0*sqrt(x1))) - 0.1260 - 0.3516*2.0*x1 + 3.0*0.2843*x1*x1 - 0.1015*4.0*x1*x1*x1 );
       // m_qdNormal[YY] = -1.0;
     // }
    //else if ((x1 == 0.0))
     // {
       // m_qdNormal[XX] = -1.0;
       // m_qdNormal[YY] = 0.0;
     // }

     // m_qdNormal /= m_qdNormal.norm2();

      computeNormalFluxAndJacob(*m_states[1], m_qdNormal, m_flux[1], m_fluxJacob[1]);

      ///Normal at local node 2 (mid node)
      m_CP2N.ComputeBNormal(m_nodes,0.5,m_qdNormal);
      //Jacobian at local node 2
      const CFreal J2 = m_qdNormal.norm2();
      m_qdNormal = -1.0/J2 * m_qdNormal;
 
      //CFreal x2 = (*m_nodes[2])[XX];
      //CFreal y2 = (*m_nodes[2])[YY];

      //if ((y2 >= 0.0) && (x2 != 0.0))
     // {
       // m_qdNormal[XX] = -0.6*(( 0.2969/(2.0*sqrt(x2))) - 0.1260 - 0.3516*2.0*x2 + 3.0*0.2843*x2*x2 - 0.1015*4.0*x2*x2*x2 );
       // m_qdNormal[YY] = 1.0;
     // }
    //else if ((y2 <= 0.0) && (x2 != 0.0))
     // {
      //  m_qdNormal[XX] = -0.6*(( 0.2969/(2.0*sqrt(x2))) - 0.1260 - 0.3516*2.0*x2 + 3.0*0.2843*x2*x2 - 0.1015*4.0*x2*x2*x2 );
       // m_qdNormal[YY] = -1.0;
     // }
    //else if ((x2 == 0.0))
     // {
       // m_qdNormal[XX] = -1.0;
       // m_qdNormal[YY] = 0.0;
     // }

//      m_qdNormal /= m_qdNormal.norm2();

      computeNormalFluxAndJacob(*m_states[2], m_qdNormal, m_flux[2], m_fluxJacob[2]);

    // distribute contributions to the three nodes
    for (CFuint iEq = 0; iEq < nbEqs; ++iEq) {
      if (!isUpdated[m_states[0]->getLocalID()]){
        rhs(m_states[0]->getLocalID(), iEq, nbEqs) -= m_alpha*onesixth*J0*m_flux[0][iEq] + halfoEminusAlpha*onesixth*J1*m_flux[1][iEq] + halfoEminusAlpha*twothird*J2*m_flux[2][iEq];

        _flagState[m_states[0]->getLocalID()] = true;
      }

      if (!isUpdated[m_states[1]->getLocalID()]){
        rhs(m_states[1]->getLocalID(), iEq, nbEqs) -= m_alpha*J1*onesixth*m_flux[1][iEq] + halfoEminusAlpha*onesixth*J0*m_flux[0][iEq] + halfoEminusAlpha*twothird*J2*m_flux[2][iEq];

        _flagState[m_states[1]->getLocalID()] = true;
      }

      if (!isUpdated[m_states[2]->getLocalID()]){
        rhs(m_states[2]->getLocalID(), iEq, nbEqs) -= m_alpha*twothird*J2*m_flux[2][iEq] + halfoEminusAlpha*onesixth*J1*m_flux[1][iEq] + halfoEminusAlpha*onesixth*J0*m_flux[0][iEq];

        _flagState[m_states[2]->getLocalID()] = true;
      }

    }

    // compute the jacobian contributions for the involved states

    for (CFuint iEq = 0; iEq < nbEqs; ++iEq) {
      for (CFuint jEq = 0; jEq < nbEqs; ++jEq) {
        if (!isUpdated[m_states[0]->getLocalID()]){
          acc->setValue(0, 0, iEq, jEq,            m_alpha*onesixth*J0*m_fluxJacob[0](iEq, jEq));
          acc->setValue(0, 1, iEq, jEq, halfoEminusAlpha*onesixth*J1*m_fluxJacob[1](iEq, jEq));
          acc->setValue(0, 2, iEq, jEq, halfoEminusAlpha*twothird*J2*m_fluxJacob[2](iEq, jEq));
        }

        if (!isUpdated[m_states[1]->getLocalID()]){
          acc->setValue(1, 0, iEq, jEq, halfoEminusAlpha*onesixth*J0*m_fluxJacob[0](iEq, jEq));
          acc->setValue(1, 1, iEq, jEq,            m_alpha*onesixth*J1*m_fluxJacob[1](iEq, jEq));
          acc->setValue(1, 2, iEq, jEq, halfoEminusAlpha*twothird*J2*m_fluxJacob[2](iEq, jEq));

        }

        if (!isUpdated[m_states[2]->getLocalID()]){
          acc->setValue(2, 0, iEq, jEq, halfoEminusAlpha*onesixth*J0*m_fluxJacob[0](iEq, jEq));
          acc->setValue(2, 1, iEq, jEq, halfoEminusAlpha*onesixth*J1*m_fluxJacob[1](iEq, jEq));
          acc->setValue(2, 2, iEq, jEq,            m_alpha*twothird*J2*m_fluxJacob[2](iEq, jEq));
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

void WeakSlipWall2DHOImplIsoP2::setFaceNormal
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
