#include "FluctSplit/FluctSplitNavierStokes.hh"
#include "FluctSplit/HONavierStokes/WeakSlipWallEulerHO2DIsoP2.hh"
#include "FluctSplit/InwardNormalsData.hh"
#include "Framework/CFL.hh"
#include "Framework/PhysicalModel.hh"
#include "Framework/MethodCommandProvider.hh"
#include "NavierStokes/Euler2DVarSet.hh"
#include "Framework/MeshData.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Physics::NavierStokes;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<WeakSlipWallEulerHO2DIsoP2,
		      FluctuationSplitData,
		      FluctSplitNavierStokesModule>
weakSlipWallEulerho2DIsoP2Provider("WeakSlipWallEulerHO2DIsoP2");

//////////////////////////////////////////////////////////////////////////////

void WeakSlipWallEulerHO2DIsoP2::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOption< CFreal >("alpha","distribution coefficient");
}

//////////////////////////////////////////////////////////////////////////////

WeakSlipWallEulerHO2DIsoP2::WeakSlipWallEulerHO2DIsoP2(const std::string& name) :
  FluctuationSplitCom(name),
  socket_rhs("rhs"),
  socket_states("states"),
  socket_isUpdated("isUpdated"),
  socket_normals("normals"),
  socket_faceNeighCell("faceNeighCell"),
  _flagState(),
  _varSet(CFNULL),
  _physicalData(),
  m_weights(3), 
  m_qpPos(3)

{

   addConfigOptionsTo(this);
  _alpha = 1.0;
   setParameter("alpha",&_alpha);
}

//////////////////////////////////////////////////////////////////////////////

WeakSlipWallEulerHO2DIsoP2::~WeakSlipWallEulerHO2DIsoP2()
{
  delete m_qdState;
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
WeakSlipWallEulerHO2DIsoP2::needsSockets()
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

void WeakSlipWallEulerHO2DIsoP2::setup()
{

  _varSet = getMethodData().getUpdateVar().d_castTo<Euler2DVarSet>();

  // get the data handle for the states
  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();
  _flagState.resize(states.size());
  _flagState = false;

  // set up the physical data
  _varSet->getModel()->resizePhysicalData(_physicalData);

 CFuint max_nb_states = 0;
 Common::SafePtr<GeometricEntityPool<StdTrsGeoBuilder> >
  geoBuilder = getMethodData().getStdTrsGeoBuilder();

  StdTrsGeoBuilder::GeoData& geoData = geoBuilder->getDataGE();

  std::vector< Common::SafePtr<TopologicalRegionSet> >& trss = getTrsList();

  std::vector< SafePtr<TopologicalRegionSet> >::iterator itrs = trss.begin();


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

    // vector of fluxes of each state
    m_flux.resize(max_nb_states);
    const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();

    for ( CFuint iState = 0; iState < max_nb_states; ++iState)
      m_flux[iState].resize(nbEqs);



    // release the face
    geoBuilder->releaseGE();
    }
  }

  // vector of fluxes of each state
   m_flux.resize(max_nb_states);
   const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
   for ( CFuint iState = 0; iState < max_nb_states; ++iState)
      m_flux[iState].resize(nbEqs);

  // vector of states of the face
    m_states.resize(max_nb_states);

      //Weights at Gauss points
      m_weights[0] = 5./36.;
      m_weights[1] = 8./36.;
      m_weights[2] = 5./36.;

      //Position of quadrature points 
      //Integration domain <-0.25;0.25> considered

      const CFreal s = std::sqrt(0.6);

      m_qpPos[0] = -0.25*s;
      m_qpPos[1] =  0.0;
      m_qpPos[2] =  0.25*s;

      m_qdState = new State();

}

//////////////////////////////////////////////////////////////////////////////

void WeakSlipWallEulerHO2DIsoP2::executeOnTrs()
{
  CFAUTOTRACE;
  DataHandle<bool> isUpdated = socket_isUpdated.getDataHandle();
  DataHandle<CFreal> rhs = socket_rhs.getDataHandle();

// // dump the rhs to a file
//   ofstream file_in ("weakslip_cp2p2.in");
//   for (CFuint i = 0; i < rhs.size() ; ++i) file_in << rhs[i] << "\n";
//   file_in.close();
//   for (CFuint i = 0; i < rhs.size() ; ++i) rhs[i] = 0.0;




  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  const CFuint dim = PhysicalModelStack::getActive()->getDim();


  m_totFlux0.resize(nbEqs);
  m_totFlux1.resize(nbEqs);
  m_totFlux.resize(nbEqs);
  m_qdFlux.resize(nbEqs);
  m_qdNormal.resize(dim);



  const CFreal oEminusAlpha = 1. - _alpha;

 Common::SafePtr<GeometricEntityPool<StdTrsGeoBuilder> >
    geoBuilder = getMethodData().getStdTrsGeoBuilder();

  StdTrsGeoBuilder::GeoData& geoData = geoBuilder->getDataGE();
  geoData.trs = getCurrentTRS();

  const CFuint nbFaces = getCurrentTRS()->getLocalNbGeoEnts();


  DataHandle< InwardNormalsData*> normals = socket_normals.getDataHandle();
  DataHandle<Common::Trio<CFuint, CFuint, Common::SafePtr<Framework::TopologicalRegionSet> > >
  faceNeighCell = socket_faceNeighCell.getDataHandle();


  for (CFuint iFace = 0; iFace < nbFaces; ++iFace) {

    geoData.idx = iFace;
    GeometricEntity& face = *geoBuilder->buildGE();
//     vector<Node*>& fnodes = *face.getNodes();

    // set number of states on one face
//     const CFuint nb_node = face.nbNodes();
//     const CFuint nb_state = face.nbStates();

//     CFout << "Nr of nodes = " << nb_node << ", nr of states = " << nb_state << "\n";

    const CFuint faceID = face.getID();

//     const CFuint cellTrsID = faceNeighCell[faceID].first;
//     const CFuint iFaceLocal = faceNeighCell[faceID].second;


    Common::SafePtr<TopologicalRegionSet> cellTrs = faceNeighCell[faceID].third;
//     const CFuint cellLocalID = cellTrs->getLocalGeoID(cellTrsID);



  m_nodes = *face.getNodes();
  m_states = *face.getStates();



   //###############################################################################################
   //                  INTEGRATE THE FLUX ACROSS THE HO FACE
   //###############################################################################################


    for(CFuint iEq = 0; iEq < nbEqs; ++iEq) m_totFlux0[iEq] = 0.0;
    for(CFuint iEq = 0; iEq < nbEqs; ++iEq) m_totFlux1[iEq] = 0.0;


    for(CFuint iQd = 0; iQd < 3; ++iQd) {
      const CFreal xiQd = m_qpPos[iQd] + 0.25;
      m_CP2N.ComputeBNormal(m_nodes,xiQd,m_qdNormal);

      m_qdNormal = -1.0 * m_qdNormal;

      (*m_qdState) = (1-3*xiQd+2*xiQd*xiQd)*(*m_states[0]) + (2*xiQd*xiQd-xiQd)*(*m_states[1]) + 4*xiQd*(1-xiQd)*(*m_states[2]);

      const CFreal jacob = m_qdNormal.norm2();
      const CFreal invJacob = 1.0/jacob;
      computeNormalFlux((*m_qdState), invJacob * m_qdNormal, m_qdFlux);
      m_totFlux0 = m_totFlux0 + jacob * m_weights[iQd] * m_qdFlux;
    }



    for(CFuint iQd = 0; iQd < 3; ++iQd) {
      const CFreal xiQd = m_qpPos[iQd] + 0.75;
      m_CP2N.ComputeBNormal(m_nodes,xiQd,m_qdNormal);

      m_qdNormal = -1.0 * m_qdNormal;

      (*m_qdState) = (1-3*xiQd+2*xiQd*xiQd)*(*m_states[0]) + (2*xiQd*xiQd-xiQd)*(*m_states[1]) + 4*xiQd*(1-xiQd)*(*m_states[2]);

      const CFreal jacob = m_qdNormal.norm2();
      const CFreal invJacob = 1.0/jacob;
      computeNormalFlux((*m_qdState), invJacob * m_qdNormal, m_qdFlux);
      m_totFlux1 = m_totFlux1 + jacob * m_weights[iQd] * m_qdFlux;
    }


    m_totFlux = m_totFlux0 + m_totFlux1;


//#################################################################################################
//                             THE ACTUAL BC - MARTIN
//#################################################################################################




//    for (CFuint iEq = 0; iEq < nbEqs; ++iEq) {
// 
//       if (!isUpdated[m_states[0]->getLocalID()]){
//         rhs(m_states[0]->getLocalID(), iEq, nbEqs) -= 1./6.0 * m_totFlux[iEq];
//         _flagState[m_states[0]->getLocalID()] = true;
//       }
// 
// 
//       if (!isUpdated[m_states[1]->getLocalID()]){
//         rhs(m_states[1]->getLocalID(), iEq, nbEqs) -=  1./6.0 * m_totFlux[iEq];
// 
//         _flagState[m_states[1]->getLocalID()] = true;
//       }
// 
// 
//       if (!isUpdated[m_states[2]->getLocalID()]){
//         rhs(m_states[2]->getLocalID(), iEq, nbEqs) -= 2./3.0 * m_totFlux[iEq];
// 
//         _flagState[m_states[2]->getLocalID()] = true;
//       }
// 
//     }


//#################################################################################################
//                                 //MARTIN
//#################################################################################################


//#################################################################################################
//                                 MARTIN - EDWIN LIKE MODIFICATION
//#################################################################################################


//     for (CFuint iEq = 0; iEq < nbEqs; ++iEq) {
// 
//       if (!isUpdated[m_states[0]->getLocalID()]){
//          rhs(m_states[0]->getLocalID(), iEq, nbEqs) -= _alpha*m_totFlux0[iEq] + oEminusAlpha*m_totFlux0[iEq]; 
//         _flagState[m_states[0]->getLocalID()] = true;
//       }
// 
//       if (!isUpdated[m_states[2]->getLocalID()]){
//          rhs(m_states[2]->getLocalID(), iEq, nbEqs) -= oEminusAlpha*m_totFlux0[iEq] + _alpha*m_totFlux0[iEq];
//          rhs(m_states[2]->getLocalID(), iEq, nbEqs) -= _alpha*m_totFlux1[iEq] + oEminusAlpha*m_totFlux1[iEq];
//         _flagState[m_states[2]->getLocalID()] = true;
//       }
// 
//       if (!isUpdated[m_states[1]->getLocalID()]){
//          rhs(m_states[1]->getLocalID(), iEq, nbEqs) -= oEminusAlpha*m_totFlux1[iEq] + _alpha*m_totFlux1[iEq]; 
//         _flagState[m_states[1]->getLocalID()] = true;
//       }
// 
//     }

//#################################################################################################
//                                 //MARTIN - EDWIN LIKE MODIFICATION
//#################################################################################################



//#################################################################################################
//                                 EDWIN's BC
//#################################################################################################


      //The normal is scaled by the length of the total face (2 subfaces)

      // Normal at local node 0 (first end node)
//       m_CP2N.ComputeBNormal(m_nodes,0.0,m_qdNormal);
//       //Jacobian at local node 0
//       const CFreal J0 = m_qdNormal.norm2();
//       m_qdNormal = -1.0/J0 * m_qdNormal;
// 
//       computeNormalFlux(*m_states[0], m_qdNormal, m_flux[0]);
// 
//       ///Normal at local node 1 (second end node)
//       m_CP2N.ComputeBNormal(m_nodes,1.0,m_qdNormal);
//       //Jacobian at local node 1
//       const CFreal J1 = m_qdNormal.norm2();
//       m_qdNormal = -1.0/J1 * m_qdNormal;
// 
//       computeNormalFlux(*m_states[1], m_qdNormal, m_flux[1]);
// 
//       ///Normal at local node 2 (mid node)
//       m_CP2N.ComputeBNormal(m_nodes,0.5,m_qdNormal);
//       //Jacobian at local node 2
//       const CFreal J2 = m_qdNormal.norm2();
//       m_qdNormal = -1.0/J2 * m_qdNormal;
// 
//       computeNormalFlux(*m_states[2], m_qdNormal, m_flux[2]);
// 
// 
// 
// 
//    for (CFuint iEq = 0; iEq < nbEqs; ++iEq) {
// 
//       if (!isUpdated[m_states[0]->getLocalID()]){
//          rhs(m_states[0]->getLocalID(), iEq, nbEqs) -= _alpha*0.5*J0*m_flux[0][iEq] + oEminusAlpha*0.5*J2*m_flux[2][iEq]; 
//         _flagState[m_states[0]->getLocalID()] = true;
//       }
// 
//       if (!isUpdated[m_states[2]->getLocalID()]){
//          rhs(m_states[2]->getLocalID(), iEq, nbEqs) -= oEminusAlpha*0.5*J0*m_flux[0][iEq] + _alpha*0.5*J2*m_flux[2][iEq];
//          rhs(m_states[2]->getLocalID(), iEq, nbEqs) -= _alpha*0.5*J2*m_flux[2][iEq] + oEminusAlpha*0.5*J1*m_flux[1][iEq];
//         _flagState[m_states[2]->getLocalID()] = true;
//       }
// 
//       if (!isUpdated[m_states[1]->getLocalID()]){
//          rhs(m_states[1]->getLocalID(), iEq, nbEqs) -= oEminusAlpha*0.5*J2*m_flux[2][iEq] + _alpha*0.5*J1*m_flux[1][iEq]; 
//         _flagState[m_states[1]->getLocalID()] = true;
//       }
// 
// 
//     }


//#################################################################################################
//                                 // EDWIN's BC
//#################################################################################################



//#################################################################################################
//                                 NADEGE
//#################################################################################################


      ///The normals are rescaled to unit length!

      ///Normal at local node 0 (first end node)
      m_CP2N.ComputeBNormal(m_nodes,0.0,m_qdNormal);
      //Jacobian at local node 0
      const CFreal J0 = m_qdNormal.norm2();
      m_qdNormal = -1.0/J0 * m_qdNormal;

      computeNormalFlux(*m_states[0], m_qdNormal, m_flux[0]);

      ///Normal at local node 1 (second end node)
      m_CP2N.ComputeBNormal(m_nodes,1.0,m_qdNormal);
      //Jacobian at local node 1
      const CFreal J1 = m_qdNormal.norm2();
      m_qdNormal = -1.0/J1 * m_qdNormal;

      computeNormalFlux(*m_states[1], m_qdNormal, m_flux[1]);

      ///Normal at local node 2 (mid node)
      m_CP2N.ComputeBNormal(m_nodes,0.5,m_qdNormal);
      //Jacobian at local node 2
      const CFreal J2 = m_qdNormal.norm2();
      m_qdNormal = -1.0/J2 * m_qdNormal;

      computeNormalFlux(*m_states[2], m_qdNormal, m_flux[2]);




    /// distribute contributions to the THREE nodes
    for (CFuint iEq = 0; iEq < nbEqs; ++iEq) {

      if (!isUpdated[m_states[0]->getLocalID()]){
        rhs(m_states[0]->getLocalID(), iEq, nbEqs) -= _alpha*(1.0/6.0)*J0*m_flux[0][iEq] + (oEminusAlpha/2.0)*(1.0/6.0)*J1*m_flux[1][iEq] + (oEminusAlpha/2.0)*(2.0/3.0)*J2*m_flux[2][iEq];

        _flagState[m_states[0]->getLocalID()] = true;
      }


      if (!isUpdated[m_states[1]->getLocalID()]){
        rhs(m_states[1]->getLocalID(), iEq, nbEqs) -=  _alpha*(1.0/6.0)*J1*m_flux[1][iEq] + (oEminusAlpha/2.0)*(1.0/6.0)*J0*m_flux[0][iEq] + (oEminusAlpha/2.0)*(2.0/3.0)*J2*m_flux[2][iEq];

        _flagState[m_states[1]->getLocalID()] = true;
      }


      if (!isUpdated[m_states[2]->getLocalID()]){
        rhs(m_states[2]->getLocalID(), iEq, nbEqs) -= _alpha*(2.0/3.0)*J2*m_flux[2][iEq] + (oEminusAlpha/2.0)*(1.0/6.0)*J1*m_flux[1][iEq] + (oEminusAlpha/2.0)*(1.0/6.0)*J0*m_flux[0][iEq];

        _flagState[m_states[2]->getLocalID()] = true;
      }

    }



//#################################################################################################
//                                 // NADEGE
//#################################################################################################




   geoBuilder->releaseGE();
  }

  ///@todo this could be done in a more efficient way ...
  const CFuint nbStates = _flagState.size();
  for (CFuint i = 0; i < nbStates; ++i) {
    if (_flagState[i] == true) {
      isUpdated[i] = true;
    }
  }


// dump the rhs to a file
// ofstream file_out ("weakslip_cp2p2.out");
// for (CFuint i = 0; i < rhs.size() ; ++i) file_out << rhs[i] << "\n";
// file_out.close();
// CF_DEBUG_EXIT;


}

//////////////////////////////////////////////////////////////////////////////

void WeakSlipWallEulerHO2DIsoP2::computeNormalFlux(const State& state,
                                            const RealVector& normal,
                                            RealVector& flux)
{
  State& ss = *(const_cast<State*>(&state));
  _varSet->computePhysicalData(ss, _physicalData);

  const CFreal un = _physicalData[EulerTerm::VX]*normal[XX] +  _physicalData[EulerTerm::VY]*normal[YY];

  const CFreal rho = _physicalData[EulerTerm::RHO];

  //if (!getMethodData().isResidualTransformationNeeded()) {
  flux[0] = rho*un;
  flux[1] = un*rho*_physicalData[EulerTerm::VX];
  flux[2] = un*rho*_physicalData[EulerTerm::VY];
  flux[3] = un*rho*_physicalData[EulerTerm::H];
  //}
  //   else {
  //     _tmpFlux[0] = rho*un;
  //     _tmpFlux[1] = un*rho*_physicalData[EulerTerm::VX];
  //     _tmpFlux[2] = un*rho*_physicalData[EulerTerm::VY];
  //     _tmpFlux[3] = un*rho*_physicalData[EulerTerm::H];

  //     // set the transformation from update to solution in update
  //     SafePtr<VarSetMatrixTransformer> solToUpdateInUpdate =
  //       getMethodData().getSolToUpdateInUpdateMatTrans();

  //     solToUpdateInUpdate->setMatrix(state);
  //     const RealMatrix& tMatrix = *solToUpdateInUpdate->getMatrix();
  //     flux = tMatrix*_tmpFlux;
  //   }
}

//////////////////////////////////////////////////////////////////////////////

void WeakSlipWallEulerHO2DIsoP2::setFaceNormal(const CFuint faceID,
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

//   CFout << "Cell local ID = " << cellLocalID+1 << "\n";
//   CFout << "\tFace ID = " << iFaceLocal << "\n";

  normal[XX] = -normals[cellLocalID]->getFaceNormComp(iFaceLocal, XX);
  normal[YY] = -normals[cellLocalID]->getFaceNormComp(iFaceLocal, YY);
}

//////////////////////////////////////////////////////////////////////////////

} // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
