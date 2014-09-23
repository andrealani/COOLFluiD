#include "FluctSplit/FluctSplitNavierStokes.hh"
#include "FluctSplit/HONavierStokes/WeakSlipWallEulerHO2DIsoP3.hh"
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

MethodCommandProvider<WeakSlipWallEulerHO2DIsoP3,
		      FluctuationSplitData,
		      FluctSplitNavierStokesModule>
weakSlipWallEulerho2DIsoP3Provider("WeakSlipWallEulerHO2DIsoP3");

//////////////////////////////////////////////////////////////////////////////

void WeakSlipWallEulerHO2DIsoP3::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOption< CFreal >("alpha","distribution coefficient");
}

//////////////////////////////////////////////////////////////////////////////

WeakSlipWallEulerHO2DIsoP3::WeakSlipWallEulerHO2DIsoP3(const std::string& name) :
  FluctuationSplitCom(name),
  socket_rhs("rhs"),
  socket_states("states"),
  socket_isUpdated("isUpdated"),
  socket_normals("normals"),
  socket_faceNeighCell("faceNeighCell"),
  _flagState(),
  _varSet(CFNULL),
  _physicalData(),
  m_weights(nbQdPt), 
  m_qpPos(nbQdPt)

{

   addConfigOptionsTo(this);
  _alpha = 1.0;
   setParameter("alpha",&_alpha);
}

//////////////////////////////////////////////////////////////////////////////

WeakSlipWallEulerHO2DIsoP3::~WeakSlipWallEulerHO2DIsoP3()
{
  delete m_qdState;
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
WeakSlipWallEulerHO2DIsoP3::needsSockets()
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

void WeakSlipWallEulerHO2DIsoP3::setup()
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
    m_nodeFlux.resize(max_nb_states);
    const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();

    for ( CFuint iState = 0; iState < max_nb_states; ++iState)
      m_nodeFlux[iState].resize(nbEqs);



    // release the face
    geoBuilder->releaseGE();
    }
  }

  // vector of fluxes of each state
   m_nodeFlux.resize(max_nb_states);
   const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
   for ( CFuint iState = 0; iState < max_nb_states; ++iState)
      m_nodeFlux[iState].resize(nbEqs);

  // vector of states of the face
    m_face_states.resize(max_nb_states);

      //Weights at Gauss points
      m_weights[0] = 5./54.;
      m_weights[1] = 8./54.;
      m_weights[2] = 5./54.;

  //Position of quadrature points in reference space < -1.0/6.0 : 1.0/6.0 >

      const CFreal s = std::sqrt(0.6);

      m_qpPos[0] = -1.0/6.0*s;
      m_qpPos[1] =  0.0;
      m_qpPos[2] =  1.0/6.0*s;

      m_qdState = new State();

}

//////////////////////////////////////////////////////////////////////////////

void WeakSlipWallEulerHO2DIsoP3::executeOnTrs()
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


  m_subFaceFlux.resize(3); ///P3 face has 3 subfaces
  m_subFaceFlux[0].resize(nbEqs);
  m_subFaceFlux[1].resize(nbEqs);
  m_subFaceFlux[2].resize(nbEqs);



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



  m_face_nodes = *face.getNodes();
  m_face_states = *face.getStates();



   //###############################################################################################
   //                  INTEGRATE THE FLUX ACROSS THE HO FACE
   //###############################################################################################


//     for(CFuint iEq = 0; iEq < nbEqs; ++iEq) m_subFaceFlux[0][iEq] = 0.0;
//     for(CFuint iEq = 0; iEq < nbEqs; ++iEq) m_subFaceFlux[1][iEq] = 0.0;
//     for(CFuint iEq = 0; iEq < nbEqs; ++iEq) m_subFaceFlux[2][iEq] = 0.0;

    m_subFaceFlux[0] = 0.0;
    m_subFaceFlux[1] = 0.0;
    m_subFaceFlux[2] = 0.0;

    ///Integrate flux across first subface:
    for(CFuint iQd = 0; iQd < nbQdPt; ++iQd) {
      const CFreal xiQd = m_qpPos[iQd] + 1.0/6.0; // 1.0/6.0 is the center of the first subface in reference frame <0,1>
      m_CP3N.ComputeBNormal(m_face_nodes,xiQd,m_qdNormal);

      m_qdNormal = -1.0 * m_qdNormal;

      const CFreal L0 = 1.0-xiQd;
      const CFreal L1 = xiQd;

      (*m_qdState) = 0.5*(3.0*L0-1.0)*(3.0*L0-2.0)*L0*(*m_face_states[0]) + 0.5*(3.0*L1-1.0)*(3.0*L1-2.0)*L1*(*m_face_states[1]) + \
                    4.5*L0*L1*(3.0*L0-1.0)*(*m_face_states[2]) + 4.5*L0*L1*(3.0*L1-1.0)*(*m_face_states[3]);

      const CFreal jacob = m_qdNormal.norm2();
      computeNormalFlux((*m_qdState), 1.0/jacob * m_qdNormal, m_qdFlux);
      m_subFaceFlux[0] = m_subFaceFlux[0] + jacob * m_weights[iQd] * m_qdFlux;
    }


    ///Integrate flux across middle subface:
    for(CFuint iQd = 0; iQd < nbQdPt; ++iQd) {
      const CFreal xiQd = m_qpPos[iQd] + 0.5;
      m_CP3N.ComputeBNormal(m_face_nodes,xiQd,m_qdNormal);

      m_qdNormal = -1.0 * m_qdNormal;

      const CFreal L0 = 1.0-xiQd;
      const CFreal L1 = xiQd;

      (*m_qdState) = 0.5*(3.0*L0-1.0)*(3.0*L0-2.0)*L0*(*m_face_states[0]) + 0.5*(3.0*L1-1.0)*(3.0*L1-2.0)*L1*(*m_face_states[1]) + \
                    4.5*L0*L1*(3.0*L0-1.0)*(*m_face_states[2]) + 4.5*L0*L1*(3.0*L1-1.0)*(*m_face_states[3]);

      const CFreal jacob = m_qdNormal.norm2();
      computeNormalFlux((*m_qdState), 1.0/jacob * m_qdNormal, m_qdFlux);
      m_subFaceFlux[1] = m_subFaceFlux[1] + jacob * m_weights[iQd] * m_qdFlux;
    }


    ///Integrate flux across the last subface:
    for(CFuint iQd = 0; iQd < nbQdPt; ++iQd) {
      const CFreal xiQd = m_qpPos[iQd] + 5.0/6.0;
      m_CP3N.ComputeBNormal(m_face_nodes,xiQd,m_qdNormal);

      m_qdNormal = -1.0 * m_qdNormal;

      const CFreal L0 = 1.0-xiQd;
      const CFreal L1 = xiQd;

      (*m_qdState) = 0.5*(3.0*L0-1.0)*(3.0*L0-2.0)*L0*(*m_face_states[0]) + 0.5*(3.0*L1-1.0)*(3.0*L1-2.0)*L1*(*m_face_states[1]) + \
                    4.5*L0*L1*(3.0*L0-1.0)*(*m_face_states[2]) + 4.5*L0*L1*(3.0*L1-1.0)*(*m_face_states[3]);

      const CFreal jacob = m_qdNormal.norm2();
      computeNormalFlux((*m_qdState), 1.0/jacob * m_qdNormal, m_qdFlux);
      m_subFaceFlux[2] = m_subFaceFlux[2] + jacob * m_weights[iQd] * m_qdFlux;
    }


//#################################################################################################
//                                 NADEGE
//#################################################################################################


      ///The normals are rescaled to unit length!

      ///Normal at local node 0 (first end node)
      m_CP3N.ComputeBNormal(m_face_nodes,0.0,m_qdNormal);
      //Jacobian at local node 0
      const CFreal J0 = m_qdNormal.norm2();
      m_qdNormal = -1.0/J0 * m_qdNormal;

      computeNormalFlux(*m_face_states[0], m_qdNormal, m_nodeFlux[0]);

      ///Normal at local node 1 (second end node)
      m_CP3N.ComputeBNormal(m_face_nodes,1.0,m_qdNormal);
      //Jacobian at local node 1
      const CFreal J1 = m_qdNormal.norm2();
      m_qdNormal = -1.0/J1 * m_qdNormal;

      computeNormalFlux(*m_face_states[1], m_qdNormal, m_nodeFlux[1]);

      ///Normal at local node 2 (first inner node)
      m_CP3N.ComputeBNormal(m_face_nodes,1.0/3.0,m_qdNormal);
      //Jacobian at local node 2
      const CFreal J2 = m_qdNormal.norm2();
      m_qdNormal = -1.0/J2 * m_qdNormal;

      computeNormalFlux(*m_face_states[2], m_qdNormal, m_nodeFlux[2]);

      ///Normal at local node 3 (second inner node)
      m_CP3N.ComputeBNormal(m_face_nodes,2.0/3.0,m_qdNormal);
      //Jacobian at local node 3
      const CFreal J3 = m_qdNormal.norm2();
      m_qdNormal = -1.0/J3 * m_qdNormal;

      computeNormalFlux(*m_face_states[3], m_qdNormal, m_nodeFlux[3]);



    /// distribute contributions to the THREE nodes
    for (CFuint iEq = 0; iEq < nbEqs; ++iEq) {

//       if (!isUpdated[m_face_states[0]->getLocalID()]){
        rhs(m_face_states[0]->getLocalID(), iEq, nbEqs) -= _alpha*(1.0/8.0)*J0*m_nodeFlux[0][iEq] + (oEminusAlpha/3.0)*(1.0/8.0)*J1*m_nodeFlux[1][iEq] + (oEminusAlpha/3.0)*(3.0/8.0)*J2*m_nodeFlux[2][iEq]+(oEminusAlpha/3.0)*(3.0/8.0)*J3*m_nodeFlux[3][iEq];

//         _flagState[m_face_states[0]->getLocalID()] = true;
//       }

//       if (!isUpdated[m_face_states[1]->getLocalID()]){
        rhs(m_face_states[1]->getLocalID(), iEq, nbEqs) -=  _alpha*(1.0/8.0)*J1*m_nodeFlux[1][iEq] + (oEminusAlpha/3.0)*(1.0/8.0)*J0*m_nodeFlux[0][iEq] + (oEminusAlpha/3.0)*(3.0/8.0)*J2*m_nodeFlux[2][iEq]+(oEminusAlpha/3.0)*(3.0/8.0)*J3*m_nodeFlux[3][iEq];

//         _flagState[m_face_states[1]->getLocalID()] = true;
//       }

//       if (!isUpdated[m_face_states[2]->getLocalID()]){
        rhs(m_face_states[2]->getLocalID(), iEq, nbEqs) -= _alpha*(3.0/8.0)*J2*m_nodeFlux[2][iEq] + (oEminusAlpha/3.0)*(1.0/8.0)*J0*m_nodeFlux[0][iEq] + (oEminusAlpha/3.0)*(1.0/8.0)*J1*m_nodeFlux[1][iEq]+(oEminusAlpha/3.0)*(1.0/8.0)*J3*m_nodeFlux[3][iEq];

//         _flagState[m_face_states[2]->getLocalID()] = true;
//       }

//       if (!isUpdated[m_face_states[3]->getLocalID()]){
        rhs(m_face_states[3]->getLocalID(), iEq, nbEqs) -= _alpha*(3.0/8.0)*J3*m_nodeFlux[3][iEq] + (oEminusAlpha/3.0)*(1.0/8.0)*J0*m_nodeFlux[0][iEq] + (oEminusAlpha/3.0)*(1.0/8.0)*J1*m_nodeFlux[1][iEq]+(oEminusAlpha/3.0)*(1.0/8.0)*J2*m_nodeFlux[2][iEq];

//         _flagState[m_face_states[3]->getLocalID()] = true;
//       }



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

void WeakSlipWallEulerHO2DIsoP3::computeNormalFlux(const State& state,
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

void WeakSlipWallEulerHO2DIsoP3::setFaceNormal(const CFuint faceID,
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
