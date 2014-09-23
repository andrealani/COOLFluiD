#include "Framework/CFL.hh"
#include "Framework/PhysicalModel.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/MeshData.hh"

#include "NavierStokes/Euler2DVarSet.hh"

#include "FluctSplit/FluctuationSplitStrategy.hh"
#include "FluctSplit/HONavierStokes/FluctSplitHONavierStokes.hh"
#include "FluctSplit/HONavierStokes/WeakSlipWallEulerHO2DIsoP2Upwind.hh"
#include "FluctSplit/InwardNormalsData.hh"

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

MethodCommandProvider<WeakSlipWallEulerHO2DIsoP2Upwind,
		      FluctuationSplitData,
		      FluctSplitHONavierStokesModule>
weakSlipWallEulerHO2DIsoP2UpwindProvider("WeakSlipWallEulerHO2DIsoP2Upwind");

//////////////////////////////////////////////////////////////////////////////

  ///List of faces of the three subtriangles that are adjacent to the wall boundary:
  CFuint WeakSlipWallEulerHO2DIsoP2Upwind::m_subFaceIdx0[6] = {2, 0, 1, 5, 3, 4};
  CFuint WeakSlipWallEulerHO2DIsoP2Upwind::m_subFaceIdx1[6] = {3, 4, 5, 6, 7, 8};
  CFuint WeakSlipWallEulerHO2DIsoP2Upwind::m_subFaceIdx2[6] = {7, 8, 6, 1, 2, 0};

  CFuint WeakSlipWallEulerHO2DIsoP2Upwind::m_subElemIdx0[2] = {0, 1};
  CFuint WeakSlipWallEulerHO2DIsoP2Upwind::m_subElemIdx1[2] = {1, 2};
  CFuint WeakSlipWallEulerHO2DIsoP2Upwind::m_subElemIdx2[2] = {2, 0};

  CFuint WeakSlipWallEulerHO2DIsoP2Upwind::m_subTri0[6] = { 0, 3, 5, 3, 1, 4 };
  CFuint WeakSlipWallEulerHO2DIsoP2Upwind::m_subTri1[6] = { 1, 4, 3, 4, 2, 5 };
  CFuint WeakSlipWallEulerHO2DIsoP2Upwind::m_subTri2[6] = { 2, 5, 4, 5, 0, 3 };
//////////////////////////////////////////////////////////////////////////////

void WeakSlipWallEulerHO2DIsoP2Upwind::defineConfigOptions(Config::OptionList& options)
{
}

//////////////////////////////////////////////////////////////////////////////

WeakSlipWallEulerHO2DIsoP2Upwind::WeakSlipWallEulerHO2DIsoP2Upwind(const std::string& name) :
  FluctuationSplitCom(name),
  socket_rhs("rhs"),
  socket_states("states"),
  socket_isUpdated("isUpdated"),
  socket_normals("normals"),
  socket_faceNeighCell("faceNeighCell"),
  _flagState(),
  _varSet(CFNULL),
  _physicalData(),
  m_weights(nbQdPts),
  m_qpPos(nbQdPts),
  m_subface_center(2),
  matrix_face_norms(3,2),
  matrix_node_norms(3,2),
  vector_face_areas(3),
  vector_node_areas(3),
  m_linearStates(CFNULL)

{

   addConfigOptionsTo(this);

   const CFuint dim = DIM_2D;
   const CFuint nbfaces = 3;
   const CFuint nbnodes = 3;

   /// @todo this is a memory leak, they are not destroyed in the destructor of InwardNormalsData
   CFreal * faceNormals = new CFreal[nbfaces*dim];
   CFreal * faceAreas   = new CFreal[nbfaces];
   CFreal * nodeNormals = new CFreal[nbnodes*dim];
   CFreal * nodeAreas   = new CFreal[nbnodes];
   // place them in the inwardnormals
   m_subcell_normals =
     new InwardNormalsData(faceNormals,faceAreas,nodeNormals,nodeAreas,nbfaces,0);


}

//////////////////////////////////////////////////////////////////////////////

WeakSlipWallEulerHO2DIsoP2Upwind::~WeakSlipWallEulerHO2DIsoP2Upwind()
{
  delete m_qdState;
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
WeakSlipWallEulerHO2DIsoP2Upwind::needsSockets()
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

void WeakSlipWallEulerHO2DIsoP2Upwind::setup()
{
  _varSet = getMethodData().getUpdateVar().d_castTo<Euler2DVarSet>();

  // get the data handle for the states
  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();
  _flagState.resize(states.size());
  _flagState = false;

  // set up the physical data
  _varSet->getModel()->resizePhysicalData(_physicalData);

 CFuint max_nb_nodes = 0;
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

      const CFuint nb_node = currFace.nbNodes();
      const CFuint nb_state = currFace.nbStates();

      max_nb_nodes = std::max(max_nb_nodes, nb_node);
      max_nb_states = std::max(max_nb_states, nb_state);

    // release the face
    geoBuilder->releaseGE();
    }
  }

  // vector of fluxes of each state
   m_flux.resize(max_nb_states);


   const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();

   for ( CFuint iState = 0; iState < max_nb_states; ++iState)
    m_flux[iState].resize(nbEqs);

  m_subStates.resize(3);   // 3 states in each sub element
  m_subresidual.resize(3); // 3 residuals in each sub element

  m_subresidual[0].resize(nbEqs);
  m_subresidual[1].resize(nbEqs);
  m_subresidual[2].resize(nbEqs);

  // vector of nodes of the face
  m_face_nodes.resize(max_nb_nodes);
  // vector of states of the face
  m_face_states.resize(max_nb_states);

  m_subface_center[0] = 0.25;
  m_subface_center[1] = 0.75;

  //Weights at Gauss points
//   m_weights[0] = 0.25;
//   m_weights[1] = 0.25;

  m_weights[0] = 5./36.;
  m_weights[1] = 8./36.;
  m_weights[2] = 5./36.;

//   m_weights[0] = 0.2369268850*0.25;
//   m_weights[1] = 0.4786286705*0.25;
//   m_weights[2] = 0.5688888889*0.25;
//   m_weights[3] = 0.4786286705*0.25;
//   m_weights[4] = 0.2369268850*0.25;


  //Position of quadrature points in reference space < -1/4 : 1/4 >


//     m_qpPos[0] = -0.2;
//     m_qpPos[1] = 0.2;

  const CFreal s = std::sqrt(0.6);
// 
  m_qpPos[0] = -0.25*s;
  m_qpPos[1] =  0.0;
  m_qpPos[2] =  0.25*s;


//   m_qpPos[0] = -0.9061798459*0.25;
//   m_qpPos[1] = -0.5384693101*0.25;
//   m_qpPos[2] = 0.0;
//   m_qpPos[3] = 0.5384693101*0.25;
//   m_qpPos[4] = 0.9061798459*0.25;

  m_qdState = new State();

  ///NEW: splitter setup
  m_splitter    = getMethodData().getSplitter();
  m_solutionVar = getMethodData().getSolutionVar();
  m_updateVar   = getMethodData().getUpdateVar();

  // tell the splitters to compute the betas
  getMethodData().getDistributionData().computeBetas = true;

  m_betasInSubTriag.resize(nbEqs,nbEqs);


//---------------------------------------------------------------------------------------------------

  // allocating data for the temporary local residual
  const CFuint maxNbStatesInCell = MeshDataStack::getActive()->Statistics().getMaxNbStatesInCell();
  m_residual.resize(maxNbStatesInCell);
  for (CFuint i = 0; i < maxNbStatesInCell; ++i)
  {
    m_residual[i].resize(nbEqs);
    m_residual[i] = 0.0;
  }

//---------------------------------------------------------------------------------------------------


}

//////////////////////////////////////////////////////////////////////////////

void WeakSlipWallEulerHO2DIsoP2Upwind::executeOnTrs()
{
  CFAUTOTRACE;
  DataHandle<bool> isUpdated = socket_isUpdated.getDataHandle();
  DataHandle<CFreal> rhs = socket_rhs.getDataHandle();


  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  const CFuint dim = PhysicalModelStack::getActive()->getDim();

  DistributionData& distdata = getMethodData().getDistributionData();

  m_subFaceFlux.resize(nbEqs);

  m_qdFlux.resize(nbEqs);
  m_qdNormal.resize(dim);

 Common::SafePtr<GeometricEntityPool<StdTrsGeoBuilder> >
    geoBuilder = getMethodData().getStdTrsGeoBuilder();

  StdTrsGeoBuilder::GeoData& geoData = geoBuilder->getDataGE();
  geoData.trs = getCurrentTRS();

  const CFuint nbFaces = getCurrentTRS()->getLocalNbGeoEnts();


  DataHandle< InwardNormalsData*> normals = socket_normals.getDataHandle();
  DataHandle<Common::Trio<CFuint, CFuint, Common::SafePtr<Framework::TopologicalRegionSet> > >
  faceNeighCell = socket_faceNeighCell.getDataHandle();

/*
  ///NEW:GET NODAL CONNECTIVITY:
  Common::SafePtr<MeshData::ConnTable> cellNodes = MeshDataStack::getActive()->getConnectivity("cellNodes_InnerCells");

  ///NEW:GET NODAL COORDINATES
//   DataHandle < Framework::Node*, Framework::GLOBAL > nodes = getCFmeshData().getNodesHandle();
  DataHandle < Framework::Node*, Framework::GLOBAL > nodes = MeshDataStack::getActive()->getNodeDataSocketSink().getDataHandle();

  ///NEW:GET STATES CONNECTIVITY:
  Common::SafePtr<MeshData::ConnTable> cellStates = MeshDataStack::getActive()->getConnectivity("cellStates_InnerCells");

  ///NEW:GET STATES
  DataHandle < Framework::State*, Framework::GLOBAL > states = MeshDataStack::getActive()->getStateDataSocketSink().getDataHandle();
*/

//---------------------------------------------------------------------------------------------------

  // create a builder for neighbour cells
  Framework::GeometricEntityPool<Framework::StdTrsGeoBuilder> geoBuilderCell;
  geoBuilderCell.setup();
  StdTrsGeoBuilder::GeoData& geoDataCell = geoBuilderCell.getDataGE();

//   DistributionData& ddata = getMethodData().getDistributionData();

//---------------------------------------------------------------------------------------------------

 for (CFuint iFace = 0; iFace < nbFaces; ++iFace) {

  geoData.idx = iFace;
  GeometricEntity& face = *geoBuilder->buildGE();

  const CFuint faceID = face.getID();

  m_face_nodes = *face.getNodes();
  m_face_states = *face.getStates();


//#################################################################################################
//                                 Distribute the fluxes to boundary nodes
//#################################################################################################

  ///Find out to which cell this face belongs:
  const CFuint cellTrsID = faceNeighCell[faceID].first;
//     const CFuint iFaceLocal = faceNeighCell[faceID].second;
  Common::SafePtr<TopologicalRegionSet> cellTrs = faceNeighCell[faceID].third;
  const CFuint cellLocalID = cellTrs->getLocalGeoID(cellTrsID);

//---------------------------------------------------------------------------------------------------

  // build the cell
  geoDataCell.trs = faceNeighCell[faceID].third;
  geoDataCell.idx = cellLocalID;
  GeometricEntity& cell = *geoBuilderCell.buildGE();
  
//   vector<Node*> *const cell_nodes = cell.getNodes(); ///We actually do not need the nodes
  vector<State*> *const cell_states = cell.getStates();

  distdata.cell   = &cell;
  distdata.cellID = cell.getID();
  distdata.states = cell_states;

//   getMethodData().getFluctSplitStrategy()->computeFluctuation(m_residual);
//---------------------------------------------------------------------------------------------------

  const CFuint faceIdx0 = m_face_nodes[0]->getGlobalID();
  const CFuint faceIdx1 = m_face_nodes[1]->getGlobalID();

  const CFuint cellIdx0 = (*cell_states)[0]->getGlobalID();
  const CFuint cellIdx1 = (*cell_states)[1]->getGlobalID();
  const CFuint cellIdx2 = (*cell_states)[2]->getGlobalID();


  //Remap the local coordinates so that the local the local indices 0 3 4 1 of P3P3 triangle are
  //adjacent to wall boundary. This is just shuffling of indexes which are predefined in m_localTriIDx[0|1|2]
  m_subFaceIdx = m_subFaceIdx0;
  m_subElemIdx = m_subElemIdx0;
  m_subTri = m_subTri0;

  if(faceIdx0 == cellIdx1 && faceIdx1 == cellIdx2) {
    m_subFaceIdx = m_subFaceIdx1;
    m_subElemIdx = m_subElemIdx1;
    m_subTri = m_subTri1;
  }
  else if(faceIdx0 == cellIdx2 && faceIdx1 == cellIdx0) {
    m_subFaceIdx = m_subFaceIdx2;
    m_subElemIdx = m_subElemIdx2;
    m_subTri = m_subTri2;
  }
   else if (!(faceIdx0 == cellIdx0 && faceIdx1 == cellIdx1)) {
  CFout << "Boundary cell: " << cellIdx0 << "," << cellIdx1 << "," << cellIdx2 << "\n";
  CFout << "Boundary face: " << faceIdx0 << "," << faceIdx1 << "," << m_face_nodes[1]->getGlobalID() << "\n";

     CF_DEBUG_POINT;
     CF_DEBUG_EXIT;

   }


  InwardNormalsData& cellnormals = (*socket_normals.getDataHandle() [cellLocalID]);

  CFuint ibdrysubface = 0;

///Loop over subtriangles adjacent to wall boundary
for(CFuint isubelem = 0; isubelem < 2; ++isubelem, ++ibdrysubface) { 

    m_subFaceFlux = 0.0;

    //INTEGRATE THE FLUX ACROSS THE HO FACE
    //Compute the flux through boundary face - this is the 'corrective residual'
    for(CFuint iQd = 0; iQd < nbQdPts; ++iQd) {
      const CFreal xiQd = m_qpPos[iQd] + m_subface_center[ibdrysubface];
      m_CP2N.ComputeBNormal(m_face_nodes,xiQd,m_qdNormal);

      m_qdNormal = -1.0 * m_qdNormal;

      const CFreal p0 = (1.0-3.0*xiQd+2.0*xiQd*xiQd);
      const CFreal p1 = xiQd*(2.0*xiQd-1.0);
      const CFreal p2 = 4.0*xiQd*(1.0-xiQd);

      (*m_qdState) = p0*(*m_face_states[0]) + p1*(*m_face_states[1]) + p2*(*m_face_states[2]);

      const CFreal jacob = m_qdNormal.norm2();
      const CFreal invJacob = 1.0/jacob;
      computeNormalFlux((*m_qdState), invJacob * m_qdNormal, m_qdFlux);

      m_subFaceFlux = m_subFaceFlux + jacob * m_weights[iQd] * m_qdFlux;
    }


  ///set indexes of subtriangle: vertices ...
  const CFuint v0 = m_subTri[3*isubelem];
  const CFuint v1 = m_subTri[3*isubelem+1];
  const CFuint v2 = m_subTri[3*isubelem+2];

  m_subStates[0] = (*cell_states)[v0];
  m_subStates[1] = (*cell_states)[v1];
  m_subStates[2] = (*cell_states)[v2];

  /// ... and faces

  const CFuint f0 = m_subFaceIdx[3*isubelem];
  const CFuint f1 = m_subFaceIdx[3*isubelem+1];
  const CFuint f2 = m_subFaceIdx[3*isubelem+2];

  // stupidily put the normals data into a matrix to put it in the normals data
  matrix_face_norms (0,XX) = cellnormals.getFaceNormComp(f0,XX);
  matrix_face_norms (1,XX) = cellnormals.getFaceNormComp(f1,XX);
  matrix_face_norms (2,XX) = cellnormals.getFaceNormComp(f2,XX);

  matrix_face_norms (0,YY) = cellnormals.getFaceNormComp(f0,YY);
  matrix_face_norms (1,YY) = cellnormals.getFaceNormComp(f1,YY);
  matrix_face_norms (2,YY) = cellnormals.getFaceNormComp(f2,YY);

  vector_face_areas[0] = cellnormals.getAreaFace(f0);
  vector_face_areas[1] = cellnormals.getAreaFace(f1);
  vector_face_areas[2] = cellnormals.getAreaFace(f2);

  // this is always the same because triangle face normals are opposite to nodal normals
  matrix_node_norms (0,XX) = - matrix_face_norms (1,XX);
  matrix_node_norms (1,XX) = - matrix_face_norms (2,XX);
  matrix_node_norms (2,XX) = - matrix_face_norms (0,XX);

  matrix_node_norms (0,YY) = - matrix_face_norms (1,YY);
  matrix_node_norms (1,YY) = - matrix_face_norms (2,YY);
  matrix_node_norms (2,YY) = - matrix_face_norms (0,YY);

  vector_node_areas[0] = vector_face_areas[1];
  vector_node_areas[1] = vector_face_areas[2];
  vector_node_areas[2] = vector_face_areas[0];

  // set them in the temporary InwardNormalsData to pass to splitter
  m_subcell_normals->setFaceNormals ( matrix_face_norms );
  m_subcell_normals->setNodalNormals( matrix_node_norms );

  m_subcell_normals->setFaceAreas ( vector_face_areas );
  m_subcell_normals->setNodalAreas( vector_node_areas );

  //pass the substates corresponding to vertices of subtriangle to distribution data 
//   distdata.subStates = &m_subStates;

  // compute the residual and the upwind parameters k in this cell
  distdata.tStates = computeConsistentStates(&m_subStates);
  distdata.subStates = &m_subStates;

  m_splitter->computeK(m_subStates, m_subcell_normals);

  // transform fluxes of subelement to distribution variables
  SafePtr<RealVector> phi = &distdata.phi;
  *phi = *getMethodData().getSolutionToDistribMatTrans()->transformFromRef(&m_subFaceFlux);

  SafePtr<vector<RealMatrix> >& currBetaMatrix = distdata.currBetaMat;
  currBetaMatrix = &distdata.betaMats[m_subElemIdx[isubelem]];
  cf_assert (currBetaMatrix.isNotNull());


  //Splitter computes the beta matrices
  m_splitter->distribute(m_subresidual);

//   CFout << "Residuals computed by splitter:\n";
//   CFout << m_subresidual[0] << "\n";
//   CFout << m_subresidual[1] << "\n";
//   CFout << m_subresidual[2] << "\n";
// 
// 
//   m_subresidual[0] = 0.0;
//   m_subresidual[1] = 0.0;
//   m_subresidual[2] = 0.0;
// 
//   //Let's try to distribute manually
//   m_subresidual[0] = ((*currBetaMatrix)[0]) * m_subFaceFlux;
//   m_subresidual[1] = ((*currBetaMatrix)[1]) * m_subFaceFlux;
//   m_subresidual[2] = ((*currBetaMatrix)[2]) * m_subFaceFlux;
// 
//   CFout << "Residuals computed by multiplication:\n";
//   CFout << m_subresidual[0] << "\n";
//   CFout << m_subresidual[1] << "\n";
//   CFout << m_subresidual[2] << "\n";
// 
//    std::cin.get();


    //Update node 0 of subtriangle.
    for (CFuint iEq = 0; iEq < nbEqs; ++iEq) {
      if (!isUpdated[(*cell_states)[v0]->getLocalID()]){
//         _flagState[(*cell_states)[v0]->getLocalID()] = true;
        rhs((*cell_states)[v0]->getLocalID(), iEq, nbEqs) -= m_subresidual[0][iEq];
      }
    }

    for (CFuint iEq = 0; iEq < nbEqs; ++iEq) {
      if (!isUpdated[(*cell_states)[v1]->getLocalID()]){
//         _flagState[(*cell_states)[v1]->getLocalID()] = true;
        rhs((*cell_states)[v1]->getLocalID(), iEq, nbEqs) -= m_subresidual[1][iEq];

      }
    }

    for (CFuint iEq = 0; iEq < nbEqs; ++iEq) {
      if (!isUpdated[(*cell_states)[v2]->getLocalID()]){
        rhs((*cell_states)[v2]->getLocalID(), iEq, nbEqs) -= m_subresidual[2][iEq];
//         _flagState[(*cell_states)[v2]->getLocalID()] = true;
      }
    }

} ///Loop over subtriangles adjacent to wall boundary

  geoBuilderCell.releaseGE();
  geoBuilder->releaseGE();

  }



  ///@todo this could be done in a more efficient way ...
//   const CFuint nbStates = _flagState.size();
//   for (CFuint i = 0; i < nbStates; ++i) {
//     if (_flagState[i] == true) {
//       isUpdated[i] = true;
//     }
//   }


// dump the rhs to a file
// ofstream file_out ("weakslip_cp2p2.out");
// for (CFuint i = 0; i < rhs.size() ; ++i) file_out << rhs[i] << "\n";
// file_out.close();
// CF_DEBUG_EXIT;


}

//////////////////////////////////////////////////////////////////////////////

void WeakSlipWallEulerHO2DIsoP2Upwind::computeNormalFlux(const State& state,
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


} // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
