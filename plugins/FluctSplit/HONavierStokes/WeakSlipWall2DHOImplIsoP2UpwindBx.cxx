#include "FluctSplit/HONavierStokes/WeakSlipWall2DHOImplIsoP2UpwindBx.hh"
#include "FluctSplit/InwardNormalsData.hh"
#include "Framework/PhysicalModel.hh"
#include "Framework/BlockAccumulator.hh"
#include "Framework/LSSMatrix.hh"
#include "Framework/MeshData.hh"

#include "FluctSplit/HONavierStokes/FluctSplitHONavierStokes.hh"
#include "MathTools/MatrixInverter.hh"

#include "NavierStokes/EulerTerm.hh"
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

  ///List of faces of the three subtriangles that are adjacent to the wall boundary:
  CFuint WeakSlipWall2DHOImplIsoP2UpwindBx::m_subFaceIdx0[6] = {2, 0, 1, 5, 3, 4};
  CFuint WeakSlipWall2DHOImplIsoP2UpwindBx::m_subFaceIdx1[6] = {3, 4, 5, 6, 7, 8};
  CFuint WeakSlipWall2DHOImplIsoP2UpwindBx::m_subFaceIdx2[6] = {7, 8, 6, 1, 2, 0};

  CFuint WeakSlipWall2DHOImplIsoP2UpwindBx::m_subElemIdx0[2] = {0, 1};
  CFuint WeakSlipWall2DHOImplIsoP2UpwindBx::m_subElemIdx1[2] = {1, 2};
  CFuint WeakSlipWall2DHOImplIsoP2UpwindBx::m_subElemIdx2[2] = {2, 0};

  CFuint WeakSlipWall2DHOImplIsoP2UpwindBx::m_subTri0[6] = { 0, 3, 5, 3, 1, 4 };
  CFuint WeakSlipWall2DHOImplIsoP2UpwindBx::m_subTri1[6] = { 1, 4, 3, 4, 2, 5 };
  CFuint WeakSlipWall2DHOImplIsoP2UpwindBx::m_subTri2[6] = { 2, 5, 4, 5, 0, 3 };

//////////////////////////////////////////////////////////////////////////////

void WeakSlipWall2DHOImplIsoP2UpwindBx::defineConfigOptions(Config::OptionList& options)
{
}

//////////////////////////////////////////////////////////////////////////////

WeakSlipWall2DHOImplIsoP2UpwindBx::WeakSlipWall2DHOImplIsoP2UpwindBx(const std::string& name) :
  FluctuationSplitCom(name),
  socket_rhs("rhs"),
  socket_states("states"),
  socket_nodes("nodes"),
  socket_isUpdated("isUpdated"),
  socket_normals("normals"),
  socket_updateCoeff("updateCoeff"),
  socket_faceNeighCell("faceNeighCell"),
  _physicalData(),
  m_subFaceJacob(),
  _im0(),
  _in0(),
  _im1(),
  _in1(),
  _flagState(),
  _tJacob(),
  m_subface_center(2),
  m_weights(nbQdPts),
  m_qpPos(nbQdPts),
  matrix_face_norms(3,2),
  matrix_node_norms(3,2),
  vector_face_areas(3),
  vector_node_areas(3),
  m_linearStates(CFNULL),
  m_phiN(0),
  m_phiLDA(0),
  m_sumKplusU(),
  m_sumKU(),
  m_sumKplus(),
  m_uTemp(),
  m_uMin(),
  m_tmp(),
  m_kPlus(0),
  m_k(0),
  m_kMin(0),
  m_eValues(0),
  m_theta(),
  m_inverter(CFNULL),
  m_invK(),
  m_uInflow(),
  _pData()
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

WeakSlipWall2DHOImplIsoP2UpwindBx::~WeakSlipWall2DHOImplIsoP2UpwindBx()
{
  deletePtr( m_qdState );
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
WeakSlipWall2DHOImplIsoP2UpwindBx::needsSockets()
{
//   CFout << "Socket allocation\n";

  std::vector<Common::SafePtr<BaseDataSocketSink> > result;

  result.push_back(&socket_rhs);
  result.push_back(&socket_states);
  result.push_back(&socket_isUpdated);
  result.push_back(&socket_normals);
  result.push_back(&socket_faceNeighCell);
  result.push_back(&socket_updateCoeff);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

void WeakSlipWall2DHOImplIsoP2UpwindBx::setup()
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


  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();

  ///Resize the matrices that hold the product beta * flux_jacobian
  m_distJacobFrom0.resize(nbEqs,nbEqs);
  m_distJacobFrom1.resize(nbEqs,nbEqs);
  m_distJacobFrom3.resize(nbEqs,nbEqs);
  m_distJacobToState.resize(nbEqs,nbEqs);

  m_qdJacob.resize(nbEqs,nbEqs);

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
  m_weights[0] = 5./36.;
  m_weights[1] = 8./36.;
  m_weights[2] = 5./36.;

//   m_weights[0] = 0.2369268850*0.25;
//   m_weights[1] = 0.4786286705*0.25;
//   m_weights[2] = 0.5688888889*0.25;
//   m_weights[3] = 0.4786286705*0.25;
//   m_weights[4] = 0.2369268850*0.25;


  //Position of quadrature points in reference space < -1.0/6.0 : 1.0/6.0 >

  const CFreal s = std::sqrt(0.6);

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


//   const CFuint dim = PhysicalModelStack::getActive()->getDim();
//   m_qdNormal.resize(dim);


  m_phiN.resize(3);

  m_phiN[0].resize(nbEqs);
  m_phiN[1].resize(nbEqs);
  m_phiN[2].resize(nbEqs);

  m_phiLDA.resize(3);

  m_phiLDA[0].resize(nbEqs);
  m_phiLDA[1].resize(nbEqs);
  m_phiLDA[2].resize(nbEqs);
 m_maxNbStatesInCell = MeshDataStack::getActive()->Statistics().getMaxNbStatesInCell();

  m_kPlus.resize(m_maxNbStatesInCell);
  m_k.resize(m_maxNbStatesInCell);
  m_kMin.resize(m_maxNbStatesInCell);
  m_eValues.resize(m_maxNbStatesInCell);

  for (CFuint i = 0; i < m_maxNbStatesInCell; ++i) {
    m_kPlus[i] = new RealMatrix(nbEqs, nbEqs);
    m_kMin[i] = new RealMatrix(nbEqs, nbEqs);
    m_k[i] = new RealMatrix(nbEqs, nbEqs);
    m_eValues[i] = new RealVector(nbEqs);
  }
  m_adimNormal.resize(2);

     m_sumKplusU.resize(nbEqs);
  m_sumKU.resize(nbEqs);
  m_sumKplus.resize(nbEqs,nbEqs);
  m_uTemp.resize(nbEqs);
  m_uMin.resize(nbEqs);
  m_tmp.resize(nbEqs,nbEqs);
  m_invK.resize(nbEqs, nbEqs);
  m_uInflow.resize(nbEqs);

  m_inverter = MatrixInverter::create(nbEqs, false);
  _cterm = Framework::PhysicalModelStack::getActive()->getImplementor()->
    getConvectiveTerm().d_castTo<Physics::NavierStokes::EulerTerm>();

  _cterm->resizePhysicalData(_pData);

  _grad.resize(Framework::PhysicalModelStack::getActive()->getDim());
// set _choiceVar to pressure if default value is required
  
    _varID = Physics::NavierStokes::EulerTerm::P;
 

  // resize the beta's storages in the MethodData
  const CFuint nbSubTriangles = 4;
  const CFuint nbNodesInTriangle = 3;

  DistributionData::BetaMatrices& betaMats =
    getMethodData().getDistributionData().betaMats;
  betaMats.resize(nbSubTriangles);

  for (CFuint t = 0; t < nbSubTriangles; ++t) {
    betaMats[t].resize(nbNodesInTriangle);

    for (CFuint i = 0; i < nbNodesInTriangle; ++i) {
      betaMats[t][i].resize(nbEqs, nbEqs);
    }
  }

  // tell the splitters to compute the betas
  getMethodData().getDistributionData().computeBetas = true;

 }

//////////////////////////////////////////////////////////////////////////////

void WeakSlipWall2DHOImplIsoP2UpwindBx::executeOnTrs()
{
  CFAUTOTRACE;


//   CF_DEBUG_POINT;
  DataHandle<bool> isUpdated = socket_isUpdated.getDataHandle();
  DataHandle<CFreal> rhs = socket_rhs.getDataHandle();

  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  const CFuint dim = PhysicalModelStack::getActive()->getDim();

  DistributionData& distdata = getMethodData().getDistributionData();

  m_subFaceFlux.resize(nbEqs);

  m_subFaceJacob.resize(nbQdPts);

  for(CFuint iQd = 0; iQd < nbQdPts; ++iQd)
    m_subFaceJacob[iQd].resize(nbEqs,nbEqs);


  m_qdFlux.resize(nbEqs);
  m_qdNormal.resize(dim);

  m_p0.resize(nbQdPts);
  m_p1.resize(nbQdPts);
  m_p3.resize(nbQdPts);


//   RealVector normal(dim);


  SafePtr<LinearSystemSolver> lss =
    getMethodData().getLinearSystemSolver()[0];

  SafePtr<LSSMatrix> jacobMatrix = lss->getMatrix();

  // block accumulator 2*2


  Common::SafePtr<GeometricEntityPool<StdTrsGeoBuilder> >
    geoBuilder = getMethodData().getStdTrsGeoBuilder();

  StdTrsGeoBuilder::GeoData& geoData = geoBuilder->getDataGE();
  geoData.trs = getCurrentTRS();

  const CFuint nbFaces = getCurrentTRS()->getLocalNbGeoEnts();


  DataHandle< InwardNormalsData*> normals = socket_normals.getDataHandle();
  DataHandle<Common::Trio<CFuint, CFuint, Common::SafePtr<Framework::TopologicalRegionSet> > >
  faceNeighCell = socket_faceNeighCell.getDataHandle();

//---------------------------------------------------------------------------------------------------

  // create a builder for neighbour cells
  Framework::GeometricEntityPool<Framework::StdTrsGeoBuilder> geoBuilderCell;
  geoBuilderCell.setup();
  StdTrsGeoBuilder::GeoData& geoDataCell = geoBuilderCell.getDataGE();

//   DistributionData& ddata = getMethodData().getDistributionData();

//---------------------------------------------------------------------------------------------------


  for (CFuint iFace = 0; iFace < nbFaces; ++iFace) {

    geoData.idx = iFace;
    GeometricEntity& currFace = *geoBuilder->buildGE();

    const CFuint faceID = currFace.getID();

    m_face_nodes  = *currFace.getNodes();
    m_face_states = *currFace.getStates();

    const CFuint nb_state = currFace.nbStates();

    // set the face normal
  //  setFaceNormal(currFace->getID(), normal);

//  auto_ptr<BlockAccumulator> acc(lss->createBlockAccumulator(nb_state, nb_state, nbEqs));

 for ( CFuint iState = 0; iState < nb_state; ++iState){
      m_face_states[iState] = currFace.getState(iState);
      m_face_nodes[iState] = currFace.getNode(iState);
//       acc->setRowColIndex(iState, m_face_states[iState]->getLocalID());
    }




//#################################################################################################
//                                 Distribute the fluxes to boundary nodes
//#################################################################################################
///*COMMENT*


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

  //std::cout<<"CellID"<<distdata.cellID<<std::endl;
  const CFuint nb_cell_states  = cell_states->size();
  auto_ptr<BlockAccumulator> acc(lss->createBlockAccumulator(nb_cell_states, nb_cell_states, nbEqs));

  for ( CFuint istate = 0; istate <  nb_cell_states; ++istate){
      acc->setRowColIndex(istate, (*cell_states)[istate]->getLocalID());
  }

  //Initialize the values in the accumulator:
  acc->setValue(0.0);




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

   }


  InwardNormalsData& cellnormals = (*socket_normals.getDataHandle() [cellLocalID]);

  CFuint ibdrysubface = 0;

///Loop over subtriangles adjacent to wall boundary
for(CFuint isubelem = 0; isubelem < 2; ++isubelem, ++ibdrysubface) { 

   //INTEGRATE THE FLUX ACROSS THE HO FACE

    m_subFaceFlux = 0.0;
    for(CFuint iQd=0; iQd<nbQdPts; ++iQd) {
      m_subFaceJacob[iQd] = 0.0;
    }

    for(CFuint iQd = 0; iQd < nbQdPts; ++iQd) {
      const CFreal xiQd = m_qpPos[iQd] + m_subface_center[ibdrysubface];
      m_CP2N.ComputeBNormal(m_face_nodes,xiQd,m_qdNormal);

      m_qdNormal = -1.0 * m_qdNormal;

      m_p0[iQd] = (1-3*xiQd+2*xiQd*xiQd);
      m_p1[iQd] = (2*xiQd*xiQd-xiQd);
      m_p3[iQd] = 4*xiQd*(1-xiQd);
      (*m_qdState) = m_p0[iQd]*(*m_face_states[0]) + m_p1[iQd]*(*m_face_states[1]) + m_p3[iQd]*(*m_face_states[2]);

      const CFreal jacob = m_qdNormal.norm2();
      const CFreal invJacob = 1.0/jacob;
      computeNormalFluxAndJacob((*m_qdState), invJacob * m_qdNormal, m_qdFlux,m_qdJacob);
      m_subFaceFlux = m_subFaceFlux + jacob * m_weights[iQd] * m_qdFlux;
      m_subFaceJacob[iQd] = jacob * m_weights[iQd] * m_qdJacob;
    }





  ///set indexes of subtriangle: vertices ...
  const CFuint v0 = m_subTri[3*isubelem];
  const CFuint v1 = m_subTri[3*isubelem+1];
  const CFuint v2 = m_subTri[3*isubelem+2];

  m_subStates[0] = (*cell_states)[v0];
  m_subStates[1] = (*cell_states)[v1];
  m_subStates[2] = (*cell_states)[v2];
distdata.tStates = computeConsistentStates(&m_subStates);
  /// ... and faces

  const CFuint f0 = m_subFaceIdx[3*isubelem];
  const CFuint f1 = m_subFaceIdx[3*isubelem+1];
  const CFuint f2 = m_subFaceIdx[3*isubelem+2];

  ///Set of local indexes of the large (p2p2) face 
  ///which is adjacent to domain boundary
  const CFuint ib0 = m_subTri[0];
  const CFuint ib1 = m_subTri[4];
  const CFuint ib3 = m_subTri[1];


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

  // compute the residual and the upwind parameters k in this cell
  distdata.tStates = computeConsistentStates(&m_subStates);
  //pass the substates corresponding to vertices of subtriangle to distribution data 
  distdata.subStates = &m_subStates;

//  m_splitter->computeK(m_subStates, m_subcell_normals);

   computeK(m_subStates,m_kPlus);          
//  std::cout<<"kplus"<<m_kPlus[0]<<std::endl;
 
  // transform fluxes of subelement to distribution variables
  SafePtr<RealVector> phi = &distdata.phi;
  *phi = *getMethodData().getSolutionToDistribMatTrans()->transformFromRef(&m_subFaceFlux);


  SafePtr<std::vector<RealMatrix> >& currBetaMatrix = distdata.currBetaMat;
  currBetaMatrix = &distdata.betaMats[m_subElemIdx[isubelem]];
  cf_assert (currBetaMatrix.isNotNull());
  


   distributeN(m_kPlus,m_phiN);
   computeBlendingCoeff(m_subStates,m_theta);
   distributeLDA(m_kPlus,m_phiLDA);

   m_subresidual[0] = m_theta*m_phiN[0] + (1.0 - m_theta)*m_phiLDA[0];
   m_subresidual[1] = m_theta*m_phiN[1] + (1.0 - m_theta)*m_phiLDA[1];
   m_subresidual[2] = m_theta*m_phiN[2] + (1.0 - m_theta)*m_phiLDA[2]; 
//                 std::cout<<"kplus"<<m_subresidual[1]<<std::endl;
//                 std::cout<<"kplus"<<m_subresidual[2]<<std::endl;
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


  ///Calculate contributions to node 0 from nodes 0, 1, and 3
  m_distJacobFrom0 = 0.0;
  for(CFuint iQd = 0; iQd < nbQdPts; ++iQd) {
    m_distJacobFrom0 = m_distJacobFrom0 + m_p0[iQd]*m_subFaceJacob[iQd] ;
  }
 
  m_distJacobFrom1 = 0.0;
  for(CFuint iQd = 0; iQd < nbQdPts; ++iQd) {
    m_distJacobFrom1 = m_distJacobFrom1 + m_p1[iQd]*m_subFaceJacob[iQd];
  }

  m_distJacobFrom3 = 0.0;
  for(CFuint iQd = 0; iQd < nbQdPts; ++iQd) {
    m_distJacobFrom3 = m_distJacobFrom3 + m_p3[iQd]*m_subFaceJacob[iQd];
  }

  ///UPDATE OF SUBELEMENT STATE 0
  ///Contribution to sub-element state 0 from large p2p2 element state 0
  m_distJacobToState=((*currBetaMatrix)[0]) * m_distJacobFrom0 /*+ m_theta*(*m_kPlus[0])*m_invK*((*m_kPlus[1])+(*m_kPlus[2]))*/;

    for (CFuint iEq = 0; iEq < nbEqs; ++iEq) {
      for (CFuint jEq = 0; jEq < nbEqs; ++jEq) {
        if (!isUpdated[(*cell_states)[v0]->getLocalID()]){
          acc->addValue(v0, ib0, iEq, jEq,m_distJacobToState(iEq,jEq));
        }
      }
    }

  ///Contribution to sub-element state 0 from large p2p2 element state 1
  m_distJacobToState=((*currBetaMatrix)[0]) * m_distJacobFrom1 ;

    for (CFuint iEq = 0; iEq < nbEqs; ++iEq) {
      for (CFuint jEq = 0; jEq < nbEqs; ++jEq) {
        if (!isUpdated[(*cell_states)[v0]->getLocalID()]){
          acc->addValue(v0, ib1, iEq, jEq,m_distJacobToState(iEq,jEq));
        }
      }
    }

  ///Contribution to sub-element state 0 from large p2p2 element state 3
  m_distJacobToState=((*currBetaMatrix)[0]) * m_distJacobFrom3 /*- m_theta*(*m_kPlus[0])*m_invK*((*m_kPlus[1]))*/;

    for (CFuint iEq = 0; iEq < nbEqs; ++iEq) {
      for (CFuint jEq = 0; jEq < nbEqs; ++jEq) {
        if (!isUpdated[(*cell_states)[v0]->getLocalID()]){
          acc->addValue(v0, ib3, iEq, jEq,m_distJacobToState(iEq,jEq));
        }
      }
    }
/*
  m_distJacobToState= 0.0 - m_theta*(*m_kPlus[0])*m_invK*((*m_kPlus[2]));
    for (CFuint iEq = 0; iEq < nbEqs; ++iEq) {
      for (CFuint jEq = 0; jEq < nbEqs; ++jEq) {
        if (!isUpdated[(*cell_states)[v0]->getLocalID()]){
          acc->addValue(v0, v2, iEq, jEq,m_distJacobToState(iEq,jEq));
        }
      }
    }*/

  ///UPDATE OF SUBELEMENT STATE 1
  ///Contribution to sub-element state 1 from large p2p2 element state 0
  m_distJacobToState=((*currBetaMatrix)[1]) * m_distJacobFrom0 /*- m_theta*(*m_kPlus[1])*m_invK*((*m_kPlus[0]))*/;

    for (CFuint iEq = 0; iEq < nbEqs; ++iEq) {
      for (CFuint jEq = 0; jEq < nbEqs; ++jEq) {
        if (!isUpdated[(*cell_states)[v1]->getLocalID()]){
          acc->addValue(v1, ib0, iEq, jEq,m_distJacobToState(iEq,jEq));
        }
      }
    }

  ///Contribution to sub-element state 1 from large p2p2 element state 1
  m_distJacobToState=((*currBetaMatrix)[1]) * m_distJacobFrom1 ;

    for (CFuint iEq = 0; iEq < nbEqs; ++iEq) {
      for (CFuint jEq = 0; jEq < nbEqs; ++jEq) {
        if (!isUpdated[(*cell_states)[v1]->getLocalID()]){
          acc->addValue(v1, ib1, iEq, jEq,m_distJacobToState(iEq,jEq));
        }
      }
    }

   

  ///Contribution to sub-element state 1 from large p2p2 element state 3
  m_distJacobToState=((*currBetaMatrix)[1]) * m_distJacobFrom3 /*+ m_theta*(*m_kPlus[1])*m_invK*((*m_kPlus[0])+(*m_kPlus[2]))*/;

    for (CFuint iEq = 0; iEq < nbEqs; ++iEq) {
      for (CFuint jEq = 0; jEq < nbEqs; ++jEq) {
        if (!isUpdated[(*cell_states)[v1]->getLocalID()]){
          acc->addValue(v1, ib3, iEq, jEq,m_distJacobToState(iEq,jEq));
        }
      }
    }

///Contribution to sub-element state 1 from large p2p2 element state 5
//   m_distJacobToState= 0.0 - m_theta*(*m_kPlus[1])*m_invK*((*m_kPlus[2]));
// 
//     for (CFuint iEq = 0; iEq < nbEqs; ++iEq) {
//       for (CFuint jEq = 0; jEq < nbEqs; ++jEq) {
//         if (!isUpdated[(*cell_states)[v1]->getLocalID()]){
//           acc->addValue(v1, v2, iEq, jEq,m_distJacobToState(iEq,jEq));
//         }
//       }
//     }

  ///UPDATE OF SUBELEMENT STATE 2
  ///Contribution to sub-element state 2 from large p2p2 element state 0
  m_distJacobToState=((*currBetaMatrix)[2]) * m_distJacobFrom0 /*- m_theta*(*m_kPlus[2])*m_invK*((*m_kPlus[0])) */;

    for (CFuint iEq = 0; iEq < nbEqs; ++iEq) {
      for (CFuint jEq = 0; jEq < nbEqs; ++jEq) {
        if (!isUpdated[(*cell_states)[v2]->getLocalID()]){
          acc->addValue(v2, ib0, iEq, jEq,m_distJacobToState(iEq,jEq));
        }
      }
    }

  ///Contribution to sub-element state 2 from large p2p2 element state 1
  m_distJacobToState=((*currBetaMatrix)[2]) * m_distJacobFrom1;

    for (CFuint iEq = 0; iEq < nbEqs; ++iEq) {
      for (CFuint jEq = 0; jEq < nbEqs; ++jEq) {
        if (!isUpdated[(*cell_states)[v2]->getLocalID()]){
          acc->addValue(v2, ib1, iEq, jEq,m_distJacobToState(iEq,jEq));
        }
      }
    }

  ///Contribution to sub-element state 2 from large p2p2 element state 3
  m_distJacobToState=((*currBetaMatrix)[2]) * m_distJacobFrom3 /*- m_theta*(*m_kPlus[2])*m_invK*((*m_kPlus[1]))*/;

    for (CFuint iEq = 0; iEq < nbEqs; ++iEq) {
      for (CFuint jEq = 0; jEq < nbEqs; ++jEq) {
        if (!isUpdated[(*cell_states)[v2]->getLocalID()]){
          acc->addValue(v2, ib3, iEq, jEq,m_distJacobToState(iEq,jEq));
        }
      }
    }

///Contribution to sub-element state 2 from large p2p2 element state 3
  m_distJacobToState=((*currBetaMatrix)[2]) * m_distJacobFrom3 /*+ m_theta*(*m_kPlus[2])*m_invK*((*m_kPlus[0])+(*m_kPlus[1]))*/;

    for (CFuint iEq = 0; iEq < nbEqs; ++iEq) {
      for (CFuint jEq = 0; jEq < nbEqs; ++jEq) {
        if (!isUpdated[(*cell_states)[v2]->getLocalID()]){
          acc->addValue(v2, v2, iEq, jEq,m_distJacobToState(iEq,jEq));
        }
      }
    }


} ///Loop over subtriangles adjacent to wall boundary



    // add the contributions to the jacobian
    jacobMatrix->addValues(*acc);

    acc->reset();

    // release the cell
    geoBuilderCell.releaseGE();
    // release the face
    geoBuilder->releaseGE();
  }
}
//////////////////////////////////////////////////////////////////////////////////////
void WeakSlipWall2DHOImplIsoP2UpwindBx::computeK(const std::vector<Framework::State*>& states,
                std::vector<RealMatrix*>& m_kPlus){
  
  DataHandle< CFreal> updateCoeff = socket_updateCoeff.getDataHandle();
  const CFuint m_nbStatesInCell = states.size();

  const CFreal m_invDim = 1.0/2.0;
  // The transformation of the normal is needed if the coordinate system is rotated
  for (CFuint iState = 0; iState < m_nbStatesInCell; ++iState) {
    // The transformation of the normal is needed if the coordinate system
    // is rotated. The normal is adimensionalized, so there is need to multiply
    // by the nodeArea when computing the k parameter

    for (CFuint iDim = 0; iDim < 2; ++iDim) {
        m_adimNormal[iDim] = m_subcell_normals->getNodalNormComp(iState, iDim);
    }
    m_adimNormal *= 1. / m_subcell_normals->getAreaNode(iState);
   
    //std::cout<<m_subcell_normals <<std::endl;

    getMethodData().getDistribVar()->splitJacobian(*m_kPlus[iState],
                 *m_kMin[iState],
                 *m_eValues[iState],
                 m_adimNormal);

    CFreal m_nodeArea = m_subcell_normals->getAreaNode(iState);

    *m_kPlus[iState] *= m_invDim * m_nodeArea;
    *m_kMin[iState]  *= m_invDim * m_nodeArea;

    if (!getMethodData().getDistributionData().isPerturb) {
      const CFreal maxEigenValue = std::max(0.0, m_eValues[iState]->max());

      updateCoeff[states[iState]->getLocalID()] += m_invDim*m_nodeArea*maxEigenValue;

    }
  }



}
//////////////////////////////////////////////////////////////////////////////
void WeakSlipWall2DHOImplIsoP2UpwindBx::distributeN(std::vector<RealMatrix*> & m_kPlus,vector<RealVector> & phiN){

  const vector<State*>& tStates = *getMethodData().getDistributionData().tStates;
  const RealVector& phiT = getMethodData().getDistributionData().phi;


  m_tmp = (*m_kPlus[0])  + (*m_kMin[0]);
/*
   for (CFuint i=0; i < tmp.size(); ++i) {
      sumAbsK[i] = std::abs(tmp[i]);
    }*/

  m_sumKU = m_tmp*(*tStates[0]);

  m_sumKplusU = (*m_kPlus[0]) * (*tStates[0]);

  m_sumKplus = *m_kPlus[0];


  for (CFuint iState = 1; iState < 3; ++iState) {

    m_sumKplusU += (*m_kPlus[iState])*(*tStates[iState]);
    m_sumKplus += *m_kPlus[iState];

    m_tmp = (*m_kPlus[iState])  + (*m_kMin[iState]);
    m_sumKU += m_tmp*(*tStates[iState]);

  }
  m_inverter->invert(m_sumKplus, m_invK);

  CFLogDebugMax( "invK = " << "\n" <<m_invK << "\n");

  m_uInflow = m_invK * (m_sumKplusU - m_sumKU);

  CFLogDebugMax( "uInflow = " << "\n" << m_uInflow << "\n");


  for (CFuint iState = 0; iState < 3; ++iState) {

    phiN[iState] = (*m_kPlus[iState])*(*tStates[iState] - m_uInflow);
    m_tmp = (*m_kPlus[iState])*m_invK;
    phiN[iState] -= m_tmp*(m_sumKU - phiT);

}


}

//////////////////////////////////////////////////////////////////////////////
void WeakSlipWall2DHOImplIsoP2UpwindBx::distributeLDA(std::vector<RealMatrix*> & m_kPlus,vector<RealVector> & phiLDA){

 const RealVector& phiT = getMethodData().getDistributionData().phi;

  m_sumKplus = *m_kPlus[0];
  for (CFuint iState = 1; iState < 3; ++iState) {
    m_sumKplus  += *m_kPlus[iState];
  }
  m_inverter->invert(m_sumKplus, m_invK);
  m_uTemp = m_invK*phiT;

  for (CFuint iState = 0; iState < 3; ++iState) {
    phiLDA[iState] = (*m_kPlus[iState])*m_uTemp;
  if (getMethodData().getDistributionData().computeBetas) {
      (*getMethodData().getDistributionData().currBetaMat)[iState] =
  (*m_kPlus[iState])*m_invK;
    }
   }


}
//////////////////////////////////////////////////////////////////////////////
void WeakSlipWall2DHOImplIsoP2UpwindBx::computeBlendingCoeff(const std::vector<Framework::State*>& states,CFreal & result)
{

  CFreal vol = getMethodData().getDistributionData().cell->computeVolume();

   vol /= 4.0;


  m_updateVar   = getMethodData().getUpdateVar();
  const RealVector& lData =  _cterm->getPhysicalData();

  const CFuint nbStates = states.size();
  const CFuint dim = PhysicalModelStack::getActive()->getDim();

  _grad = 0.0;
  for (CFuint i = 0; i < nbStates; ++i) {
    m_updateVar->computePhysicalData(*states[i], _pData);
    for (CFuint iDim = 0; iDim < dim; ++iDim) {
      _grad[iDim] += _pData[_varID]*m_subcell_normals->getNodalNormComp(i,iDim);
    }
  }

  cf_assert(dim == DIM_2D);
  const CFreal h = 2.0*std::sqrt(static_cast<CFreal>(vol/MathTools::MathConsts::CFrealPi()));
  CFreal theta = 0.0;

     theta = _grad[XX] * lData[EulerTerm::VX] + _grad[YY] * lData[EulerTerm::VY];

     theta *= std::pow(    (std::pow(lData[EulerTerm::VX],2) + std::pow(lData[EulerTerm::VY],2)), 0.5   );

//      updateScalingValues();
    _length = 1.0;
    _deltaP = 43000.0;
    _speed = 280.0;



     theta *= _length/(vol*dim*_deltaP * _speed * _speed);   //modified!!

  // AL: this was missing
  theta = max(0.,theta);

  theta = min(1.0,theta*theta*h);
  result = theta;


//   CFreal theta = _grad[XX]*lData[EulerTerm::VX] + _grad[YY]*lData[EulerTerm::VY];
//   if (dim == DIM_3D) {
//     theta += _grad[ZZ]*lData[EulerTerm::VZ];
//   }
//   //   cout << "grad*V  = " << theta << endl;
//   //   cout << "_length = " << _length << endl;
//   //   cout << "_deltaP = " << _deltaP << endl;
//   //   cout << "_speed = " << _speed << endl;
// 
//   theta *= _length/(vol*dim*_deltaP*_speed);
// 
//   result = min(1.0,theta*theta*h);

  // this->_alpha = min(1.0,theta*h);
  // cout << "theta = " << this->_alpha << endl;
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
