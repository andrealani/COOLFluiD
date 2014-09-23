#include "Common/BadValueException.hh"

#include "Framework/MethodStrategyProvider.hh"
#include "Framework/ContourIntegrator.hh"
#include "Framework/MeshData.hh"
#include "FluctSplit/FluctSplitHO.hh"
#include "Framework/SubSystemStatus.hh"
#include "Framework/BaseTerm.hh"

#include "FluctSplit/HOCRD_SplitStrategyIsoP2.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::MathTools;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

MethodStrategyProvider<HOCRD_SplitStrategyIsoP2,
                       FluctuationSplitData,
                       FluctuationSplitStrategy,
                       FluctSplitHOModule>
aHOCRDFluctSplitIsoP2StrategyProvider("HOCRDIsoP2");

//////////////////////////////////////////////////////////////////////////////

HOCRD_SplitStrategyIsoP2::HOCRD_SplitStrategyIsoP2(const std::string& name) :
  FluctuationSplitStrategy(name),
  socket_updateCoeff("updateCoeff"),
  m_solutionVar(CFNULL),
  m_updateVar(CFNULL),
  m_unitFaceNormals(),
  m_scaledFaceNormals(),
  matrix_face_norms(3,2),
  matrix_node_norms(3,2),
  vector_face_areas(3),
  vector_node_areas(3)
{
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

HOCRD_SplitStrategyIsoP2::~HOCRD_SplitStrategyIsoP2()
{
  if (isSetup()) unsetup();
}

//////////////////////////////////////////////////////////////////////////////

void HOCRD_SplitStrategyIsoP2::unsetup()
{
  deletePtr(m_subcell_normals);

  // AL: this could be changed ...
  // where are deleted the qdstates ???
  for (CFuint i = 0; i < m_qdExtraVars.size(); ++i) {
    deletePtr(m_qdExtraVars[i]);
  }
}

//////////////////////////////////////////////////////////////////////////////

void HOCRD_SplitStrategyIsoP2::setup()
{
  CFAUTOTRACE;

  // first call parent method
  FluctuationSplitStrategy::setup();

  DataHandle< CFreal> updateCoeff = socket_updateCoeff.getDataHandle();

  if (getMethodData().isMultipleSplitter())
    throw Common::BadValueException(FromHere(), "Cannot use HOCRD on curved elements with multiple splitters");

  m_splitter    = getMethodData().getSplitter();
  m_solutionVar = getMethodData().getSolutionVar();
  m_updateVar   = getMethodData().getUpdateVar();

  m_unitFaceNormals.resize(MeshDataStack::getActive()->Statistics().getMaxNbFacesInCell());
  for (CFuint  i = 0; i < m_unitFaceNormals.size(); ++i)
  {
    m_unitFaceNormals[i].resize(PhysicalModelStack::getActive()->getDim());
  }

    // number of quadrature points used to compute the fluctuation
  //const CFuint nbQdPts = 3;


  m_scaledFaceNormals.resize(nbQdPts);
  for (CFuint  i = 0; i < m_scaledFaceNormals.size(); ++i)
  {
    m_scaledFaceNormals[i].resize(PhysicalModelStack::getActive()->getDim());
  }
  
  // physical data evaluated in the quadrature points
  m_pdata.resize(nbQdPts);
  for (CFuint  i = 0; i < nbQdPts; ++i) {
    PhysicalModelStack::getActive()->getImplementor()->getConvectiveTerm()->
      resizePhysicalData(m_pdata[i]);
  }  
  
  // this is of fundamental importance
  // sets the maximum number of quadrature points in contour integration
  
  // sub elemt table
  substates.resize(3);   // 3 states in each sub element
  subresidual.resize(3); // 3 residuals in each sub element

  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();

  subresidual[0].resize(nbEqs);
  subresidual[1].resize(nbEqs);
  subresidual[2].resize(nbEqs);

  m_phisubT.resize(4); // P2 triangle has 4 sub triangles

  m_phisubT[0] = new RealVector(nbEqs);
  m_phisubT[1] = new RealVector(nbEqs);
  m_phisubT[2] = new RealVector(nbEqs);
  m_phisubT[3] = new RealVector(nbEqs);


  // sub elemt table : contain the faces of each sub-element
  subelemtable.resize(4,3); // 4 sub elems with 3 faces each
  subelemfacedir.resize(4,3);

  subelemtable(0,0) = 0;
  subelemtable(0,1) = 1;
  subelemtable(0,2) = 2;

  subelemtable(1,0) = 3;
  subelemtable(1,1) = 4;
  subelemtable(1,2) = 5;

  subelemtable(2,0) = 6;
  subelemtable(2,1) = 7;
  subelemtable(2,2) = 8;

  subelemtable(3,0) = 0; // faces have negative orientation
  subelemtable(3,1) = 4;
  subelemtable(3,2) = 8;

// sub elemt table : contain the orientation of the faces of each sub-element
  subelemfacedir(0,0) = 1.;
  subelemfacedir(0,1) = 1.;
  subelemfacedir(0,2) = 1.;

  subelemfacedir(1,0) = 1.;
  subelemfacedir(1,1) = 1.;
  subelemfacedir(1,2) = 1.;

  subelemfacedir(2,0) = 1.;
  subelemfacedir(2,1) = 1.;
  subelemfacedir(2,2) = 1.;

  subelemfacedir(3,0) = -1.; // faces have negative orientation
  subelemfacedir(3,1) = -1.;
  subelemfacedir(3,2) = -1.;


  // sub face table : contain the node that contain each face
  subfacetable.resize(9,2); // 9 sub faces with 2 states each

  subfacetable(0,0) = 3;
  subfacetable(0,1) = 5;

  subfacetable(1,0) = 5;
  subfacetable(1,1) = 0;

  subfacetable(2,0) = 0;
  subfacetable(2,1) = 3;

  subfacetable(3,0) = 1;
  subfacetable(3,1) = 4;

  subfacetable(4,0) = 4;
  subfacetable(4,1) = 3;

  subfacetable(5,0) = 3;
  subfacetable(5,1) = 1;

  subfacetable(6,0) = 4;
  subfacetable(6,1) = 2;

  subfacetable(7,0) = 2;
  subfacetable(7,1) = 5;

  subfacetable(8,0) = 5;
  subfacetable(8,1) = 4;

  faceflux.resize(subfacetable.nbRows()); // one flux per sub face
  for (CFuint i = 0; i < faceflux.size(); ++i)
    faceflux[i].resize(nbEqs);

  xi_ref.resize(6);
  eta_ref.resize(6);

  xi_ref[0] = 0.0;
  xi_ref[1] = 1.0;
  xi_ref[2] = 0.0;
  xi_ref[3] = 0.5;
  xi_ref[4] = 0.5;
  xi_ref[5] = 0.0;

  eta_ref[0] = 0.0;
  eta_ref[1] = 0.0;
  eta_ref[2] = 1.0;
  eta_ref[3] = 0.0;
  eta_ref[4] = 0.5;
  eta_ref[5] = 0.5;

  qd0.resize(nbQdPts); // quadrature points per face
  qd1.resize(nbQdPts); // quadrature points per face


  ///Quadrature with 3 points
 const CFreal s  = std::sqrt( 0.6 );
 const CFreal a0 = ( 1.0 - s )*0.5;
 const CFreal a1 = ( 1.0 + s )*0.5;

 qd0[0] = a0;  qd1[0] = a1;
 qd0[1] = a1;  qd1[1] = a0;
 qd0[2] = .5;  qd1[2] = .5;



  ///Quadrature with 5 points
//   const CFreal a0 = ( 1.0 - 0.9061798459 ) * 0.5;
//   const CFreal a1 = ( 1.0 + 0.9061798459 ) * 0.5;
//   const CFreal b0 = ( 1.0 - 0.5384693101 ) * 0.5;
//   const CFreal b1 = ( 1.0 + 0.5384693101 ) * 0.5;


//   qd0[0] = a0;  qd1[0] = a1;
//   qd0[1] = b0;  qd1[1] = b1;
//   qd0[2] = 0.5; qd1[2] = 0.5;
//   qd0[3] = b1;  qd1[3] = b0;
//   qd0[4] = a1;  qd1[4] = a0;


  wqd.resize(nbQdPts); // 3 or 5 or ... quadrature points per face

///Define the weights: 3 quadrature points
 wqd[0] = 5.0/36.0;
 wqd[1] = 5.0/36.0;
 wqd[2] = 8.0/36.0;

///Define the weights: 5 quadrature points
//   wqd[0] = 0.2369268850*0.25;
//   wqd[1] = 0.4786286705*0.25;
//   wqd[2] = 0.5688888889*0.25;
//   wqd[3] = 0.4786286705*0.25;
//   wqd[4] = 0.2369268850*0.25;




  qdstates.resize(nbQdPts); // 3 quadrature points per face
//  qdstates[0] = new State();
//  qdstates[1] = new State();
//  qdstates[2] = new State();

  for(CFuint iState=0;iState < nbQdPts; ++iState) {
	qdstates[iState] = new State();
  }


  const CFuint extra_var_size = getMethodData().getUpdateVar()->getExtraPhysicalVarsSize();
  m_qdExtraVars.resize(qdstates.size());
  for (CFuint i = 0; i < m_qdExtraVars.size(); ++i) {
    m_qdExtraVars[i] = new RealVector(extra_var_size);
  }

  qdnodes.resize(nbQdPts); // 3 quadrature points per face

  for(CFuint iNode = 0;iNode<nbQdPts;++iNode) {
	qdnodes[iNode] = new Node();
	qdstates[iNode]->setSpaceCoordinates(qdnodes[iNode]);
  }

/*
  qdnodes[0] = new Node();
  qdstates[0]->setSpaceCoordinates(qdnodes[0]);
  qdnodes[1] = new Node();
  qdstates[1]->setSpaceCoordinates(qdnodes[1]);
  qdnodes[2] = new Node();
  qdstates[2]->setSpaceCoordinates(qdnodes[2]);
*/



  facenormal.resize(2); // only 2D

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

void HOCRD_SplitStrategyIsoP2::computeFluctuation(vector<RealVector>& residual)
{
  DistributionData& distdata = getMethodData().getDistributionData();

  vector<State*>& states = *distdata.states;


  cf_assert (states.size() == 6); // P2 triangles for solution space
  cf_assert (residual.size() >= states.size()); // state residual

  // reset the residual because we will accumulate the sub element contributions
  for (CFuint i = 0; i < residual.size(); ++i) { residual[i] = 0.0; }

  InwardNormalsData& cellnormals = (*socket_normals.getDataHandle() [distdata.cellID]);

  cf_assert (cellnormals.nbFaces()  == 9); // P2 triangles have 9 normals as for now


    const CFreal nx0 = cellnormals.getFaceNormComp(0,XX);
    const CFreal ny0 = cellnormals.getFaceNormComp(0,YY);

    const CFreal nx1 = cellnormals.getFaceNormComp(1,XX);
    const CFreal ny1 = cellnormals.getFaceNormComp(1,YY);

    const CFreal nx2 = cellnormals.getFaceNormComp(2,XX);
    const CFreal ny2 = cellnormals.getFaceNormComp(2,YY);

    const CFreal nx3 = cellnormals.getFaceNormComp(3,XX);
    const CFreal ny3 = cellnormals.getFaceNormComp(3,YY);

    const CFreal nx4 = cellnormals.getFaceNormComp(4,XX);
    const CFreal ny4 = cellnormals.getFaceNormComp(4,YY);

    const CFreal nx5 = cellnormals.getFaceNormComp(5,XX);
    const CFreal ny5 = cellnormals.getFaceNormComp(5,YY);

    const CFreal nx6 = cellnormals.getFaceNormComp(6,XX);
    const CFreal ny6 = cellnormals.getFaceNormComp(6,YY);

    const CFreal nx7 = cellnormals.getFaceNormComp(7,XX);
    const CFreal ny7 = cellnormals.getFaceNormComp(7,YY);

    const CFreal nx8 = cellnormals.getFaceNormComp(8,XX);
    const CFreal ny8 = cellnormals.getFaceNormComp(8,YY);

  (*m_phisubT[0]) = 0.;
  (*m_phisubT[1]) = 0.;
  (*m_phisubT[2]) = 0.;
  (*m_phisubT[3]) = 0.;


  computeHOCurvedFluctuation();



  /*****         Triangle 0-3-5          *****/
  // faces 0-1-2

  substates[0] = states[0];
  substates[1] = states[3];
  substates[2] = states[5];

  // stupidily put the normals data into a matrix to put it in the normals data
  matrix_face_norms (0,XX) = nx2;
  matrix_face_norms (1,XX) = nx0;
  matrix_face_norms (2,XX) = nx1;

  matrix_face_norms (0,YY) = ny2;
  matrix_face_norms (1,YY) = ny0;
  matrix_face_norms (2,YY) = ny1;


  vector_face_areas[0] = cellnormals.getAreaFace(2);
  vector_face_areas[1] = cellnormals.getAreaFace(0);
  vector_face_areas[2] = cellnormals.getAreaFace(1);

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

  distdata.tStates = computeConsistentStates(&substates);
  distdata.subStates = &substates;

  m_splitter->computeK(substates, m_subcell_normals);

  // transform fluxes of subelement to distribution variables
  SafePtr<RealVector> phi = &distdata.phi;
  *phi = *getMethodData().getSolutionToDistribMatTrans()->transformFromRef(m_phisubT[0]);

  SafePtr<vector<RealMatrix> >& currBetaMatrix = distdata.currBetaMat;

  currBetaMatrix = &distdata.betaMats[0];

  cf_assert (currBetaMatrix.isNotNull());
  m_splitter->distribute(subresidual);

  residual[0] += subresidual[0];
  residual[3] += subresidual[1];
  residual[5] += subresidual[2];

  /*****         Triangle 3-1-4          *****/
  // faces 4-5-3

  substates[0] = states[3];
  substates[1] = states[1];
  substates[2] = states[4];

  // stupidily put the normals data into a matrix to put it in the normals data
  matrix_face_norms (0,XX) = nx5;
  matrix_face_norms (1,XX) = nx3;
  matrix_face_norms (2,XX) = nx4;

  matrix_face_norms (0,YY) = ny5;
  matrix_face_norms (1,YY) = ny3;
  matrix_face_norms (2,YY) = ny4;

  vector_face_areas[0] = cellnormals.getAreaFace(5);
  vector_face_areas[1] = cellnormals.getAreaFace(3);
  vector_face_areas[2] = cellnormals.getAreaFace(4);

  // this is always the same because triangles face normals are opposite to nodal normals
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

  distdata.tStates = computeConsistentStates(&substates);
  *distdata.subStates = substates;
  // compute the residual and the upwind parameters k in this cell

  m_splitter->computeK(substates, m_subcell_normals);

  // transform fluxes of subelement to distribution variables
  *phi = *getMethodData().getSolutionToDistribMatTrans()->transformFromRef(m_phisubT[1]);

  currBetaMatrix = &distdata.betaMats[1];
  cf_assert (currBetaMatrix.isNotNull());
  m_splitter->distribute(subresidual);

  residual[3] += subresidual[0];
  residual[1] += subresidual[1];
  residual[4] += subresidual[2];

  /*****         Triangle 5-4-2          *****/
  // faces 7-8-6

  substates[0] = states[5];
  substates[1] = states[4];
  substates[2] = states[2];

  // stupidily put the normals data into a matrix to put it in the normals data
  matrix_face_norms (0,XX) = nx8;
  matrix_face_norms (1,XX) = nx6;
  matrix_face_norms (2,XX) = nx7;

  matrix_face_norms (0,YY) = ny8;
  matrix_face_norms (1,YY) = ny6;
  matrix_face_norms (2,YY) = ny7;

  vector_face_areas[0] = cellnormals.getAreaFace(8);
  vector_face_areas[1] = cellnormals.getAreaFace(6);
  vector_face_areas[2] = cellnormals.getAreaFace(7);

  // this is always the same because triangles face normals are opposite to nodal normals
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

  distdata.tStates = computeConsistentStates(&substates);
  *distdata.subStates = substates;
  // compute the residual and the upwind parameters k in this cell

  m_splitter->computeK(substates, m_subcell_normals);

  // transform fluxes of subelement to distribution variables
  *phi = *getMethodData().getSolutionToDistribMatTrans()->transformFromRef(m_phisubT[2]);

  currBetaMatrix = &distdata.betaMats[2];
  cf_assert (currBetaMatrix.isNotNull());
  m_splitter->distribute(subresidual);

  residual[5] += subresidual[0];
  residual[4] += subresidual[1];
  residual[2] += subresidual[2];

  /*****         Triangle 4-5-3          *****/
  // faces 0-5-6 NEGATIVE!!!

  substates[0] = states[4];
  substates[1] = states[5];
  substates[2] = states[3];

  // stupidily put the normals data into a matrix to put it in the normals data
  matrix_face_norms (0,XX) = - nx8;
  matrix_face_norms (1,XX) = - nx0;
  matrix_face_norms (2,XX) = - nx4;

  matrix_face_norms (0,YY) = - ny8;
  matrix_face_norms (1,YY) = - ny0;
  matrix_face_norms (2,YY) = - ny4;

  vector_face_areas[0] = cellnormals.getAreaFace(8);
  vector_face_areas[1] = cellnormals.getAreaFace(0);
  vector_face_areas[2] = cellnormals.getAreaFace(4);

  // this is always the same because triangles face normals are opposite to nodal normals
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

  distdata.tStates = computeConsistentStates(&substates);
  *distdata.subStates = substates;

  m_splitter->computeK(substates, m_subcell_normals);

  // transform fluxes of subelement to distribution variables
  *phi = *getMethodData().getSolutionToDistribMatTrans()->transformFromRef(m_phisubT[3]);

  currBetaMatrix = &distdata.betaMats[3];
  cf_assert (currBetaMatrix.isNotNull());
  m_splitter->distribute(subresidual);

  residual[4] += subresidual[0];
  residual[5] += subresidual[1];
  residual[3] += subresidual[2];


//   if (output) { CF_DEBUG_EXIT; }
}

//////////////////////////////////////////////////////////////////////////////

void HOCRD_SplitStrategyIsoP2::computeHOCurvedFluctuation()
{

  cf_assert (qdstates.size() == nbQdPts); // only triags so three quadrature points per face

  vector<State*>& states = *getMethodData().getDistributionData().states;

  State& state0 = *(states[0]);
  State& state1 = *(states[1]);
  State& state2 = *(states[2]);
  State& state3 = *(states[3]);
  State& state4 = *(states[4]);
  State& state5 = *(states[5]);

  vector<Node*>& nodes = *getMethodData().getDistributionData().cell->getNodes();


  cf_assert (nodes.size()  == 6); // P2 triangles for geometry space

  const CFreal x0 = (*nodes[0])[XX];
  const CFreal x1 = (*nodes[1])[XX];
  const CFreal x2 = (*nodes[2])[XX];

  const CFreal x3 = (*nodes[3])[XX];
  const CFreal x4 = (*nodes[4])[XX];
  const CFreal x5 = (*nodes[5])[XX];

  const CFreal y0 = (*nodes[0])[YY];
  const CFreal y1 = (*nodes[1])[YY];
  const CFreal y2 = (*nodes[2])[YY];

  const CFreal y3 = (*nodes[3])[YY];
  const CFreal y4 = (*nodes[4])[YY];
  const CFreal y5 = (*nodes[5])[YY];



  InwardNormalsData& cellnormals = (*socket_normals.getDataHandle() [getMethodData().getDistributionData().cellID]);

 // Computation  of the fluctuation of each faces of the cell
for (CFuint iFace = 0; iFace < subfacetable.nbRows(); ++iFace)
{
    facenormal[XX] = - cellnormals.getFaceNormComp(iFace,XX);
    facenormal[YY] = - cellnormals.getFaceNormComp(iFace,YY);


    for (CFuint iQd = 0; iQd < nbQdPts; ++iQd)
    {
      /// Position of integration point in reference space
      const CFreal xi  = qd0[iQd] * xi_ref[subfacetable(iFace,0)] + qd1[iQd] * xi_ref[subfacetable(iFace,1)];
      const CFreal eta = qd0[iQd] * eta_ref[subfacetable(iFace,0)] + qd1[iQd] * eta_ref[subfacetable(iFace,1)];

      /// Values of shape functions in reference space:
      L1 = 1-xi-eta;
      L2 = xi;
      L3 = eta;

      const CFreal SF0 = L1*( 2.0*L1 - 1.0 );
      const CFreal SF1 = L2*( 2.0*L2 - 1.0 );
      const CFreal SF2 = L3*( 2.0*L3 - 1.0 );
      const CFreal SF3 = 4.0*L1*L2;
      const CFreal SF4 = 4.0*L3*L2;
      const CFreal SF5 = 4.0*L1*L3;


      /// Compute the position of the integration point in 'physical' space
      (*qdnodes[iQd])[XX] =  SF0 * x0 + SF1 * x1 + SF2 * x2 + SF3 * x3 + SF4 * x4 + SF5 * x5;
      (*qdnodes[iQd])[YY] =  SF0 * y0 + SF1 * y1 + SF2 * y2 + SF3 * y3 + SF4 * y4 + SF5 * y5;

      /// Values of shape functions in physical space:
      (*qdstates[iQd]) = SF0 * state0 +  SF1 * state1 + SF2 * state2 + SF3* state3 + SF4 * state4 + SF5 * state5;

      /// Compute normals at the integration points
      m_CP2N.ComputeNormal(nodes,iFace,xi,eta,m_scaledFaceNormals[iQd]);

      m_scaledFaceNormals[iQd] *= -1.0;


    }
    
    computeStatesData(nbQdPts, m_updateVar, qdstates, m_pdata, m_qdExtraVars);// three quadrature points per face
    
    faceflux[iFace] = 0.;
    for (CFuint iQd = 0; iQd < nbQdPts; ++iQd)
    {
//       faceflux[iFace] += wqd[iQd] * m_updateVar->getFlux()((*qdstates[iQd]),facenormal);
//       faceflux[iFace] += wqd[iQd] * m_updateVar->getFlux()((*qdstates[iQd]),m_scaledFaceNormals[iQd]);


        const CFreal jacob = m_scaledFaceNormals[iQd].norm2();
        const CFreal invJacob = 1.0/jacob;

        faceflux[iFace] += jacob * wqd[iQd] * m_updateVar->getFlux()(m_pdata[iQd],invJacob*m_scaledFaceNormals[iQd]);
    }



  }

  cf_assert (m_phisubT.size() == subelemtable.nbRows());

 // recomposition of the fluctuation of each sub cell (taking care of the orientation)
  for (CFuint iRow = 0; iRow < subelemtable.nbRows(); ++iRow)
  {
    RealVector& phi = (*m_phisubT[iRow]);
    phi = 0.;
    for (CFuint jCol = 0; jCol < subelemtable.nbCols(); ++jCol)
    {
      phi -= subelemfacedir(iRow,jCol) * faceflux[subelemtable(iRow,jCol)];
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<Framework::BaseDataSocketSink> >
HOCRD_SplitStrategyIsoP2::needsSockets()
{
  std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > result
    = FluctuationSplitStrategy::needsSockets();
  result.push_back(&socket_updateCoeff);
  return result;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD
