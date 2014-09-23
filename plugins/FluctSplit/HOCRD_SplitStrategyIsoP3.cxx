#include "FluctSplit/HOCRD_SplitStrategyIsoP3.hh"

#include "Common/BadValueException.hh"
#include "Framework/MethodStrategyProvider.hh"
#include "Framework/ContourIntegrator.hh"
#include "Framework/MeshData.hh"
#include "Framework/BaseTerm.hh"
#include "FluctSplit/FluctSplitHO.hh"
#include "Framework/SubSystemStatus.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::MathTools;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

MethodStrategyProvider<HOCRD_SplitStrategyIsoP3,
                       FluctuationSplitData,
                       FluctuationSplitStrategy,
                       FluctSplitHOModule>
aHOCRDFluctSplitIsoP3StrategyProvider("HOCRDIsoP3");

//////////////////////////////////////////////////////////////////////////////

HOCRD_SplitStrategyIsoP3::HOCRD_SplitStrategyIsoP3(const std::string& name) :
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

HOCRD_SplitStrategyIsoP3::~HOCRD_SplitStrategyIsoP3()
{
  if (isSetup()) unsetup();
}

//////////////////////////////////////////////////////////////////////////////

void HOCRD_SplitStrategyIsoP3::unsetup()
{
  deletePtr(m_subcell_normals);

  // AL: this could be changed ...
  // where are deleted the qdstates ???
  for (CFuint i = 0; i < m_qdExtraVars.size(); ++i) {
    deletePtr(m_qdExtraVars[i]);
  }


  for(CFuint i=0; i<nbQdPts; ++i) delete qdstates[i];
}

//////////////////////////////////////////////////////////////////////////////

void HOCRD_SplitStrategyIsoP3::setup()
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
  
  // sub elemt table
  substates.resize(3);   // 3 states in each sub element
  subresidual.resize(3); // 3 residuals in each sub element

  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();

  subresidual[0].resize(nbEqs);
  subresidual[1].resize(nbEqs);
  subresidual[2].resize(nbEqs);

  m_phisubT.resize(9); // P3 triangle has 9 sub triangles

  m_phisubT[0] = new RealVector(nbEqs);
  m_phisubT[1] = new RealVector(nbEqs);
  m_phisubT[2] = new RealVector(nbEqs);
  m_phisubT[3] = new RealVector(nbEqs);
  m_phisubT[4] = new RealVector(nbEqs);
  m_phisubT[5] = new RealVector(nbEqs);
  m_phisubT[6] = new RealVector(nbEqs);
  m_phisubT[7] = new RealVector(nbEqs);
  m_phisubT[8] = new RealVector(nbEqs);

  // sub elemt table : contain the faces of each sub-element
  subelemtable.resize(9,3); // 9 sub elems with 3 faces each
  subelemfacedir.resize(9,3);

  subelemtable(0,0) = 0;
  subelemtable(0,1) = 1;
  subelemtable(0,2) = 2;

  subelemtable(1,0) = 3;
  subelemtable(1,1) = 4;
  subelemtable(1,2) = 5;

  subelemtable(2,0) = 6;
  subelemtable(2,1) = 7;
  subelemtable(2,2) = 8;

  subelemtable(3,0) = 9;
  subelemtable(3,1) = 10;
  subelemtable(3,2) = 11;

  subelemtable(4,0) = 12;
  subelemtable(4,1) = 13;
  subelemtable(4,2) = 14;

  subelemtable(5,0) = 15;
  subelemtable(5,1) = 16;
  subelemtable(5,2) = 17;

  subelemtable(6,0) = 0;  // faces have negative orientation
  subelemtable(6,1) = 16;
  subelemtable(6,2) = 14;

  subelemtable(7,0) = 15; // faces have negative orientation
  subelemtable(7,1) = 4;
  subelemtable(7,2) = 11;

  subelemtable(8,0) = 12; // faces have negative orientation
  subelemtable(8,1) = 10;
  subelemtable(8,2) = 8;

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

  subelemfacedir(3,0) = 1.;
  subelemfacedir(3,1) = 1.;
  subelemfacedir(3,2) = 1.;

  subelemfacedir(4,0) = 1.;
  subelemfacedir(4,1) = 1.;
  subelemfacedir(4,2) = 1.;

  subelemfacedir(5,0) = 1.;
  subelemfacedir(5,1) = 1.;
  subelemfacedir(5,2) = 1.;

  subelemfacedir(6,0) = -1.;
  subelemfacedir(6,1) = -1.;
  subelemfacedir(6,2) = -1.;

  subelemfacedir(7,0) = -1.;
  subelemfacedir(7,1) = -1.;
  subelemfacedir(7,2) = -1.;

  subelemfacedir(8,0) = -1.;
  subelemfacedir(8,1) = -1.;
  subelemfacedir(8,2) = -1.;

  // sub face table : contain the node that contain each face
  subfacetable.resize(18,2); // 18 sub faces with 2 states each

  subfacetable(0,0) = 3;
  subfacetable(0,1) = 8;

  subfacetable(1,0) = 8;
  subfacetable(1,1) = 0;

  subfacetable(2,0) = 0;
  subfacetable(2,1) = 3;

  subfacetable(3,0) = 1;
  subfacetable(3,1) = 5;

  subfacetable(4,0) = 5;
  subfacetable(4,1) = 4;

  subfacetable(5,0) = 4;
  subfacetable(5,1) = 1;

  subfacetable(6,0) = 6;
  subfacetable(6,1) = 2;

  subfacetable(7,0) = 2;
  subfacetable(7,1) = 7;

  subfacetable(8,0) = 7;
  subfacetable(8,1) = 6;

  subfacetable(9,0) = 5;
  subfacetable(9,1) = 6;

  subfacetable(10,0) = 6;
  subfacetable(10,1) = 9;

  subfacetable(11,0) = 9;
  subfacetable(11,1) = 5;

  subfacetable(12,0) = 9;
  subfacetable(12,1) = 7;

  subfacetable(13,0) = 7;
  subfacetable(13,1) = 8;

  subfacetable(14,0) = 8;
  subfacetable(14,1) = 9;

  subfacetable(15,0) = 4;
  subfacetable(15,1) = 9;

  subfacetable(16,0) = 9;
  subfacetable(16,1) = 3;

  subfacetable(17,0) = 3;
  subfacetable(17,1) = 4;

  faceflux.resize(subfacetable.nbRows()); // one flux per sub face

  for (CFuint i = 0; i < faceflux.size(); ++i)
    faceflux[i].resize(nbEqs);


//   const CFreal s  = std::sqrt( 0.6 );
//   const CFreal a0 = ( 1.0 - s )*0.5;
//   const CFreal a1 = ( 1.0 + s )*0.5;

  ///change
  const CFreal a0 = 0.9061798459;
  const CFreal a1 = 0.5384693101;

  /// /&/

  xi_ref.resize(10);
  eta_ref.resize(10);

  xi_ref[0] = 0.0;
  xi_ref[1] = 1.0;
  xi_ref[2] = 0.0;
  xi_ref[3] = 1.0/3.0;
  xi_ref[4] = 2.0/3.0;
  xi_ref[5] = 2.0/3.0;
  xi_ref[6] = 1.0/3.0;
  xi_ref[7] = 0.0;
  xi_ref[8] = 0.0;
  xi_ref[9] = 1.0/3.0;

  eta_ref[0] = 0.0;
  eta_ref[1] = 0.0;
  eta_ref[2] = 1.0;
  eta_ref[3] = 0.0;
  eta_ref[4] = 0.0;
  eta_ref[5] = 1.0/3.0;
  eta_ref[6] = 2.0/3.0;
  eta_ref[7] = 2.0/3.0;
  eta_ref[8] = 1.0/3.0;
  eta_ref[9] = 1.0/3.0;


//   qd0.resize(nbQdPts); // 3 quadrature points per face
//   qd1.resize(nbQdPts); // 3 quadrature points per face
// 
//   qd0[0] = a0;  qd1[0] = a1;
//   qd0[1] = a1;  qd1[1] = a0;
//   qd0[2] = .5;  qd1[2] = .5;
// 
//   wqd.resize(nbQdPts); // 3 quadrature points per face
//   //The weights in reference interval <-1,1> are: 5/9, 8/9 and 5/9
//   //We integrate in interval of length 1/3 (= length of one subface of P3P3 triangle)
//   //Therefore we have to divide the weights by a factor 2/(1/3) = 6
//   //Finally, the numerical integral should be multiplied by the length of reference interval
//   //Including this factor in the integration weights, we have:
// 
//   wqd[0] = 5.0/54.0;
//   wqd[1] = 5.0/54.0;
//   wqd[2] = 8.0/54.0;
// 
//   qdstates.resize(3); // 3 quadrature points per face
//   qdstates[0] = new State();
//   qdstates[1] = new State();
//   qdstates[2] = new State();


  ///change
  qd0.resize(nbQdPts); // 5 quadrature points per face
  qd1.resize(nbQdPts); // 5 quadrature points per face

  qd0[0] = (1.0+a0)*0.5;  qd1[0] = (1.0-a0)*0.5;
  qd0[1] = (1.0+a1)*0.5;  qd1[1] = (1.0-a1)*0.5;
  qd0[2] = 0.5;           qd1[2] = 0.5;
  qd0[3] = (1.0-a1)*0.5;  qd1[3] = (1.0+a1)*0.5;
  qd0[4] = (1.0-a0)*0.5;  qd1[4] = (1.0+a0)*0.5;

  wqd.resize(nbQdPts); 

  const CFreal oneSixth = 1.0/6.0;

  wqd[0] = 0.2369268850*oneSixth;
  wqd[1] = 0.4786286705*oneSixth;
  wqd[2] = 0.5688888889*oneSixth;
  wqd[3] = 0.4786286705*oneSixth;
  wqd[4] = 0.2369268850*oneSixth;



  qdstates.resize(nbQdPts); // 5 quadrature points per face
  qdstates[0] = new State();
  qdstates[1] = new State();
  qdstates[2] = new State();
  qdstates[3] = new State();
  qdstates[4] = new State();

  /// /&/


  const CFuint extra_var_size = getMethodData().getUpdateVar()->getExtraPhysicalVarsSize();
  m_qdExtraVars.resize(qdstates.size());
  for (CFuint i = 0; i < m_qdExtraVars.size(); ++i) {
    m_qdExtraVars[i] = new RealVector(extra_var_size);
  }

//   qdnodes.resize(3); // 3 quadrature points per face
//   qdnodes[0] = new Node();
//   qdstates[0]->setSpaceCoordinates(qdnodes[0]);
//   qdnodes[1] = new Node();
//   qdstates[1]->setSpaceCoordinates(qdnodes[1]);
//   qdnodes[2] = new Node();
//   qdstates[2]->setSpaceCoordinates(qdnodes[2]);

  ///change
  qdnodes.resize(nbQdPts); // 5 quadrature points per face

  ///Allocate variables for storing quadrature data
  for(CFuint iQdPt=0;iQdPt<nbQdPts;++iQdPt) {
    qdnodes[iQdPt] = new Node();
    qdstates[iQdPt]->setSpaceCoordinates(qdnodes[iQdPt]);
  }

//   qdnodes[0] = new Node();
//   qdstates[0]->setSpaceCoordinates(qdnodes[0]);
//   qdnodes[1] = new Node();
//   qdstates[1]->setSpaceCoordinates(qdnodes[1]);
//   qdnodes[2] = new Node();
//   qdstates[2]->setSpaceCoordinates(qdnodes[2]);
//   qdnodes[3] = new Node();
//   qdstates[3]->setSpaceCoordinates(qdnodes[3]);
//   qdnodes[4] = new Node();
//   qdstates[4]->setSpaceCoordinates(qdnodes[4]);

  /// /&/


  facenormal.resize(2); // only 2D

  // resize the beta's storages in the MethodData
  const CFuint nbSubTriangles = 9;
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

void HOCRD_SplitStrategyIsoP3::computeFluctuation(vector<RealVector>& residual)
{ 
  
  DistributionData& distdata = getMethodData().getDistributionData();

  vector<State*>& states = *distdata.states;


  cf_assert (states.size() == 10); // P2 triangles for solution space
  cf_assert (residual.size() >= states.size()); // state residual

  // reset the residual because we will accumulate the sub element contributions
  for (CFuint i = 0; i < residual.size(); ++i) { residual[i] = 0.0; }

  InwardNormalsData& cellnormals = (*socket_normals.getDataHandle() [distdata.cellID]);

  cf_assert (cellnormals.nbFaces()  == 18); // P3 triangles have 18 normals as for now


  (*m_phisubT[0]) = 0.;
  (*m_phisubT[1]) = 0.;
  (*m_phisubT[2]) = 0.;
  (*m_phisubT[3]) = 0.;
  (*m_phisubT[4]) = 0.;
  (*m_phisubT[5]) = 0.;
  (*m_phisubT[6]) = 0.;
  (*m_phisubT[7]) = 0.;
  (*m_phisubT[8]) = 0.;

  computeHOCurvedFluctuation();

  /*****         Triangle 0-3-8          *****/
  // faces 2-0-1

  for(CFuint isubelem=0;isubelem<9;++isubelem) {

  //Face indexes of this subtriangle:
  const CFuint f0 = subelemtable(isubelem,2);
  const CFuint f1 = subelemtable(isubelem,0);
  const CFuint f2 = subelemtable(isubelem,1);

  //Subface vertex: read the face as it is stored in subfacetable, 
  //or in reverse order (starting with local node 1) in case of negative face orientation?
  if(subelemfacedir(isubelem,2) > 0) m_subfaceVertex = 0; else m_subfaceVertex = 1;

  //Vertices of this subtriangle:
  const CFuint v0 = subfacetable(f0,m_subfaceVertex);
  const CFuint v1 = subfacetable(f1,m_subfaceVertex);
  const CFuint v2 = subfacetable(f2,m_subfaceVertex);

//   CFout << "Subelement " << isubelem << "\n";
//   CFout << "Vertices: " << v0 << " , " << v1 << " , " << v2 << "\n";
//   CFout << "Faces: " << f0 << " , " << f1 << " , " << f2 << "\n\n";

  //States belonging to subelement:
  substates[0] = states[v0];
  substates[1] = states[v1];
  substates[2] = states[v2];

  // stupidily put the normals data into a matrix to put it in the normals data
  matrix_face_norms (0,XX) = subelemfacedir(isubelem,2)*cellnormals.getFaceNormComp(f0,XX);
  matrix_face_norms (1,XX) = subelemfacedir(isubelem,0)*cellnormals.getFaceNormComp(f1,XX);
  matrix_face_norms (2,XX) = subelemfacedir(isubelem,1)*cellnormals.getFaceNormComp(f2,XX);

  matrix_face_norms (0,YY) = subelemfacedir(isubelem,2)*cellnormals.getFaceNormComp(f0,YY);
  matrix_face_norms (1,YY) = subelemfacedir(isubelem,0)*cellnormals.getFaceNormComp(f1,YY);;
  matrix_face_norms (2,YY) = subelemfacedir(isubelem,1)*cellnormals.getFaceNormComp(f2,YY);


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

  distdata.tStates = computeConsistentStates(&substates);
  distdata.subStates = &substates;

 
  m_splitter->computeK(substates, m_subcell_normals);

  // transform fluxes of subelement to distribution variables
  SafePtr<RealVector> phi = &distdata.phi;
  *phi = *getMethodData().getSolutionToDistribMatTrans()->transformFromRef(m_phisubT[isubelem]);

  SafePtr<vector<RealMatrix> >& currBetaMatrix = distdata.currBetaMat;

  currBetaMatrix = &distdata.betaMats[isubelem];

  cf_assert (currBetaMatrix.isNotNull());
  m_splitter->distribute(subresidual);

  residual[v0] += subresidual[0];
  residual[v1] += subresidual[1];
  residual[v2] += subresidual[2];


  } //Loop over subelements

//   CF_DEBUG_EXIT;

}

//////////////////////////////////////////////////////////////////////////////

void HOCRD_SplitStrategyIsoP3::computeHOCurvedFluctuation()
{

  cf_assert (qdstates.size() == nbQdPts); // only triags so five quadrature points per face
 
  vector<State*>& states = *getMethodData().getDistributionData().states;

  State& state0 = *(states[0]);
  State& state1 = *(states[1]);
  State& state2 = *(states[2]);
  State& state3 = *(states[3]);
  State& state4 = *(states[4]);
  State& state5 = *(states[5]);
  State& state6 = *(states[6]);
  State& state7 = *(states[7]);
  State& state8 = *(states[8]);
  State& state9 = *(states[9]);

  vector<Node*>& nodes = *getMethodData().getDistributionData().cell->getNodes();


  cf_assert (nodes.size()  == 10); // P3 triangles for geometry space

  const CFreal x0 = (*nodes[0])[XX];
  const CFreal x1 = (*nodes[1])[XX];
  const CFreal x2 = (*nodes[2])[XX];

  const CFreal x3 = (*nodes[3])[XX];
  const CFreal x4 = (*nodes[4])[XX];
  const CFreal x5 = (*nodes[5])[XX];
  const CFreal x6 = (*nodes[6])[XX];
  const CFreal x7 = (*nodes[7])[XX];
  const CFreal x8 = (*nodes[8])[XX];

  const CFreal x9 = (*nodes[9])[XX];


  const CFreal y0 = (*nodes[0])[YY];
  const CFreal y1 = (*nodes[1])[YY];
  const CFreal y2 = (*nodes[2])[YY];

  const CFreal y3 = (*nodes[3])[YY];
  const CFreal y4 = (*nodes[4])[YY];
  const CFreal y5 = (*nodes[5])[YY];
  const CFreal y6 = (*nodes[6])[YY];
  const CFreal y7 = (*nodes[7])[YY];
  const CFreal y8 = (*nodes[8])[YY];

  const CFreal y9 = (*nodes[9])[YY];

  InwardNormalsData& cellnormals = (*socket_normals.getDataHandle() [getMethodData().getDistributionData().cellID]);

 // Computation  of the fluctuation of each faces of the cell
for (CFuint iFace = 0; iFace < subfacetable.nbRows(); ++iFace)
{
    facenormal[XX] = - cellnormals.getFaceNormComp(iFace,XX);
    facenormal[YY] = - cellnormals.getFaceNormComp(iFace,YY);

    ///change
    for (CFuint iQd = 0; iQd < nbQdPts; ++iQd)
    {
      /// Position of integration point in reference space
      const CFreal xi  = qd0[iQd] * xi_ref[subfacetable(iFace,0)] + qd1[iQd] * xi_ref[subfacetable(iFace,1)];
      const CFreal eta = qd0[iQd] * eta_ref[subfacetable(iFace,0)] + qd1[iQd] * eta_ref[subfacetable(iFace,1)];

      /// Values of shape functions in reference space:
      const CFreal L0 = 1.0-xi-eta;
      const CFreal L1 = xi;
      const CFreal L2 = eta;

      const CFreal SF0 = 0.5 * (3.0*L0-1.0) * (3.0*L0-2.0) * L0;
      const CFreal SF1 = 0.5 * (3.0*L1-1.0) * (3.0*L1-2.0) * L1;
      const CFreal SF2 = 0.5 * (3.0*L2-1.0) * (3.0*L2-2.0) * L2;
      const CFreal SF3 = 4.5 * L0 * L1 * (3.0*L0-1.0);
      const CFreal SF4 = 4.5 * L0 * L1 * (3.0*L1-1.0);
      const CFreal SF5 = 4.5 * L1 * L2 * (3.0*L1-1.0);
      const CFreal SF6 = 4.5 * L1 * L2 * (3.0*L2-1.0);
      const CFreal SF7 = 4.5 * L0 * L2 * (3.0*L2-1.0);
      const CFreal SF8 = 4.5 * L0 * L2 * (3.0*L0-1.0);
      const CFreal SF9 = 27.0* L0 * L1 * L2;

      /// Compute the position of the integration point in 'physical' space
      (*qdnodes[iQd])[XX] =  SF0 * x0 + SF1 * x1 + SF2 * x2 + SF3 * x3 + SF4 * x4 + SF5 * x5 + SF6 * x6 + SF7 * x7 + SF8 * x8 + SF9 * x9;
      (*qdnodes[iQd])[YY] =  SF0 * y0 + SF1 * y1 + SF2 * y2 + SF3 * y3 + SF4 * y4 + SF5 * y5 + SF6 * y6 + SF7 * y7 + SF8 * y8 + SF9 * y9;

      /// Values of shape functions in physical space:
      (*qdstates[iQd]) = SF0 * state0 + SF1 * state1 + SF2 * state2 + SF3 * state3 + SF4 * state4 + SF5 * state5 + \
                         SF6 * state6 + SF7 * state7 + SF8 * state8 + SF9 * state9;

      /// Compute normals at the integration points
      m_CP3N.ComputeNormal(nodes,iFace,xi,eta,m_scaledFaceNormals[iQd]);
      m_scaledFaceNormals[iQd] *= -1.0;

    }

    //change
    computeStatesData(nbQdPts, m_updateVar, qdstates, m_pdata, m_qdExtraVars);
    
    faceflux[iFace] = 0.;

    ///change
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
HOCRD_SplitStrategyIsoP3::needsSockets()
{
  std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > result
    = FluctuationSplitStrategy::needsSockets();
  result.push_back(&socket_updateCoeff);
  return result;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD
