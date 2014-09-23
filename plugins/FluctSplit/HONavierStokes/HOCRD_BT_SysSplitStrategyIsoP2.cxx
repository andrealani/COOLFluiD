#include "FluctSplit/HONavierStokes/HOCRD_BT_SysSplitStrategyIsoP2.hh"

#include "Common/BadValueException.hh"
#include "Framework/MethodStrategyProvider.hh"
#include "Framework/ContourIntegrator.hh"
#include "Framework/MeshData.hh"
#include "Framework/BaseTerm.hh"
#include "FluctSplit/HONavierStokes/FluctSplitHONavierStokes.hh"
#include "MathTools/MatrixInverter.hh"

#include "NavierStokes/EulerTerm.hh"
#include "Framework/SubSystemStatus.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Physics::NavierStokes;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

MethodStrategyProvider<HOCRD_BT_SysSplitStrategyIsoP2,
                       FluctuationSplitData,
                       FluctuationSplitStrategy,
                       FluctSplitHONavierStokesModule>
aHOCRDBTFluctSplitIsoP2StrategyProvider("HOCRD_BT_IsoP2");

//////////////////////////////////////////////////////////////////////////////

void HOCRD_BT_SysSplitStrategyIsoP2::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< CFreal>("Delta","Delta of variable.");
  options.addConfigOption< CFreal>("Length","Reference Length.");
  options.addConfigOption< CFreal>("Speed","Reference Speed.");
  options.addConfigOption< std::string>("VarName","Variable name.");
  options.addConfigOption< bool >("HO","High order discretization");
  options.addConfigOption< CFuint >("MaxNbSubElems","Maximum number of sub-elements for high-order computations.");
  options.addConfigOption< bool >  ("StoreThetas","Store the thetas for visualization");
  options.addConfigOption< bool >  ("UmaxTheta","Use the maximum of the thetas");
  options.addConfigOption< std::string>("Shockdetector","Which shock detetecto to use");
}

//////////////////////////////////////////////////////////////////////////////

HOCRD_BT_SysSplitStrategyIsoP2::HOCRD_BT_SysSplitStrategyIsoP2(const std::string& name) :
  FluctuationSplitStrategy(name),
  socket_updateCoeff("updateCoeff"),
  socket_thetas("thetas",false),
  m_solutionVar(CFNULL),
  m_updateVar(CFNULL),
  m_unitFaceNormals(),
  m_scaledFaceNormals(),
   m_phisubT(),
  m_phiN1(0),
  m_phiN2(0),
  m_phiN3(0),
  m_phiN4(0),
  m_phi(0),
  subelemfacedir(0,0),
  subelemtable(0,0),
  subfacetable(0,0),
  substates(),
  subresidual(),
  faceflux(),
  qd0(),
  qd1(),
  wqd(),
  qdstates(),
  facenormal(),
  m_sumKplusU(),
  m_sumKU(),
  m_sumKplus(),
  m_uTemp(),
  m_uMin(),
  m_tmp(),
  m_k1Plus(0),
  m_k2Plus(0),
  m_k3Plus(0),
  m_k4Plus(0),
  m_k(0),
  m_kMin(0),
  m_eValues(0),
  m_adimNormal(),
  m_theta1(),
  m_theta2(),
  m_theta3(),
  m_theta4(),
  m_theta(),
  m_inverter(CFNULL),
  m_invK(),
  m_uInflow(),
  matrix_face_norms(3,2),
  matrix_node_norms(3,2),
  vector_face_areas(3),
  vector_node_areas(3),
  _pData()
{
  addConfigOptionsTo(this);

_deltaP = 0.0;
  setParameter("Delta",&_deltaP);

  _length = 1.0;
  setParameter("Length",&_length);

  _speed = 0.0;
  setParameter("Speed",&_speed);

  _varName = "p";
  setParameter("VarName",&_varName);

  m_max_nbsubcells = 1;
  setParameter("MaxNbSubElems",&m_max_nbsubcells);

  m_store_thetas = false;
  setParameter("StoreThetas",&m_store_thetas);

  m_use_max_theta = true;
  setParameter("UmaxTheta",&m_use_max_theta);

  _isHO = false;
  this->setParameter("HO",&_isHO);

   
  _sh_detector = "Jirka";
  // Choice between Jirka shock capturing and 
  // The improved one of Antonino, the one of Antonino
  // is implemented only in 2D
  this->setParameter("Shockdetector",&_sh_detector);

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

HOCRD_BT_SysSplitStrategyIsoP2::~HOCRD_BT_SysSplitStrategyIsoP2()
{
  if (isSetup()) unsetup();
}

//////////////////////////////////////////////////////////////////////////////

void HOCRD_BT_SysSplitStrategyIsoP2::unsetup()
{
  for (CFuint i = 0; i < m_k1Plus.size(); ++i)
  {
    deletePtr(m_k1Plus[i]);
    deletePtr(m_k2Plus[i]);
    deletePtr(m_k3Plus[i]);
    deletePtr(m_k4Plus[i]);
    deletePtr(m_k[i]);
    deletePtr(m_kMin[i]);
  }
  deletePtr(m_inverter);

  deletePtr(m_subcell_normals);

  // AL: this could be changed ...
  // where are deleted the qdstates ???
  for (CFuint i = 0; i < m_qdExtraVars.size(); ++i) {
    deletePtr(m_qdExtraVars[i]);
  }
}

//////////////////////////////////////////////////////////////////////////////

void HOCRD_BT_SysSplitStrategyIsoP2::setup()
{
  CFAUTOTRACE;

  // first call parent method
  FluctuationSplitStrategy::setup();

  DataHandle< CFreal> updateCoeff = socket_updateCoeff.getDataHandle();

  if (getMethodData().isMultipleSplitter())
    throw Common::BadValueException(FromHere(), "Cannot use HOCRD on curved elements with multiple splitters");
  m_solutionVar = getMethodData().getSolutionVar();
  m_updateVar   = getMethodData().getUpdateVar();

  m_unitFaceNormals.resize(MeshDataStack::getActive()->Statistics().getMaxNbFacesInCell());
  for (CFuint  i = 0; i < m_unitFaceNormals.size(); ++i) {
    m_unitFaceNormals[i].resize(PhysicalModelStack::getActive()->getDim());
  }

   // number of quadrature points used to compute the fluctuation
  const CFuint nbQdPts = 3;


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

  m_phisubT.resize(4); // P2 triangle has 4 sub triangles

  m_phisubT[0] = new RealVector(nbEqs);
  m_phisubT[1] = new RealVector(nbEqs);
  m_phisubT[2] = new RealVector(nbEqs);
  m_phisubT[3] = new RealVector(nbEqs);

  m_phiN1.resize(3);

  m_phiN1[0].resize(nbEqs);
  m_phiN1[1].resize(nbEqs);
  m_phiN1[2].resize(nbEqs);


  m_phiN2.resize(3);

  m_phiN2[0].resize(nbEqs);
  m_phiN2[1].resize(nbEqs);
  m_phiN2[2].resize(nbEqs);

  m_phiN3.resize(3);

  m_phiN3[0].resize(nbEqs);
  m_phiN3[1].resize(nbEqs);
  m_phiN3[2].resize(nbEqs);

  m_phiN4.resize(3);

  m_phiN4[0].resize(nbEqs);
  m_phiN4[1].resize(nbEqs);
  m_phiN4[2].resize(nbEqs);


 m_phi.resize(3);

  m_phi[0].resize(nbEqs);
  m_phi[1].resize(nbEqs);
  m_phi[2].resize(nbEqs);


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

  const CFreal s  = std::sqrt( 0.6 );
  const CFreal a0 = ( 1.0 - s )*0.5;
  const CFreal a1 = ( 1.0 + s )*0.5;

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

  qd0[0] = a0;  qd1[0] = a1;
  qd0[1] = a1;  qd1[1] = a0;
  qd0[2] = .5;  qd1[2] = .5;

  wqd.resize(nbQdPts); // 3 quadrature points per face
  wqd[0] = 5.0/18.0;
  wqd[1] = 5.0/18.0;
  wqd[2] = 8.0/18.0;

  qdstates.resize(3); // 3 quadrature points per face
  qdstates[0] = new State();
  qdstates[1] = new State();
  qdstates[2] = new State();

  const CFuint extra_var_size = getMethodData().getUpdateVar()->getExtraPhysicalVarsSize();
  m_qdExtraVars.resize(qdstates.size());
  for (CFuint i = 0; i < m_qdExtraVars.size(); ++i) {
    m_qdExtraVars[i] = new RealVector(extra_var_size);
  }

  qdnodes.resize(3); // 3 quadrature points per face
  qdnodes[0] = new Node();
  qdstates[0]->setSpaceCoordinates(qdnodes[0]);
  qdnodes[1] = new Node();
  qdstates[1]->setSpaceCoordinates(qdnodes[1]);
  qdnodes[2] = new Node();
  qdstates[2]->setSpaceCoordinates(qdnodes[2]);

  facenormal.resize(2); // only 2D

  m_maxNbStatesInCell = MeshDataStack::getActive()->Statistics().getMaxNbStatesInCell();

  m_k1Plus.resize(m_maxNbStatesInCell);
  m_k2Plus.resize(m_maxNbStatesInCell);
  m_k3Plus.resize(m_maxNbStatesInCell);
  m_k4Plus.resize(m_maxNbStatesInCell);
  m_k.resize(m_maxNbStatesInCell);
  m_kMin.resize(m_maxNbStatesInCell);
  m_eValues.resize(m_maxNbStatesInCell);

  for (CFuint i = 0; i < m_maxNbStatesInCell; ++i) {
    m_k1Plus[i] = new RealMatrix(nbEqs, nbEqs);
    m_k2Plus[i] = new RealMatrix(nbEqs, nbEqs);
    m_k3Plus[i] = new RealMatrix(nbEqs,nbEqs );
    m_k4Plus[i] = new RealMatrix(nbEqs, nbEqs);
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
  if (_varName == "rho") {
    _varID = Physics::NavierStokes::EulerTerm::RHO;
  }
  else {
    _varID = Physics::NavierStokes::EulerTerm::P;
  }

  if ( m_store_thetas && !socket_thetas.isConnected())
    throw Common::BadValueException (FromHere(),"User required storing of thetas but socket is not connected");



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

void HOCRD_BT_SysSplitStrategyIsoP2::computeFluctuation(vector<RealVector>& residual)
{
  DistributionData& distdata = getMethodData().getDistributionData();

  vector<State*>& states = *distdata.states;


   const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();

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

  //First we compute the residual of N scheme and the thetas on each sub element

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

//  m_splitter->computeK(substates, m_subcell_normals);
// compute the residual and the upwind parameters k in this cell
  computeK(substates,m_k1Plus);
  // transform fluxes of subelement to distribution variables
  SafePtr<RealVector> phi = &distdata.phi;
  *phi = *getMethodData().getSolutionToDistribMatTrans()->transformFromRef(m_phisubT[0]);

  distributeN(m_k1Plus,m_phiN1);

  computeBlendingCoeff( m_theta1);

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

 // m_splitter->computeK(substates, m_subcell_normals);

  // compute the residual and the upwind parameters k in this cell
  computeK(substates,m_k2Plus);

  // transform fluxes of subelement to distribution variables
  *phi = *getMethodData().getSolutionToDistribMatTrans()->transformFromRef(m_phisubT[1]);

  distributeN(m_k2Plus,m_phiN2);
computeBlendingCoeff(m_theta2);

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

 // m_splitter->computeK(substates, m_subcell_normals);
 // compute the residual and the upwind parameters k in this cell
  computeK(substates,m_k3Plus);

  // transform fluxes of subelement to distribution variables
  *phi = *getMethodData().getSolutionToDistribMatTrans()->transformFromRef(m_phisubT[2]);

  distributeN(m_k3Plus,m_phiN3);
computeBlendingCoeff(m_theta3);

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

//  m_splitter->computeK(substates, m_subcell_normals);
// compute the residual and the upwind parameters k in this cell
  computeK(substates,m_k4Plus);

  // transform fluxes of subelement to distribution variables
  *phi = *getMethodData().getSolutionToDistribMatTrans()->transformFromRef(m_phisubT[3]);

  distributeN(m_k4Plus,m_phiN4);
  computeBlendingCoeff(m_theta4);
// Then we take the max of the thetas
  // this is outside because is used for the artificial viscosity
//   m_theta = max( max( max(m_theta1, m_theta2), m_theta3), m_theta4);
//   if (m_use_max_theta)
//   {
//     m_theta1 = max(m_min_theta, m_theta);
//     m_theta2 = m_theta1;
//     m_theta3 = m_theta1;
//     m_theta4 = m_theta1;
//   }
//   else
//   {
//     m_theta1 = max(m_min_theta, m_theta1);
//     m_theta2 = max(m_min_theta, m_theta2);
//     m_theta3 = max(m_min_theta, m_theta3);
//     m_theta4 = max(m_min_theta, m_theta4);
//   }

  if (m_store_thetas && !distdata.isPerturb)
  {

      Framework::DataHandle< CFreal > thetas = socket_thetas.getDataHandle();
      for (CFuint iEq = 0; iEq < nbEqs; ++iEq)
      {
        thetas( distdata.cellID*m_max_nbsubcells + 0, iEq, nbEqs) = m_theta1;
        thetas( distdata.cellID*m_max_nbsubcells + 1, iEq, nbEqs) = m_theta2;
        thetas( distdata.cellID*m_max_nbsubcells + 2, iEq, nbEqs) = m_theta3;
        thetas( distdata.cellID*m_max_nbsubcells + 3, iEq, nbEqs) = m_theta4;
      }
    
  }

//Then we use this theta to build the blending


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

  //m_splitter->computeK(substates, m_subcell_normals);

  // transform fluxes of subelement to distribution variables
   *phi = *getMethodData().getSolutionToDistribMatTrans()->transformFromRef(m_phisubT[0]);

  SafePtr<vector<RealMatrix> >& currBetaMatrix = distdata.currBetaMat;

  currBetaMatrix = &distdata.betaMats[0];

  cf_assert (currBetaMatrix.isNotNull());
  distributeLDA(m_k1Plus,m_phi);

  residual[0] += m_theta1*m_phiN1[0] + (1.0 - m_theta1)*m_phi[0];
  residual[3] += m_theta1*m_phiN1[1] + (1.0 - m_theta1)*m_phi[1];
  residual[5] += m_theta1*m_phiN1[2] + (1.0 - m_theta1)*m_phi[2];

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

 // m_splitter->computeK(substates, m_subcell_normals);

  // transform fluxes of subelement to distribution variables
  *phi = *getMethodData().getSolutionToDistribMatTrans()->transformFromRef(m_phisubT[1]);

  currBetaMatrix = &distdata.betaMats[1];
  cf_assert (currBetaMatrix.isNotNull());

  distributeLDA(m_k2Plus,m_phi);


  residual[3] += m_theta2*m_phiN2[0] + (1.0 - m_theta2)*m_phi[0];
  residual[1] += m_theta2*m_phiN2[1] + (1.0 - m_theta2)*m_phi[1];
  residual[4] += m_theta2*m_phiN2[2] + (1.0 - m_theta2)*m_phi[2];

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

  //m_splitter->computeK(substates, m_subcell_normals);

  // transform fluxes of subelement to distribution variables
  *phi = *getMethodData().getSolutionToDistribMatTrans()->transformFromRef(m_phisubT[2]);

  currBetaMatrix = &distdata.betaMats[2];
  cf_assert (currBetaMatrix.isNotNull());
  distributeLDA(m_k3Plus,m_phi);

  residual[5] += m_theta3*m_phiN3[0] + (1.0 - m_theta3)*m_phi[0];
  residual[4] += m_theta3*m_phiN3[1] + (1.0 - m_theta3)*m_phi[1];
  residual[2] += m_theta3*m_phiN3[2] + (1.0 - m_theta3)*m_phi[2];

   /*****         Triangle 4-5-3          *****/

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

  //m_splitter->computeK(substates, m_subcell_normals);

  // transform fluxes of subelement to distribution variables
  *phi = *getMethodData().getSolutionToDistribMatTrans()->transformFromRef(m_phisubT[3]);

  currBetaMatrix = &distdata.betaMats[3];
  cf_assert (currBetaMatrix.isNotNull());

  distributeLDA(m_k4Plus,m_phi);

  residual[4] += m_theta4*m_phiN4[0] + (1.0 - m_theta4)*m_phi[0];
  residual[5] += m_theta4*m_phiN4[1] + (1.0 - m_theta4)*m_phi[1];
  residual[3] += m_theta4*m_phiN4[2] + (1.0 - m_theta4)*m_phi[2]; 


}

//////////////////////////////////////////////////////////////////////////////

void HOCRD_BT_SysSplitStrategyIsoP2::computeHOCurvedFluctuation()
{

  cf_assert (qdstates.size() == 3); // only triags so three quadrature points per face

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


    for (CFuint iQd = 0; iQd < 3; ++iQd)
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
    
    computeStatesData(3, m_updateVar, qdstates, m_pdata, m_qdExtraVars); // three quadrature points per face

    faceflux[iFace] = 0.;
    for (CFuint iQd = 0; iQd < 3; ++iQd)
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
//////////////////////////////////////////////////////////////////////////////////////
void HOCRD_BT_SysSplitStrategyIsoP2::computeK(const std::vector<Framework::State*>& states,
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
void HOCRD_BT_SysSplitStrategyIsoP2::distributeN(std::vector<RealMatrix*> & m_kPlus,vector<RealVector> & phiN){

  const vector<State*>& tStates = *getMethodData().getDistributionData().tStates;
  const RealVector& phiT = getMethodData().getDistributionData().phi;


  m_tmp = (*m_kPlus[0])  + (*m_kMin[0]);

  //  for (CFuint i=0; i < tmp.size(); ++i) {
  //     sumAbsK[i] = std::abs(tmp[i]);
  //   }

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
void HOCRD_BT_SysSplitStrategyIsoP2::distributeLDA(std::vector<RealMatrix*> & m_kPlus,vector<RealVector> & phiLDA){

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
void HOCRD_BT_SysSplitStrategyIsoP2::computeBlendingCoeff(CFreal & result)
{
  vector<State*> states = substates;

  CFreal vol = getMethodData().getDistributionData().cell->computeVolume();

  if (_isHO) vol /= 4.0;


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
  if (_sh_detector == "Jirka"){
    theta = _grad[XX]*lData[EulerTerm::VX] + _grad[YY]*lData[EulerTerm::VY];
    if (dim == DIM_3D) {
      theta += _grad[ZZ]*lData[EulerTerm::VZ];
    } 
  
//    updateScalingValues();
  
    theta *= _length/(vol*dim*_deltaP*_speed);
  }

  else if (_sh_detector == "Anton"){
     theta = _grad[XX] * lData[EulerTerm::VX] + _grad[YY] * lData[EulerTerm::VY];

     theta *= std::pow(    (std::pow(lData[EulerTerm::VX],2) + std::pow(lData[EulerTerm::VY],2)), 0.5   );

//      updateScalingValues();

     theta *= _length/(vol*dim*_deltaP * _speed * _speed);   //modified!!


  }


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

std::vector<Common::SafePtr<Framework::BaseDataSocketSink> >
HOCRD_BT_SysSplitStrategyIsoP2::needsSockets()
{
  std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > result
    = FluctuationSplitStrategy::needsSockets();
   result.push_back(&socket_updateCoeff);
   result.push_back(&socket_thetas);
   return result;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD
