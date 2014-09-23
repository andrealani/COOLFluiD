#include "Common/BadValueException.hh"
#include "Framework/MethodStrategyProvider.hh"
#include "Framework/ContourIntegrator.hh"
#include "Framework/MeshData.hh"
#include "Framework/BaseTerm.hh"
#include "FluctSplit/HONavierStokes/FluctSplitHONavierStokes.hh"
#include "MathTools/MatrixInverter.hh"

#include "NavierStokes/EulerTerm.hh"

#include "FluctSplit/HONavierStokes/HOCRD_BT_SysSplitStrategy.hh"

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

MethodStrategyProvider<HOCRD_BT_SysSplitStrategy,
                       FluctuationSplitData,
                       FluctuationSplitStrategy,
                       FluctSplitHONavierStokesModule>
hocrdbtsysFluctSplitStrategyProvider("HOCRD_SysBT");

//////////////////////////////////////////////////////////////////////////////

void HOCRD_BT_SysSplitStrategy::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< CFreal>("Delta","Delta of variable.");
  options.addConfigOption< CFreal>("Length","Reference Length.");
  options.addConfigOption< CFreal>("Speed","Reference Speed.");
  options.addConfigOption< std::string>("VarName","Variable name.");
  options.addConfigOption< bool >("HO","High order discretization");
  options.addConfigOption< CFuint >("MaxNbSubElems","Maximum number of sub-elements for high-order computations.");
  options.addConfigOption< bool >  ("StoreThetas","Store the thetas for visualization");
  options.addConfigOption< bool >  ("UmaxTheta","Use the maximum of the thetas");
  options.addConfigOption< bool >  ("AddArtificialViscosity","Add artificial viscosity for shock monotonicity");
  options.addConfigOption< CFreal >  ("d0","Intial viscosity");
  options.addConfigOption< CFreal, Config::DynamicOption<> >("MinTheta","Minimum theta, used to keep a minimum diffusion");
  options.addConfigOption< std::string>("Shockdetector","Which shock detetecto to use");
}


//////////////////////////////////////////////////////////////////////////////

HOCRD_BT_SysSplitStrategy::HOCRD_BT_SysSplitStrategy(const std::string& name) :
  FluctuationSplitStrategy(name),
  socket_updateCoeff("updateCoeff"),
  socket_thetas("thetas",false),
  m_solutionVar(CFNULL),
  m_updateVar(CFNULL),
  m_unitFaceNormals(),
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
  _pData(),
  _grad(),
  m_min_theta(0.),
   states(),
  _values(),
  m_res_art_visc(0),
  m_d0(),
  _phi_diff_bub(0),
  _phi_diff_gal(0),
  _phi_diff_bub_split(0),
  kappa(0,0),
  F1(),
  F2(),
  F3(),
  m_cellVolume(),
  artvisc_qd0(),
  artvisc_qd1(),
  artvisc_qd2(),
  artvisc_wqd(),
  _normal(),
  m_sc(0.0)
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

  m_add_artvisc = false;
  setParameter("AddArtificialViscosity",&m_add_artvisc);

  m_d0 = 1.0;
  setParameter("d0",&m_d0);

  m_min_theta = 0.;
  setParameter("MinTheta",&m_min_theta);


  _isHO = false;
  this->setParameter("HO",&_isHO);

 // Choice between Jirka shock capturing and 
  // The improved one of Antonino, the one of Antonino
  // is implemented only in 2D
  _sh_detector = "Jirka";
  this->setParameter("Shockdetector",&_sh_detector);
}

//////////////////////////////////////////////////////////////////////////////

HOCRD_BT_SysSplitStrategy::~HOCRD_BT_SysSplitStrategy()
{
  if (isSetup()) unsetup();
}

//////////////////////////////////////////////////////////////////////////////

void HOCRD_BT_SysSplitStrategy::unsetup()
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

  // AL: this could be changed ...
  // where are deleted the qdstates ???
  for (CFuint i = 0; i < m_qdExtraVars.size(); ++i) {
   deletePtr(m_qdExtraVars[i]);
  }
}

//////////////////////////////////////////////////////////////////////////////

void HOCRD_BT_SysSplitStrategy::setup()
{
  CFAUTOTRACE;

  // first call parent method
  FluctuationSplitStrategy::setup();

  DataHandle< CFreal> updateCoeff = socket_updateCoeff.getDataHandle();

  if (getMethodData().isMultipleSplitter())
    throw BadValueException (FromHere(),"Cannot use HOCRD with multiple splitters");

  m_solutionVar = getMethodData().getSolutionVar();
  m_updateVar   = getMethodData().getUpdateVar();

  m_unitFaceNormals.resize(MeshDataStack::getActive()->Statistics().getMaxNbFacesInCell());
  for (CFuint  i = 0; i < m_unitFaceNormals.size(); ++i) {
    m_unitFaceNormals[i].resize(PhysicalModelStack::getActive()->getDim());
  }

   // number of quadrature point used to compute the fluctuatuion
  const CFuint nbQdPts = 3;
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

  m_res_art_visc.resize(6);
  m_res_art_visc[0].resize(nbEqs);
  m_res_art_visc[1].resize(nbEqs);
  m_res_art_visc[2].resize(nbEqs);
  m_res_art_visc[3].resize(nbEqs);
  m_res_art_visc[4].resize(nbEqs);
  m_res_art_visc[5].resize(nbEqs);


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

  //Set up for the artificial viscosity

  _phi_diff_bub = new RealVector(nbEqs);
  _phi_diff_bub_split.resize(3); // 3 residuals in each sub element

  _phi_diff_bub_split[0] = new RealVector(nbEqs);
  _phi_diff_bub_split[1] = new RealVector(nbEqs);
  _phi_diff_bub_split[2] = new RealVector(nbEqs);

  _phi_diff_gal.resize(6); // one galerkin function per state

  for (CFuint iState = 0; iState < 6; ++ iState)
     _phi_diff_gal[iState] = new RealVector(nbEqs);

  // kappa is used to "distribute" the part of the fluctuation with the bubble function
   kappa.resize(4,6); // 4 sub-elements 6 nodes
  CFreal one_fourth = 1.0/4.0;

  kappa(0,0) = 3.0*one_fourth;
  kappa(0,1) = -one_fourth;
  kappa(0,2) = -one_fourth;
  kappa(0,3) = 5.0*one_fourth;
  kappa(0,4) = one_fourth;
  kappa(0,5) = 5.0*one_fourth;

  kappa(1,0) = -one_fourth;
  kappa(1,1) = 3.0*one_fourth;
  kappa(1,2) = -one_fourth;
  kappa(1,3) = 5.0*one_fourth;
  kappa(1,4) = 5.0*one_fourth;
  kappa(1,5) = one_fourth;

  kappa(2,0) = -one_fourth;
  kappa(2,1) = -one_fourth;
  kappa(2,2) = 3.0*one_fourth;
  kappa(2,3) = one_fourth;
  kappa(2,4) = 5.0*one_fourth;
  kappa(2,5) = 5.0*one_fourth;

  kappa(3,0) = -one_fourth;
  kappa(3,1) = -one_fourth;
  kappa(3,2) = -one_fourth;
  kappa(3,3) = 5.0*one_fourth;
  kappa(3,4) = 5.0*one_fourth;
  kappa(3,5) = 5.0*one_fourth;
  // tell the splitters to compute the betas
  getMethodData().getDistributionData().computeBetas = true;

  F1.resize(nbEqs);
  F2.resize(nbEqs);
  F3.resize(nbEqs);


  // Setup : should create a function to do all this setup

  artvisc_qd0.resize(4); // quadrature points per face
  artvisc_qd1.resize(4); // quadrature points per face
  artvisc_qd2.resize(4);

  artvisc_qd0[0] = 1.0/3.0;  artvisc_qd1[0] = 1.0/3.0; artvisc_qd2[0] = 1.0/3.0;
  artvisc_qd0[1] = 0.6;  artvisc_qd1[1] = 0.2; artvisc_qd2[1] = 0.2;
  artvisc_qd0[2] = 0.2;  artvisc_qd1[2] = 0.6; artvisc_qd2[2] = 0.2;
  artvisc_qd0[3] = 0.2;  artvisc_qd1[3] = 0.2; artvisc_qd2[3] = 0.6;

  artvisc_wqd.resize(4); // 4 quadrature points per surface

  artvisc_wqd[0] = -27.0/48.0;
  artvisc_wqd[1] = 25.0/48.0;
  artvisc_wqd[2] = 25.0/48.0;
  artvisc_wqd[3] = 25.0/48.0;
 _normal.resize(PhysicalModelStack::getActive()->getDim());
}

//////////////////////////////////////////////////////////////////////////////

void HOCRD_BT_SysSplitStrategy::computeFluctuation(vector<RealVector>& residual)
{

   const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();

  states = *getMethodData().getDistributionData().states;
  cf_assert(states.size() == 6); // P2 triangles for solution space
  cf_assert(residual.size() >= states.size()); // state residual

  // reset the residual because we will accumulate the sub element contributions
  for (CFuint i = 0; i < residual.size(); ++i)
  {
    residual[i] = 0.0;
  }

  InwardNormalsData& cellnormals = (*socket_normals.getDataHandle()
				    [getMethodData().getDistributionData().cellID]);

  cf_assert(cellnormals.nbFaces()  == 3); // triangles have 3 faces

  // normals are half scale
  cellnormals.scale(0.5);
  computeHOFluctuation();
  //First we compute the residual of N scheme and the thetas on each sub element

  /*****         Triangle 0-3-5          *****/

  substates[0] = states[0];
  substates[1] = states[3];
  substates[2] = states[5];

  getMethodData().getDistributionData().tStates =
    computeConsistentStates(&substates);

  // compute the residual and the upwind parameters k in this cell
  computeK(substates,m_k1Plus);
    // transform fluxes of subelement to distribution variables
  SafePtr<RealVector> phi = &getMethodData().getDistributionData().phi;
  *phi = *getMethodData().getSolutionToDistribMatTrans()->transformFromRef(m_phisubT[0]);

  distributeN(m_k1Plus,m_phiN1);


  computeBlendingCoeff( m_theta1);

   /*****         Triangle 3-1-4          *****/

  substates[0] = states[3];
  substates[1] = states[1];
  substates[2] = states[4];

  getMethodData().getDistributionData().tStates =
    computeConsistentStates(&substates);

  // compute the residual and the upwind parameters k in this cell
  computeK(substates,m_k2Plus);

   // transform fluxes of subelement to distribution variables

  *phi = *getMethodData().getSolutionToDistribMatTrans()->transformFromRef(m_phisubT[1]);

  distributeN(m_k2Plus,m_phiN2);
computeBlendingCoeff(m_theta2);

  /*****         Triangle 5-4-2          *****/

  substates[0] = states[5];
  substates[1] = states[4];
  substates[2] = states[2];

  getMethodData().getDistributionData().tStates =
    computeConsistentStates(&substates);

  // compute the residual and the upwind parameters k in this cell
  computeK(substates,m_k3Plus);

   // transform fluxes of subelement to distribution variables

  *phi = *getMethodData().getSolutionToDistribMatTrans()->transformFromRef(m_phisubT[2]);

  distributeN(m_k3Plus,m_phiN3);
computeBlendingCoeff(m_theta3);


 /*****         Triangle 4-5-3          *****/

  substates[0] = states[4];
  substates[1] = states[5];
  substates[2] = states[3];

  getMethodData().getDistributionData().tStates = computeConsistentStates(&substates);

  // compute the residual and the upwind parameters k in this cell
  // the oriantation of the normal in this sub-element is oposite to the one of the element
  cellnormals.scale(-0.5);

 // compute the residual and the upwind parameters k in this cell
  computeK(substates,m_k4Plus);
   // transform fluxes of subelement to distribution variables

  *phi = *getMethodData().getSolutionToDistribMatTrans()->transformFromRef(m_phisubT[3]);

  distributeN(m_k4Plus,m_phiN4);
  computeBlendingCoeff(m_theta4);

  // Then we take the max of the thetas
  // this is outside because is used for the artificial viscosity
  m_theta = max( max( max(m_theta1, m_theta2), m_theta3), m_theta4);
  if (m_use_max_theta)
  {
    m_theta1 = max(m_min_theta, m_theta);
    m_theta2 = m_theta1;
    m_theta3 = m_theta1;
    m_theta4 = m_theta1;
  }
  else
  {
    m_theta1 = max(m_min_theta, m_theta1);
    m_theta2 = max(m_min_theta, m_theta2);
    m_theta3 = max(m_min_theta, m_theta3);
    m_theta4 = max(m_min_theta, m_theta4);
  }

  DistributionData& distdata = getMethodData().getDistributionData();
  if (m_store_thetas && !distdata.isPerturb)
  {
    for (CFuint iSubCell = 0; iSubCell < 4; ++iSubCell)
    {

      Framework::DataHandle< CFreal > thetas = socket_thetas.getDataHandle();
      for (CFuint iEq = 0; iEq < nbEqs; ++iEq)
      {
        thetas( distdata.cellID*m_max_nbsubcells + iSubCell, iEq, nbEqs) = m_theta;
      }
    }
  }

//Then we use this theta to build the blending


  /*****         Triangle 0-3-5          *****/

  substates[0] = states[0];
  substates[1] = states[3];
  substates[2] = states[5];

  getMethodData().getDistributionData().tStates =
    computeConsistentStates(&substates);

  *phi = *getMethodData().getSolutionToDistribMatTrans()->transformFromRef(m_phisubT[0]);


  SafePtr<vector<RealMatrix> >& currBetaMatrix =
    getMethodData().getDistributionData().currBetaMat;

  currBetaMatrix = &getMethodData().getDistributionData().betaMats[0];

  distributeLDA(m_k1Plus,m_phi);

  residual[0] += m_theta1*m_phiN1[0] + (1.0 - m_theta1)*m_phi[0];
  residual[3] += m_theta1*m_phiN1[1] + (1.0 - m_theta1)*m_phi[1];
  residual[5] += m_theta1*m_phiN1[2] + (1.0 - m_theta1)*m_phi[2];

  /*****         Triangle 3-1-4          *****/

  substates[0] = states[3];
  substates[1] = states[1];
  substates[2] = states[4];

  getMethodData().getDistributionData().tStates =
    computeConsistentStates(&substates);

  *phi = *getMethodData().getSolutionToDistribMatTrans()->transformFromRef(m_phisubT[1]);

  currBetaMatrix = &getMethodData().getDistributionData().betaMats[1];
  cf_assert (currBetaMatrix.isNotNull());

  distributeLDA(m_k2Plus,m_phi);


  residual[3] += m_theta2*m_phiN2[0] + (1.0 - m_theta2)*m_phi[0];
  residual[1] += m_theta2*m_phiN2[1] + (1.0 - m_theta2)*m_phi[1];
  residual[4] += m_theta2*m_phiN2[2] + (1.0 - m_theta2)*m_phi[2];

  /*****         Triangle 5-4-2          *****/

  substates[0] = states[5];
  substates[1] = states[4];
  substates[2] = states[2];

  getMethodData().getDistributionData().tStates =
    computeConsistentStates(&substates);

  *phi = *getMethodData().getSolutionToDistribMatTrans()->transformFromRef(m_phisubT[2]);


  currBetaMatrix = &getMethodData().getDistributionData().betaMats[2];
  cf_assert (currBetaMatrix.isNotNull());

  distributeLDA(m_k3Plus,m_phi);

  residual[5] += m_theta3*m_phiN3[0] + (1.0 - m_theta3)*m_phi[0];
  residual[4] += m_theta3*m_phiN3[1] + (1.0 - m_theta3)*m_phi[1];
  residual[2] += m_theta3*m_phiN3[2] + (1.0 - m_theta3)*m_phi[2];

   /*****         Triangle 4-5-3          *****/

  substates[0] = states[4];
  substates[1] = states[5];
  substates[2] = states[3];

  getMethodData().getDistributionData().tStates =
    computeConsistentStates(&substates);

  *phi = *getMethodData().getSolutionToDistribMatTrans()->transformFromRef(m_phisubT[3]);

  currBetaMatrix = &getMethodData().getDistributionData().betaMats[3];
  cf_assert (currBetaMatrix.isNotNull());

  distributeLDA(m_k4Plus,m_phi);

  residual[4] += m_theta4*m_phiN4[0] + (1.0 - m_theta4)*m_phi[0];
  residual[5] += m_theta4*m_phiN4[1] + (1.0 - m_theta4)*m_phi[1];
  residual[3] += m_theta4*m_phiN4[2] + (1.0 - m_theta4)*m_phi[2];

  cellnormals.unscale();

 if (m_add_artvisc)
  {
  computeArtificialViscosity(m_res_art_visc);

  residual[0] += m_res_art_visc[0];
  residual[1] += m_res_art_visc[1];
  residual[2] += m_res_art_visc[2];
  residual[3] += m_res_art_visc[3];
  residual[4] += m_res_art_visc[4];
  residual[5] += m_res_art_visc[5];
}
}

//////////////////////////////////////////////////////////////////////////////

void HOCRD_BT_SysSplitStrategy::computeArtificialViscosity(std::vector<RealVector>& result)
{

  DataHandle< InwardNormalsData*> normals = socket_normals.getDataHandle();
  DataHandle< CFreal> updateCoeff = socket_updateCoeff.getDataHandle();

  const CFuint nbCellStates = states.size();
  CFuint m_cellID = getMethodData().getDistributionData().cellID;
  GeometricEntity *const geo = getMethodData().getDistributionData().cell;

 m_cellVolume = geo->computeVolume();

  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();

  for (CFuint iState = 0; iState < nbCellStates; ++iState)
    for (CFuint iEq = 0; iEq < nbEqs; ++iEq)
        result[iState][iEq] = 0.0;

const CFreal mu = (m_d0*sqrt(m_cellVolume)/2.0)*sin(MathTools::MathConsts::CFrealPi()*0.5*m_theta);

// unused // const CFreal ovDimCoeff = 1./(dimCoeff);
  // Triangle 0 : nodes 0-3-5

  CFuint i1 = 0;
  CFuint i2 = 3;
  CFuint i3 = 5;

  // Here we compute the two residual of the diffusive part
  // first one is the integrate of diffusion times basis function (corresponding to a galerkin)
  //the second is the integrate of diffusion times a bubble function.
  // This two integrals are done on each sub-elemement.
  // @todo : possible to compute fluct_diff_galerkin only one time for the whole triangle

  fluctuation_diff_galerkin(i1, i2, i3);
  fluctuation_diff_bubble(i1, i2, i3);


  // call beta of the sub-triangles, betas are always the one of the LDA schemes
  vector<RealMatrix> betasInTriag;
  betasInTriag.resize(3);
  for (CFuint i = 0 ; i< 3; ++i)
    betasInTriag[i].resize(4,4);
  betasInTriag = getMethodData().getDistributionData().betaMats[0];

  // The "bubble" part should be distributed using the Beta of LDA, this for consistency
  for ( CFuint iNode = 0; iNode < 3; ++iNode) {
     (*_phi_diff_bub_split[iNode]) = betasInTriag[iNode]*(*_phi_diff_bub);
  }
  // Diffusive resildual is distributed to each node of the element

  for (CFuint iEq = 0; iEq < nbEqs; ++iEq) {
      result[0][iEq] += (*_phi_diff_gal[0])[iEq] + 3.0*(*_phi_diff_bub_split[0])[iEq] - kappa(0,0)*(*_phi_diff_bub)[iEq];
      result[1][iEq] += (*_phi_diff_gal[1])[iEq] - kappa(0,1)*(*_phi_diff_bub)[iEq];
      result[2][iEq] += (*_phi_diff_gal[2])[iEq] - kappa(0,2)*(*_phi_diff_bub)[iEq];
      result[3][iEq] += (*_phi_diff_gal[3])[iEq] + 3.0*(*_phi_diff_bub_split[1])[iEq] - kappa(0,3)*(*_phi_diff_bub)[iEq];
      result[4][iEq] += (*_phi_diff_gal[4])[iEq] - kappa(0,4)*(*_phi_diff_bub)[iEq];
      result[5][iEq] += (*_phi_diff_gal[5])[iEq] + 3.0*(*_phi_diff_bub_split[2])[iEq] - kappa(0,5)*(*_phi_diff_bub)[iEq];
  }


  // Triangle 2 : nodes 3-1-4

  i1 = 3;
  i2 = 1;
  i3 = 4;

  fluctuation_diff_galerkin(i1, i2, i3);
  fluctuation_diff_bubble(i1, i2, i3);

  betasInTriag = getMethodData().getDistributionData().betaMats[1];
  for ( CFuint iNode = 0; iNode < 3; ++iNode) {
     (*_phi_diff_bub_split[iNode]) = betasInTriag[iNode]*(*_phi_diff_bub);
  }

  for (CFuint iEq = 0; iEq < nbEqs; ++iEq) {
     result[0][iEq] += (*_phi_diff_gal[0])[iEq] - kappa(1,0)*(*_phi_diff_bub)[iEq];
     result[1][iEq] += (*_phi_diff_gal[1])[iEq] + 3.0*(*_phi_diff_bub_split[1])[iEq] - kappa(1,1)*(*_phi_diff_bub)[iEq];
     result[2][iEq] += (*_phi_diff_gal[2])[iEq] - kappa(1,2)*(*_phi_diff_bub)[iEq];
     result[3][iEq] += (*_phi_diff_gal[3])[iEq] + 3.0*(*_phi_diff_bub_split[0])[iEq] - kappa(1,3)*(*_phi_diff_bub)[iEq];
     result[4][iEq] += (*_phi_diff_gal[4])[iEq] + 3.0*(*_phi_diff_bub_split[2])[iEq] - kappa(1,4)*(*_phi_diff_bub)[iEq];
     result[5][iEq] += (*_phi_diff_gal[5])[iEq] - kappa(1,5)*(*_phi_diff_bub)[iEq];
  }


  // Triangle 3 : nodes 5-4-2
  i1 = 5;
  i2 = 4;
  i3 = 2;

  fluctuation_diff_galerkin(i1, i2, i3);
  fluctuation_diff_bubble(i1, i2, i3);


  betasInTriag = getMethodData().getDistributionData().betaMats[2];

  for ( CFuint iNode = 0; iNode < 3; ++iNode) {
     (*_phi_diff_bub_split[iNode]) = betasInTriag[iNode]*(*_phi_diff_bub);
  }

  for (CFuint iEq = 0; iEq < nbEqs; ++iEq) {
      result[0][iEq] += (*_phi_diff_gal[0])[iEq] - kappa(2,0)*(*_phi_diff_bub)[iEq];
      result[1][iEq] += (*_phi_diff_gal[1])[iEq] - kappa(2,1)*(*_phi_diff_bub)[iEq];
      result[2][iEq] += (*_phi_diff_gal[2])[iEq] + 3.0*(*_phi_diff_bub_split[2])[iEq] - kappa(2,2)*(*_phi_diff_bub)[iEq];
      result[3][iEq] += (*_phi_diff_gal[3])[iEq] - kappa(2,3)*(*_phi_diff_bub)[iEq];
      result[4][iEq] += (*_phi_diff_gal[4])[iEq] + 3.0*(*_phi_diff_bub_split[1])[iEq] - kappa(2,4)*(*_phi_diff_bub)[iEq];
      result[5][iEq] += (*_phi_diff_gal[5])[iEq] + 3.0*(*_phi_diff_bub_split[0])[iEq] - kappa(2,5)*(*_phi_diff_bub)[iEq];

  }

  // Triangle 4 : nodes 4-5-3
  i1 = 4;
  i2 = 5;
  i3 = 3;

  fluctuation_diff_galerkin(i1, i2, i3);
  fluctuation_diff_bubble(i1, i2, i3);

  betasInTriag = getMethodData().getDistributionData().betaMats[3];

  for ( CFuint iNode = 0; iNode < 3; ++iNode) {
     (*_phi_diff_bub_split[iNode]) = betasInTriag[iNode]*(*_phi_diff_bub);
  }

  for (CFuint iEq = 0; iEq < nbEqs; ++iEq) {
     result[0][iEq] += (*_phi_diff_gal[0])[iEq] - kappa(3,0)*(*_phi_diff_bub)[iEq];
     result[1][iEq] += (*_phi_diff_gal[1])[iEq] - kappa(3,1)*(*_phi_diff_bub)[iEq];
     result[2][iEq] += (*_phi_diff_gal[2])[iEq] - kappa(3,2)*(*_phi_diff_bub)[iEq];
     result[3][iEq] += (*_phi_diff_gal[3])[iEq] + 3.0*(*_phi_diff_bub_split[2])[iEq] - kappa(3,3)*(*_phi_diff_bub)[iEq];
     result[4][iEq] += (*_phi_diff_gal[4])[iEq] + 3.0*(*_phi_diff_bub_split[0])[iEq] - kappa(3,4)*(*_phi_diff_bub)[iEq];
     result[5][iEq] += (*_phi_diff_gal[5])[iEq] + 3.0*(*_phi_diff_bub_split[1])[iEq] - kappa(3,5)*(*_phi_diff_bub)[iEq];

  }

    //update coefficient

    const CFreal diffCoeff = mu ;
    CFreal coeff = diffCoeff/(m_cellVolume*4.0);
    if (!getMethodData().getDistributionData().isPerturb) {
      for (CFuint i = 0 ; i < 3; ++i) {
        const CFreal faceArea = (normals[m_cellID]->getAreaNode(i));

        updateCoeff[geo->getState(i)->getLocalID()] += coeff*faceArea*faceArea;
      }

      coeff = 2.0*diffCoeff/(3.0*m_cellVolume);
      CFreal dot_prod = (normals[m_cellID]->getNodalNormComp(0,XX))*(normals[m_cellID]->getNodalNormComp(1,XX)) +
                        (normals[m_cellID]->getNodalNormComp(0,YY))*(normals[m_cellID]->getNodalNormComp(1,YY));

      CFreal faceArea1 = (normals[m_cellID]->getAreaNode(0));
      CFreal faceArea2 = (normals[m_cellID]->getAreaNode(1));
      updateCoeff[geo->getState(3)->getLocalID()] += coeff*(faceArea1*faceArea1 + faceArea2*faceArea2 + dot_prod);
                                                     ;

      faceArea1 = (normals[m_cellID]->getAreaNode(1));
      faceArea2 = (normals[m_cellID]->getAreaNode(2));
      dot_prod = (normals[m_cellID]->getNodalNormComp(1,XX))*(normals[m_cellID]->getNodalNormComp(2,XX)) +
                 (normals[m_cellID]->getNodalNormComp(1,YY))*(normals[m_cellID]->getNodalNormComp(2,YY));
      updateCoeff[geo->getState(4)->getLocalID()] += coeff*(faceArea1*faceArea1 + faceArea2*faceArea2 + dot_prod);

      faceArea1 = (normals[m_cellID]->getAreaNode(2));
      faceArea2 = (normals[m_cellID]->getAreaNode(0));
      dot_prod = (normals[m_cellID]->getNodalNormComp(0,XX))*(normals[m_cellID]->getNodalNormComp(2,XX)) +
                 (normals[m_cellID]->getNodalNormComp(0,YY))*(normals[m_cellID]->getNodalNormComp(2,YY));
      updateCoeff[geo->getState(5)->getLocalID()] += coeff*(faceArea1*faceArea1 + faceArea2*faceArea2 + dot_prod);
}

  for (CFuint iEq = 0; iEq < nbEqs; ++iEq) {
     result[0][iEq] = mu*result[0][iEq];
     result[1][iEq] = mu*result[1][iEq];
     result[2][iEq] = mu*result[2][iEq];
     result[3][iEq] = mu*result[3][iEq];
     result[4][iEq] = mu*result[4][iEq];
     result[5][iEq] = mu*result[5][iEq];

  }
}
//////////////////////////////////////////////////////////////////////////////

void HOCRD_BT_SysSplitStrategy::fluctuation_diff_galerkin(CFuint& i1, CFuint& i2, CFuint& i3){


  // Coordinate of the vertex of the element
  const Node& node0_vert = states[0]->getCoordinates();
  const Node& node1_vert = states[1]->getCoordinates();
  const Node& node2_vert = states[2]->getCoordinates();

  const CFreal x1_vert = node0_vert[XX];
  const CFreal x2_vert = node1_vert[XX];
  const CFreal x3_vert = node2_vert[XX];

  const CFreal y1_vert = node0_vert[YY];
  const CFreal y2_vert = node1_vert[YY];
  const CFreal y3_vert = node2_vert[YY];

  //Coordinates of the node defining the surface where the integral is computed
  //  This correspond to vertex of the sub-element on which we integrate
  Node& node0 = states[i1]->getCoordinates();
  Node& node1 = states[i2]->getCoordinates();
  Node& node2 = states[i3]->getCoordinates();

  const CFuint nbStates = states.size();
  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();

  InwardNormalsData& cellnormals = (*socket_normals.getDataHandle()[getMethodData().getDistributionData().cellID]);

  const CFreal nx1 = cellnormals.getNodalNormComp(0,XX);
  const CFreal nx2 = cellnormals.getNodalNormComp(1,XX);
  const CFreal nx3 = cellnormals.getNodalNormComp(2,XX);

  const CFreal ny1 = cellnormals.getNodalNormComp(0,YY);
  const CFreal ny2 = cellnormals.getNodalNormComp(1,YY);
  const CFreal ny3 = cellnormals.getNodalNormComp(2,YY);

  const CFuint dim = PhysicalModelStack::getActive()->getDim();
  const CFreal dimCoeff = 1./dim;
  const CFreal coeffGrad = dimCoeff/m_cellVolume;
  const CFreal inv_volume = 1.0/m_cellVolume;
  const CFreal one_eighth = 1.0/8.0;


  for (CFuint iStates = 0; iStates < nbStates; ++ iStates ){
     for (CFuint iEq = 0; iEq < nbEqs; ++ iEq)
           (*_phi_diff_gal[iStates])[iEq] = 0.0;
  }

  for (CFuint iQd = 0; iQd < 4; ++iQd) {
     //point od quadrature
     const CFreal x = artvisc_qd0[iQd] * node0[XX] + artvisc_qd1[iQd] * node1[XX] + artvisc_qd2[iQd] * node2[XX] ;
     const CFreal y = artvisc_qd0[iQd] * node0[YY] + artvisc_qd1[iQd] * node1[YY] + artvisc_qd2[iQd] * node2[YY] ;

     // Linear basis function
     CFreal L1 = 1.0 + 0.5*( ( x - x1_vert )*nx1 + ( y - y1_vert )*ny1 )*inv_volume ;
     CFreal L2 = 1.0 + 0.5*( ( x - x2_vert )*nx2 + ( y - y2_vert )*ny2 )*inv_volume ;
     CFreal L3 = 1.0 + 0.5*( ( x - x3_vert )*nx3 + ( y - y3_vert )*ny3 )*inv_volume ;


     RealVector grad_state_x;
     grad_state_x.resize(4);
     RealVector grad_state_y;
     grad_state_y.resize(4);

     for (CFuint i = 0; i<4; ++i){
         grad_state_x[i] = (nx1*(4.0*L1 - 1.0)*(*states[0])[i] +
		            nx2*(4.0*L2 - 1.0)*(*states[1])[i] +
		            nx3*(4.0*L3 - 1.0)*(*states[2])[i] +
		            4.0*(nx1*L2 + nx2*L1)*(*states[3])[i] +
		            4.0*(nx2*L3 + nx3*L2)*(*states[4])[i] +
		            4.0*(nx3*L1 + nx1*L3)*(*states[5])[i])*coeffGrad;

        grad_state_y[i] = (ny1*(4.0*L1 - 1.0)*(*states[0])[i] +
		           ny2*(4.0*L2 - 1.0)*(*states[1])[i] +
		           ny3*(4.0*L3 - 1.0)*(*states[2])[i] +
		           4.0*(ny1*L2 + ny2*L1)*(*states[3])[i] +
		           4.0*(ny2*L3 + ny3*L2)*(*states[4])[i] +
		           4.0*(ny3*L1 + ny1*L3)*(*states[5])[i])*coeffGrad;
        }

  for (CFuint i = 0 ; i <nbEqs; ++i)
	{
          F1 = grad_state_x[i]*nx1 + grad_state_y[i] *ny1;
          F2 = grad_state_x[i]*nx2 + grad_state_y[i] *ny2;
          F3 = grad_state_x[i]*nx3 + grad_state_y[i] *ny3;
}

      for (CFuint i = 1 ; i < 3; ++i)
	{
          (*_phi_diff_gal[0])[i] += F1[i]*(4.0*L1 - 1.0)*one_eighth*artvisc_wqd[iQd];
          (*_phi_diff_gal[1])[i] += F2[i]*(4.0*L2 - 1.0)*one_eighth*artvisc_wqd[iQd];
          (*_phi_diff_gal[2])[i] += F3[i]*(4.0*L3 - 1.0)*one_eighth*artvisc_wqd[iQd];
          (*_phi_diff_gal[3])[i] += 4.0*(F1[i]*L2 + F2[i]*L1)*one_eighth*artvisc_wqd[iQd];
          (*_phi_diff_gal[4])[i] += 4.0*(F2[i]*L3 + F3[i]*L2)*one_eighth*artvisc_wqd[iQd];
          (*_phi_diff_gal[5])[i] += 4.0*(F1[i]*L3 + F3[i]*L1)*one_eighth*artvisc_wqd[iQd];
       }

  }

}

//////////////////////////////////////////////////////////////////////////////

void HOCRD_BT_SysSplitStrategy::fluctuation_diff_bubble(CFuint& i1, CFuint& i2, CFuint& i3){


  Node& node0_vert = states[0]->getCoordinates();
  Node& node1_vert = states[1]->getCoordinates();
  Node& node2_vert = states[2]->getCoordinates();

  // Coordinate of the vertex
  const CFreal x1_vert = node0_vert[XX];
  const CFreal x2_vert = node1_vert[XX];
  const CFreal x3_vert = node2_vert[XX];

  const CFreal y1_vert = node0_vert[YY];
  const CFreal y2_vert = node1_vert[YY];
  const CFreal y3_vert = node2_vert[YY];


  Node& node0 = states[i1]->getCoordinates();
  Node& node1 = states[i2]->getCoordinates();
  Node& node2 = states[i3]->getCoordinates();

  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();

  InwardNormalsData& cellnormals = (*socket_normals.getDataHandle()[getMethodData().getDistributionData().cellID]);

  const CFreal nx1 = cellnormals.getNodalNormComp(0,XX);
  const CFreal nx2 = cellnormals.getNodalNormComp(1,XX);
  const CFreal nx3 = cellnormals.getNodalNormComp(2,XX);

  const CFreal ny1 = cellnormals.getNodalNormComp(0,YY);
  const CFreal ny2 = cellnormals.getNodalNormComp(1,YY);
  const CFreal ny3 = cellnormals.getNodalNormComp(2,YY);

  const CFuint dim = PhysicalModelStack::getActive()->getDim();
  const CFreal dimCoeff = 1./dim;
  const CFreal coeffGrad = dimCoeff/m_cellVolume;
  const CFreal inv_volume = 1.0/m_cellVolume;
  const CFreal one_third = 1.0/3.0;

  // coordinates of the gravity center of the sub-element
  const CFreal xg = (node0[XX] + node1[XX] + node2[XX])*one_third;
  const CFreal yg = (node0[YY] + node1[YY] + node2[YY])*one_third;


  // The buble function is linear peacewise linear. Which means that its gradient is not continuous
  // Then, we just need one node per constant part. i choose to take one point at the gravity center of the 3 triangkes
  // composed by nodes and gravity center of the sub-element.
  for (CFuint iEq = 0; iEq < nbEqs; ++ iEq)
     (*_phi_diff_bub)[iEq] = 0.0;

  /**************************************************/
  /***          1st small triangle i1-i3-ng       ***/
  /**************************************************/

  double x = (xg + node0[XX] + node2[XX])*one_third ;
  double y = (yg + node0[YY] + node2[YY])*one_third ;

  CFreal L1 = 1.0 + 0.5*( ( x - x1_vert )*nx1 + ( y - y1_vert )*ny1 )*inv_volume ;
  CFreal L2 = 1.0 + 0.5*( ( x - x2_vert )*nx2 + ( y - y2_vert )*ny2 )*inv_volume ;
  CFreal L3 = 1.0 + 0.5*( ( x - x3_vert )*nx3 + ( y - y3_vert )*ny3 )*inv_volume ;

  RealVector grad_state_x;
  grad_state_x.resize(nbEqs);
  RealVector grad_state_y;
  grad_state_y.resize(nbEqs);

  for (CFuint i = 0; i<4; ++i){
     grad_state_x[i] = (nx1*(4.0*L1 - 1.0)*(*states[0])[i] +
		        nx2*(4.0*L2 - 1.0)*(*states[1])[i] +
		        nx3*(4.0*L3 - 1.0)*(*states[2])[i] +
		        4.0*(nx1*L2 + nx2*L1)*(*states[3])[i] +
		        4.0*(nx2*L3 + nx3*L2)*(*states[4])[i] +
		        4.0*(nx3*L1 + nx1*L3)*(*states[5])[i])*coeffGrad;

     grad_state_y[i] = (ny1*(4.0*L1 - 1.0)*(*states[0])[i] +
		        ny2*(4.0*L2 - 1.0)*(*states[1])[i] +
		        ny3*(4.0*L3 - 1.0)*(*states[2])[i] +
		        4.0*(ny1*L2 + ny2*L1)*(*states[3])[i] +
		        4.0*(ny2*L3 + ny3*L2)*(*states[4])[i] +
		        4.0*(ny3*L1 + ny1*L3)*(*states[5])[i])*coeffGrad;
  }

  _normal[XX] = (node0[YY] - node2[YY]);
  _normal[YY] = (node2[XX] - node0[XX]);

 for (CFuint i = 0 ; i <nbEqs; ++i)
	{
          F1[i] = grad_state_x[i]*_normal[XX] + grad_state_y[i] *_normal[YY];

}


  for (CFuint i = 0 ; i < nbEqs; ++i){
     (*_phi_diff_bub)[i] = F1[i]*m_cellVolume/
     (12.0*((node2[XX]*yg - xg*node2[YY]) - node0[XX]*(yg - node2[YY]) + node0[YY]*(xg - node2[XX])));
     }

  /**************************************************/
  /***         2nd small triangle  i3-i2-g     ***/
  /**************************************************/

  x = (xg + node1[XX] + node2[XX])*one_third ;
  y = (yg + node1[YY] + node2[YY])*one_third ;

  L1 = 1.0 + 0.5*( ( x - x1_vert )*nx1 + ( y - y1_vert )*ny1 )*inv_volume ;
  L2 = 1.0 + 0.5*( ( x - x2_vert )*nx2 + ( y - y2_vert )*ny2 )*inv_volume ;
  L3 = 1.0 + 0.5*( ( x - x3_vert )*nx3 + ( y - y3_vert )*ny3 )*inv_volume ;
for (CFuint i = 0; i<4; ++i){
  grad_state_x[i] = (nx1*(4.0*L1 - 1.0)*(*states[0])[i] +
		  nx2*(4.0*L2 - 1.0)*(*states[1])[i] +
		  nx3*(4.0*L3 - 1.0)*(*states[2])[i] +
		  4.0*(nx1*L2 + nx2*L1)*(*states[3])[i] +
		  4.0*(nx2*L3 + nx3*L2)*(*states[4])[i] +
		  4.0*(nx3*L1 + nx1*L3)*(*states[5])[i])*coeffGrad;

  grad_state_y[i] = (ny1*(4.0*L1 - 1.0)*(*states[0])[i] +
		  ny2*(4.0*L2 - 1.0)*(*states[1])[i] +
		  ny3*(4.0*L3 - 1.0)*(*states[2])[i] +
		  4.0*(ny1*L2 + ny2*L1)*(*states[3])[i] +
		  4.0*(ny2*L3 + ny3*L2)*(*states[4])[i] +
		  4.0*(ny3*L1 + ny1*L3)*(*states[5])[i])*coeffGrad;
}
  _normal[XX] = (node1[YY] - node2[YY]);
  _normal[YY] = (node2[XX] - node1[XX]);

 for (CFuint i = 0 ; i <nbEqs; ++i)
	{
          F1[i] = grad_state_x[i]*_normal[XX] + grad_state_y[i] *_normal[YY];

}
  for (CFuint i = 1 ; i < 3; ++i){
     (*_phi_diff_bub)[i] += F1[i]*m_cellVolume/
	                    (12.0*((node2[XX]*yg - xg*node2[YY]) - node1[XX]*(yg - node2[YY]) + node1[YY]*(xg - node2[XX])));
      }

  /**************************************************/
  /***         3rd small triangle  i1-i3-ig       ***/
  /**************************************************/

  x = (xg + node0[XX] + node1[XX])*one_third ;
  y = (yg + node0[YY] + node1[YY])*one_third ;

  L1 = 1.0 + 0.5*( ( x - x1_vert )*nx1 + ( y - y1_vert )*ny1 )*inv_volume ;
  L2 = 1.0 + 0.5*( ( x - x2_vert )*nx2 + ( y - y2_vert )*ny2 )*inv_volume ;
  L3 = 1.0 + 0.5*( ( x - x3_vert )*nx3 + ( y - y3_vert )*ny3 )*inv_volume ;
for (CFuint i = 0; i<4; ++i){
  grad_state_x[i] = (nx1*(4.0*L1 - 1.0)*(*states[0])[i] +
		  nx2*(4.0*L2 - 1.0)*(*states[1])[i] +
		  nx3*(4.0*L3 - 1.0)*(*states[2])[i] +
		  4.0*(nx1*L2 + nx2*L1)*(*states[3])[i] +
		  4.0*(nx2*L3 + nx3*L2)*(*states[4])[i] +
		  4.0*(nx3*L1 + nx1*L3)*(*states[5])[i])*coeffGrad;

  grad_state_y[i] = (ny1*(4.0*L1 - 1.0)*(*states[0])[i] +
		  ny2*(4.0*L2 - 1.0)*(*states[1])[i] +
		  ny3*(4.0*L3 - 1.0)*(*states[2])[i] +
		  4.0*(ny1*L2 + ny2*L1)*(*states[3])[i] +
		  4.0*(ny2*L3 + ny3*L2)*(*states[4])[i] +
		  4.0*(ny3*L1 + ny1*L3)*(*states[5])[i])*coeffGrad;
}
  _normal[XX] = (node0[YY] - node1[YY]);
  _normal[YY] = (node1[XX] - node0[XX]);


 for (CFuint i = 0 ; i <nbEqs; ++i)
	{
          F1[i] = grad_state_x[i]*_normal[XX] + grad_state_y[i] *_normal[YY];

}

  for (CFuint i = 0 ; i < nbEqs; ++i){
      (*_phi_diff_bub)[i] += F1[i]*m_cellVolume/
	                     (12.0*((node1[XX]*yg - xg*node1[YY]) - node0[XX]*(yg - node1[YY]) + node0[YY]*(xg - node1[XX])));

      }


}
//////////////////////////////////////////////////////////////////////////////

void HOCRD_BT_SysSplitStrategy::computeHOFluctuation()
{

  cf_assert(qdstates.size() == 3); // only triags so three quadrature points per face

  vector<State*>& states = *getMethodData().getDistributionData().states;

  const State& state0 = *(states[0]);
  const State& state1 = *(states[1]);
  const State& state2 = *(states[2]);
  const State& state3 = *(states[3]);
  const State& state4 = *(states[4]);
  const State& state5 = *(states[5]);

  vector<Node*>& nodes = *getMethodData().getDistributionData().cell->getNodes();

  cf_assert(nodes.size()  == 3); // P1 triangles for geometry space

  const CFreal x1 = (*nodes[0])[XX];
  const CFreal x2 = (*nodes[1])[XX];
  const CFreal x3 = (*nodes[2])[XX];

  const CFreal y1 = (*nodes[0])[YY];
  const CFreal y2 = (*nodes[1])[YY];
  const CFreal y3 = (*nodes[2])[YY];

  InwardNormalsData& cellnormals = (*socket_normals.getDataHandle()
				    [getMethodData().getDistributionData().cellID]);
  // we multiply by two because the normals scale has been changed in computeFluctuationAndUpdateCoeff
  // and is 0.5
  const CFreal nx1 = 2. * cellnormals.getNodalNormComp(0,XX);
  const CFreal nx2 = 2. * cellnormals.getNodalNormComp(1,XX);
  const CFreal nx3 = 2. * cellnormals.getNodalNormComp(2,XX);

  const CFreal ny1 = 2. * cellnormals.getNodalNormComp(0,YY);
  const CFreal ny2 = 2. * cellnormals.getNodalNormComp(1,YY);
  const CFreal ny3 = 2. * cellnormals.getNodalNormComp(2,YY);

  const CFreal inv_volume = 1.0 / getMethodData().getDistributionData().cell->computeVolume();

 // Computation  of the fluctuation of each faces of the cell
for (CFuint iFace = 0; iFace < subfacetable.nbRows(); ++iFace)
  {
    facenormal[XX] = cellnormals.getNodalNormComp(iFace%3,XX);
    facenormal[YY] = cellnormals.getNodalNormComp(iFace%3,YY);


    for (CFuint iQd = 0; iQd < 3; ++iQd)
    {

      const Node& node0 = states[subfacetable(iFace,0)]->getCoordinates();
      const Node& node1 = states[subfacetable(iFace,1)]->getCoordinates();
      const CFreal x = qd0[iQd] * node0[XX] + qd1[iQd] * node1[XX];
      const CFreal y = qd0[iQd] * node0[YY] + qd1[iQd] * node1[YY];
      const CFreal L1 = 1.0 + 0.5*( ( x - x1 )*nx1 + ( y - y1 )*ny1 ) * inv_volume ;
      const CFreal L2 = 1.0 + 0.5*( ( x - x2 )*nx2 + ( y - y2 )*ny2 ) * inv_volume ;
      const CFreal L3 = 1.0 + 0.5*( ( x - x3 )*nx3 + ( y - y3 )*ny3 ) * inv_volume ;

      (*qdstates[iQd]) = (L1*( 2.0*L1 - 1.0 ) * state0) +
                         (L2*( 2.0*L2 - 1.0 ) * state1) +
                         (L3*( 2.0*L3 - 1.0 ) * state2) +
                         (4.0*L1*L2           * state3) +
                         (4.0*L3*L2           * state4) +
                         (4.0*L1*L3           * state5);
    }
    
    computeStatesData(3, m_updateVar, qdstates, m_pdata, m_qdExtraVars); // three quadrature points per face
    
    faceflux[iFace] = 0.;
    for (CFuint iQd = 0; iQd < 3; ++iQd)
    {

      faceflux[iFace] += wqd[iQd] * m_updateVar->getFlux()(m_pdata[iQd],facenormal);
    }
  }

  cf_assert(m_phisubT.size() == subelemtable.nbRows());

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
void HOCRD_BT_SysSplitStrategy::computeK(const std::vector<Framework::State*>& states,
                std::vector<RealMatrix*>& m_kPlus){
  InwardNormalsData& cellnormals = (*socket_normals.getDataHandle()
				    [getMethodData().getDistributionData().cellID]);

  DataHandle< CFreal> updateCoeff = socket_updateCoeff.getDataHandle();
  const CFuint m_nbStatesInCell = states.size();

  const CFreal m_invDim = 1.0/2.0;
  // The transformation of the normal is needed if the coordinate system is rotated
  for (CFuint iState = 0; iState < m_nbStatesInCell; ++iState) {
    // The transformation of the normal is needed if the coordinate system
    // is rotated. The normal is adimensionalized, so there is need to multiply
    // by the nodeArea when computing the k parameter

    for (CFuint iDim = 0; iDim < 2; ++iDim) {
        m_adimNormal[iDim] = cellnormals.getNodalNormComp(iState, iDim);
    }
    m_adimNormal *= 1. / cellnormals.getAreaNode(iState);

    getMethodData().getDistribVar()->splitJacobian(*m_kPlus[iState],
					       *m_kMin[iState],
					       *m_eValues[iState],
					       m_adimNormal);

    CFreal m_nodeArea = cellnormals.getAreaNode(iState);

    *m_kPlus[iState] *= m_invDim * m_nodeArea;
    *m_kMin[iState]  *= m_invDim * m_nodeArea;

    if (!getMethodData().getDistributionData().isPerturb) {
      const CFreal maxEigenValue = std::max(0.0, m_eValues[iState]->max());

      updateCoeff[states[iState]->getLocalID()] += m_invDim*m_nodeArea*maxEigenValue;

    }
  }



}
//////////////////////////////////////////////////////////////////////////////
void HOCRD_BT_SysSplitStrategy::distributeN(std::vector<RealMatrix*> & m_kPlus,vector<RealVector> & phiN){

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
void HOCRD_BT_SysSplitStrategy::distributeLDA(std::vector<RealMatrix*> & m_kPlus,vector<RealVector> & phiLDA){

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
void HOCRD_BT_SysSplitStrategy::computeBlendingCoeff(CFreal & result)
{
  vector<State*> states = substates;

  InwardNormalsData& cellnormals = (*socket_normals.getDataHandle()
				    [getMethodData().getDistributionData().cellID]);
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
      _grad[iDim] += _pData[_varID]*cellnormals.getNodalNormComp(i,iDim);
    }
  }
  cf_assert(dim == DIM_2D);
  const CFreal h = 2.0*std::sqrt(static_cast<CFreal>(vol/MathTools::MathConsts::CFrealPi()));

  CFreal sc;
  if (_sh_detector == "Jirka"){
  sc = _grad[XX]*lData[EulerTerm::VX] + _grad[YY]*lData[EulerTerm::VY];
  if (dim == DIM_3D) {
    sc += _grad[ZZ]*lData[EulerTerm::VZ];
  }
  //   cout << "grad*V  = " << theta << endl;
  //   cout << "_length = " << _length << endl;
  //   cout << "_deltaP = " << _deltaP << endl;
  //   cout << "_speed = " << _speed << endl;

  sc *= _length/(vol*dim*_deltaP*_speed);
  m_sc = sc;
  }
  else if (_sh_detector == "Anton"){
  sc = _grad[XX]*lData[EulerTerm::VX] + _grad[YY]*lData[EulerTerm::VY];
  if (dim == DIM_3D) {
    sc += _grad[ZZ]*lData[EulerTerm::VZ];
  }
  sc *= std::pow(    (std::pow(lData[EulerTerm::VX],2) + std::pow(lData[EulerTerm::VY],2)), 0.5   );
  //   cout << "grad*V  = " << theta << endl;
  //   cout << "_length = " << _length << endl;
  //   cout << "_deltaP = " << _deltaP << endl;
  //   cout << "_speed = " << _speed << endl;

  sc *= _length/(vol*dim*_deltaP*_speed * _speed);
  m_sc = sc;
}
  //second test
  result = min(1.0,(max(0.0, m_sc))*(max(0.0, m_sc))*h);

  // this->_alpha = min(1.0,theta*h);
  // cout << "theta = " << this->_alpha << endl;
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<Framework::BaseDataSocketSink> >
HOCRD_BT_SysSplitStrategy::needsSockets()
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
