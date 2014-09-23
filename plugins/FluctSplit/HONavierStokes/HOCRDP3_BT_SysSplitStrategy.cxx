#include "FluctSplit/HONavierStokes/HOCRDP3_BT_SysSplitStrategy.hh"

#include "Common/BadValueException.hh"
#include "Framework/MethodStrategyProvider.hh"
#include "Framework/ContourIntegrator.hh"
#include "Framework/MeshData.hh"
#include "Framework/BaseTerm.hh"
#include "FluctSplit/FluctSplitHO.hh"

#include "MathTools/MatrixInverter.hh"

#include "NavierStokes/EulerTerm.hh"

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

MethodStrategyProvider<HOCRDP3_BT_SysSplitStrategy,
                       FluctuationSplitData,
                       FluctuationSplitStrategy,
                       FluctSplitHOModule>
hocrp3dbtFluctSplitStrategyProvider("HOCRDP3_BT");
//////////////////////////////////////////////////////////////////////////////

void HOCRDP3_BT_SysSplitStrategy::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< CFreal>("Delta","Delta of variable.");
  options.addConfigOption< CFreal>("Length","Reference Length.");
  options.addConfigOption< CFreal>("Speed","Reference Speed.");
  options.addConfigOption< std::string>("VarName","Variable name.");
  options.addConfigOption< bool >("HO","High order discretization");
  options.addConfigOption< CFuint >("MaxNbSubElems","Maximum number of sub-elements for high-order computations.");
  options.addConfigOption< bool >  ("StoreThetas","Store the thetas for visualization");
  options.addConfigOption< bool >  ("UmaxTheta","Use the maximum of the thetas");
  options.addConfigOption< CFreal, Config::DynamicOption<> >("MinTheta","Minimum theta, used to keep a minimum diffusion");
  options.addConfigOption< std::string>("Shockdetector","Which shock detetecto to use");
}
//////////////////////////////////////////////////////////////////////////////

HOCRDP3_BT_SysSplitStrategy::HOCRDP3_BT_SysSplitStrategy(const std::string& name) :
  FluctuationSplitStrategy(name),
  socket_updateCoeff("updateCoeff"),
  socket_thetas("thetas",false),
  m_solutionVar(CFNULL),
  m_updateVar(CFNULL),
  m_unitFaceNormals(),
  m_phiN1(0),
  m_phiN2(0),
  m_phiN3(0),
  m_phiN4(0),
  m_phiN5(0),
  m_phiN6(0),
  m_phiN7(0),
  m_phiN8(0),
  m_phiN9(0),
  m_phi(0),
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
  m_k5Plus(0),
  m_k6Plus(0),
  m_k7Plus(0),
  m_k8Plus(0),
  m_k9Plus(0),
  m_k(0),
  m_kMin(0),
  m_eValues(0),
  m_theta1(),
  m_theta2(),
  m_theta3(),
  m_theta4(),
  m_theta5(),
  m_theta6(),
  m_theta7(),
  m_theta8(),
  m_theta9(),
  m_theta(),
  m_inverter(CFNULL),
  m_invK(),
  m_uInflow(),
  _pData(),
  _grad(),
  m_min_theta(0.),
  m_adimNormal()
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

HOCRDP3_BT_SysSplitStrategy::~HOCRDP3_BT_SysSplitStrategy()
{
  if (isSetup()) unsetup();
}

//////////////////////////////////////////////////////////////////////////////

void HOCRDP3_BT_SysSplitStrategy::unsetup()
{
  // AL: this could be changed ...
  // where are deleted the qdstates ???
  for (CFuint i = 0; i < m_qdExtraVars.size(); ++i) {
    deletePtr(m_qdExtraVars[i]);
  }

  for (CFuint i = 0; i < m_k1Plus.size(); ++i)
  {
    deletePtr(m_k1Plus[i]);
    deletePtr(m_k2Plus[i]);
    deletePtr(m_k3Plus[i]);
    deletePtr(m_k4Plus[i]);
    deletePtr(m_k5Plus[i]);
    deletePtr(m_k6Plus[i]);
    deletePtr(m_k7Plus[i]);
    deletePtr(m_k8Plus[i]);
    deletePtr(m_k9Plus[i]);
    deletePtr(m_k[i]);
    deletePtr(m_kMin[i]);
  }
  deletePtr(m_inverter);
}

//////////////////////////////////////////////////////////////////////////////

void HOCRDP3_BT_SysSplitStrategy::setup()
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
  const CFuint nbQdPts = 5;

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

  for (CFuint isub_elem = 0; isub_elem < m_phisubT.size(); ++ isub_elem)
    m_phisubT[isub_elem] = new RealVector(nbEqs);

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

  const CFreal s1  = (1.0/3.0)*sqrt(5.0-2.0*sqrt(10.0/7.0));
  const CFreal a01 = ( 1.0 - s1 )*0.5;
  const CFreal a11 = ( 1.0 + s1 )*0.5;

  const CFreal s2  = (1.0/3.0)*sqrt(5.0+2.0*sqrt(10.0/7.0));
  const CFreal a02 = ( 1.0 - s2 )*0.5;
  const CFreal a12 = ( 1.0 + s2 )*0.5;

  qd0.resize(nbQdPts); // quadrature points per face
  qd1.resize(nbQdPts); // quadrature points per face

  qd0[0] = a01;  qd1[0] = a11;
  qd0[1] = a11;  qd1[1] = a01;
  qd0[2] = a02;  qd1[2] = a12;
  qd0[3] = a12;  qd1[3] = a02;
  qd0[4] = .5;   qd1[4] = .5;

  wqd.resize(nbQdPts); // 3 quadrature points per face
  wqd[0] = (322.0+13.0*sqrt(70.0))/1800.0;
  wqd[1] = (322.0+13.0*sqrt(70.0))/1800.0;
  wqd[2] = (322.0-13.0*sqrt(70.0))/1800.0;
  wqd[3] = (322.0-13.0*sqrt(70.0))/1800.0;
  wqd[4] = 128.0/450.0;

  qdstates.resize(nbQdPts); // 3 quadrature points per face
  qdstates[0] = new State();
  qdstates[1] = new State();
  qdstates[2] = new State();
  qdstates[3] = new State();
  qdstates[4] = new State();


  const CFuint extra_var_size = getMethodData().getUpdateVar()->getExtraPhysicalVarsSize();
  m_qdExtraVars.resize(qdstates.size());
  for (CFuint i = 0; i < m_qdExtraVars.size(); ++i) {
    m_qdExtraVars[i] = new RealVector(extra_var_size);
  }

  qdnodes.resize(nbQdPts); // 3 quadrature points per face
  qdnodes[0] = new Node();
  qdstates[0]->setSpaceCoordinates(qdnodes[0]);
  qdnodes[1] = new Node();
  qdstates[1]->setSpaceCoordinates(qdnodes[1]);
  qdnodes[2] = new Node();
  qdstates[2]->setSpaceCoordinates(qdnodes[2]);
  qdnodes[3] = new Node();
  qdstates[3]->setSpaceCoordinates(qdnodes[3]);
  qdnodes[4] = new Node();
  qdstates[4]->setSpaceCoordinates(qdnodes[4]);

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

  m_phiN5.resize(3);

  m_phiN5[0].resize(nbEqs);
  m_phiN5[1].resize(nbEqs);
  m_phiN5[2].resize(nbEqs);

  m_phiN6.resize(3);

  m_phiN6[0].resize(nbEqs);
  m_phiN6[1].resize(nbEqs);
  m_phiN6[2].resize(nbEqs);


  m_phiN7.resize(3);

  m_phiN7[0].resize(nbEqs);
  m_phiN7[1].resize(nbEqs);
  m_phiN7[2].resize(nbEqs);

  m_phiN8.resize(3);

  m_phiN8[0].resize(nbEqs);
  m_phiN8[1].resize(nbEqs);
  m_phiN8[2].resize(nbEqs);

  m_phiN9.resize(3);

  m_phiN9[0].resize(nbEqs);
  m_phiN9[1].resize(nbEqs);
  m_phiN9[2].resize(nbEqs);

  m_phi.resize(3);

  m_phi[0].resize(nbEqs);
  m_phi[1].resize(nbEqs);
  m_phi[2].resize(nbEqs);

  m_maxNbStatesInCell = MeshDataStack::getActive()->Statistics().getMaxNbStatesInCell();

  m_k1Plus.resize(m_maxNbStatesInCell);
  m_k2Plus.resize(m_maxNbStatesInCell);
  m_k3Plus.resize(m_maxNbStatesInCell);
  m_k4Plus.resize(m_maxNbStatesInCell);
  m_k5Plus.resize(m_maxNbStatesInCell);
  m_k6Plus.resize(m_maxNbStatesInCell);
  m_k7Plus.resize(m_maxNbStatesInCell);
  m_k8Plus.resize(m_maxNbStatesInCell);
  m_k9Plus.resize(m_maxNbStatesInCell);
  m_k.resize(m_maxNbStatesInCell);
  m_kMin.resize(m_maxNbStatesInCell);
  m_eValues.resize(m_maxNbStatesInCell);

  for (CFuint i = 0; i < m_maxNbStatesInCell; ++i) {
    m_k1Plus[i] = new RealMatrix(nbEqs, nbEqs);
    m_k2Plus[i] = new RealMatrix(nbEqs, nbEqs);
    m_k3Plus[i] = new RealMatrix(nbEqs,nbEqs );
    m_k4Plus[i] = new RealMatrix(nbEqs, nbEqs);
    m_k5Plus[i] = new RealMatrix(nbEqs, nbEqs);
    m_k6Plus[i] = new RealMatrix(nbEqs, nbEqs);
    m_k7Plus[i] = new RealMatrix(nbEqs,nbEqs );
    m_k8Plus[i] = new RealMatrix(nbEqs, nbEqs);
    m_k9Plus[i] = new RealMatrix(nbEqs, nbEqs);
    m_kMin[i] = new RealMatrix(nbEqs, nbEqs);
    m_k[i] = new RealMatrix(nbEqs, nbEqs);
    m_eValues[i] = new RealVector(nbEqs);
  }

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

  m_adimNormal.resize(2);
}

//////////////////////////////////////////////////////////////////////////////

void HOCRDP3_BT_SysSplitStrategy::computeFluctuation(vector<RealVector>& residual)
{
  vector<State*>& states = *getMethodData().getDistributionData().states;
  cf_assert(states.size() == 10); // P3 triangles for solution space
  cf_assert(residual.size() >= states.size()); // state residual
  DistributionData& ddata = getMethodData().getDistributionData();
  // reset the residual because we will accumulate the sub element contributions
  for (CFuint i = 0; i < residual.size(); ++i)
  {
    residual[i] = 0.0;
  }
  InwardNormalsData& cellnormals = (*socket_normals.getDataHandle()
				    [getMethodData().getDistributionData().cellID]);

  cf_assert(cellnormals.nbFaces()  == 3); // triangles have 3 faces

  CFreal onethird = 1.0/3.0;

  // normals are half scale
  cellnormals.scale(onethird);
  computeHOP3Fluctuation();

  /*****         Triangle 0-3-8          *****/

  substates[0] = states[0];
  substates[1] = states[3];
  substates[2] = states[8];

  getMethodData().getDistributionData().tStates =
    computeConsistentStates(&substates);
  
  ddata.subStates = &substates;
  // compute the residual and the upwind parameters k in this cell
  computeK(substates,m_k1Plus);

  // transform fluxes of subelement to distribution variables
  SafePtr<RealVector> phi = &getMethodData().getDistributionData().phi;
  *phi = *getMethodData().getSolutionToDistribMatTrans()->transformFromRef(m_phisubT[0]);
 
  distributeN(m_k1Plus,m_phiN1);


  computeBlendingCoeff( m_theta1);

  // cout << "*currBetaMatrix0 = " << (*currBetaMatrix)[0] << endl;
  // cout << "*currBetaMatrix1 = " << (*currBetaMatrix)[1] << endl;
  // cout << "*currBetaMatrix2 = " << (*currBetaMatrix)[2] << endl;

  /*****         Triangle 4-1-5          *****/

  substates[0] = states[4];
  substates[1] = states[1];
  substates[2] = states[5];

  getMethodData().getDistributionData().tStates = computeConsistentStates(&substates);
  

  // compute the residual and the upwind parameters k in this cell
  computeK(substates,m_k2Plus);

  // transform fluxes of subelement to distribution variables
  *phi = *getMethodData().getSolutionToDistribMatTrans()->transformFromRef(m_phisubT[1]);


  distributeN(m_k2Plus,m_phiN2);
  computeBlendingCoeff(m_theta2);


  /*****         Triangle 7-6-2          *****/
  substates[0] = states[7];
  substates[1] = states[6];
  substates[2] = states[2];

  getMethodData().getDistributionData().tStates = computeConsistentStates(&substates);

  ddata.subStates = &substates;

  // compute the residual and the upwind parameters k in this cell
  computeK(substates,m_k3Plus);


  // transform fluxes of subelement to distribution variables
  *phi = *getMethodData().getSolutionToDistribMatTrans()->transformFromRef(m_phisubT[2]);


  distributeN(m_k3Plus,m_phiN3);
  computeBlendingCoeff(m_theta3);


  /*****         Triangle 9-5-6          *****/

  substates[0] = states[9];
  substates[1] = states[5];
  substates[2] = states[6];

  getMethodData().getDistributionData().tStates = computeConsistentStates(&substates);

   // compute the residual and the upwind parameters k in this cell
  computeK(substates,m_k4Plus);

  // transform fluxes of subelement to distribution variables
   *phi = *getMethodData().getSolutionToDistribMatTrans()->transformFromRef(m_phisubT[3]);


  distributeN(m_k4Plus,m_phiN4);
  computeBlendingCoeff(m_theta4);

 /*****         Triangle 8-9-7          *****/

  substates[0] = states[8];
  substates[1] = states[9];
  substates[2] = states[7];

  getMethodData().getDistributionData().tStates = computeConsistentStates(&substates);
;
 // compute the residual and the upwind parameters k in this cell
  computeK(substates,m_k5Plus);

  // transform fluxes of subelement to distribution variables
  *phi = *getMethodData().getSolutionToDistribMatTrans()->transformFromRef(m_phisubT[4]);

  
  distributeN(m_k5Plus,m_phiN5);

  computeBlendingCoeff(m_theta5);


/*****         Triangle 3-4-9          *****/
  substates[0] = states[3];
  substates[1] = states[4];
  substates[2] = states[9];

   getMethodData().getDistributionData().tStates = computeConsistentStates(&substates);
 
  // compute the residual and the upwind parameters k in this cell
  computeK(substates,m_k6Plus);

  // transform fluxes of subelement to distribution variables
  *phi = *getMethodData().getSolutionToDistribMatTrans()->transformFromRef(m_phisubT[5]);

  distributeN(m_k6Plus,m_phiN6);
  computeBlendingCoeff(m_theta6);


  // compute the residual and the upwind parameters k in these cells
  // the oriantation of the normal in this sub-element is oposite to the one of the element
  cellnormals.scale(-onethird);

/*****         Triangle 9-8-3          *****/
  substates[0] = states[9];
  substates[1] = states[8];
  substates[2] = states[3];

  getMethodData().getDistributionData().tStates = computeConsistentStates(&substates);

  // compute the residual and the upwind parameters k in this cell
  computeK(substates,m_k7Plus);
  
  // transform fluxes of subelement to distribution variables
  *phi = *getMethodData().getSolutionToDistribMatTrans()->transformFromRef(m_phisubT[6]);

  distributeN(m_k7Plus,m_phiN7);
  computeBlendingCoeff(m_theta7);


/*****         Triangle 5-9-4          *****/
  substates[0] = states[5];
  substates[1] = states[9];
  substates[2] = states[4];

  getMethodData().getDistributionData().tStates = computeConsistentStates(&substates);

  // compute the residual and the upwind parameters k in this cell
  computeK(substates,m_k8Plus);

  // transform fluxes of subelement to distribution variables
  *phi = *getMethodData().getSolutionToDistribMatTrans()->transformFromRef(m_phisubT[7]);

  distributeN(m_k8Plus,m_phiN8);
  computeBlendingCoeff(m_theta8);

/*****         Triangle 6-7-9          *****/
  substates[0] = states[6];
  substates[1] = states[7];
  substates[2] = states[9];

  getMethodData().getDistributionData().tStates = computeConsistentStates(&substates);

  // compute the residual and the upwind parameters k in this cell
  computeK(substates,m_k9Plus);

  // transform fluxes of subelement to distribution variables
  *phi = *getMethodData().getSolutionToDistribMatTrans()->transformFromRef(m_phisubT[8]);


  distributeN(m_k9Plus,m_phiN9);
  computeBlendingCoeff(m_theta9);
  cellnormals.unscale();

  // Then we take the max of the thetas
  // this is outside because is used for the artificial viscosity
  m_theta = max( max( max(m_theta1, m_theta2), m_theta3), m_theta4);
  m_theta = max(max( max( m_theta, m_theta5), m_theta6), m_theta7);
  m_theta = max( max( m_theta, m_theta8), m_theta9);
  if (m_use_max_theta)
  {
    m_theta1 = max(m_min_theta, m_theta);
    m_theta2 = m_theta1;
    m_theta3 = m_theta1;
    m_theta4 = m_theta1;
    m_theta5 = m_theta1;
    m_theta6 = m_theta1;
    m_theta7 = m_theta1;
    m_theta8 = m_theta1;
    m_theta9 = m_theta1;

  }
  else
  {
    m_theta1 = max(m_min_theta, m_theta1);
    m_theta2 = max(m_min_theta, m_theta2);
    m_theta3 = max(m_min_theta, m_theta3);
    m_theta4 = max(m_min_theta, m_theta4);
    m_theta5 = max(m_min_theta, m_theta5);
    m_theta6 = max(m_min_theta, m_theta6);
    m_theta7 = max(m_min_theta, m_theta7);
    m_theta8 = max(m_min_theta, m_theta8);
    m_theta9 = max(m_min_theta, m_theta9);
    
  }

  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  DistributionData& distdata = getMethodData().getDistributionData();
  if (m_store_thetas && !distdata.isPerturb)
  {
    for (CFuint iSubCell = 0; iSubCell < 9; ++iSubCell)
    {

      Framework::DataHandle< CFreal > thetas = socket_thetas.getDataHandle();
      for (CFuint iEq = 0; iEq < nbEqs; ++iEq)
      {
        thetas( distdata.cellID*m_max_nbsubcells + iSubCell, iEq, nbEqs) = m_theta;
      }
    }
  }

//Then we use this theta to build the blending

  /*****         Triangle 0-3-8          *****/

  substates[0] = states[0];
  substates[1] = states[3];
  substates[2] = states[8];

  getMethodData().getDistributionData().tStates =
    computeConsistentStates(&substates);
  
  ddata.subStates = &substates;

  *phi = *getMethodData().getSolutionToDistribMatTrans()->transformFromRef(m_phisubT[0]);
 
  SafePtr<vector<RealMatrix> >& currBetaMatrix =
    getMethodData().getDistributionData().currBetaMat;


  // cout <<"getMethodData().getDistributionData().betaMats.size() = " <<
  //  getMethodData().getDistributionData().betaMats.size()  << endl;

  currBetaMatrix = &getMethodData().getDistributionData().betaMats[0];

  distributeLDA(m_k1Plus,m_phi);

  residual[0] += m_theta1*m_phiN1[0] + (1.0 - m_theta1)*m_phi[0];
  residual[3] += m_theta1*m_phiN1[1] + (1.0 - m_theta1)*m_phi[1];
  residual[8] += m_theta1*m_phiN1[2] + (1.0 - m_theta1)*m_phi[2];


  /*****         Triangle 4-1-5          *****/
  substates[0] = states[4];
  substates[1] = states[1];
  substates[2] = states[5];

  getMethodData().getDistributionData().tStates = computeConsistentStates(&substates);
  
  ddata.subStates = &substates;
  cellnormals.scale(onethird);


  // transform fluxes of subelement to distribution variables
  *phi = *getMethodData().getSolutionToDistribMatTrans()->transformFromRef(m_phisubT[1]);

  currBetaMatrix = &getMethodData().getDistributionData().betaMats[1];
  cf_assert(currBetaMatrix.isNotNull());

  distributeLDA(m_k2Plus,m_phi);


  residual[4] += m_theta2*m_phiN2[0] + (1.0 - m_theta2)*m_phi[0];
  residual[1] += m_theta2*m_phiN2[1] + (1.0 - m_theta2)*m_phi[1];
  residual[5] += m_theta2*m_phiN2[2] + (1.0 - m_theta2)*m_phi[2];

  /*****         Triangle 7-6-2          *****/
  substates[0] = states[7];
  substates[1] = states[6];
  substates[2] = states[2];

  getMethodData().getDistributionData().tStates = computeConsistentStates(&substates);

  ddata.subStates = &substates;

  // transform fluxes of subelement to distribution variables
  *phi = *getMethodData().getSolutionToDistribMatTrans()->transformFromRef(m_phisubT[2]);

  currBetaMatrix = &getMethodData().getDistributionData().betaMats[2];
  cf_assert(currBetaMatrix.isNotNull());

  distributeLDA(m_k3Plus,m_phi);

  residual[7] += m_theta3*m_phiN3[0] + (1.0 - m_theta3)*m_phi[0];
  residual[6] += m_theta3*m_phiN3[1] + (1.0 - m_theta3)*m_phi[1];
  residual[2] += m_theta3*m_phiN3[2] + (1.0 - m_theta3)*m_phi[2];

  /*****         Triangle 9-5-6          *****/

  substates[0] = states[9];
  substates[1] = states[5];
  substates[2] = states[6];

  getMethodData().getDistributionData().tStates = computeConsistentStates(&substates);

  *phi = *getMethodData().getSolutionToDistribMatTrans()->transformFromRef(m_phisubT[3]);
  currBetaMatrix = &getMethodData().getDistributionData().betaMats[3];
  cf_assert(currBetaMatrix.isNotNull());

  distributeLDA(m_k4Plus,m_phi);

  residual[9] += m_theta4*m_phiN4[0] + (1.0 - m_theta4)*m_phi[0];
  residual[5] += m_theta4*m_phiN4[1] + (1.0 - m_theta4)*m_phi[1];
  residual[6] += m_theta4*m_phiN4[2] + (1.0 - m_theta4)*m_phi[2];

 /*****         Triangle 8-9-7          *****/

  substates[0] = states[8];
  substates[1] = states[9];
  substates[2] = states[7];

  getMethodData().getDistributionData().tStates = computeConsistentStates(&substates);

  ddata.subStates = &substates;

  // transform fluxes of subelement to distribution variables
  *phi = *getMethodData().getSolutionToDistribMatTrans()->transformFromRef(m_phisubT[4]);

  currBetaMatrix = &getMethodData().getDistributionData().betaMats[4];
  cf_assert(currBetaMatrix.isNotNull());

  distributeLDA(m_k5Plus,m_phi);

  residual[8] += m_theta5*m_phiN5[0] + (1.0 - m_theta5)*m_phi[0];
  residual[9] += m_theta5*m_phiN5[1] + (1.0 - m_theta5)*m_phi[1];
  residual[7] += m_theta5*m_phiN5[2] + (1.0 - m_theta5)*m_phi[2];


/*****         Triangle 3-4-9          *****/
  substates[0] = states[3];
  substates[1] = states[4];
  substates[2] = states[9];

  getMethodData().getDistributionData().tStates = computeConsistentStates(&substates);

  ddata.subStates = &substates;

  // transform fluxes of subelement to distribution variables
  *phi = *getMethodData().getSolutionToDistribMatTrans()->transformFromRef(m_phisubT[5]);

  currBetaMatrix = &getMethodData().getDistributionData().betaMats[5];
  cf_assert(currBetaMatrix.isNotNull());

  distributeLDA(m_k6Plus,m_phi);

  residual[3] += m_theta6*m_phiN6[0] + (1.0 - m_theta6)*m_phi[0];
  residual[4] += m_theta6*m_phiN6[1] + (1.0 - m_theta6)*m_phi[1];
  residual[9] += m_theta6*m_phiN6[2] + (1.0 - m_theta6)*m_phi[2];


/*****         Triangle 9-8-3          *****/
  substates[0] = states[9];
  substates[1] = states[8];
  substates[2] = states[3];

  getMethodData().getDistributionData().tStates = computeConsistentStates(&substates);

  ddata.subStates = &substates;

  // transform fluxes of subelement to distribution variables
  *phi = *getMethodData().getSolutionToDistribMatTrans()->transformFromRef(m_phisubT[6]);

  currBetaMatrix = &getMethodData().getDistributionData().betaMats[6];

  distributeLDA(m_k7Plus,m_phi);

  residual[9] += m_theta7*m_phiN7[0] + (1.0 - m_theta7)*m_phi[0];
  residual[8] += m_theta7*m_phiN7[1] + (1.0 - m_theta7)*m_phi[1];
  residual[3] += m_theta7*m_phiN7[2] + (1.0 - m_theta7)*m_phi[2];

/*****         Triangle 5-9-4          *****/
  substates[0] = states[5];
  substates[1] = states[9];
  substates[2] = states[4];

  getMethodData().getDistributionData().tStates = computeConsistentStates(&substates);

  ddata.subStates = &substates;

  // transform fluxes of subelement to distribution variables
  *phi = *getMethodData().getSolutionToDistribMatTrans()->transformFromRef(m_phisubT[7]);

  currBetaMatrix = &getMethodData().getDistributionData().betaMats[7];
  cf_assert(currBetaMatrix.isNotNull());

  distributeLDA(m_k8Plus,m_phi);

  residual[5] += m_theta8*m_phiN8[0] + (1.0 - m_theta8)*m_phi[0];
  residual[9] += m_theta8*m_phiN8[1] + (1.0 - m_theta8)*m_phi[1];
  residual[4] += m_theta8*m_phiN8[2] + (1.0 - m_theta8)*m_phi[2];

/*****         Triangle 6-7-9          *****/
  substates[0] = states[6];
  substates[1] = states[7];
  substates[2] = states[9];

  getMethodData().getDistributionData().tStates = computeConsistentStates(&substates);

  ddata.subStates = &substates;

  // transform fluxes of subelement to distribution variables
  *phi = *getMethodData().getSolutionToDistribMatTrans()->transformFromRef(m_phisubT[8]);

  currBetaMatrix = &getMethodData().getDistributionData().betaMats[8];
  cf_assert(currBetaMatrix.isNotNull());

  distributeLDA(m_k9Plus,m_phi);

  residual[6] += m_theta9*m_phiN9[0] + (1.0 - m_theta9)*m_phi[0];
  residual[7] += m_theta9*m_phiN9[1] + (1.0 - m_theta9)*m_phi[1];
  residual[9] += m_theta9*m_phiN9[2] + (1.0 - m_theta9)*m_phi[2];

  cellnormals.unscale();

//Then we use this theta to build the blending
}

//////////////////////////////////////////////////////////////////////////////

void HOCRDP3_BT_SysSplitStrategy::computeHOP3Fluctuation()
{

  cf_assert(qdstates.size() == 5); // only triags so three quadrature points per face

  vector<State*>& states = *getMethodData().getDistributionData().states;

  const State& state0 = *(states[0]);
  const State& state1 = *(states[1]);
  const State& state2 = *(states[2]);
  const State& state3 = *(states[3]);
  const State& state4 = *(states[4]);
  const State& state5 = *(states[5]);
  const State& state6 = *(states[6]);
  const State& state7 = *(states[7]);
  const State& state8 = *(states[8]);
  const State& state9 = *(states[9]);


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
  const CFreal nx1 = 3. * cellnormals.getNodalNormComp(0,XX);
  const CFreal nx2 = 3. * cellnormals.getNodalNormComp(1,XX);
  const CFreal nx3 = 3. * cellnormals.getNodalNormComp(2,XX);

  const CFreal ny1 = 3. * cellnormals.getNodalNormComp(0,YY);
  const CFreal ny2 = 3. * cellnormals.getNodalNormComp(1,YY);
  const CFreal ny3 = 3. * cellnormals.getNodalNormComp(2,YY);


  const CFreal inv_volume = 1.0 / getMethodData().getDistributionData().cell->computeVolume();
// CF_DEBUG_POINT;

 // Computation  of the fluctuation of each faces of the cell
for (CFuint iFace = 0; iFace < subfacetable.nbRows(); ++iFace)
  {
    facenormal[XX] = cellnormals.getNodalNormComp(iFace%3,XX);
    facenormal[YY] = cellnormals.getNodalNormComp(iFace%3,YY);


    for (CFuint iQd = 0; iQd < 5; ++iQd)
    {

      const Node& node0 = states[subfacetable(iFace,0)]->getCoordinates();
      const Node& node1 = states[subfacetable(iFace,1)]->getCoordinates();
      const CFreal x = qd0[iQd] * node0[XX] + qd1[iQd] * node1[XX];
      const CFreal y = qd0[iQd] * node0[YY] + qd1[iQd] * node1[YY];
      const CFreal L1 = 1.0 + 0.5*( ( x - x1 )*nx1 + ( y - y1 )*ny1 ) * inv_volume ;
      const CFreal L2 = 1.0 + 0.5*( ( x - x2 )*nx2 + ( y - y2 )*ny2 ) * inv_volume ;
      const CFreal L3 = 1.0 + 0.5*( ( x - x3 )*nx3 + ( y - y3 )*ny3 ) * inv_volume ;

      const CFreal a0 = 0.5*( 3.0*L1 - 1.0 )*( 3.0*L1 - 2.0 )*L1 ;
      const CFreal a1 = 0.5*( 3.0*L2 - 1.0 )*( 3.0*L2 - 2.0 )*L2 ;
      const CFreal a2 = 0.5*( 3.0*L3 - 1.0 )*( 3.0*L3 - 2.0 )*L3 ;

      const CFreal a3 = 4.5*L1*L2*( 3.0*L1 - 1.0 ) ;
      const CFreal a4 = 4.5*L1*L2*( 3.0*L2 - 1.0 ) ;

      const CFreal a5 = 4.5*L2*L3*( 3.0*L2 - 1.0 ) ;
      const CFreal a6 = 4.5*L2*L3*( 3.0*L3 - 1.0 ) ;

      const CFreal a7 = 4.5*L3*L1*( 3.0*L3 - 1.0 ) ;
      const CFreal a8 = 4.5*L3*L1*( 3.0*L1 - 1.0 ) ;

      const CFreal a9 = 27.0*L1*L2*L3;

      (*qdstates[iQd]) = a0*state0 + a1*state1 + a2*state2 + a3*state3 +
                         a4*state4 + a5*state5 + a6*state6 + a7*state7 +
                         a8*state8 + a9*state9;

      (*qdnodes[iQd])[XX] = x;
      (*qdnodes[iQd])[YY] = y;
    }
    
    computeStatesData(5, m_updateVar, qdstates, m_pdata, m_qdExtraVars); // three quadrature points per face

    faceflux[iFace] = 0.;
    for (CFuint iQd = 0; iQd < 5; ++iQd)
    {

      faceflux[iFace] += wqd[iQd] * m_updateVar->getFlux()((*qdstates[iQd]),facenormal);
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

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<Framework::BaseDataSocketSink> >
HOCRDP3_BT_SysSplitStrategy::needsSockets()
{
  std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > result
    = FluctuationSplitStrategy::needsSockets();
  result.push_back(&socket_updateCoeff);
   result.push_back(&socket_thetas);
  return result;
}

//////////////////////////////////////////////////////////////////////////////////////
void HOCRDP3_BT_SysSplitStrategy::computeK(const std::vector<Framework::State*>& states,
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
void HOCRDP3_BT_SysSplitStrategy::distributeN(std::vector<RealMatrix*> & m_kPlus,vector<RealVector> & phiN){

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
void HOCRDP3_BT_SysSplitStrategy::distributeLDA(std::vector<RealMatrix*> & m_kPlus,vector<RealVector> & phiLDA){

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
void HOCRDP3_BT_SysSplitStrategy::computeBlendingCoeff(CFreal & result)
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
  
  CFreal sc = 0.;
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
}
  //second test
  result = min(1.0,(max(0.0, sc))*(max(0.0, sc))*h);

  // this->_alpha = min(1.0,theta*h);
  // cout << "theta = " << this->_alpha << endl;
}
//////////////////////////////////////////////////////////////////////////////
    }//end namespace FluctSplit

}// end namespace COOLFluiD
