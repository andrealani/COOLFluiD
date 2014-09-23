
#include "Framework/MethodStrategyProvider.hh"
#include "Framework/ContourIntegrator.hh"
#include "Framework/MeshData.hh"
#include "Framework/SubSystemStatus.hh"

#include "FluctSplit/FluctSplit.hh"
#include "FluctSplit/UnstP2SUPG_ArtDiffStrategy.hh"

#include "MathTools/MatrixInverter.hh"

#include "NavierStokes/EulerVarSet.hh"
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

MethodStrategyProvider<UnstP2SUPG_ArtDiffStrategy,
                       FluctuationSplitData,
                       ArtificialDiffusionStrategy,
                       FluctSplitModule>
                       unstP2supgArtDiffStrategyProvider("UnstP2SUPG");

//////////////////////////////////////////////////////////////////////////////

void UnstP2SUPG_ArtDiffStrategy::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< bool >("WithShockDetection","Add shock detection to the SUPG artificial diffusion.");
  options.addConfigOption< bool >("SimpScaling","Use a simplified scaling.");
}

//////////////////////////////////////////////////////////////////////////////

UnstP2SUPG_ArtDiffStrategy::UnstP2SUPG_ArtDiffStrategy(const std::string& name) :
  ArtificialDiffusionStrategy(name),
  socket_pastStates("pastStates"),
  socket_interStates("interStates"),
  m_pastStates(0),
  m_states(0),
  m_interStates(0),
  m_linstate(0),
  m_nabla_u0N1(0),
  m_nabla_u1N1(0),
  m_nabla_u2N1(0),
  m_nabla_u0n(0),
  m_nabla_u1n(0),
  m_nabla_u2n(0),
  m_nabla_u0n1(0),
  m_nabla_u1n1(0),
  m_nabla_u2n1(0),
  m_Res0(0),
  m_Res1(0),
  m_Res2(0),
  m_tempmat(0,0),
  m_kPlus(0),
  m_inverter(CFNULL),
  m_sumKplus(0,0),
  m_tau(0,0),
  m_adimNormal(0),
  m_eValues(0),
  m_rightEv(0,0),
  m_leftEv(0,0),
  m_eValuesP(0),
  _cterm(CFNULL),
  m_theta(0),
  m_min_states(0),
  m_max_states(0),
  m_PastRes0(0),
  m_PastRes1(0),
  m_PastRes2(0),
  m_nbEqs(0)
{
  CFAUTOTRACE;
  addConfigOptionsTo(this);

  m_with_shock_detect = true;
  setParameter("WithShockDetection",&m_with_shock_detect);
  // if it is on false then we use tau=(sumkplus)^-1 otherwise we use tau=alpha*Id
  // which is less computing demanding
  m_simplified_scaling = false;
  setParameter("SimpScaling", &m_simplified_scaling);

}


//////////////////////////////////////////////////////////////////////////////

UnstP2SUPG_ArtDiffStrategy::~UnstP2SUPG_ArtDiffStrategy()
{
}

//////////////////////////////////////////////////////////////////////////////

void UnstP2SUPG_ArtDiffStrategy::setup()
{
  CFAUTOTRACE;

  // first call parent method
  ArtificialDiffusionStrategy::setup();

  m_nbEqs = PhysicalModelStack::getActive()->getNbEq();
  CFuint dim   = PhysicalModelStack::getActive()->getDim();

 const CFuint maxNbStatesInCell = MeshDataStack::getActive()->Statistics().getMaxNbStatesInCell();

  m_interStates.resize(maxNbStatesInCell);
  m_states.resize(maxNbStatesInCell);
  m_pastStates.resize(maxNbStatesInCell);
  // Resizing all the states of the past and temporary states
  for (CFuint i = 0; i < maxNbStatesInCell; ++i) {
    m_interStates[i] = new State();
    m_states[i] = new State();
    m_pastStates[i] = new State();
   }


  m_linstate.resize(maxNbStatesInCell);
  for (CFuint i = 0; i < maxNbStatesInCell; ++i) {
    m_linstate[i] = new State();
  }

  m_nabla_u0N1.resize(dim);
  m_nabla_u1N1.resize(dim);
  m_nabla_u2N1.resize(dim);
  m_nabla_u0n.resize(dim);
  m_nabla_u1n.resize(dim);
  m_nabla_u2n.resize(dim);
  m_nabla_u0n1.resize(dim);
  m_nabla_u1n1.resize(dim);
  m_nabla_u2n1.resize(dim);

  for (CFuint i = 0; i < dim; ++i) {
             m_nabla_u0N1[i].resize(m_nbEqs);
             m_nabla_u1N1[i].resize(m_nbEqs);
             m_nabla_u2N1[i].resize(m_nbEqs);
             m_nabla_u0n[i].resize(m_nbEqs);
             m_nabla_u1n[i].resize(m_nbEqs);
             m_nabla_u2n[i].resize(m_nbEqs);
             m_nabla_u0n1[i].resize(m_nbEqs);
             m_nabla_u1n1[i].resize(m_nbEqs);
             m_nabla_u2n1[i].resize(m_nbEqs);
  }

  m_Res0.resize(m_nbEqs);
  m_Res1.resize(m_nbEqs);
  m_Res2.resize(m_nbEqs);

  m_tempmat.resize(m_nbEqs,m_nbEqs);

  m_kPlus.resize(3);

  for (CFuint i = 0; i < 3; ++i) {
    m_kPlus[i] = new RealMatrix(m_nbEqs, m_nbEqs);

  }


  m_inverter = MatrixInverter::create(m_nbEqs, false);

   m_sumKplus.resize(m_nbEqs,m_nbEqs);

  m_tau.resize(m_nbEqs,m_nbEqs);

  m_adimNormal.resize(dim);

  m_eValues.resize(m_nbEqs);
  m_rightEv.resize(m_nbEqs,m_nbEqs);
  m_leftEv.resize(m_nbEqs,m_nbEqs);
  m_eValuesP.resize(m_nbEqs);

_cterm = Framework::PhysicalModelStack::getActive()->getImplementor()->
    getConvectiveTerm().d_castTo<Physics::NavierStokes::EulerTerm>();

  m_theta.resize(m_nbEqs);
  m_max_states.resize(m_nbEqs);
  m_min_states.resize(m_nbEqs);

  const CFuint nbCells =
   MeshDataStack::getActive()->getTrs("InnerCells")->getLocalNbGeoEnts();

  // Resizing the storage for the full past Residuals
  m_PastRes0.resize(nbCells);
  m_PastRes1.resize(nbCells);
  m_PastRes2.resize(nbCells);


  for(CFuint iCell = 0; iCell < nbCells; ++iCell) {

    m_PastRes0[iCell].resize(m_nbEqs);
    m_PastRes1[iCell].resize(m_nbEqs);
    m_PastRes2[iCell].resize(m_nbEqs);

  }
 }

//////////////////////////////////////////////////////////////////////////////

void UnstP2SUPG_ArtDiffStrategy::unsetup()
{
  CFAUTOTRACE;
  for (CFuint i = 0; i < m_kPlus.size(); ++i)
  {
    deletePtr(m_kPlus[i]);
  }
  deletePtr(m_inverter);

 for (CFuint i = 0; i < m_interStates.size(); ++i)
  {
    deletePtr(m_interStates[i]);
    deletePtr(m_states[i]);
    deletePtr(m_pastStates[i]);
  }

}

//////////////////////////////////////////////////////////////////////////////

void UnstP2SUPG_ArtDiffStrategy::addArtificialDiff(std::vector<RealVector>& residual)
{

// The first iteration needs a special function because at this point
  // only two layer are available althought we need three.

  CFuint nbIter = SubSystemStatusStack::getActive()->getNbIter();/*
if(nbIter != 1) doaddArtificialDiff(residual);*/

  if(nbIter == 1) doFirstaddArtificialDiff(residual);
  else doaddArtificialDiff(residual);

}

//////////////////////////////////////////////////////////////////////////////
void UnstP2SUPG_ArtDiffStrategy::doaddArtificialDiff(std::vector<RealVector>& residual)
{
  DistributionData& ddata = getMethodData().getDistributionData();
  vector<State*>& states = *ddata.states;
  const CFuint nbStatesInCell = states.size();
  const CFreal volume = ddata.cell->computeVolume();
  const CFreal gradcoeff = 1.0 / (2.0*volume);
  const CFreal onethirdArea = volume/3.0;
  const RealVector& lData =  _cterm->getPhysicalData();
  const CFuint CellID = ddata.cellID;
  InwardNormalsData& cellnormals = (*socket_normals.getDataHandle()
				    [ddata.cellID]);
  const CFreal nx1 = cellnormals.getNodalNormComp(0,XX);
  const CFreal nx2 = cellnormals.getNodalNormComp(1,XX);
  const CFreal nx3 = cellnormals.getNodalNormComp(2,XX);

  const CFreal ny1 = cellnormals.getNodalNormComp(0,YY);
  const CFreal ny2 = cellnormals.getNodalNormComp(1,YY);
  const CFreal ny3 = cellnormals.getNodalNormComp(2,YY);

  const CFreal dt1 = SubSystemStatusStack::getActive()->getDT();
  const CFreal dt0 =SubSystemStatusStack::getActive()->getPreviousDT();
  const CFreal q = dt1/dt0;
  const CFreal Jpast  = -dt1*(q*q)/(6.0*(1.0 + q));
  const CFreal Jinter  = dt1*(3.0 + q)/6.0;
  const CFreal Jpresent  = dt1*(3.0 + 2.0*q)/(6.0*(1.0 + q));

  // Every thing that concern intermediate and past states is done only once
  if (SubSystemStatusStack::getActive()->isFirstStep()){

    // Get the past states
    DataHandle<State*> pastStatesStorage = socket_pastStates.getDataHandle();
    for (CFuint i = 0; i < nbStatesInCell; ++i) {
      const CFuint stateID = states[i]->getLocalID();
      *m_pastStates[i] = *pastStatesStorage[stateID];
    }


    // Get the intemediate states
    DataHandle<State*> interStatesStorage = socket_interStates.getDataHandle();
    for (CFuint i = 0; i < nbStatesInCell; ++i) {
      const CFuint stateID = states[i]->getLocalID();
      *m_interStates[i] = *interStatesStorage[stateID];
    }

    // Compute the Jacobian at node0 time n
    for (CFuint i = 0; i < nbStatesInCell; ++i) {
      (*m_linstate[i]) = (*m_interStates[0]);
    }


    vector<State*> * toto = computeConsistentStates(&m_linstate);
    _distribVar->computeJacobians();
    const std::vector<RealMatrix> jacob0_n = *PhysicalModelStack::getActive()->getImplementor()->getJacobians();

    // Compute the Jacobian at node1 time n
    for (CFuint i = 0; i < nbStatesInCell; ++i) {
      (*m_linstate[i]) = (*m_interStates[1]);
    }
    toto = computeConsistentStates(&m_linstate);
    _distribVar->computeJacobians();
    const std::vector<RealMatrix> jacob1_n = *PhysicalModelStack::getActive()->getImplementor()->getJacobians();


    // Compute the Jacobian at node2 time n
    for (CFuint i = 0; i < nbStatesInCell; ++i) {
      (*m_linstate[i]) = (*m_interStates[2]);
    }
    toto = computeConsistentStates(&m_linstate);
    _distribVar->computeJacobians();
    const std::vector<RealMatrix> jacob2_n = *PhysicalModelStack::getActive()->getImplementor()->getJacobians();

    // Compute the Jacobian at node0 time n-1
    for (CFuint i = 0; i < nbStatesInCell; ++i) {
      (*m_linstate[i]) = (*m_pastStates[0]);
    }

    toto = computeConsistentStates(&m_linstate);
    _distribVar->computeJacobians();
    const std::vector<RealMatrix> jacob0_n1 = *PhysicalModelStack::getActive()->getImplementor()->getJacobians();

    // Compute the Jacobian at node1 time n-1
    for (CFuint i = 0; i < nbStatesInCell; ++i) {
      (*m_linstate[i]) = (*m_pastStates[1]);
    }

    toto = computeConsistentStates(&m_linstate);
    _distribVar->computeJacobians();
    const std::vector<RealMatrix> jacob1_n1 = *PhysicalModelStack::getActive()->getImplementor()->getJacobians();

    // Compute the Jacobian at node2 time n-1
    for (CFuint i = 0; i < nbStatesInCell; ++i) {
      (*m_linstate[i]) = (*m_pastStates[2]);
    }
    toto = computeConsistentStates(&m_linstate);
    _distribVar->computeJacobians();
    const std::vector<RealMatrix> jacob2_n1 = *PhysicalModelStack::getActive()->getImplementor()->getJacobians();

    // Compute the tStates of each time level
    vector<State*> *const tPastStates = computeConsistentStates(&m_pastStates);
    vector<State*> *const tInterStates = computeConsistentStates(&m_interStates);

    vector<State*> & m_tPastStates = *tPastStates;
    vector<State*> & m_tInterStates = *tInterStates;

    // Compute the gradient of u at node0 time n
    m_nabla_u0n[XX] = (3.0*nx1*(*m_tInterStates[0]) -     nx2*(*m_tInterStates[1]) - nx3*(*m_tInterStates[2]) +
                       4.0*nx2*(*m_tInterStates[3]) + 4.0*nx3*(*m_tInterStates[5]))*gradcoeff;
    m_nabla_u0n[YY] = (3.0*ny1*(*m_tInterStates[0]) -     ny2*(*m_tInterStates[1]) - ny3*(*m_tInterStates[2]) +
                       4.0*ny2*(*m_tInterStates[3]) + 4.0*ny3*(*m_tInterStates[5]))*gradcoeff;


   // Compute the gradient of u at node0 time n-1
   m_nabla_u0n1[XX] = (3.0*nx1*(*m_tPastStates[0]) -     nx2*(*m_tPastStates[1]) - nx3*(*m_tPastStates[2]) +
                       4.0*nx2*(*m_tPastStates[3]) + 4.0*nx3*(*m_tPastStates[5]))*gradcoeff;
   m_nabla_u0n1[YY] = (3.0*ny1*(*m_tPastStates[0]) -     ny2*(*m_tPastStates[1]) - ny3*(*m_tPastStates[2]) +
                       4.0*ny2*(*m_tPastStates[3]) + 4.0*ny3*(*m_tPastStates[5]))*gradcoeff;

 // Compute the gradient of u at node1 time n
   m_nabla_u1n[XX] =     (-nx1*(*m_tInterStates[0]) + 3.0*nx2*(*m_tInterStates[1]) - nx3*(*m_tInterStates[2]) +
                      4.0* nx1*(*m_tInterStates[3]) + 4.0*nx3*(*m_tInterStates[4]))*gradcoeff;
   m_nabla_u1n[YY] =     (-ny1*(*m_tInterStates[0]) + 3.0*ny2*(*m_tInterStates[1]) - ny3*(*m_tInterStates[2]) +
                      4.0* ny1*(*m_tInterStates[3]) + 4.0*ny3*(*m_tInterStates[4]))*gradcoeff;


   // Compute the gradient of u at node1 time n-1
   m_nabla_u1n1[XX] =    (-nx1*(*m_tPastStates[0]) + 3.0*nx2*(*m_tPastStates[1]) - nx3*(*m_tPastStates[2]) +
                      4.0*nx1*(*m_tPastStates[3]) + 4.0*nx3*(*m_tPastStates[4]))*gradcoeff;
   m_nabla_u1n1[YY] =    (-ny1*(*m_tPastStates[0]) + 3.0*ny2*(*m_tPastStates[1]) - ny3*(*m_tPastStates[2]) +
                      4.0*ny1*(*m_tPastStates[3]) + 4.0*ny3*(*m_tPastStates[4]))*gradcoeff;

  // Compute the gradient of u at node2 time n
   m_nabla_u2n[XX] =     (-nx1*(*m_tInterStates[0]) -     nx2*(*m_tInterStates[1]) + 3.0*nx3*(*m_tInterStates[2]) +
                      4.0*nx2*(*m_tInterStates[4]) + 4.0*nx1*(*m_tInterStates[5]))*gradcoeff;
   m_nabla_u2n[YY] =     (-ny1*(*m_tInterStates[0]) -     ny2*(*m_tInterStates[1]) + 3.0*ny3*(*m_tInterStates[2]) +
                      4.0*ny2*(*m_tInterStates[4]) + 4.0*ny1*(*m_tInterStates[5]))*gradcoeff;


   // Compute the gradient of u at node2 time n-1
   m_nabla_u2n1[XX] =    (-nx1*(*m_tPastStates[0]) -     nx2*(*m_tPastStates[1]) + 3.0*nx3*(*m_tPastStates[2]) +
                      4.0*nx2*(*m_tPastStates[4]) + 4.0*nx1*(*m_tPastStates[5]))*gradcoeff;
   m_nabla_u2n1[YY] =    (-ny1*(*m_tPastStates[0]) -     ny2*(*m_tPastStates[1]) + 3.0*ny3*(*m_tPastStates[2]) +
                      4.0*ny2*(*m_tPastStates[4]) + 4.0*ny1*(*m_tPastStates[5]))*gradcoeff;

    // Compute the residual in node0, node1 and node2
    m_Res0 = (-1.0)*(*m_interStates[0]) +
    Jpast    * (jacob0_n1[XX] * m_nabla_u0n1[XX] + jacob0_n1[YY] * m_nabla_u0n1[YY]) +
    Jinter   * (jacob0_n[XX]  * m_nabla_u0n[XX]  + jacob0_n[YY]  * m_nabla_u0n[YY]) ;

    m_Res1 = (-1.0)*(*m_interStates[1]) +
    Jpast    * (jacob1_n1[XX] * m_nabla_u1n1[XX] + jacob1_n1[YY] * m_nabla_u1n1[YY]) +
    Jinter   * (jacob1_n[XX]  * m_nabla_u1n[XX]  + jacob1_n[YY]  * m_nabla_u1n[YY]);

    m_Res2 = (-1.0)*(*m_interStates[2]) +
    Jpast    * (jacob2_n1[XX] * m_nabla_u2n1[XX] + jacob2_n1[YY] * m_nabla_u2n1[YY]) +
    Jinter   * (jacob2_n[XX]  * m_nabla_u2n[XX]  + jacob2_n[YY]  * m_nabla_u2n[YY]) ;

    for (CFuint iEq = 0; iEq < m_nbEqs; ++iEq){
        m_PastRes0[CellID][iEq] = m_Res0[iEq];
        m_PastRes1[CellID][iEq] = m_Res1[iEq];
        m_PastRes2[CellID][iEq] = m_Res2[iEq];
    }

  }

  if (!SubSystemStatusStack::getActive()->isFirstStep()){
    for (CFuint iEq = 0; iEq < m_nbEqs; ++iEq){
      m_Res0[iEq] = m_PastRes0[CellID][iEq];
      m_Res1[iEq] = m_PastRes1[CellID][iEq];
      m_Res2[iEq] = m_PastRes2[CellID][iEq];
    }
  }

  // Compute the Jacobian at node0 time n+1
  for (CFuint i = 0; i < nbStatesInCell; ++i) {
    (*m_linstate[i]) = (*states[0]);
  }
  vector<State*> * toto = computeConsistentStates(&m_linstate);
  _distribVar->computeJacobians();
  const std::vector<RealMatrix> jacob0_N1 = *PhysicalModelStack::getActive()->getImplementor()->getJacobians();

  // Compute the Jacobian at node1 time n+1
  for (CFuint i = 0; i < nbStatesInCell; ++i) {
    (*m_linstate[i]) = (*states[1]);
  }
  toto = computeConsistentStates(&m_linstate);
  _distribVar->computeJacobians();
  const std::vector<RealMatrix> jacob1_N1 = *PhysicalModelStack::getActive()->getImplementor()->getJacobians();

  // Compute the Jacobian at node2 time n+1
  for (CFuint i = 0; i < nbStatesInCell; ++i) {
    (*m_linstate[i]) = (*states[2]);
  }
  toto = computeConsistentStates(&m_linstate);
  _distribVar->computeJacobians();
  const std::vector<RealMatrix> jacob2_N1 = *PhysicalModelStack::getActive()->getImplementor()->getJacobians();

  // Compute the tStates of each time level
  vector<State*> *const tStates = computeConsistentStates(&states);

  vector<State*> & m_tStates = *tStates;

  const CFreal rho = lData[EulerTerm::RHO];
  const CFreal P = lData[EulerTerm::P];
  const CFreal vel2 = lData[EulerTerm::V];
  const CFreal gamma = lData[EulerTerm::GAMMA];
  const CFreal c = sqrt(gamma*P/rho);

  //Compute the scaling coefficient
  if(m_simplified_scaling){

    CFreal m_alpha = 0.0;
    for (CFuint iState = 0; iState < 3; ++iState){
      for (CFuint iDim = 0; iDim < 2; ++iDim) {
        m_adimNormal[iDim] = cellnormals.getNodalNormComp(iState, iDim);
      }
    CFreal m_nodeArea = cellnormals.getAreaNode(iState);
    m_adimNormal *= 1. / m_nodeArea;

    getMethodData().getDistribVar()->computeEigenValuesVectors(m_rightEv, m_leftEv, m_eValues,m_adimNormal);

    for(CFuint iEq=0; iEq<m_nbEqs; ++iEq)
      m_alpha = max(m_alpha,m_eValues[iEq]);
    }

    m_tau = 0.0;
    for(CFuint iEq=0; iEq<m_nbEqs; ++iEq)
      m_tau(iEq,iEq) = 1.0/m_alpha;

    //m_tau = m_identity/m_alpha;
  }
  else{
    for (CFuint iState = 0; iState < 3; ++iState){
      for (CFuint iDim = 0; iDim < 2; ++iDim) {
        m_adimNormal[iDim] = cellnormals.getNodalNormComp(iState, iDim);
      }
      CFreal m_nodeArea = cellnormals.getAreaNode(iState);
      m_adimNormal *= 1. / m_nodeArea;

      getMethodData().getDistribVar()->computeEigenValuesVectors(m_rightEv, m_leftEv, m_eValues,m_adimNormal);


      CFreal epsilon = 0.01*(vel2+c);

      for(CFuint iEq=0; iEq<m_nbEqs; ++iEq){
        if (m_eValues[iEq] > epsilon){
          m_eValuesP[iEq] = max(0.,m_eValues[iEq]);
        }
        else {
          m_eValuesP[iEq] = (m_eValues[iEq] + epsilon)*(m_eValues[iEq] + epsilon)/(4.0*epsilon);
        }
      }

      *m_kPlus[iState] = m_rightEv*(m_eValuesP*m_leftEv);
    }

    m_sumKplus = *m_kPlus[0];
    for (CFuint iState = 1; iState < 3; ++iState){
      m_sumKplus += *m_kPlus[iState];
    }

    m_inverter->invert(m_sumKplus, m_tau);
  }
  // Compute the gradient of u at node0 time n+1
  m_nabla_u0N1[XX] = (3.0*nx1*(*m_tStates[0]) - nx2*(*m_tStates[1]) - nx3*(*m_tStates[2]) +
                      4.0*nx2*(*m_tStates[3]) + 4.0*nx3*(*m_tStates[5]))*gradcoeff;
  m_nabla_u0N1[YY] = (3.0*ny1*(*m_tStates[0]) - ny2*(*m_tStates[1]) - ny3*(*m_tStates[2]) +
                      4.0*ny2*(*m_tStates[3]) + 4.0*ny3*(*m_tStates[5]))*gradcoeff;


  // Compute the gradient of u at node1 time n+1
  m_nabla_u1N1[XX] =    (-nx1*(*m_tStates[0]) + 3.0*nx2*(*m_tStates[1]) - nx3*(*m_tStates[2]) +
                      4.0*nx1*(*m_tStates[3]) + 4.0*nx3*(*m_tStates[4]))*gradcoeff;
  m_nabla_u1N1[YY] =    (-ny1*(*m_tStates[0]) + 3.0*ny2*(*m_tStates[1]) - ny3*(*m_tStates[2]) +
                      4.0*ny1*(*m_tStates[3]) + 4.0*ny3*(*m_tStates[4]))*gradcoeff;



  // Compute the gradient of u at node2 time n+1
  m_nabla_u2N1[XX] =    (-nx1*(*m_tStates[0]) -     nx2*(*m_tStates[1]) + 3.0*nx3*(*m_tStates[2]) +
                      4.0*nx2*(*m_tStates[4]) + 4.0*nx1*(*m_tStates[5]))*gradcoeff;
  m_nabla_u2N1[YY] =    (-ny1*(*m_tStates[0]) -     ny2*(*m_tStates[1]) + 3.0*ny3*(*m_tStates[2]) +
                      4.0*ny2*(*m_tStates[4]) + 4.0*ny1*(*m_tStates[5]))*gradcoeff;



  // Compute the residual in node0, node1 and node2
  m_Res0 += (*states[0]) +
            Jpresent * (jacob0_N1[XX] * m_nabla_u0N1[XX] + jacob0_N1[YY] * m_nabla_u0N1[YY]) ;

  m_Res1 += (*states[1]) +
            Jpresent * (jacob1_N1[XX] * m_nabla_u1N1[XX] + jacob1_N1[YY] * m_nabla_u1N1[YY]) ;

  m_Res2 += (*states[2]) +
            Jpresent * (jacob2_N1[XX] * m_nabla_u2N1[XX] + jacob2_N1[YY] * m_nabla_u2N1[YY]) ;

  //Dissipation for node 0
  m_tempmat = sqrt(volume)*gradcoeff*( 3.0)*(nx1 * jacob0_N1[XX] + ny1 * jacob0_N1[YY])*m_tau;
  residual[0] = m_tempmat*m_Res0*onethirdArea;
  m_tempmat = sqrt(volume)*gradcoeff*(-1.0)*(nx1 * jacob1_N1[XX] + ny1 * jacob1_N1[YY])*m_tau;
  residual[0] += m_tempmat*m_Res1*onethirdArea;
  m_tempmat = sqrt(volume)*gradcoeff*(-1.0)*(nx1 * jacob2_N1[XX] + ny1 * jacob2_N1[YY])*m_tau;
  residual[0] += m_tempmat*m_Res2*onethirdArea;

  //Dissipation for node 1
  m_tempmat = sqrt(volume)*gradcoeff*(-1.0)*(nx2 * jacob0_N1[XX] + ny2 * jacob0_N1[YY])*m_tau;
  residual[1] = m_tempmat*m_Res0*onethirdArea;
  m_tempmat = sqrt(volume)*gradcoeff*( 3.0)*(nx2 * jacob1_N1[XX] + ny2 * jacob1_N1[YY])*m_tau;
  residual[1] += m_tempmat*m_Res1*onethirdArea;
  m_tempmat = sqrt(volume)*gradcoeff*(-1.0)*(nx2 * jacob2_N1[XX] + ny2 * jacob2_N1[YY])*m_tau;
  residual[1] += m_tempmat*m_Res2*onethirdArea;

  //Dissipation for node 2
  m_tempmat = sqrt(volume)*gradcoeff*(-1.0)*(nx3 * jacob0_N1[XX] + ny3 * jacob0_N1[YY])*m_tau;
  residual[2] = m_tempmat*m_Res0*onethirdArea;
  m_tempmat = sqrt(volume)*gradcoeff*(-1.0)*(nx3 * jacob1_N1[XX] + ny3 * jacob1_N1[YY])*m_tau;
  residual[2] += m_tempmat*m_Res1*onethirdArea;
  m_tempmat = sqrt(volume)*gradcoeff*( 3.0)*(nx3 * jacob2_N1[XX] + ny3 * jacob2_N1[YY])*m_tau;
  residual[2] += m_tempmat*m_Res2*onethirdArea;

  //Dissipation for node 3
  m_tempmat = sqrt(volume)*gradcoeff*( 4.0)*(nx2 * jacob0_N1[XX] + ny2 * jacob0_N1[YY])*m_tau;
  residual[3] = m_tempmat*m_Res0*onethirdArea;
  m_tempmat = sqrt(volume)*gradcoeff*( 4.0)*(nx1 * jacob1_N1[XX] + ny1 * jacob1_N1[YY])*m_tau;
  residual[3] += m_tempmat*m_Res1*onethirdArea;

  //Dissipation for node 4
  m_tempmat = sqrt(volume)*gradcoeff*( 4.0)*(nx3 * jacob1_N1[XX] + ny3 * jacob1_N1[YY])*m_tau;
  residual[4] = m_tempmat*m_Res1*onethirdArea;
  m_tempmat = sqrt(volume)*gradcoeff*( 4.0)*(nx2 * jacob2_N1[XX] + ny2 * jacob2_N1[YY])*m_tau;
  residual[4] += m_tempmat*m_Res2*onethirdArea;

  //Dissipation for node 5
  m_tempmat = sqrt(volume)*gradcoeff*( 4.0)*(nx3 * jacob0_N1[XX] + ny3 * jacob0_N1[YY])*m_tau;
  residual[5] = m_tempmat*m_Res0*onethirdArea;
  m_tempmat = sqrt(volume)*gradcoeff*( 4.0)*(nx1 * jacob2_N1[XX] + ny1 * jacob2_N1[YY])*m_tau;
  residual[5] += m_tempmat*m_Res2*onethirdArea;


  CFuint nbStates = states.size();

  // by default is one, meaning no shock detection
  m_theta = 1.0;
  if (m_with_shock_detect)
  {

    m_min_states =   MathConsts::CFrealMax();
    m_max_states = - MathConsts::CFrealMax();

    for(CFuint iEq = 0; iEq<m_nbEqs; ++iEq){
      for (CFuint iState=0;iState<nbStates; ++iState){
        m_min_states[iEq] = std::min((*states[iState])[iEq],m_min_states[iEq]);
        m_max_states[iEq] = std::max((*states[iState])[iEq],m_max_states[iEq]);
      }

      m_theta[iEq] -= std::abs(m_max_states[iEq] - m_min_states[iEq])/(std::abs(m_max_states[iEq]) + std::abs(m_min_states[iEq])+10.0E-10);
    }
  }

  CFreal m_thetamin = m_theta.min();
  for (CFuint iState=0;iState<nbStates; ++iState){
    residual[iState] *= m_thetamin;

  }
}

//////////////////////////////////////////////////////////////////////////////
void UnstP2SUPG_ArtDiffStrategy::doFirstaddArtificialDiff(std::vector<RealVector>& residual)
{
  DistributionData& ddata = getMethodData().getDistributionData();
  vector<State*>& states = *ddata.states;
  const CFuint nbStatesInCell = states.size();
  const CFreal volume = ddata.cell->computeVolume();
  const CFreal gradcoeff = 1.0 / (2.0*volume);
  const CFreal onethirdArea = volume/3.0;
  const CFreal cellID = ddata.cellID;
  const RealVector& lData =  _cterm->getPhysicalData();
  InwardNormalsData& cellnormals = (*socket_normals.getDataHandle()
				    [ddata.cellID]);
  const CFreal nx1 = cellnormals.getNodalNormComp(0,XX);
  const CFreal nx2 = cellnormals.getNodalNormComp(1,XX);
  const CFreal nx3 = cellnormals.getNodalNormComp(2,XX);

  const CFreal ny1 = cellnormals.getNodalNormComp(0,YY);
  const CFreal ny2 = cellnormals.getNodalNormComp(1,YY);
  const CFreal ny3 = cellnormals.getNodalNormComp(2,YY);

  const CFreal dt = SubSystemStatusStack::getActive()->getDT();
  const CFreal J  = dt/2.0;

  if (SubSystemStatusStack::getActive()->isFirstStep()){
    // Get the intemediate states
    DataHandle<State*> interStatesStorage = socket_interStates.getDataHandle();
    for (CFuint i = 0; i < nbStatesInCell; ++i) {
      const CFuint stateID = (*ddata.states)[i]->getLocalID();
      *m_interStates[i] = *interStatesStorage[stateID];
    }


    // Compute the Jacobian at node0 time n
    for (CFuint i = 0; i < nbStatesInCell; ++i) {
      (*m_linstate[i]) = (*m_interStates[0]);
    }
    vector<State*> * toto = computeConsistentStates(&m_linstate);
    _distribVar->computeJacobians();
    const std::vector<RealMatrix> jacob0_n = *PhysicalModelStack::getActive()->getImplementor()->getJacobians();


    // Compute the Jacobian at node1 time n
    for (CFuint i = 0; i < nbStatesInCell; ++i) {
      (*m_linstate[i]) = (*m_interStates[1]);
    }
    toto = computeConsistentStates(&m_linstate);
    _distribVar->computeJacobians();
    const std::vector<RealMatrix> jacob1_n = *PhysicalModelStack::getActive()->getImplementor()->getJacobians();


    // Compute the Jacobian at node2 time n
    for (CFuint i = 0; i < nbStatesInCell; ++i) {
      (*m_linstate[i]) = (*m_interStates[2]);
    }
    toto = computeConsistentStates(&m_linstate);
    _distribVar->computeJacobians();
    const std::vector<RealMatrix> jacob2_n = *PhysicalModelStack::getActive()->getImplementor()->getJacobians();

    vector<State*> *const tInterStates = computeConsistentStates(&m_interStates);

    vector<State*> & m_tInterStates = *tInterStates;

    // Compute the gradient of u at node0 time n
    m_nabla_u0n[XX] = (3.0*nx1*(*m_tInterStates[0]) -     nx2*(*m_tInterStates[1]) - nx3*(*m_tInterStates[2]) +
                       4.0*nx2*(*m_tInterStates[3]) + 4.0*nx3*(*m_tInterStates[5]))*gradcoeff;
    m_nabla_u0n[YY] = (3.0*ny1*(*m_tInterStates[0]) -     ny2*(*m_tInterStates[1]) - ny3*(*m_tInterStates[2]) +
                       4.0*ny2*(*m_tInterStates[3]) + 4.0*ny3*(*m_tInterStates[5]))*gradcoeff;

    // Compute the gradient of u at node1 time n
    m_nabla_u1n[XX] =   (-nx1*(*m_tInterStates[0]) + 3.0*nx2*(*m_tInterStates[1]) - nx3*(*m_tInterStates[2]) +
                      4.0*nx1*(*m_tInterStates[3]) + 4.0*nx3*(*m_tInterStates[4]))*gradcoeff;
    m_nabla_u1n[YY] =   (-ny1*(*m_tInterStates[0]) + 3.0*ny2*(*m_tInterStates[1]) - ny3*(*m_tInterStates[2]) +
                      4.0*ny1*(*m_tInterStates[3]) + 4.0*ny3*(*m_tInterStates[4]))*gradcoeff;

    // Compute the gradient of u at node2 time n
    m_nabla_u2n[XX] =   (-nx1*(*m_tInterStates[0]) -     nx2*(*m_tInterStates[1]) + 3.0*nx3*(*m_tInterStates[2]) +
                      4.0*nx2*(*m_tInterStates[4]) + 4.0*nx1*(*m_tInterStates[5]))*gradcoeff;
    m_nabla_u2n[YY] =   (-ny1*(*m_tInterStates[0]) -     ny2*(*m_tInterStates[1]) + 3.0*ny3*(*m_tInterStates[2]) +
                      4.0*ny2*(*m_tInterStates[4]) + 4.0*ny1*(*m_tInterStates[5]))*gradcoeff;

    m_Res0 = (-1.0)*(*m_interStates[0]) +
             J   * (jacob0_n[XX]  * m_nabla_u0n[XX]  + jacob0_n[YY]  * m_nabla_u0n[YY]) ;

    m_Res1 = (-1.0)*(*m_interStates[1]) +
             J   * (jacob1_n[XX]  * m_nabla_u1n[XX]  + jacob1_n[YY]  * m_nabla_u1n[YY]) ;

    m_Res2 = (-1.0)*(*m_interStates[2]) +
             J   * (jacob2_n[XX]  * m_nabla_u2n[XX]  + jacob2_n[YY]  * m_nabla_u2n[YY]) ;

    for (CFuint iEq = 0; iEq < m_nbEqs; ++iEq){
      m_PastRes0[cellID][iEq] = m_Res0[iEq];
      m_PastRes1[cellID][iEq] = m_Res1[iEq];
      m_PastRes2[cellID][iEq] = m_Res2[iEq];
    }

  }

  if (!SubSystemStatusStack::getActive()->isFirstStep()){
    for (CFuint iEq = 0; iEq < m_nbEqs; ++iEq){
        m_Res0[iEq] = m_PastRes0[cellID][iEq];
        m_Res1[iEq] = m_PastRes1[cellID][iEq];
        m_Res2[iEq] = m_PastRes2[cellID][iEq];
    }
  }
  // Compute the Jacobian at node0 time n+1
  for (CFuint i = 0; i < nbStatesInCell; ++i) {
    (*m_linstate[i]) = (*states[0]);
  }

  vector<State*> * toto = computeConsistentStates(&m_linstate);
  _distribVar->computeJacobians();
  const std::vector<RealMatrix> jacob0_N1 = *PhysicalModelStack::getActive()->getImplementor()->getJacobians();


  // Compute the Jacobian at node1 time n+1
  for (CFuint i = 0; i < nbStatesInCell; ++i) {
    (*m_linstate[i]) = (*states[1]);
  }
  toto = computeConsistentStates(&m_linstate);
  _distribVar->computeJacobians();
  const std::vector<RealMatrix> jacob1_N1 = *PhysicalModelStack::getActive()->getImplementor()->getJacobians();


  // Compute the Jacobian at node2 time n+1
  for (CFuint i = 0; i < nbStatesInCell; ++i) {
    (*m_linstate[i]) = (*states[2]);
  }
  toto = computeConsistentStates(&m_linstate);
  _distribVar->computeJacobians();
  const std::vector<RealMatrix> jacob2_N1 = *PhysicalModelStack::getActive()->getImplementor()->getJacobians();

  // Compute the tStates of each time level
  vector<State*> *const tStates = computeConsistentStates(&states);

  vector<State*> & m_tStates = *tStates;


  const CFreal rho = lData[EulerTerm::RHO];
  const CFreal P = lData[EulerTerm::P];
  const CFreal vel2 = lData[EulerTerm::V];
  const CFreal gamma = lData[EulerTerm::GAMMA];
  const CFreal c = sqrt(gamma*P/rho);
  //Compute the scaling coefficient
  if(m_simplified_scaling){
    CFreal m_alpha = 0.0;
    for (CFuint iState = 0; iState < 3; ++iState){
      for (CFuint iDim = 0; iDim < 2; ++iDim) {
        m_adimNormal[iDim] = cellnormals.getNodalNormComp(iState, iDim);
      }
    CFreal m_nodeArea = cellnormals.getAreaNode(iState);
    m_adimNormal *= 1. / m_nodeArea;

    getMethodData().getDistribVar()->computeEigenValuesVectors(m_rightEv, m_leftEv, m_eValues,m_adimNormal);

    for(CFuint iEq=0; iEq<m_nbEqs; ++iEq)
      m_alpha = max(m_alpha,m_eValues[iEq]);
    }
    m_tau = 0.0;
    for(CFuint iEq=0; iEq<m_nbEqs; ++iEq)
      m_tau(iEq,iEq) = 1.0/m_alpha;
  }
  else{
    for (CFuint iState = 0; iState < 3; ++iState){
      for (CFuint iDim = 0; iDim < 2; ++iDim) {
        m_adimNormal[iDim] = cellnormals.getNodalNormComp(iState, iDim);
      }
      CFreal m_nodeArea = cellnormals.getAreaNode(iState);
      m_adimNormal *= 1. / m_nodeArea;

      getMethodData().getDistribVar()->computeEigenValuesVectors(m_rightEv, m_leftEv, m_eValues,m_adimNormal);


      CFreal epsilon = 0.01*(vel2+c);

      for(CFuint iEq=0; iEq<m_nbEqs; ++iEq)
        if (m_eValues[iEq] > epsilon){
          m_eValuesP[iEq] = max(0.,m_eValues[iEq]);
        }
        else {
          m_eValuesP[iEq] = (m_eValues[iEq] + epsilon)*(m_eValues[iEq] + epsilon)/(4.0*epsilon);
        }


      *m_kPlus[iState] = m_rightEv*(m_eValuesP*m_leftEv);
    }

    m_sumKplus = *m_kPlus[0];
    for (CFuint iState = 1; iState < 3; ++iState){
      m_sumKplus += *m_kPlus[iState];
    }


    m_inverter->invert(m_sumKplus, m_tau);

  }



  // Compute the gradient of u at node0 time n+1
  m_nabla_u0N1[XX] = (3.0*nx1*(*m_tStates[0]) -     nx2*(*m_tStates[1]) - nx3*(*m_tStates[2]) +
                      4.0*nx2*(*m_tStates[3]) + 4.0*nx3*(*m_tStates[5]))*gradcoeff;
  m_nabla_u0N1[YY] = (3.0*ny1*(*m_tStates[0]) -     ny2*(*m_tStates[1]) - ny3*(*m_tStates[2]) +
                      4.0*ny2*(*m_tStates[3]) + 4.0*ny3*(*m_tStates[5]))*gradcoeff;


  // Compute the gradient of u at node1 time n+1
  m_nabla_u1N1[XX] =    (-nx1*(*m_tStates[0]) + 3.0*nx2*(*m_tStates[1]) - nx3*(*m_tStates[2]) +
                      4.0*nx1*(*m_tStates[3]) + 4.0*nx3*(*m_tStates[4]))*gradcoeff;
  m_nabla_u1N1[YY] =    (-ny1*(*m_tStates[0]) + 3.0*ny2*(*m_tStates[1]) - ny3*(*m_tStates[2]) +
                      4.0*ny1*(*m_tStates[3]) + 4.0*ny3*(*m_tStates[4]))*gradcoeff;


  // Compute the gradient of u at node2 time n+1
  m_nabla_u2N1[XX] =    (-nx1*(*m_tStates[0])-      nx2*(*m_tStates[1]) + 3.0*nx3*(*m_tStates[2]) +
                      4.0*nx2*(*m_tStates[4]) + 4.0*nx1*(*m_tStates[5]))*gradcoeff;
  m_nabla_u2N1[YY] =    (-ny1*(*m_tStates[0]) -     ny2*(*m_tStates[1]) + 3.0*ny3*(*m_tStates[2]) +
                      4.0*ny2*(*m_tStates[4]) + 4.0*ny1*(*m_tStates[5]))*gradcoeff;

  // Compute the residual in node0, node1 and node2
  m_Res0 += (*states[0]) +
            J * (jacob0_N1[XX] * m_nabla_u0N1[XX] + jacob0_N1[YY] * m_nabla_u0N1[YY]) ;

  m_Res1 += (*states[1]) +
            J * (jacob1_N1[XX] * m_nabla_u1N1[XX] + jacob1_N1[YY] * m_nabla_u1N1[YY]) ;

  m_Res2 += (*states[2]) +
            J * (jacob2_N1[XX] * m_nabla_u2N1[XX] + jacob2_N1[YY] * m_nabla_u2N1[YY]) ;

  //Dissipation for node 0
  m_tempmat = sqrt(volume)*gradcoeff*( 3.0)*(nx1 * jacob0_N1[XX] + ny1 * jacob0_N1[YY])*m_tau;
  residual[0] = m_tempmat*m_Res0*onethirdArea;
  m_tempmat = sqrt(volume)*gradcoeff*(-1.0)*(nx1 * jacob1_N1[XX] + ny1 * jacob1_N1[YY])*m_tau;
  residual[0] += m_tempmat*m_Res1*onethirdArea;
  m_tempmat = sqrt(volume)*gradcoeff*(-1.0)*(nx1 * jacob2_N1[XX] + ny1 * jacob2_N1[YY])*m_tau;
  residual[0] += m_tempmat*m_Res2*onethirdArea;

  //Dissipation for node 1
  m_tempmat = sqrt(volume)*gradcoeff*(-1.0)*(nx2 * jacob0_N1[XX] + ny2 * jacob0_N1[YY])*m_tau;
  residual[1] = m_tempmat*m_Res0*onethirdArea;
  m_tempmat = sqrt(volume)*gradcoeff*( 3.0)*(nx2 * jacob1_N1[XX] + ny2 * jacob1_N1[YY])*m_tau;
  residual[1] += m_tempmat*m_Res1*onethirdArea;
  m_tempmat = sqrt(volume)*gradcoeff*(-1.0)*(nx2 * jacob2_N1[XX] + ny2 * jacob2_N1[YY])*m_tau;
  residual[1] += m_tempmat*m_Res2*onethirdArea;

  //Dissipation for node 2
  m_tempmat = sqrt(volume)*gradcoeff*(-1.0)*(nx3 * jacob0_N1[XX] + ny3 * jacob0_N1[YY])*m_tau;
  residual[2] = m_tempmat*m_Res0*onethirdArea;
  m_tempmat = sqrt(volume)*gradcoeff*(-1.0)*(nx3 * jacob1_N1[XX] + ny3 * jacob1_N1[YY])*m_tau;
  residual[2] += m_tempmat*m_Res1*onethirdArea;
  m_tempmat = sqrt(volume)*gradcoeff*( 3.0)*(nx3 * jacob2_N1[XX] + ny3 * jacob2_N1[YY])*m_tau;
  residual[2] += m_tempmat*m_Res2*onethirdArea;

  //Dissipation for node 3
  m_tempmat = sqrt(volume)*gradcoeff*( 4.0)*(nx2 * jacob0_N1[XX] + ny2 * jacob0_N1[YY])*m_tau;
  residual[3] = m_tempmat*m_Res0*onethirdArea;
  m_tempmat = sqrt(volume)*gradcoeff*( 4.0)*(nx1 * jacob1_N1[XX] + ny1 * jacob1_N1[YY])*m_tau;
  residual[3] += m_tempmat*m_Res1*onethirdArea;

  //Dissipation for node 4
  m_tempmat = sqrt(volume)*gradcoeff*( 4.0)*(nx3 * jacob1_N1[XX] + ny3 * jacob1_N1[YY])*m_tau;
  residual[4] = m_tempmat*m_Res1*onethirdArea;
  m_tempmat = sqrt(volume)*gradcoeff*( 4.0)*(nx2 * jacob2_N1[XX] + ny2 * jacob2_N1[YY])*m_tau;
  residual[4] += m_tempmat*m_Res2*onethirdArea;

//   //Dissipation for node 5
  m_tempmat = sqrt(volume)*gradcoeff*( 4.0)*(nx3 * jacob0_N1[XX] + ny3 * jacob0_N1[YY])*m_tau;
  residual[5] = m_tempmat*m_Res0*onethirdArea;
  m_tempmat = sqrt(volume)*gradcoeff*( 4.0)*(nx1 * jacob2_N1[XX] + ny1 * jacob2_N1[YY])*m_tau;
  residual[5] += m_tempmat*m_Res2*onethirdArea;

  CFuint nbStates = states.size();

   // by default is one, meaning no shock detection
   m_theta = 1.0;
   if (m_with_shock_detect){

    m_min_states =   MathConsts::CFrealMax();
    m_max_states = - MathConsts::CFrealMax();

    for(CFuint iEq = 0; iEq<m_nbEqs; ++iEq){
      for (CFuint iState=0;iState<nbStates; ++iState){
        m_min_states[iEq] = std::min((*states[iState])[iEq],m_min_states[iEq]);
        m_max_states[iEq] = std::max((*states[iState])[iEq],m_max_states[iEq]);
      }

      m_theta[iEq] -= std::abs(m_max_states[iEq] - m_min_states[iEq])/(std::abs(m_max_states[iEq]) + std::abs(m_min_states[iEq])+10.0E-10);
    }
  }

  CFreal m_thetamin = m_theta.min();
  for (CFuint iState=0;iState<nbStates; ++iState)
    {
      residual[iState] *= m_thetamin;

    }
}

//////////////////////////////////////////////////////////////////////////////

    }// End namespace FluctSplit

}// End namespace COOLFluiD
