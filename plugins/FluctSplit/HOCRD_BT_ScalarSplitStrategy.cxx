#include "Common/BadValueException.hh"
#include "Framework/MethodStrategyProvider.hh"
#include "Framework/ContourIntegrator.hh"
#include "Framework/MeshData.hh"
#include "Framework/BaseTerm.hh"

#include "FluctSplit/FluctSplitHO.hh"
#include "FluctSplit/HOCRD_BT_ScalarSplitStrategy.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::MathTools;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

MethodStrategyProvider<HOCRD_BT_ScalarSplitStrategy,
                       FluctuationSplitData,
                       FluctuationSplitStrategy,
                       FluctSplitHOModule>
hocrdbtscalarFluctSplitStrategyProvider("HOCRD_ScalarBT");

//////////////////////////////////////////////////////////////////////////////

HOCRD_BT_ScalarSplitStrategy::HOCRD_BT_ScalarSplitStrategy(const std::string& name) :
  FluctuationSplitStrategy(name),
  socket_updateCoeff("updateCoeff"),
  m_solutionVar(CFNULL),
  m_updateVar(CFNULL),
  m_unitFaceNormals(),
  m_phi(0),
  m_sumKplusU(),
  m_sumKplus(),
  m_uTemp(),
  m_uMin(),
  m_temp(),
  m_k1Plus(0),
  m_k1Min(0),
  m_k2Plus(0),
  m_k2Min(0),
  m_k3Plus(0),
  m_k3Min(0),
  m_k4Plus(0),
  m_k4Min(0),
  m_k(0),
  m_adimNormal(),
  m_theta1(0),
  m_theta2(0),
  m_theta3(0),
  m_theta4(0)
{
}

//////////////////////////////////////////////////////////////////////////////

HOCRD_BT_ScalarSplitStrategy::~HOCRD_BT_ScalarSplitStrategy()
{
}

//////////////////////////////////////////////////////////////////////////////

void HOCRD_BT_ScalarSplitStrategy::setup()
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
  m_k1Min.resize(m_maxNbStatesInCell);
  m_k2Plus.resize(m_maxNbStatesInCell);
  m_k2Min.resize(m_maxNbStatesInCell);
  m_k3Plus.resize(m_maxNbStatesInCell);
  m_k3Min.resize(m_maxNbStatesInCell);
  m_k4Plus.resize(m_maxNbStatesInCell);
  m_k4Min.resize(m_maxNbStatesInCell);
  m_k.resize(m_maxNbStatesInCell);

  for (CFuint i = 0; i < m_maxNbStatesInCell; ++i) {
    m_k1Plus[i].resize(nbEqs);
    m_k1Min[i].resize(nbEqs);
    m_k2Plus[i].resize(nbEqs);
    m_k2Min[i].resize(nbEqs);
    m_k3Plus[i].resize(nbEqs);
    m_k3Min[i].resize(nbEqs);
    m_k4Plus[i].resize(nbEqs);
    m_k4Min[i].resize(nbEqs);
    m_k[i].resize(nbEqs);
  }
  m_adimNormal.resize(2);

  m_theta1.resize(nbEqs);
  m_theta2.resize(nbEqs);
  m_theta3.resize(nbEqs);
  m_theta4.resize(nbEqs);

  m_theta.resize(nbEqs);

  m_sumKplusU.resize(nbEqs);
  m_sumKplus.resize(nbEqs);
  m_uTemp.resize(nbEqs);
  m_uMin.resize(nbEqs);
  m_temp.resize(nbEqs);
}

//////////////////////////////////////////////////////////////////////////////

void HOCRD_BT_ScalarSplitStrategy::computeFluctuation(vector<RealVector>& residual)
{

//   std::cout << " +++ begin compute fluctuation +++\n" << std::endl;

   const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();

 vector<State*>& states = *getMethodData().getDistributionData().states;
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
  computeK(substates,m_k1Plus, m_k1Min);

   // transform fluxes of subelement to distribution variables
  SafePtr<RealVector> phi = &getMethodData().getDistributionData().phi;
  *phi = *getMethodData().getSolutionToDistribMatTrans()->transformFromRef(m_phisubT[0]);

  distributeN(m_k1Plus,m_phiN1);


  for (CFuint iEq = 0; iEq < nbEqs; ++iEq){
    m_phitot = std::abs(m_phiN1[0][iEq]) + std::abs(m_phiN1[1][iEq]) + std::abs(m_phiN1[2][iEq]) ;
    m_theta1[iEq] = std::abs((*phi)[iEq])/max(MathTools::MathConsts::CFrealEps(), m_phitot);
  }


   /*****         Triangle 3-1-4          *****/

  substates[0] = states[3];
  substates[1] = states[1];
  substates[2] = states[4];

  getMethodData().getDistributionData().tStates =
    computeConsistentStates(&substates);

  // compute the residual and the upwind parameters k in this cell
  computeK(substates,m_k2Plus, m_k2Min);

   // transform fluxes of subelement to distribution variables

  *phi = *getMethodData().getSolutionToDistribMatTrans()->transformFromRef(m_phisubT[1]);

  distributeN(m_k2Plus,m_phiN2);

  for (CFuint iEq = 0; iEq < nbEqs; ++iEq){
    m_phitot = std::abs(m_phiN2[0][iEq]) + std::abs(m_phiN2[1][iEq]) + std::abs(m_phiN2[2][iEq]) ;
    m_theta2[iEq] = std::abs((*phi)[iEq])/max(MathTools::MathConsts::CFrealEps(), m_phitot);
  }

  /*****         Triangle 5-4-2          *****/

  substates[0] = states[5];
  substates[1] = states[4];
  substates[2] = states[2];

  getMethodData().getDistributionData().tStates =
    computeConsistentStates(&substates);

  // compute the residual and the upwind parameters k in this cell
  computeK(substates,m_k3Plus, m_k3Min);

   // transform fluxes of subelement to distribution variables

  *phi = *getMethodData().getSolutionToDistribMatTrans()->transformFromRef(m_phisubT[2]);

  distributeN(m_k3Plus,m_phiN3);

  for (CFuint iEq = 0; iEq < nbEqs; ++iEq){
    m_phitot = std::abs(m_phiN3[0][iEq]) + std::abs(m_phiN3[1][iEq]) + std::abs(m_phiN3[2][iEq]) ;
    m_theta3[iEq] = std::abs((*phi)[iEq])/max(MathTools::MathConsts::CFrealEps(), m_phitot);
  }

 /*****         Triangle 4-5-3          *****/

  substates[0] = states[4];
  substates[1] = states[5];
  substates[2] = states[3];

  getMethodData().getDistributionData().tStates = computeConsistentStates(&substates);

  // compute the residual and the upwind parameters k in this cell
  // the oriantation of the normal in this sub-element is oposite to the one of the element
  cellnormals.scale(-0.5);

 // compute the residual and the upwind parameters k in this cell
  computeK(substates,m_k4Plus, m_k4Min);
   // transform fluxes of subelement to distribution variables

  *phi = *getMethodData().getSolutionToDistribMatTrans()->transformFromRef(m_phisubT[3]);

  distributeN(m_k4Plus,m_phiN4);

  for (CFuint iEq = 0; iEq < nbEqs; ++iEq){

    m_phitot = std::abs(m_phiN4[0][iEq]) + std::abs(m_phiN4[1][iEq]) + std::abs(m_phiN4[2][iEq]) ;

    m_theta4[iEq] = std::abs((*phi)[iEq])/max(MathTools::MathConsts::CFrealEps(), m_phitot);
  }

// Then we take the max of the thetas

   for (CFuint iEq = 0; iEq < nbEqs; ++iEq){
    m_theta[iEq] = max(m_theta1[iEq],m_theta2[iEq]);
    m_theta[iEq] = max(m_theta[iEq],m_theta3[iEq]);
    m_theta[iEq] = max(m_theta[iEq],m_theta4[iEq]);
    }

//Then we use this theta to build the blending


  /*****         Triangle 0-3-5          *****/

  substates[0] = states[0];
  substates[1] = states[3];
  substates[2] = states[5];

  getMethodData().getDistributionData().tStates =
    computeConsistentStates(&substates);

  *phi = *getMethodData().getSolutionToDistribMatTrans()->transformFromRef(m_phisubT[0]);
  distributeLDA(m_k1Plus, m_phi);

  residual[0] += m_theta*m_phiN1[0] + (1.0 - m_theta1)*m_phi[0];
  residual[3] += m_theta*m_phiN1[1] + (1.0 - m_theta1)*m_phi[1];
  residual[5] += m_theta*m_phiN1[2] + (1.0 - m_theta1)*m_phi[2];

  /*****         Triangle 3-1-4          *****/

  substates[0] = states[3];
  substates[1] = states[1];
  substates[2] = states[4];

  getMethodData().getDistributionData().tStates =
    computeConsistentStates(&substates);

  *phi = *getMethodData().getSolutionToDistribMatTrans()->transformFromRef(m_phisubT[1]);

  distributeLDA(m_k2Plus,m_phi);

  residual[3] += m_theta*m_phiN2[0] + (1.0 - m_theta2)*m_phi[0];
  residual[1] += m_theta*m_phiN2[1] + (1.0 - m_theta2)*m_phi[1];
  residual[4] += m_theta*m_phiN2[2] + (1.0 - m_theta2)*m_phi[2];

  /*****         Triangle 5-4-2          *****/

  substates[0] = states[5];
  substates[1] = states[4];
  substates[2] = states[2];

  getMethodData().getDistributionData().tStates =
    computeConsistentStates(&substates);

  *phi = *getMethodData().getSolutionToDistribMatTrans()->transformFromRef(m_phisubT[2]);

  distributeLDA(m_k3Plus,m_phi);

  residual[5] += m_theta*m_phiN3[0] + (1.0 - m_theta3)*m_phi[0];
  residual[4] += m_theta*m_phiN3[1] + (1.0 - m_theta3)*m_phi[1];
  residual[2] += m_theta*m_phiN3[2] + (1.0 - m_theta3)*m_phi[2];

   /*****         Triangle 4-5-3          *****/

  substates[0] = states[4];
  substates[1] = states[5];
  substates[2] = states[3];

  getMethodData().getDistributionData().tStates =
    computeConsistentStates(&substates);

  *phi = *getMethodData().getSolutionToDistribMatTrans()->transformFromRef(m_phisubT[3]);

  distributeLDA(m_k4Plus,m_phi);

  residual[4] += m_theta*m_phiN4[0] + (1.0 - m_theta4)*m_phi[0];
  residual[5] += m_theta*m_phiN4[1] + (1.0 - m_theta4)*m_phi[1];
  residual[3] += m_theta*m_phiN4[2] + (1.0 - m_theta4)*m_phi[2];

}

//////////////////////////////////////////////////////////////////////////////

void HOCRD_BT_ScalarSplitStrategy::computeHOFluctuation()
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
void HOCRD_BT_ScalarSplitStrategy::computeK(const std::vector<Framework::State*>& states,
                std::vector<RealVector>& m_kPlus,
                std::vector<RealVector>& m_kMin){
  InwardNormalsData& cellnormals = (*socket_normals.getDataHandle()
				    [getMethodData().getDistributionData().cellID]);

  DataHandle< CFreal> updateCoeff = socket_updateCoeff.getDataHandle();
  const CFuint m_nbStatesInCell = states.size();
  const CFreal m_invDim = 1.0/2.0;
  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  // The transformation of the normal is needed if the coordinate system is rotated
  for (CFuint iState = 0; iState < m_nbStatesInCell; ++iState) {
    // The transformation of the normal is needed if the coordinate system
    // is rotated. The normal is adimensionalized, so there is need to multiply
    // by the nodeArea when computing the k parameter

    for (CFuint iDim = 0; iDim < 2; ++iDim) {
        m_adimNormal[iDim] = cellnormals.getNodalNormComp(iState, iDim);
    }
    m_adimNormal *= 1. / cellnormals.getAreaNode(iState);

    // uses the linearized state to compute the K vector
    getMethodData().getDistribVar()->computeScalarJacobian(m_adimNormal, m_k[iState]);

    CFreal nodeArea = cellnormals.getAreaNode(iState);
    nodeArea   *= m_invDim;
    m_k[iState] *= nodeArea;

    // calculate Kplus and Kminus
    for (CFuint iEq = 0; iEq < nbEqs; ++iEq)
    {
     m_kPlus[iState][iEq] = max(0.,m_k[iState][iEq]);
     m_kMin[iState][iEq]  = min(0.,m_k[iState][iEq]);

  }
  if ( !getMethodData().getDistributionData().isPerturb) {
      updateCoeff[states[iState]->getLocalID()] += m_kPlus[iState].max();
    }
  }


}
//////////////////////////////////////////////////////////////////////////////
void HOCRD_BT_ScalarSplitStrategy::distributeN(std::vector<RealVector> & m_kPlus,vector<RealVector> & phiN){

  const vector<State*>& tStates = *getMethodData().getDistributionData().tStates;
  const RealVector& phiT = getMethodData().getDistributionData().phi;

  m_sumKplusU = m_kPlus[0] * (*tStates[0]);

  m_sumKplus  = m_kPlus[0];

  for (CFuint iState = 1; iState < 3; ++iState) {
    m_sumKplusU += m_kPlus[iState] * (*tStates[iState]);
    m_sumKplus  += m_kPlus[iState];
  }
  m_temp = ( m_sumKplusU - phiT );
  m_uTemp = m_temp / m_sumKplus;

  for (CFuint iState = 0; iState < 3; ++iState) {
    m_uMin = *tStates[iState] - m_uTemp;
    phiN[iState] = m_kPlus[iState]*m_uMin;

  }
}

//////////////////////////////////////////////////////////////////////////////

void HOCRD_BT_ScalarSplitStrategy::distributeLDA ( std::vector<RealVector> & m_kPlus, vector<RealVector> & phiLDA )
{
//   std::cout << " +++ begin distriuteLDA +++\n" << std::endl;
  // unused //  const vector<State*>& tStates = *getMethodData().getDistributionData().tStates;
  const RealVector& phiT = getMethodData().getDistributionData().phi;

  m_sumKplus = m_kPlus[0];
  for (CFuint iState = 1; iState < 3; ++iState) {
    m_sumKplus += m_kPlus[iState];
  }
  m_uTemp = phiT / m_sumKplus;

  for (CFuint iState = 0; iState < 3; ++iState) {
//     CFlog <<  iState  << " state \n"  << CFendl;
//     CFlog <<  phiLDA[iState]  << "\n"  << CFendl;
//     CFlog <<  m_kPlus[iState] << "\n" << CFendl;
    phiLDA[iState] = m_kPlus[iState]*m_uTemp;
  }
//   CF_DEBUG_EXIT;
}


//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<Framework::BaseDataSocketSink> >
HOCRD_BT_ScalarSplitStrategy::needsSockets()
{
   std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > result
     = FluctuationSplitStrategy::needsSockets();
   result.push_back(&socket_updateCoeff);
   return result;
}

//////////////////////////////////////////////////////////////////////////////

void HOCRD_BT_ScalarSplitStrategy::unsetup()
{
  // AL: this could be changed ...
  // where are deleted the qdstates ???
  for (CFuint i = 0; i < m_qdExtraVars.size(); ++i)
    deletePtr(m_qdExtraVars[i]);
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit

} // namespace COOLFluiD
