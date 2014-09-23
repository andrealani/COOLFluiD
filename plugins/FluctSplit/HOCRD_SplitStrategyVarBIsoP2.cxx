#include "Common/BadValueException.hh"

#include "Framework/MethodStrategyProvider.hh"
#include "Framework/ContourIntegrator.hh"
#include "Framework/MeshData.hh"
#include "FluctSplit/FluctSplitHO.hh"
#include "Framework/SubSystemStatus.hh"
#include "Framework/BaseTerm.hh"

#include "FluctSplit/HOCRD_SplitStrategyVarBIsoP2.hh"

#include "LinearAdvSys/LinearAdvSysVarSet.hh"
#include "LinearAdv/LinearAdv2DVarSet.hh"
#include <boost/concept_check.hpp>

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Physics::LinearAdv;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

MethodStrategyProvider<HOCRD_SplitStrategyVarBIsoP2,
                       FluctuationSplitData,
                       FluctuationSplitStrategy,
                       FluctSplitHOModule>
aHOCRDFluctSplitVarBIsoP2StrategyProvider("HOCRDVarBIsoP2");

//////////////////////////////////////////////////////////////////////////////

HOCRD_SplitStrategyVarBIsoP2::HOCRD_SplitStrategyVarBIsoP2(const std::string& name) :
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

HOCRD_SplitStrategyVarBIsoP2::~HOCRD_SplitStrategyVarBIsoP2()
{
  if (isSetup()) unsetup();
}

//////////////////////////////////////////////////////////////////////////////

void HOCRD_SplitStrategyVarBIsoP2::unsetup()
{
  deletePtr(m_subcell_normals);

  // AL: this could be changed ...
  // where are deleted the qdstates ???
  for (CFuint i = 0; i < m_qdExtraVars.size(); ++i) {
    deletePtr(m_qdExtraVars[i]);
  }
}

//////////////////////////////////////////////////////////////////////////////

void HOCRD_SplitStrategyVarBIsoP2::setup()
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


  
  // physical data evaluated in the quadrature points
  m_pdata.resize(nbQdPts);
  for (CFuint  i = 0; i < nbQdPts; ++i) {
    PhysicalModelStack::getActive()->getImplementor()->getConvectiveTerm()->
      resizePhysicalData(m_pdata[i]);
  }  
  
// Iso-parametric triangles

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

  /* Nab: Quadrature with 3 points.
   * Quadratic order.
   * The quadrature points are placed in the midpoint of edges.
   * Numbered as usual */
  qd0[0] = 0.5;  qd1[0] = 0.0;
  qd0[1] = 0.5;  qd1[1] = 0.5;
  qd0[2] = 0.0;  qd1[2] = 0.5;
  
  wqd.resize(nbQdPts); // 3 or 5 or ... quadrature points per face
  wqd[0] = 1.0/6.0;
  wqd[1] = 1.0/6.0;
  wqd[2] = 1.0/6.0;


  sfValues.resize(nbQdPts,3);
  for(CFuint qdPt = 0; qdPt < nbQdPts; ++qdPt)
  {
    sfValues(qdPt,0) = 1. - qd0[qdPt] - qd1[qdPt];
    sfValues(qdPt,1) = qd0[qdPt];
    sfValues(qdPt,2) = qd1[qdPt];
  }
  
  sfDerivativesXi.resize(nbQdPts,3);
  sfDerivativesEta.resize(nbQdPts,3);
  for(CFuint qdPt = 0; qdPt < nbQdPts; ++qdPt)
  {
    sfDerivativesXi(qdPt,0) = -1.0;
    sfDerivativesXi(qdPt,1) =  1.0;
    sfDerivativesXi(qdPt,2) =  0.0;

    sfDerivativesEta(qdPt,0) = -1.0;
    sfDerivativesEta(qdPt,1) =  0.0;
    sfDerivativesEta(qdPt,2) =  1.0;
  }

  sfDerivativesX.resize(nbQdPts,3);
  sfDerivativesY.resize(nbQdPts,3);

  qdstates.resize(nbQdPts); // 3 quadrature points per face

  
  for(CFuint iState=0;iState < nbQdPts; ++iState) {
	qdstates[iState] = new State();
  }

  /* Nab: No se que hace esto */
  const CFuint extra_var_size = getMethodData().getUpdateVar()->getExtraPhysicalVarsSize();
  m_qdExtraVars.resize(qdstates.size());
  for (CFuint i = 0; i < m_qdExtraVars.size(); ++i) {
    m_qdExtraVars[i] = new RealVector(extra_var_size);
  }

  qdnodes.resize(nbQdPts); // 3 quadrature points per face

  // Nab: Match coordinates to every state: Node -> State
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



  // tell the splitters to compute the betas
  getMethodData().getDistributionData().computeBetas = false;

}

//////////////////////////////////////////////////////////////////////////////

void HOCRD_SplitStrategyVarBIsoP2::computeFluctuation(vector<RealVector>& residual)
{
//   DistributionData& distdata = getMethodData().getDistributionData();
// 
//   vector<State*>& states = *distdata.states;
// 
// 
//   cf_assert (states.size() == 6); // P2 triangles for solution space
//   cf_assert (residual.size() >= states.size()); // state residual

  // reset the residual because we will accumulate the sub element contributions
  for (CFuint i = 0; i < residual.size(); ++i) { residual[i] = 0.0; }

  /// Coordinates of the nodes in the physical space
 //vector<Node*>& nodes = *getMethodData().getDistributionData().cell->getNodes();

 DistributionData& distdata = getMethodData().getDistributionData();

 vector<Node*>&  nodes  = *distdata.cell->getNodes();
 vector<State*>& states = *distdata.states;

  cf_assert (nodes.size()  == 3);

  vector<RealVector> lambda;
  RealVector jacQdp;
  vector<RealVector> betaQdp;
  RealVector nablaF;
  
  
  // It gets the advection speed in the CFcase and stores them in lambda
  getVelQuadPoints(lambda, nodes);

  /*
  for(CFuint i=0;i<lambda.size();i++)
    for(CFuint j=0;j<lambda[i].size();j++)
      cout << lambda[i][j] << "\t";
  cout << "\n";
  */
  
  // It computes the jacobian in the quadrature points and stores them in jacQdp
  computeJacQuadPoints(jacQdp, nodes);

  /*
  for(CFuint i=0;i<jacQdp.size();i++)
      cout << jacQdp[i] << "\t";
  cout << "\n";
  */
  
  // It computes the value of the shape functions at the quadrature points
  // THIS IS NOT NEEDED ANYMORE: SHAPE FUNCTION VALUES ARE COMPUTED IN SETUP
  // computesShapeFunctionQuadPoints(shapeFunctionsQdp,coordinates);

  // It computes the value of the gradient of the shape functions at the quadrature points
  computesGradShapeFunctionQuadPoints(nodes, sfDerivativesX, sfDerivativesY);
  
//   for(CFuint i = 0; i<3;++i)
//     for(CFuint j = 0; j<nbQdPts;++j)
//       for(CFuint k = 0; k<coordinates.size();++k)
// 	cout << gradShapeFunctionQdp[i][j][k] << "\t";
//   cout << "\n";
  
  // It computes the beta(i,xq)
  computesBetaQuadPoints(lambda,sfDerivativesX,sfDerivativesY,betaQdp);
  
  /*
  for(CFuint j = 0; j<nbQdPts;++j)
  {
   std::cout << "Quadrature point " << j << std::endl;
   for(CFuint i = 0; i<3;++i)
   {
     std::cout << betaQdp[i][j] << "\t";
   }
   std::cout << std::endl;
  }

  std::cout << std::endl;
  */

  // It computes nabla(lambda*u) at the quadrature points
  computesNablaF(lambda, sfDerivativesX, sfDerivativesY, nablaF);
  
  /*
  // It computes de residuals

const RealVector& jac, const RealVector& weight,
                  const std::vector<RealVector>& beta, const RealVector& nablaF,
                  std::vector<RealVector>& resid
  */
  computeRes(jacQdp, wqd, betaQdp, nablaF, residual);

  /*
  std::cout << " Residuals in cell " << getMethodData().getDistributionData().cellID << ":" << std::endl;
  for(CFuint i = 0; i < residual.size(); ++i)
  {
    std::cout << residual[i] << std::endl;
  }
  std::cout << std::endl  << std::endl;
  */

  // Set the update coefficients
  DataHandle< CFreal> updateCoeff = socket_updateCoeff.getDataHandle();
  
  for(CFuint q = 0; q < nbQdPts; ++q)
  {
    for(CFuint iState = 0; iState < states.size(); ++iState)
    {
      const CFreal eigenvalue =  lambda[XX][q] * sfDerivativesX(q,iState) + lambda[YY][q] * sfDerivativesY(q,iState) ;
      updateCoeff[states[iState]->getLocalID()] += jacQdp[q]*wqd[q]*std::max(eigenvalue,1.e-10);
    }
  }




  /* Nab: No computeK and no splitter
  // m_splitter->computeK(substates, m_subcell_normals);

  // transform fluxes of subelement to distribution variables
  *phi = *getMethodData().getSolutionToDistribMatTrans()->transformFromRef(m_phisubT[3]);

*/
}

//////////////////////////////////////////////////////////////////////////////

void HOCRD_SplitStrategyVarBIsoP2::getVelQuadPoints(vector<RealVector>& vel, const std::vector<Node*>& coord)
{
  /* Dummy code.
   * It should do a real call to the CFcase information. However the DistributionVariable
   * and the DistributionData gradients are empty
   */
  
  // The resize should be in the setup
  vel.resize(DIM_2D);
  for(CFuint i = 0; i<DIM_2D;++i)
  {
    vel[i].resize(nbQdPts);
  }

  // Constant advection velocity
  for(CFuint i = 0; i<nbQdPts; ++i)
  {
    vel[XX][i] = 0.1;
    vel[YY][i] = 1.0;
  }
  
 
//   for(CFuint i = 0; i<nbQdPts; ++i){
//     CFout << "Advection velocity in qd pt = [" << vel[0][i] << "," << vel[1][i] << "]" << std::endl; 
//   }
  
 
  
  
  /* Advection velocity */
  /* Usually this information is taken in to step, first by a calling to 
   * the splitter, then the splitter calls the an istance of ConvectiveVarSet.
   * The cast is due to getModel is a LinearAdvSysVarSet function, it doesnt
   * exist in ConvectiveVarSet.
   * In the future, this part should be coded in different files, in order
   * to be consistent with the CoolFluid way to work. */

//   First way: Copied from LinearAdv2DVarSet
/*
  const RealVector& linearData = (getMethodData().getDistribVar().d_castTo<LinearAdv2DVarSet>())->getModel()->getPhysicalData();
  // i dont know if it is the information in the node or the average over the cell
  cout << linearData[LinearAdvTerm::VX] << "\t";
  cout << linearData[LinearAdvTerm::VY] << "\n";
  for(CFuint i = 0; i<nbQdPts; ++i){
    vel[0][i] = linearData[LinearAdvTerm::VX];
    vel[1][i] = linearData[LinearAdvTerm::VY];
  }
  */
  
//   Not working
//   Second way: Copied from LinearAdvTerm 
//   vector<RealMatrix>* const jacobians = PhysicalModelStack::getActive()->getImplementor()->getJacobians();
//   for(CFuint i = 0; i<nbQdPts; ++i){
//     vel[XX][i] = (*jacobians)[XX][0];
//     vel[YY][i] = (*jacobians)[YY][0];
//     cout << vel[XX][i] << "\t" << vel[YY][i] << "\n";
//   }
 
}

void HOCRD_SplitStrategyVarBIsoP2::computeJacQuadPoints(RealVector& jac, const vector<Node*>& coord){
  // The resize should be in the setup
  jac.resize(nbQdPts);

  
  /* Code for P1 Elements: (the jacobian is constant)
   * This line computes the transformation matrix (x,y)=f(xi,eta)
   * for P1 and then the jacobian, in a row.
   */
  for(CFuint i=0; i< nbQdPts; ++i)
  {
     jac[i] =   ( (*coord[1])[XX] - (*coord[0])[XX]) * ((*coord[2])[YY] - (*coord[0])[YY])
              - ( (*coord[2])[XX] - (*coord[0])[XX]) * ((*coord[1])[YY] - (*coord[0])[YY]);
     //jac[i] = std::abs(jac[i]);

    //std::cout << "Jacobian = " << jac[i] << std::endl;
  }
}


void HOCRD_SplitStrategyVarBIsoP2::computesShapeFunctionQuadPoints(vector<RealVector>& shapeFun, std::vector<RealVector>& coord){
  // P1 shape functions. It should be implemented using CooldFluid shape functions
  // 3 is the number of nodes
  //getMethodData().getDistributionData().states->size();
  shapeFun.resize(3);
  for(CFuint i = 0; i<3;++i)
    shapeFun[i].resize(nbQdPts);
  
  // shapeFun[node][quad point]
  shapeFun[0][0] = 0.5;
  shapeFun[0][1] = 0.0;
  shapeFun[0][2] = 0.5;
  shapeFun[1][0] = 0.5;
  shapeFun[1][1] = 0.5;
  shapeFun[1][2] = 0.0;
  shapeFun[2][0] = 0.0;
  shapeFun[2][1] = 0.5;
  shapeFun[2][2] = 0.5;
}


void HOCRD_SplitStrategyVarBIsoP2::computesGradShapeFunctionQuadPoints(const std::vector<Framework::Node*>& coord,
                                                                       RealMatrix& dPhidX,
                                                                       RealMatrix& dPhidY)
{
  for(CFuint q = 0; q < nbQdPts; ++q)
  {

    const CFreal dXdXi = sfDerivativesXi(q,0) * (*coord[0])[XX] +
                         sfDerivativesXi(q,1) * (*coord[1])[XX] +
                         sfDerivativesXi(q,2) * (*coord[2])[XX];

    const CFreal dXdEta = sfDerivativesEta(q,0) * (*coord[0])[XX] +
                          sfDerivativesEta(q,1) * (*coord[1])[XX] +
                          sfDerivativesEta(q,2) * (*coord[2])[XX];

    const CFreal dYdXi = sfDerivativesXi(q,0) * (*coord[0])[YY] +
                         sfDerivativesXi(q,1) * (*coord[1])[YY] +
                         sfDerivativesXi(q,2) * (*coord[2])[YY];

    const CFreal dYdEta = sfDerivativesEta(q,0) * (*coord[0])[YY] +
                          sfDerivativesEta(q,1) * (*coord[1])[YY] +
                          sfDerivativesEta(q,2) * (*coord[2])[YY];

    const CFreal jacobian = dXdXi * dXdEta - dYdXi * dYdEta;
    const CFreal invJacobian = 1.0/std::abs(jacobian);


    for(CFuint i = 0; i < 3; ++i) // Loop over shape functions
    {
      dPhidX(q,i) = invJacobian * (   dYdEta * sfDerivativesXi(q,i) - dYdXi * sfDerivativesEta(q,i) );
      dPhidY(q,i) = invJacobian * ( - dXdEta * sfDerivativesXi(q,i) + dXdXi * sfDerivativesEta(q,i) );
    }
  }

}

void HOCRD_SplitStrategyVarBIsoP2::computesBetaQuadPoints(const vector<RealVector>& vel,
                                                          const RealMatrix& dPhidX, const RealMatrix& dPhidY,
                                                          vector<RealVector>& beta)
{
  // The resize should be in the setup
  beta.resize(3);
  for(CFuint i = 0; i<3;++i)
    beta[i].resize(nbQdPts);
  
  vector<RealVector> kplus;  /// Matrix storing K+. [node, quadtrature point]
  RealVector kplusSum;	/// Vector storing the sum of K+. [quadrature point]
  
  kplus.resize(3);
  for(CFuint i = 0; i<3;++i)
    kplus[i].resize(nbQdPts);
  
  kplusSum.resize(nbQdPts);
  for(CFuint i = 0; i<nbQdPts;++i)
    kplusSum[i]=0.0;
  
  // It computes k+ and sum(k+)
  for(CFuint i = 0; i<3;++i){
    for(CFuint j = 0; j<nbQdPts;++j){
      kplus[i][j] = vel[XX][j]*dPhidX(j,i) + vel[YY][j]*dPhidY(j,i);
      kplus[i][j] = kplus[i][j]>0.0?kplus[i][j]:0.0; // K+ must be positive
      kplusSum[j] += kplus[i][j];
    }
  }
   
  for(CFuint i = 0; i<3;++i){
    for(CFuint j = 0; j<nbQdPts;++j){
      beta[i][j] = kplus[i][j]/kplusSum[j];
//       cout << beta[i][j] << "\t";
    }
  }
//   cout << "\n";
  
}


void HOCRD_SplitStrategyVarBIsoP2::computesNablaF(const std::vector<RealVector>& vel,
                                                  const RealMatrix& dPhidX, const RealMatrix& dPhidY, RealVector& nablaF)
{
  // Scalar case
  // states is a matrix with the variable of the problem [nodes,variables,cases].
  // For the first testcases there are 4 cases. They are different initial conditions
  vector<State*>& states = *getMethodData().getDistributionData().states;
//   cout << states.size() << "\t" << states[0]->size() << "\t" << states[0][0].size() << "\n";
// 	result = 3 1 1
//   cout << (*states[0]) << "\t" << (*states[1]) << "\t" << (*states[2]) << "\n";
  nablaF.resize(nbQdPts);
  for(CFuint i = 0; i<nbQdPts;++i)
  {
    nablaF[i]=0;
    for(CFuint iState = 0; iState < states.size(); ++iState)
    {
      nablaF[i] += ( states[iState][0][0])*(vel[XX][i]*dPhidX(i,iState)+vel[YY][i]*dPhidY(i,iState) );
    }
  }
}


void HOCRD_SplitStrategyVarBIsoP2::computeRes(const RealVector& jac, const RealVector& weight,
                                              const std::vector<RealVector>& beta, const RealVector& nablaF,
                                              std::vector<RealVector>& resid)
{
  // Scalar case
  CFuint nbstates = getMethodData().getDistributionData().states->size();
//   cout << resid.size() << "\t" << resid[0].size() << "\n";
  
  for(CFuint iState = 0; iState<nbstates;++iState) // node loop
  {
    for(CFuint i = 0; i < nbQdPts; ++i) // quadrature points loop
    {
      resid[iState][0] += jac[i]*nablaF[i]*beta[iState][i]*weight[i];
//       cout << "cell: " << getMethodData().getDistributionData().cellID << "\t nabla \t" << nablaF[i] << "\t beta \t" << beta[iState][i] << "\t weight \t" << weight[i] << "\n";
    }
//    cout << resid[iState][0] << "\n";
  } 
//  cout << "\t" << getMethodData().getDistributionData().cellID << "\n";
}
    
    
//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<Framework::BaseDataSocketSink> >
HOCRD_SplitStrategyVarBIsoP2::needsSockets()
{
  std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > result
    = FluctuationSplitStrategy::needsSockets();
  result.push_back(&socket_updateCoeff);
  return result;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD
