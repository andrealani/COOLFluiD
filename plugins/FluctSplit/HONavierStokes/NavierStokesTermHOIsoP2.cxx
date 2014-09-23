#include "Environment/ObjectProvider.hh"

#include "Framework/GeometricEntity.hh"
#include "Framework/MeshData.hh"
#include "Framework/MethodStrategyProvider.hh"

#include "NavierStokes/NavierStokesVarSet.hh"
#include "NavierStokes/EulerVarSet.hh"

#include "FluctSplit/HONavierStokes/FluctSplitHONavierStokes.hh"
#include "FluctSplit/FluctuationSplitStrategy.hh"
#include "FluctSplit/HONavierStokes/NavierStokesTermHOIsoP2.hh"
#include "FluctSplit/InwardNormalsData.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Physics::NavierStokes;
using namespace COOLFluiD::MathTools;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

MethodStrategyProvider<NavierStokesTermHOIsoP2,
		       FluctuationSplitData,
		       ComputeDiffusiveTerm,
		       FluctSplitHONavierStokesModule>
navierStokesDiffusiveTermHOIsoP2Provider("NavierStokesHOIsoP2");

//////////////////////////////////////////////////////////////////////////////

  //Setup of coordinates of mother element:

  CFreal NavierStokesTermHOIsoP2::xi_ref[6] = { 0.0, 1.0, 0.0, 0.5, 0.5, 0.0 };
  CFreal NavierStokesTermHOIsoP2::eta_ref[6] = { 0.0, 0.0, 1.0, 0.0, 0.5, 0.5 };

//////////////////////////////////////////////////////////////////////////////

NavierStokesTermHOIsoP2::NavierStokesTermHOIsoP2(const std::string& name) :
  ComputeDiffusiveTerm(name),
  m_diffVar(CFNULL),
  m_updateVar(CFNULL),
  m_states(),
  m_values(),
  m_gradients(),
  m_avValues(),
  m_normal()
{
}

//////////////////////////////////////////////////////////////////////////////

NavierStokesTermHOIsoP2::~NavierStokesTermHOIsoP2()
{
  for (CFuint i = 0; i< m_gradients.size(); ++i) {
    deletePtr(m_gradients[i]);
  }
}

//////////////////////////////////////////////////////////////////////////////

void NavierStokesTermHOIsoP2::setDiffusiveVarSet(SafePtr<DiffusiveVarSet> diffVar)
{
  m_diffVar = diffVar.d_castTo<NavierStokesVarSet>();
}

//////////////////////////////////////////////////////////////////////////////

void NavierStokesTermHOIsoP2::setUpdateVarSet(SafePtr<ConvectiveVarSet> updateVar)
{
  m_updateVar = updateVar.d_castTo<EulerVarSet>();
}


//////////////////////////////////////////////////////////////////////////////

void NavierStokesTermHOIsoP2::computeDiffusiveTerm
(GeometricEntity *const geo, vector<RealVector>& result, bool updateCoeffFlag)
{
  DataHandle< InwardNormalsData*> normals = socket_normals.getDataHandle();
  DataHandle< CFreal> updateCoeff = socket_updateCoeff.getDataHandle();

  // Vector of the nodes of the cell
  m_cellNodes = geo->getNodes();

  // Vector of the states of the cell
  m_cellStates = geo->getStates();
  const CFuint nbCellStates = m_cellStates->size();
  m_cellID = getMethodData().getDistributionData().cellID;

  // store the pointers to state in another array (of RealVector*)
  for (CFuint i = 0; i < nbCellStates; ++i) {
    m_states[i] = (*m_cellStates)[i];
  }

  // compute vars that will be used to compute the gradients
  m_diffVar->setGradientVars(m_states, m_values, geo->nbStates());

  // compute the radius (axysimmetric computations)
//   CFreal radius = 0.0;

  //   if (getMethodData().isAxisymmetric()) {
  //     for (CFuint i = 0; i < nbCellStates; ++i) {
  //       const Node& node = *geo->getNode(i);
  //       radius += node[YY];
  //     }
  //     radius /= nbCellStates;
  //   }

  m_cellVolume = geo->computeVolume();
// unused //  const CFuint dim = PhysicalModelStack::getActive()->getDim();
// unused //  const CFreal dimCoeff = 1./dim;


  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();

  for (CFuint iState = 0; iState < nbCellStates; ++iState)
    for (CFuint iEq = 0; iEq < nbEqs; ++iEq)
        result[iState][iEq] = 0.0;

  NSTerm& model = m_diffVar->getModel();
const CFreal mu = model.getPhysicalData()[NSTerm::MU];
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
     (*m_phi_diff_bub_split[iNode]) = betasInTriag[iNode]*(*m_phi_diff_bub);
  }

  // Diffusive resildual is distributed to each node of the element

  for (CFuint iEq = 0; iEq < nbEqs; ++iEq) {
      result[0][iEq] += (*m_phi_diff_gal[0])[iEq] + 3.0*(*m_phi_diff_bub_split[0])[iEq] - kappa(0,0)*(*m_phi_diff_bub)[iEq];
      result[1][iEq] += (*m_phi_diff_gal[1])[iEq] - kappa(0,1)*(*m_phi_diff_bub)[iEq];
      result[2][iEq] += (*m_phi_diff_gal[2])[iEq] - kappa(0,2)*(*m_phi_diff_bub)[iEq];
      result[3][iEq] += (*m_phi_diff_gal[3])[iEq] + 3.0*(*m_phi_diff_bub_split[1])[iEq] - kappa(0,3)*(*m_phi_diff_bub)[iEq];
      result[4][iEq] += (*m_phi_diff_gal[4])[iEq] - kappa(0,4)*(*m_phi_diff_bub)[iEq];
      result[5][iEq] += (*m_phi_diff_gal[5])[iEq] + 3.0*(*m_phi_diff_bub_split[2])[iEq] - kappa(0,5)*(*m_phi_diff_bub)[iEq];
  }


  // Triangle 2 : nodes 3-1-4

  i1 = 3;
  i2 = 1;
  i3 = 4;

  fluctuation_diff_galerkin(i1, i2, i3);
  fluctuation_diff_bubble(i1, i2, i3);

  betasInTriag = getMethodData().getDistributionData().betaMats[1];

  for ( CFuint iNode = 0; iNode < 3; ++iNode) {
     (*m_phi_diff_bub_split[iNode]) = betasInTriag[iNode]*(*m_phi_diff_bub);
  }

  for (CFuint iEq = 0; iEq < nbEqs; ++iEq) {
     result[0][iEq] += (*m_phi_diff_gal[0])[iEq] - kappa(1,0)*(*m_phi_diff_bub)[iEq];
     result[1][iEq] += (*m_phi_diff_gal[1])[iEq] + 3.0*(*m_phi_diff_bub_split[1])[iEq] - kappa(1,1)*(*m_phi_diff_bub)[iEq];
     result[2][iEq] += (*m_phi_diff_gal[2])[iEq] - kappa(1,2)*(*m_phi_diff_bub)[iEq];
     result[3][iEq] += (*m_phi_diff_gal[3])[iEq] + 3.0*(*m_phi_diff_bub_split[0])[iEq] - kappa(1,3)*(*m_phi_diff_bub)[iEq];
     result[4][iEq] += (*m_phi_diff_gal[4])[iEq] + 3.0*(*m_phi_diff_bub_split[2])[iEq] - kappa(1,4)*(*m_phi_diff_bub)[iEq];
     result[5][iEq] += (*m_phi_diff_gal[5])[iEq] - kappa(1,5)*(*m_phi_diff_bub)[iEq];
  }


  // Triangle 3 : nodes 5-4-2
  i1 = 5;
  i2 = 4;
  i3 = 2;

  fluctuation_diff_galerkin(i1, i2, i3);
  fluctuation_diff_bubble(i1, i2, i3);


  betasInTriag = getMethodData().getDistributionData().betaMats[2];

  for ( CFuint iNode = 0; iNode < 3; ++iNode) {
     (*m_phi_diff_bub_split[iNode]) = betasInTriag[iNode]*(*m_phi_diff_bub);
  }


  for (CFuint iEq = 0; iEq < nbEqs; ++iEq) {
      result[0][iEq] += (*m_phi_diff_gal[0])[iEq] - kappa(2,0)*(*m_phi_diff_bub)[iEq];
      result[1][iEq] += (*m_phi_diff_gal[1])[iEq] - kappa(2,1)*(*m_phi_diff_bub)[iEq];
      result[2][iEq] += (*m_phi_diff_gal[2])[iEq] + 3.0*(*m_phi_diff_bub_split[2])[iEq] - kappa(2,2)*(*m_phi_diff_bub)[iEq];
      result[3][iEq] += (*m_phi_diff_gal[3])[iEq] - kappa(2,3)*(*m_phi_diff_bub)[iEq];
      result[4][iEq] += (*m_phi_diff_gal[4])[iEq] + 3.0*(*m_phi_diff_bub_split[1])[iEq] - kappa(2,4)*(*m_phi_diff_bub)[iEq];
      result[5][iEq] += (*m_phi_diff_gal[5])[iEq] + 3.0*(*m_phi_diff_bub_split[0])[iEq] - kappa(2,5)*(*m_phi_diff_bub)[iEq];

  }

  // Triangle 4 : nodes 4-5-3
  i1 = 4;
  i2 = 5;
  i3 = 3;

  fluctuation_diff_galerkin(i1, i2, i3);
  fluctuation_diff_bubble(i1, i2, i3);

  betasInTriag = getMethodData().getDistributionData().betaMats[3];

  for ( CFuint iNode = 0; iNode < 3; ++iNode) {
     (*m_phi_diff_bub_split[iNode]) = betasInTriag[iNode]*(*m_phi_diff_bub);
  }

  for (CFuint iEq = 0; iEq < nbEqs; ++iEq) {
     result[0][iEq] += (*m_phi_diff_gal[0])[iEq] - kappa(3,0)*(*m_phi_diff_bub)[iEq];
     result[1][iEq] += (*m_phi_diff_gal[1])[iEq] - kappa(3,1)*(*m_phi_diff_bub)[iEq];
     result[2][iEq] += (*m_phi_diff_gal[2])[iEq] - kappa(3,2)*(*m_phi_diff_bub)[iEq];
     result[3][iEq] += (*m_phi_diff_gal[3])[iEq] + 3.0*(*m_phi_diff_bub_split[2])[iEq] - kappa(3,3)*(*m_phi_diff_bub)[iEq];
     result[4][iEq] += (*m_phi_diff_gal[4])[iEq] + 3.0*(*m_phi_diff_bub_split[0])[iEq] - kappa(3,4)*(*m_phi_diff_bub)[iEq];
     result[5][iEq] += (*m_phi_diff_gal[5])[iEq] + 3.0*(*m_phi_diff_bub_split[1])[iEq] - kappa(3,5)*(*m_phi_diff_bub)[iEq];

  }

    //update coefficient
    CFreal avRho = ((*m_states[0])[0] + (*m_states[1])[0] + (*m_states[2])[0])/3.0;   
    const CFreal diffCoeff = mu / avRho;
    CFreal coeff = diffCoeff/(m_cellVolume*4.0);
    if (updateCoeffFlag) {
//       for (CFuint i = 0 ; i < 3; ++i) { 
//         const CFreal faceArea = (normals[m_cellID]->getAreaNode(i));
// 
//         updateCoeff[geo->getState(i)->getLocalID()] += coeff*faceArea*faceArea;
          const CFreal area0 = (normals[m_cellID]->getAreaFace(3)+normals[m_cellID]->getAreaFace(6));
          updateCoeff[geo->getState(0)->getLocalID()] += coeff*area0*area0;

          const CFreal area1 = (normals[m_cellID]->getAreaFace(1)+normals[m_cellID]->getAreaFace(7));
          updateCoeff[geo->getState(1)->getLocalID()] += coeff*area1*area1;

          const CFreal area2 = (normals[m_cellID]->getAreaFace(2)+normals[m_cellID]->getAreaFace(5));
          updateCoeff[geo->getState(2)->getLocalID()] += coeff*area2*area2;
//          } //for

      coeff = 2.0*diffCoeff/(3.0*m_cellVolume);
//       CFreal dot_prod = (normals[m_cellID]->getNodalNormComp(0,XX))*(normals[m_cellID]->getNodalNormComp(1,XX)) +
//                         (normals[m_cellID]->getNodalNormComp(0,YY))*(normals[m_cellID]->getNodalNormComp(1,YY));
// 
//       CFreal faceArea1 = (normals[m_cellID]->getAreaNode(0));
//       CFreal faceArea2 = (normals[m_cellID]->getAreaNode(1));
//       updateCoeff[geo->getState(3)->getLocalID()] += coeff*(faceArea1*faceArea1 + faceArea2*faceArea2 + dot_prod);

      CFreal faceArea0 = (normals[m_cellID]->getAreaFace(3)+normals[m_cellID]->getAreaFace(6));
      CFreal faceArea1 = (normals[m_cellID]->getAreaFace(1)+normals[m_cellID]->getAreaFace(7));

      m_faceNormal0[XX] = 0.5*((normals[m_cellID]->getFaceNormComp(3,XX))+(normals[m_cellID]->getFaceNormComp(6,XX)));
      m_faceNormal0[YY] = 0.5*((normals[m_cellID]->getFaceNormComp(3,YY))+(normals[m_cellID]->getFaceNormComp(6,YY)));
      m_faceNormal1[XX] = 0.5*((normals[m_cellID]->getFaceNormComp(1,XX))+(normals[m_cellID]->getFaceNormComp(7,XX)));
      m_faceNormal1[YY] = 0.5*((normals[m_cellID]->getFaceNormComp(1,YY))+(normals[m_cellID]->getFaceNormComp(7,YY)));

      CFreal dot_prod = m_faceNormal0[XX]*m_faceNormal1[XX]+m_faceNormal0[YY]*m_faceNormal1[YY];

      updateCoeff[geo->getState(3)->getLocalID()] += coeff*(faceArea0*faceArea0 + faceArea1*faceArea1 + dot_prod);

// ---------------------

      faceArea0 = (normals[m_cellID]->getAreaFace(1)+normals[m_cellID]->getAreaFace(7));
      faceArea1 = (normals[m_cellID]->getAreaFace(2)+normals[m_cellID]->getAreaFace(5));

      m_faceNormal0[XX] = 0.5*((normals[m_cellID]->getFaceNormComp(1,XX))+(normals[m_cellID]->getFaceNormComp(7,XX)));
      m_faceNormal0[YY] = 0.5*((normals[m_cellID]->getFaceNormComp(1,YY))+(normals[m_cellID]->getFaceNormComp(7,YY)));
      m_faceNormal1[XX] = 0.5*((normals[m_cellID]->getFaceNormComp(2,XX))+(normals[m_cellID]->getFaceNormComp(5,XX)));
      m_faceNormal1[YY] = 0.5*((normals[m_cellID]->getFaceNormComp(2,YY))+(normals[m_cellID]->getFaceNormComp(5,YY)));

      dot_prod = m_faceNormal0[XX]*m_faceNormal1[XX]+m_faceNormal0[YY]*m_faceNormal1[YY];

      updateCoeff[geo->getState(4)->getLocalID()] += coeff*(faceArea0*faceArea0 + faceArea1*faceArea1 + dot_prod);

// ---------------------

      faceArea0 = (normals[m_cellID]->getAreaFace(2)+normals[m_cellID]->getAreaFace(5));
      faceArea1 = (normals[m_cellID]->getAreaFace(3)+normals[m_cellID]->getAreaFace(6));

      m_faceNormal0[XX] = 0.5*((normals[m_cellID]->getFaceNormComp(2,XX))+(normals[m_cellID]->getFaceNormComp(5,XX)));
      m_faceNormal0[YY] = 0.5*((normals[m_cellID]->getFaceNormComp(2,YY))+(normals[m_cellID]->getFaceNormComp(5,YY)));
      m_faceNormal1[XX] = 0.5*((normals[m_cellID]->getFaceNormComp(3,XX))+(normals[m_cellID]->getFaceNormComp(6,XX)));
      m_faceNormal1[YY] = 0.5*((normals[m_cellID]->getFaceNormComp(3,YY))+(normals[m_cellID]->getFaceNormComp(6,YY)));

      dot_prod = m_faceNormal0[XX]*m_faceNormal1[XX]+m_faceNormal0[YY]*m_faceNormal1[YY];

      updateCoeff[geo->getState(5)->getLocalID()] += coeff*(faceArea0*faceArea0 + faceArea1*faceArea1 + dot_prod);

}

  for (CFuint iEq = 0; iEq < nbEqs; ++iEq) {
     result[0][iEq] = -result[0][iEq];
     result[1][iEq] = -result[1][iEq];
     result[2][iEq] = -result[2][iEq];
     result[3][iEq] = -result[3][iEq];
     result[4][iEq] = -result[4][iEq];
     result[5][iEq] = -result[5][iEq];

  }


}

//////////////////////////////////////////////////////////////////////////////
void NavierStokesTermHOIsoP2::fluctuation_diff_galerkin(const CFuint& i0, const CFuint& i1 , const CFuint& i2){


  // Coordinate of the vertex of the element

//   m_x[V0] = (*(*m_cellNodes)[0])[XX];
//   m_x[V1] = (*(*m_cellNodes)[1])[XX];
//   m_x[V2] = (*(*m_cellNodes)[2])[XX];
//   m_x[V3] = (*(*m_cellNodes)[3])[XX];
//   m_x[V4] = (*(*m_cellNodes)[4])[XX];
//   m_x[V5] = (*(*m_cellNodes)[5])[XX];
// 
//   m_y[V0] = (*(*m_cellNodes)[0])[YY];
//   m_y[V1] = (*(*m_cellNodes)[1])[YY];
//   m_y[V2] = (*(*m_cellNodes)[2])[YY];
//   m_y[V3] = (*(*m_cellNodes)[3])[YY];
//   m_y[V4] = (*(*m_cellNodes)[4])[YY];
//   m_y[V5] = (*(*m_cellNodes)[5])[YY];

//   InwardNormalsData& cellnormals = (*socket_normals.getDataHandle()[m_cellID]);

//   cf_assert (cellnormals.nbFaces()  == 9); // P2 triangles have 9 normals

  const CFuint nbStates = m_states.size();
  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();

//   const CFreal one_eighth = 1.0/8.0;


  for (CFuint iStates = 0; iStates < nbStates; ++ iStates ){
     for (CFuint iEq = 0; iEq < nbEqs; ++ iEq)
           (*m_phi_diff_gal[iStates])[iEq] = 0.0;
  }


  ///MAIN LOOP OVER QUADRATURE POINTS
  for (CFuint iQd = 0; iQd < 4; ++iQd) {


      /// Compute derivatives of shape functions with respect to x and y (i.e. in physical space):
      // Values of shape functions in reference space:
      const CFreal xi =  m_qd0[iQd]*xi_ref[i0]+m_qd1[iQd]*xi_ref[i1]+m_qd2[iQd]*xi_ref[i2];
      const CFreal eta = m_qd0[iQd]*eta_ref[i0]+m_qd1[iQd]*eta_ref[i1]+m_qd2[iQd]*eta_ref[i2];
      const CFreal L0 = 1.0-xi-eta;
      const CFreal L1 = xi;
      const CFreal L2 = eta;

    //Compute gradients of shape functions in physical space
    ComputeSFGradients((*m_cellNodes),xi,eta,m_J,m_dNdx,m_dNdy);

    m_gradState_x = 0.0;
    m_gradState_y = 0.0;

        for(CFuint inode = 0; inode < NNODES; ++inode) {
          m_gradState_x += m_dNdx[inode]*(*m_states[inode]);
          m_gradState_y += m_dNdy[inode]*(*m_states[inode]);
        }


     //Temperature gradient:
     const CFreal grad_T_x  = m_dNdx[V0]*m_values(3,V0)  + m_dNdx[V1]*m_values(3,V1)  +
                              m_dNdx[V2]*m_values(3,V2)  + m_dNdx[V3]*m_values(3,V3)  +
                              m_dNdx[V4]*m_values(3,V4)  + m_dNdx[V5]*m_values(3,V5);

     const CFreal grad_T_y  = m_dNdy[V0]*m_values(3,V0)  + m_dNdy[V1]*m_values(3,V1)  +
                              m_dNdy[V2]*m_values(3,V2)  + m_dNdy[V3]*m_values(3,V3)  +
                              m_dNdy[V4]*m_values(3,V4)  + m_dNdy[V5]*m_values(3,V5);


     (m_qdstates) = (L0*( 2.0*L0 - 1.0 ) * (*m_states[0])) +
                    (L1*( 2.0*L1 - 1.0 ) * (*m_states[1])) +
                    (L2*( 2.0*L2 - 1.0 ) * (*m_states[2])) +
                    (4.0*L0*L1           * (*m_states[3])) +
                    (4.0*L2*L1           * (*m_states[4])) +
                    (4.0*L0*L2           * (*m_states[5]));


     (*m_gradients[0])[XX] = m_gradState_x[0];
     (*m_gradients[0])[YY] = m_gradState_y[0];
     (*m_gradients[1])[XX] = -(m_qdstates[1]/(m_qdstates[0]*m_qdstates[0]))*m_gradState_x[0] + (1.0/m_qdstates[0])*m_gradState_x[1];
     (*m_gradients[1])[YY] = -(m_qdstates[1]/(m_qdstates[0]*m_qdstates[0]))*m_gradState_y[0] + (1.0/m_qdstates[0])*m_gradState_y[1];
     (*m_gradients[2])[XX] = -(m_qdstates[2]/(m_qdstates[0]*m_qdstates[0]))*m_gradState_x[0] + (1.0/m_qdstates[0])*m_gradState_x[2];
     (*m_gradients[2])[YY] = -(m_qdstates[2]/(m_qdstates[0]*m_qdstates[0]))*m_gradState_y[0] + (1.0/m_qdstates[0])*m_gradState_y[2];
     (*m_gradients[3])[XX] = grad_T_x;
     (*m_gradients[3])[YY] = grad_T_y;


  m_SFGrad[XX] = m_dNdx[V0];
  m_SFGrad[YY] = m_dNdy[V0];
  m_F0 = m_diffVar->getFlux(m_qdstates, m_gradients, m_SFGrad , 0.0);

  m_SFGrad[XX] = m_dNdx[V1];
  m_SFGrad[YY] = m_dNdy[V1];
  m_F1 = m_diffVar->getFlux(m_qdstates, m_gradients, m_SFGrad , 0.0);

  m_SFGrad[XX] = m_dNdx[V2];
  m_SFGrad[YY] = m_dNdy[V2];
  m_F2 = m_diffVar->getFlux(m_qdstates, m_gradients, m_SFGrad , 0.0);

  m_SFGrad[XX] = m_dNdx[V3];
  m_SFGrad[YY] = m_dNdy[V3];
  m_F3 = m_diffVar->getFlux(m_qdstates, m_gradients, m_SFGrad , 0.0);

  m_SFGrad[XX] = m_dNdx[V4];
  m_SFGrad[YY] = m_dNdy[V4];
  m_F4 = m_diffVar->getFlux(m_qdstates, m_gradients, m_SFGrad , 0.0);

  m_SFGrad[XX] = m_dNdx[V5];
  m_SFGrad[YY] = m_dNdy[V5];
  m_F5 = m_diffVar->getFlux(m_qdstates, m_gradients, m_SFGrad , 0.0);


     (*m_phi_diff_gal[0])[0] = 0.0;
     (*m_phi_diff_gal[1])[0] = 0.0;
     (*m_phi_diff_gal[2])[0] = 0.0;
     (*m_phi_diff_gal[3])[0] = 0.0;
     (*m_phi_diff_gal[4])[0] = 0.0;
     (*m_phi_diff_gal[5])[0] = 0.0;

      //Coefficient for integration: Jacobian * weight * area of reference element
      const CFreal intCoeff = m_J*m_wqd[iQd]*0.25*0.5;

      for (CFuint i = 1 ; i < 4; ++i)
	    {

          (*m_phi_diff_gal[0])[i] += m_F0[i]*intCoeff;
          (*m_phi_diff_gal[1])[i] += m_F1[i]*intCoeff;
          (*m_phi_diff_gal[2])[i] += m_F2[i]*intCoeff;
          (*m_phi_diff_gal[3])[i] += m_F3[i]*intCoeff;
          (*m_phi_diff_gal[4])[i] += m_F4[i]*intCoeff;
          (*m_phi_diff_gal[5])[i] += m_F5[i]*intCoeff;
       }

  } //Loop over quadrature points

}

//////////////////////////////////////////////////////////////////////////////
void NavierStokesTermHOIsoP2::fluctuation_diff_bubble(const CFuint& i0, const CFuint& i1, const CFuint& i2)
{

//   CFout << "\nP2P2 Vertices:\n";
//   CFout << "[" << (*(*m_cellNodes)[V0])[XX] << "," << (*(*m_cellNodes)[V0])[YY] << "]\n";
//   CFout << "[" << (*(*m_cellNodes)[V1])[XX] << "," << (*(*m_cellNodes)[V1])[YY] << "]\n";
//   CFout << "[" << (*(*m_cellNodes)[V2])[XX] << "," << (*(*m_cellNodes)[V2])[YY] << "]\n\n";

//   InwardNormalsData& cellnormals = (*socket_normals.getDataHandle()[m_cellID]);

  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  const CFreal one_third = 1.0/3.0;
  const CFreal one_fourth = 0.25;

  // coordinates of the gravity center of the sub-element
  const CFreal xig  = (xi_ref[i0]  + xi_ref[i1] + xi_ref[i2])*one_third;
  const CFreal etag = (eta_ref[i0]  + eta_ref[i1] + eta_ref[i2])*one_third;

//   m_refCoord[KSI] = xig;
//   m_refCoord[ETA] = etag;
// 
//   ComputeMappedCoord((*m_cellNodes),m_refCoord, m_physCoord);
//   const CFreal xg = m_physCoord[XX];
//   const CFreal yg = m_physCoord[YY];

//   CFout << "Center of gravity of P2P2 sub-element:\n";
//   CFout << "[xg,yg] = [" << m_physCoord << "]\n";
//   CFout << "[xg,yg] = [" << xg << "," << yg << "]\n";

  // The buble function is piecewise linear. Which means that its gradient is not continuous
  // Then, we just need one node per constant part. i choose to take one point at the gravity center of the 3 triangles
  // composed by nodes and gravity center of the sub-element.
  for (CFuint iEq = 0; iEq < nbEqs; ++ iEq)
     (*m_phi_diff_bub)[iEq] = 0.0;

  /**************************************************/
  /***          1st small triangle i0-i2-ng       ***/
  /**************************************************/

  m_refCoord[KSI] = (xig + xi_ref[i0] + xi_ref[i2])*one_third;
  m_refCoord[ETA] = (etag + eta_ref[i0] + eta_ref[i2])*one_third;

//   ComputeMappedCoord((*m_cellNodes),m_refCoord, m_physCoord);

    //Compute gradients of shape functions in physical space
    ComputeSFGradients((*m_cellNodes),m_refCoord[KSI],m_refCoord[ETA],m_J,m_dNdx,m_dNdy);

    m_gradState_x = 0.0;
    m_gradState_y = 0.0;

    //Gradients of states:
    for(CFuint inode = 0; inode < NNODES; ++inode) {
        m_gradState_x += m_dNdx[inode]*(*m_states[inode]);
        m_gradState_y += m_dNdy[inode]*(*m_states[inode]);
    }

     //Temperature gradient:
     CFreal grad_T_x  = m_dNdx[V0]*m_values(3,V0)  + m_dNdx[V1]*m_values(3,V1) +
                        m_dNdx[V2]*m_values(3,V2)  + m_dNdx[V3]*m_values(3,V3) +
                        m_dNdx[V4]*m_values(3,V4)  + m_dNdx[V5]*m_values(3,V5);

     CFreal grad_T_y  = m_dNdy[V0]*m_values(3,V0)  + m_dNdy[V1]*m_values(3,V1) +
                        m_dNdy[V2]*m_values(3,V2)  + m_dNdy[V3]*m_values(3,V3) +
                        m_dNdy[V4]*m_values(3,V4)  + m_dNdy[V5]*m_values(3,V5);

     m_L0 = 1.0 - m_refCoord[KSI] - m_refCoord[ETA];
     m_L1 = m_refCoord[KSI];
     m_L2 = m_refCoord[ETA];

     (m_qdstates) = (m_L0*( 2.0*m_L0 - 1.0 ) * (*m_states[0])) +
                    (m_L1*( 2.0*m_L1 - 1.0 ) * (*m_states[1])) +
                    (m_L2*( 2.0*m_L2 - 1.0 ) * (*m_states[2])) +
                    (4.0*m_L0*m_L1           * (*m_states[3])) +
                    (4.0*m_L2*m_L1           * (*m_states[4])) +
                    (4.0*m_L0*m_L2           * (*m_states[5]));


  (*m_gradients[0])[XX] = m_gradState_x[0];
  (*m_gradients[0])[YY] = m_gradState_y[0];
  (*m_gradients[1])[XX] = -(m_qdstates[1]/(m_qdstates[0]*m_qdstates[0]))*m_gradState_x[0] + (1.0/m_qdstates[0])*m_gradState_x[1];
  (*m_gradients[1])[YY] = -(m_qdstates[1]/(m_qdstates[0]*m_qdstates[0]))*m_gradState_y[0] + (1.0/m_qdstates[0])*m_gradState_y[1];
  (*m_gradients[2])[XX] = -(m_qdstates[2]/(m_qdstates[0]*m_qdstates[0]))*m_gradState_x[0] + (1.0/m_qdstates[0])*m_gradState_x[2];
  (*m_gradients[2])[YY] = -(m_qdstates[2]/(m_qdstates[0]*m_qdstates[0]))*m_gradState_y[0] + (1.0/m_qdstates[0])*m_gradState_y[2];
  (*m_gradients[3])[XX] = grad_T_x;
  (*m_gradients[3])[YY] = grad_T_y;


//   m_faceNormal[XX] = orient*cellnormals.getFaceNormComp(f1,XX);
//   m_faceNormal[YY] = orient*cellnormals.getFaceNormComp(f1,YY);

  ComputeBubbleGradient((*m_cellNodes),m_refCoord[KSI],m_refCoord[ETA],m_J, m_gradS);
  m_F = m_diffVar->getFlux(m_qdstates, m_gradients,m_gradS , 0.0);


  const CFreal coefT0 = 0.5*one_third*one_fourth*m_J;
  for (CFuint i = 1 ; i < nbEqs; ++i){
     (*m_phi_diff_bub)[i] = coefT0*m_F[i];
     }

  /**************************************************/
  /***         2nd small triangle  i2-i1-g     ***/
  /**************************************************/

  m_refCoord[KSI] = (xig + xi_ref[i1] + xi_ref[i2])*one_third;
  m_refCoord[ETA] = (etag + eta_ref[i1] + eta_ref[i2])*one_third;

    //Compute gradients of shape functions in physical space
    ComputeSFGradients((*m_cellNodes),m_refCoord[KSI],m_refCoord[ETA],m_J,m_dNdx,m_dNdy);

    m_gradState_x = 0.0;
    m_gradState_y = 0.0;

    //Gradients of states:
    for(CFuint inode = 0; inode < NNODES; ++inode) {
        m_gradState_x += m_dNdx[inode]*(*m_states[inode]);
        m_gradState_y += m_dNdy[inode]*(*m_states[inode]);
    }

     //Temperature gradient:
     grad_T_x  = m_dNdx[V0]*m_values(3,V0)  + m_dNdx[V1]*m_values(3,V1)  +
                 m_dNdx[V2]*m_values(3,V2)  + m_dNdx[V3]*m_values(3,V3)  +
                 m_dNdx[V4]*m_values(3,V4)  + m_dNdx[V5]*m_values(3,V5);

     grad_T_y  = m_dNdy[V0]*m_values(3,V0)  + m_dNdy[V1]*m_values(3,V1)  +
                 m_dNdy[V2]*m_values(3,V2)  + m_dNdy[V3]*m_values(3,V3)  +
                 m_dNdy[V4]*m_values(3,V4)  + m_dNdy[V5]*m_values(3,V5);

     m_L0 = 1.0 - m_refCoord[KSI] - m_refCoord[ETA];
     m_L1 = m_refCoord[KSI];
     m_L2 = m_refCoord[ETA];

     (m_qdstates) = (m_L0*( 2.0*m_L0 - 1.0 ) * (*m_states[0])) +
                    (m_L1*( 2.0*m_L1 - 1.0 ) * (*m_states[1])) +
                    (m_L2*( 2.0*m_L2 - 1.0 ) * (*m_states[2])) +
                    (4.0*m_L0*m_L1           * (*m_states[3])) +
                    (4.0*m_L2*m_L1           * (*m_states[4])) +
                    (4.0*m_L0*m_L2           * (*m_states[5]));



  (*m_gradients[0])[XX] = m_gradState_x[0];
  (*m_gradients[0])[YY] = m_gradState_y[0];
  (*m_gradients[1])[XX] = -(m_qdstates[1]/(m_qdstates[0]*m_qdstates[0]))*m_gradState_x[0] + (1.0/m_qdstates[0])*m_gradState_x[1];
  (*m_gradients[1])[YY] = -(m_qdstates[1]/(m_qdstates[0]*m_qdstates[0]))*m_gradState_y[0] + (1.0/m_qdstates[0])*m_gradState_y[1];
  (*m_gradients[2])[XX] = -(m_qdstates[2]/(m_qdstates[0]*m_qdstates[0]))*m_gradState_x[0] + (1.0/m_qdstates[0])*m_gradState_x[2];
  (*m_gradients[2])[YY] = -(m_qdstates[2]/(m_qdstates[0]*m_qdstates[0]))*m_gradState_y[0] + (1.0/m_qdstates[0])*m_gradState_y[2];
  (*m_gradients[3])[XX] = grad_T_x;
  (*m_gradients[3])[YY] = grad_T_y;


  ComputeBubbleGradient((*m_cellNodes),m_refCoord[KSI],m_refCoord[ETA],m_J, m_gradS);
  m_F = m_diffVar->getFlux(m_qdstates, m_gradients,m_gradS , 0.0);


  const CFreal coefT1 = 0.5*one_third*one_fourth*m_J;
  for (CFuint i = 1 ; i < nbEqs; ++i){
     (*m_phi_diff_bub)[i] += coefT1*m_F[i];
     }


  /**************************************************/
  /***         3rd small triangle  i0-i1-ig       ***/
  /**************************************************/

  m_refCoord[KSI] = (xig + xi_ref[i0] + xi_ref[i1])*one_third;
  m_refCoord[ETA] = (etag + eta_ref[i0] + eta_ref[i1])*one_third;

    //Compute gradients of shape functions in physical space
    ComputeSFGradients((*m_cellNodes),m_refCoord[KSI],m_refCoord[ETA],m_J,m_dNdx,m_dNdy);

    m_gradState_x = 0.0;
    m_gradState_y = 0.0;

    //Gradients of states:
    for(CFuint inode = 0; inode < NNODES; ++inode) {
        m_gradState_x += m_dNdx[inode]*(*m_states[inode]);
        m_gradState_y += m_dNdy[inode]*(*m_states[inode]);
    }

     //Temperature gradient:
     grad_T_x  = m_dNdx[V0]*m_values(3,V0)  + m_dNdx[V1]*m_values(3,V1)  +
                 m_dNdx[V2]*m_values(3,V2)  + m_dNdx[V3]*m_values(3,V3)  +
                 m_dNdx[V4]*m_values(3,V4)  + m_dNdx[V5]*m_values(3,V5);

     grad_T_y  = m_dNdy[V0]*m_values(3,V0)  + m_dNdy[V1]*m_values(3,V1)  +
                 m_dNdy[V2]*m_values(3,V2)  + m_dNdy[V3]*m_values(3,V3)  +
                 m_dNdy[V4]*m_values(3,V4)  + m_dNdy[V5]*m_values(3,V5);

     m_L0 = 1.0 - m_refCoord[KSI] - m_refCoord[ETA];
     m_L1 = m_refCoord[KSI];
     m_L2 = m_refCoord[ETA];

     (m_qdstates) = (m_L0*( 2.0*m_L0 - 1.0 ) * (*m_states[0])) +
                    (m_L1*( 2.0*m_L1 - 1.0 ) * (*m_states[1])) +
                    (m_L2*( 2.0*m_L2 - 1.0 ) * (*m_states[2])) +
                    (4.0*m_L0*m_L1           * (*m_states[3])) +
                    (4.0*m_L2*m_L1           * (*m_states[4])) +
                    (4.0*m_L0*m_L2           * (*m_states[5]));

  (*m_gradients[0])[XX] = m_gradState_x[0];
  (*m_gradients[0])[YY] = m_gradState_y[0];
  (*m_gradients[1])[XX] = -(m_qdstates[1]/(m_qdstates[0]*m_qdstates[0]))*m_gradState_x[0] + (1.0/m_qdstates[0])*m_gradState_x[1];
  (*m_gradients[1])[YY] = -(m_qdstates[1]/(m_qdstates[0]*m_qdstates[0]))*m_gradState_y[0] + (1.0/m_qdstates[0])*m_gradState_y[1];
  (*m_gradients[2])[XX] = -(m_qdstates[2]/(m_qdstates[0]*m_qdstates[0]))*m_gradState_x[0] + (1.0/m_qdstates[0])*m_gradState_x[2];
  (*m_gradients[2])[YY] = -(m_qdstates[2]/(m_qdstates[0]*m_qdstates[0]))*m_gradState_y[0] + (1.0/m_qdstates[0])*m_gradState_y[2];
  (*m_gradients[3])[XX] = grad_T_x;
  (*m_gradients[3])[YY] = grad_T_y;

  ComputeBubbleGradient((*m_cellNodes),m_refCoord[KSI],m_refCoord[ETA],m_J, m_gradS);
  m_F = m_diffVar->getFlux(m_qdstates, m_gradients,m_gradS , 0.0);


  const CFreal coefT2 = 0.5*one_third*one_fourth*m_J;
  for (CFuint i = 1 ; i < nbEqs; ++i){
     (*m_phi_diff_bub)[i] += coefT2*m_F[i];
     }

}
//////////////////////////////////////////////////////////////////////////////

void NavierStokesTermHOIsoP2::setup()
{
  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  const CFuint nbStatesInControlVolume =
    MeshDataStack::getActive()->Statistics().getMaxNbStatesInCell();
  const CFuint nbNodesInControlVolume =
    MeshDataStack::getActive()->Statistics().getMaxNbNodesInCell();
  
  m_states.resize(nbStatesInControlVolume);
  m_values.resize(nbEqs, nbStatesInControlVolume);
  m_gradients.resize(nbEqs);
  for (CFuint i = 0; i< nbEqs; ++i) {
    m_gradients[i] = new RealVector(PhysicalModelStack::getActive()->getDim());
  }

  m_avValues.resize(nbEqs);
  m_normal.resize(PhysicalModelStack::getActive()->getDim());

  // Setup : should create a function to do all this setup
  const CFuint nbQdPts = 4;
  m_qd0.resize(nbQdPts); // quadrature points per face
  m_qd1.resize(nbQdPts); // quadrature points per face
  m_qd2.resize(nbQdPts);

  m_qd0[0] = 1.0/3.0;  m_qd1[0] = 1.0/3.0; m_qd2[0] = 1.0/3.0;
  m_qd0[1] = 0.6;  m_qd1[1] = 0.2; m_qd2[1] = 0.2;
  m_qd0[2] = 0.2;  m_qd1[2] = 0.6; m_qd2[2] = 0.2;
  m_qd0[3] = 0.2;  m_qd1[3] = 0.2; m_qd2[3] = 0.6;

  m_wqd.resize(nbQdPts); // 4 quadrature points per surface

  m_wqd[0] = -27.0/48.0;
  m_wqd[1] = 25.0/48.0;
  m_wqd[2] = 25.0/48.0;
  m_wqd[3] = 25.0/48.0;


 ///Quadrature with 3 points on the face of triangular cell
//  const CFreal s  = std::sqrt( 0.6 );
//  const CFreal a0 = ( 1.0 - s )*0.5;
//  const CFreal a1 = ( 1.0 + s )*0.5;
// 
//   m_faceQd0.resize(3);
//   m_faceQd1.resize(3);
//   m_faceWQd.resize(3);
// 
//   m_faceQd0[0] = a0;  m_faceQd1[0] = a1;
//   m_faceQd0[1] = a1;  m_faceQd1[1] = a0;
//   m_faceQd0[2] = .5;  m_faceQd1[2] = .5;
// 
//   m_faceWQd[0] = 5.0/18.0;
//   m_faceWQd[1] = 5.0/18.0;
//   m_faceWQd[2] = 8.0/18.0;



  m_qdstates.resize(nbEqs);


  ///Coordinates of triangle, shape functions and their derivatives
//   m_x.resize(NNODES);
//   m_y.resize(NNODES);

//   m_dNdxi.resize(NNODES);
//   m_dNdeta.resize(NNODES);
  m_dNdx.resize(NNODES);
  m_dNdy.resize(NNODES);

  m_dxdxi = m_dxdeta = m_dydxi = m_dydeta = 0.0;

  m_gradState_x.resize(nbEqs);
  m_gradState_y.resize(nbEqs);

  m_F0.resize(nbEqs);
  m_F1.resize(nbEqs);
  m_F2.resize(nbEqs);
  m_F3.resize(nbEqs);
  m_F4.resize(nbEqs);
  m_F5.resize(nbEqs);
  m_F.resize(nbEqs);

  m_faceNormal0.resize(DIM_2D);
  m_faceNormal1.resize(DIM_2D);
  m_SFGrad.resize(DIM_2D);
  m_refCoord.resize(DIM_2D);
  m_physCoord.resize(DIM_2D);

  m_gradS.resize(DIM_2D);

  m_phi_diff_bub = new RealVector(nbEqs);
  m_phi_diff_bub_split.resize(nbNodesInControlVolume); // 3 residuals in each sub element

  m_phi_diff_bub_split[0] = new RealVector(nbEqs);
  m_phi_diff_bub_split[1] = new RealVector(nbEqs);
  m_phi_diff_bub_split[2] = new RealVector(nbEqs);

  m_phi_diff_gal.resize(nbStatesInControlVolume); // one galerkin function per state

  for (CFuint iState = 0; iState < nbStatesInControlVolume; ++ iState)
     m_phi_diff_gal[iState] = new RealVector(nbEqs);

  kappa.resize(4,6); // 4 sub-elements 6 nodes

  CFreal one_fourth = 1.0/4.0;

  // kappa is used to "distribute" the part of the fluctuation with the bubble function
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

  getMethodData().getDistributionData().computeBetas = true;
}

//////////////////////////////////////////////////////////////////////////////

void NavierStokesTermHOIsoP2::configure ( Config::ConfigArgs& args )
{
  ComputeDiffusiveTerm::configure(args);
}

//////////////////////////////////////////////////////////////////////////////

void NavierStokesTermHOIsoP2::ComputeMappedCoord(const std::vector<Framework::Node*> & nodes, const RealVector & refcoord, RealVector & physcoord)
{

  const CFreal L0 = 1.0 - refcoord[KSI] - refcoord[ETA];
  const CFreal L1 = refcoord[KSI];
  const CFreal L2 = refcoord[ETA];

  const CFreal N0 = L0*(2.0*L0-1.0);
  const CFreal N1 = L1*(2.0*L1-1.0);
  const CFreal N2 = L2*(2.0*L2-1.0);
  const CFreal N3 = 4.0*L0*L1;
  const CFreal N4 = 4.0*L2*L1;
  const CFreal N5 = 4.0*L0*L2;

  physcoord[XX] = (*nodes[V0])[XX]*N0 + (*nodes[V1])[XX]*N1 + (*nodes[V2])[XX]*N2 + \
                  (*nodes[V3])[XX]*N3 + (*nodes[V4])[XX]*N4 + (*nodes[V5])[XX]*N5;
  physcoord[YY] = (*nodes[V0])[YY]*N0 + (*nodes[V1])[YY]*N1 + (*nodes[V2])[YY]*N2 + \
                  (*nodes[V3])[YY]*N3 + (*nodes[V4])[YY]*N4 + (*nodes[V5])[YY]*N5;

}

//////////////////////////////////////////////////////////////////////////////

void NavierStokesTermHOIsoP2::ComputeSFGradients(const std::vector<Framework::Node*>& nodes, const CFreal xi, const CFreal eta, CFreal J, RealVector & dNdx, RealVector & dNdy)
{

      const CFreal L0 = 1.0-xi-eta;
      const CFreal L1 = xi;
      const CFreal L2 = eta;

      //Derivatives of shape functions in reference space
      const CFreal dNdxi0 = (1.0-4.0*L0);
      const CFreal dNdxi1 = 4.0*L1-1.0;
      const CFreal dNdxi2 = 0.0;
      const CFreal dNdxi3 = 4.0*(L0-L1);
      const CFreal dNdxi4 = 4.0*L2;
      const CFreal dNdxi5 = -4.0*L2;

      const CFreal dNdeta0 = (1.0-4.0*L0);
      const CFreal dNdeta1 = 0.0;
      const CFreal dNdeta2 = 4.0*L2-1.0;
      const CFreal dNdeta3 = -4.0*L1;
      const CFreal dNdeta4 = 4.0*L1;
      const CFreal dNdeta5 = 4.0*(L0-L2);

      m_dxdxi  = 0.0;
      m_dxdeta = 0.0;
      m_dydxi  = 0.0;
      m_dydeta = 0.0;

    //Components of 2x2 Jacobian of the transformation 
    // physical -> reference space

    const CFreal dxdxi = (*nodes[V0])[XX]*dNdxi0 + (*nodes[V1])[XX]*dNdxi1 + (*nodes[V2])[XX]*dNdxi2 + \
                         (*nodes[V3])[XX]*dNdxi3 + (*nodes[V4])[XX]*dNdxi4 + (*nodes[V5])[XX]*dNdxi5;

    const CFreal dxdeta = (*nodes[V0])[XX]*dNdeta0 + (*nodes[V1])[XX]*dNdeta1 + (*nodes[V2])[XX]*dNdeta2 + \
                          (*nodes[V3])[XX]*dNdeta3 + (*nodes[V4])[XX]*dNdeta4 + (*nodes[V5])[XX]*dNdeta5;

    const CFreal dydxi = (*nodes[V0])[YY]*dNdxi0 + (*nodes[V1])[YY]*dNdxi1 + (*nodes[V2])[YY]*dNdxi2 + \
                         (*nodes[V3])[YY]*dNdxi3 + (*nodes[V4])[YY]*dNdxi4 + (*nodes[V5])[YY]*dNdxi5;

    const CFreal dydeta = (*nodes[V0])[YY]*dNdeta0 + (*nodes[V1])[YY]*dNdeta1 + (*nodes[V2])[YY]*dNdeta2 + \
                          (*nodes[V3])[YY]*dNdeta3 + (*nodes[V4])[YY]*dNdeta4 + (*nodes[V5])[YY]*dNdeta5;

    J = dxdxi*dydeta-dxdeta*dydxi;
    const CFreal invDet = 1.0/J;

    dNdx[V0] = invDet*( dNdxi0*dydeta - dNdeta0*dydxi);
    dNdy[V0] = invDet*(-dNdxi0*dxdeta + dNdeta0*dxdxi); 

    dNdx[V1] = invDet*( dNdxi1*dydeta - dNdeta1*dydxi);
    dNdy[V1] = invDet*(-dNdxi1*dxdeta + dNdeta1*dxdxi); 

    dNdx[V2] = invDet*( dNdxi2*dydeta - dNdeta2*dydxi);
    dNdy[V2] = invDet*(-dNdxi2*dxdeta + dNdeta2*dxdxi); 

    dNdx[V3] = invDet*( dNdxi3*dydeta - dNdeta3*dydxi);
    dNdy[V3] = invDet*(-dNdxi3*dxdeta + dNdeta3*dxdxi); 

    dNdx[V4] = invDet*( dNdxi4*dydeta - dNdeta4*dydxi);
    dNdy[V4] = invDet*(-dNdxi4*dxdeta + dNdeta4*dxdxi); 

    dNdx[V5] = invDet*( dNdxi5*dydeta - dNdeta5*dydxi);
    dNdy[V5] = invDet*(-dNdxi5*dxdeta + dNdeta5*dxdxi); 


//     for(CFuint inode = 0; inode<NNODES;++inode) {
//       m_dxdxi  += m_dNdxi[inode]*m_x[inode];
//       m_dxdeta += m_dNdeta[inode]*m_x[inode];
//       m_dydxi  += m_dNdxi[inode]*m_y[inode];
//       m_dydeta += m_dNdeta[inode]*m_y[inode];
//     } 
//       
// //     const CFreal invDet = 1.0/(m_dxdxi*m_dydeta-m_dxdeta*m_dydxi);
//       const CFreal invDet = 1.0;
//       
//     //Derivatives of shape functions in physical space:
//     for(CFuint inode=0; inode<NNODES;++inode) {
//       dNdx[inode] = invDet*(m_dNdxi[inode]*m_dydeta - m_dNdeta[inode]*m_dydxi);
//       dNdy[inode] = invDet*(-m_dNdxi[inode]*m_dxdeta + m_dNdeta[inode]*m_dxdxi); 
//     } 
//       

}

//////////////////////////////////////////////////////////////////////////////

void NavierStokesTermHOIsoP2::ComputeBubbleGradient(const std::vector<Framework::Node*>& nodes, const CFreal xi, const CFreal eta, CFreal J, RealVector & gradS)
{

      const CFreal L0 = 1.0-xi-eta;
      const CFreal L1 = xi;
      const CFreal L2 = eta;

      //Derivatives of shape functions in reference space
      const CFreal dNdxi0 = (1.0-4.0*L0);
      const CFreal dNdxi1 = 4.0*L1-1.0;
      const CFreal dNdxi2 = 0.0;
      const CFreal dNdxi3 = 4.0*(L0-L1);
      const CFreal dNdxi4 = 4.0*L2;
      const CFreal dNdxi5 = -4.0*L2;

      const CFreal dNdeta0 = (1.0-4.0*L0);
      const CFreal dNdeta1 = 0.0;
      const CFreal dNdeta2 = 4.0*L2-1.0;
      const CFreal dNdeta3 = -4.0*L1;
      const CFreal dNdeta4 = 4.0*L1;
      const CFreal dNdeta5 = 4.0*(L0-L2);

      m_dxdxi  = 0.0;
      m_dxdeta = 0.0;
      m_dydxi  = 0.0;
      m_dydeta = 0.0;

    //Components of 2x2 Jacobian of the transformation 
    // physical -> reference space

    const CFreal dxdxi = (*nodes[V0])[XX]*dNdxi0 + (*nodes[V1])[XX]*dNdxi1 + (*nodes[V2])[XX]*dNdxi2 + \
                         (*nodes[V3])[XX]*dNdxi3 + (*nodes[V4])[XX]*dNdxi4 + (*nodes[V5])[XX]*dNdxi5;

    const CFreal dxdeta = (*nodes[V0])[XX]*dNdeta0 + (*nodes[V1])[XX]*dNdeta1 + (*nodes[V2])[XX]*dNdeta2 + \
                          (*nodes[V3])[XX]*dNdeta3 + (*nodes[V4])[XX]*dNdeta4 + (*nodes[V5])[XX]*dNdeta5;

    const CFreal dydxi = (*nodes[V0])[YY]*dNdxi0 + (*nodes[V1])[YY]*dNdxi1 + (*nodes[V2])[YY]*dNdxi2 + \
                         (*nodes[V3])[YY]*dNdxi3 + (*nodes[V4])[YY]*dNdxi4 + (*nodes[V5])[YY]*dNdxi5;

    const CFreal dydeta = (*nodes[V0])[YY]*dNdeta0 + (*nodes[V1])[YY]*dNdeta1 + (*nodes[V2])[YY]*dNdeta2 + \
                          (*nodes[V3])[YY]*dNdeta3 + (*nodes[V4])[YY]*dNdeta4 + (*nodes[V5])[YY]*dNdeta5;

    J = dxdxi*dydeta-dxdeta*dydxi;
    const CFreal invDet = 1.0/J;

    //Bubble function: S = 27 * L0 * L1 * L2
    const CFreal dSdxi = 27.0*L2*(L0-L1);
    const CFreal dSdeta = 27.0*L1*(L0-L2);

    gradS[XX] = invDet*( dSdxi*dydeta - dSdeta*dydxi);
    gradS[YY] = invDet*(-dSdxi*dxdeta + dSdeta*dxdxi); 

}


//////////////////////////////////////////////////////////////////////////////




    } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
