#include "Environment/ObjectProvider.hh"

#include "Framework/GeometricEntity.hh"
#include "Framework/MeshData.hh"
#include "Framework/MethodStrategyProvider.hh"

#include "LinearAdv/AdvectionDiffusionVarSet.hh"
#include "LinearAdv/LinearAdv2DVarSet.hh"

#include "FluctSplit/FluctSplitAdvectionDiffusion.hh"
#include "FluctSplit/FluctuationSplitStrategy.hh"
#include "FluctSplit/ScalarDiffusionTermHO.hh"
#include "FluctSplit/InwardNormalsData.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Physics::LinearAdv;
using namespace COOLFluiD::MathTools;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

MethodStrategyProvider<ScalarDiffusionTermHO,
		       FluctuationSplitData,
		       ComputeDiffusiveTerm,
		       FluctSplitAdvectionDiffusionModule>
scalardiffusionDiffusiveTermHOProvider("ScalarDiffusionHO");

//////////////////////////////////////////////////////////////////////////////

ScalarDiffusionTermHO::ScalarDiffusionTermHO(const std::string& name) :
  ComputeDiffusiveTerm(name),
  _diffVar(CFNULL),
  _updateVar(CFNULL),
  _states(),
  _values(),
  _gradients(),
  _avValues(),
  _normal()
{
}

//////////////////////////////////////////////////////////////////////////////

ScalarDiffusionTermHO::~ScalarDiffusionTermHO()
{
  for (CFuint i = 0; i< _values.size(); ++i) {
    deletePtr(_values[i]);
  }

  for (CFuint i = 0; i< _gradients.size(); ++i) {
    deletePtr(_gradients[i]);
  }
}

//////////////////////////////////////////////////////////////////////////////

void ScalarDiffusionTermHO::setDiffusiveVarSet(SafePtr<DiffusiveVarSet> diffVar)
{
  _diffVar = diffVar.d_castTo<Physics::LinearAdv::AdvectionDiffusionVarSet>();
}

//////////////////////////////////////////////////////////////////////////////

void ScalarDiffusionTermHO::setUpdateVarSet(SafePtr<ConvectiveVarSet> updateVar)
{
  _updateVar = updateVar.d_castTo<Physics::LinearAdv::LinearAdv2DVarSet>();
}


//////////////////////////////////////////////////////////////////////////////

void ScalarDiffusionTermHO::computeDiffusiveTerm
(GeometricEntity *const geo, vector<RealVector>& result, bool updateCoeffFlag)
{
  DataHandle< InwardNormalsData*> normals = socket_normals.getDataHandle();
  DataHandle< CFreal> updateCoeff = socket_updateCoeff.getDataHandle();
  
  // Vector of the states of thecell
  m_cellStates = geo->getStates();
  const CFuint nbCellStates = m_cellStates->size();
  m_cellID = getMethodData().getDistributionData().cellID;
  
  // store the pointers to state in another array (of RealVector*)
  for (CFuint i = 0; i < nbCellStates; ++i) {
    _states[i] = (*m_cellStates)[i];
  }
  
  // compute vars that will be used to compute the gradients
  _diffVar->setGradientVars(_states, _values, geo->nbStates());

  // compute the radius (axysimmetric computations)
//   CFreal radius = 0.0;

  //   if (getMethodData().isAxisymmetric()) {
  //     for (CFuint i = 0; i < nbCellStates; ++i) {
  //       const Node& node = *geo->getNode(i);
  //       radius += node[YY];
  //     }
  //     radius /= nbCellStates;
  //   }
  _cellVolume = geo->computeVolume();
// unused //  const CFuint dim = PhysicalModelStack::getActive()->getDim();
// unused //  const CFreal dimCoeff = 1./dim;


  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();

  for (CFuint iState = 0; iState < nbCellStates; ++iState)
    for (CFuint iEq = 0; iEq < nbEqs; ++iEq)
        result[iState][iEq] = 0.0;

  ADTerm& model = _diffVar->getModel();
const CFreal nu = model.getPhysicalData()[ADTerm::NU];
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
    const CFreal diffCoeff = nu;
    CFreal coeff = diffCoeff/(_cellVolume*4.0);
    if (updateCoeffFlag) {
      for (CFuint i = 0 ; i < 3; ++i) { 
        const CFreal faceArea = (normals[m_cellID]->getAreaNode(i));

        updateCoeff[geo->getState(i)->getLocalID()] += coeff*faceArea*faceArea;
      } 

      coeff = 2.0*diffCoeff/(3.0*_cellVolume);
      CFreal dot_prod = (normals[m_cellID]->getNodalNormComp(0,XX))*(normals[m_cellID]->getNodalNormComp(1,XX)) +
                        (normals[m_cellID]->getNodalNormComp(0,YY))*(normals[m_cellID]->getNodalNormComp(1,YY));
  
      CFreal faceArea1 = (normals[m_cellID]->getAreaNode(0));
      CFreal faceArea2 = (normals[m_cellID]->getAreaNode(1));
      updateCoeff[geo->getState(3)->getLocalID()] += coeff*(faceArea1*faceArea1 + faceArea2*faceArea2 + dot_prod);
        

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
     result[0][iEq] = -result[0][iEq];
     result[1][iEq] = -result[1][iEq];
     result[2][iEq] = -result[2][iEq];
     result[3][iEq] = -result[3][iEq];
     result[4][iEq] = -result[4][iEq];
     result[5][iEq] = -result[5][iEq];

  }


}

void ScalarDiffusionTermHO::fluctuation_diff_galerkin(CFuint& i1, CFuint& i2 , CFuint& i3){

  // Coordinate of the vertex of the element
  const Node& node0_vert = (*m_cellStates)[0]->getCoordinates();
  const Node& node1_vert = (*m_cellStates)[1]->getCoordinates();
  const Node& node2_vert = (*m_cellStates)[2]->getCoordinates();

  const CFreal x1_vert = node0_vert[XX];
  const CFreal x2_vert = node1_vert[XX];
  const CFreal x3_vert = node2_vert[XX];

  const CFreal y1_vert = node0_vert[YY];
  const CFreal y2_vert = node1_vert[YY];
  const CFreal y3_vert = node2_vert[YY];

  //Coordinates of the node defining the surface where the integral is computed
  //  This correspond to vertex of the sub-element on which we integrate
  Node& node0 = (*m_cellStates)[i1]->getCoordinates();
  Node& node1 = (*m_cellStates)[i2]->getCoordinates();
  Node& node2 = (*m_cellStates)[i3]->getCoordinates();

  const CFuint nbStates = _states.size();
  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();

  InwardNormalsData& cellnormals = (*socket_normals.getDataHandle()[m_cellID]);

  const CFreal nx1 = cellnormals.getNodalNormComp(0,XX);
  const CFreal nx2 = cellnormals.getNodalNormComp(1,XX);
  const CFreal nx3 = cellnormals.getNodalNormComp(2,XX);

  const CFreal ny1 = cellnormals.getNodalNormComp(0,YY);
  const CFreal ny2 = cellnormals.getNodalNormComp(1,YY);
  const CFreal ny3 = cellnormals.getNodalNormComp(2,YY);

  const CFuint dim = PhysicalModelStack::getActive()->getDim();
  const CFreal dimCoeff = 1./dim;
  const CFreal coeffGrad = dimCoeff/_cellVolume;
  const CFreal inv_volume = 1.0/_cellVolume;
  const CFreal one_eighth = 1.0/8.0;

  for (CFuint iStates = 0; iStates < nbStates; ++ iStates ){
     for (CFuint iEq = 0; iEq < nbEqs; ++ iEq)
           (*_phi_diff_gal[iStates])[iEq] = 0.0;
  }

  for (CFuint iQd = 0; iQd < 4; ++iQd) {
     //point od quadrature
     const CFreal x = qd0[iQd] * node0[XX] + qd1[iQd] * node1[XX] + qd2[iQd] * node2[XX] ;
     const CFreal y = qd0[iQd] * node0[YY] + qd1[iQd] * node1[YY] + qd2[iQd] * node2[YY] ;

     // Linear basis function
     CFreal L1 = 1.0 + 0.5*( ( x - x1_vert )*nx1 + ( y - y1_vert )*ny1 )*inv_volume ;
     CFreal L2 = 1.0 + 0.5*( ( x - x2_vert )*nx2 + ( y - y2_vert )*ny2 )*inv_volume ;
     CFreal L3 = 1.0 + 0.5*( ( x - x3_vert )*nx3 + ( y - y3_vert )*ny3 )*inv_volume ;

     for (CFuint i = 0; i< nbEqs; ++i){
         (*_gradients[i])[XX] = (nx1*(4.0*L1 - 1.0)*(*_states[0])[i] +
		            nx2*(4.0*L2 - 1.0)*(*_states[1])[i] +
		            nx3*(4.0*L3 - 1.0)*(*_states[2])[i] +
		            4.0*(nx1*L2 + nx2*L1)*(*_states[3])[i] +
		            4.0*(nx2*L3 + nx3*L2)*(*_states[4])[i] +
		            4.0*(nx3*L1 + nx1*L3)*(*_states[5])[i])*coeffGrad;

        (*_gradients[i])[YY] = (ny1*(4.0*L1 - 1.0)*(*_states[0])[i] +
		           ny2*(4.0*L2 - 1.0)*(*_states[1])[i] +
		           ny3*(4.0*L3 - 1.0)*(*_states[2])[i] +
		           4.0*(ny1*L2 + ny2*L1)*(*_states[3])[i] +
		           4.0*(ny2*L3 + ny3*L2)*(*_states[4])[i] +
		           4.0*(ny3*L1 + ny1*L3)*(*_states[5])[i])*coeffGrad;

        }

     _normal[XX] = nx1;
     _normal[YY] = ny1;

     RealVector F1;
     F1.resize(nbEqs);

     F1 = _diffVar->getFlux(qdstates, _gradients,_normal , 0.0);

     _normal[XX] = nx2;
     _normal[YY] = ny2;

     RealVector F2;
     F2.resize(nbEqs);
     F2 = _diffVar->getFlux(qdstates, _gradients,_normal , 0.0);

     _normal[XX] = nx3;
     _normal[YY] = ny3;

     RealVector F3;
     F3.resize(nbEqs);
     F3 = _diffVar->getFlux(qdstates, _gradients,_normal , 0.0);


      for (CFuint i = 0 ; i < nbEqs; ++i)
	{
          (*_phi_diff_gal[0])[i] += F1[i]*(4.0*L1 - 1.0)*one_eighth*wqd[iQd];
          (*_phi_diff_gal[1])[i] += F2[i]*(4.0*L2 - 1.0)*one_eighth*wqd[iQd];
          (*_phi_diff_gal[2])[i] += F3[i]*(4.0*L3 - 1.0)*one_eighth*wqd[iQd];
          (*_phi_diff_gal[3])[i] += 4.0*(F1[i]*L2 + F2[i]*L1)*one_eighth*wqd[iQd];
          (*_phi_diff_gal[4])[i] += 4.0*(F2[i]*L3 + F3[i]*L2)*one_eighth*wqd[iQd];
          (*_phi_diff_gal[5])[i] += 4.0*(F1[i]*L3 + F3[i]*L1)*one_eighth*wqd[iQd];
       }

  }

}

void ScalarDiffusionTermHO::fluctuation_diff_bubble(CFuint& i1, CFuint& i2, CFuint& i3)
{
  Node& node0_vert = (*m_cellStates)[0]->getCoordinates();
  Node& node1_vert = (*m_cellStates)[1]->getCoordinates();
  Node& node2_vert = (*m_cellStates)[2]->getCoordinates();

  // Coordinate of the vertex
  const CFreal x1_vert = node0_vert[XX];
  const CFreal x2_vert = node1_vert[XX];
  const CFreal x3_vert = node2_vert[XX];

  const CFreal y1_vert = node0_vert[YY];
  const CFreal y2_vert = node1_vert[YY];
  const CFreal y3_vert = node2_vert[YY];

  Node& node0 = (*m_cellStates)[i1]->getCoordinates();
  Node& node1 = (*m_cellStates)[i2]->getCoordinates();
  Node& node2 = (*m_cellStates)[i3]->getCoordinates();

  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();

  InwardNormalsData& cellnormals = (*socket_normals.getDataHandle()[m_cellID]);

  const CFreal nx1 = cellnormals.getNodalNormComp(0,XX);
  const CFreal nx2 = cellnormals.getNodalNormComp(1,XX);
  const CFreal nx3 = cellnormals.getNodalNormComp(2,XX);

  const CFreal ny1 = cellnormals.getNodalNormComp(0,YY);
  const CFreal ny2 = cellnormals.getNodalNormComp(1,YY);
  const CFreal ny3 = cellnormals.getNodalNormComp(2,YY);

  const CFuint dim = PhysicalModelStack::getActive()->getDim();
  const CFreal dimCoeff = 1./dim;
  const CFreal coeffGrad = dimCoeff/_cellVolume;
  const CFreal inv_volume = 1.0/_cellVolume;
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


  for (CFuint i = 0; i< nbEqs; ++i){
     (*_gradients[i])[XX] = (nx1*(4.0*L1 - 1.0)*(*_states[0])[i] +
		        nx2*(4.0*L2 - 1.0)*(*_states[1])[i] +
		        nx3*(4.0*L3 - 1.0)*(*_states[2])[i] +
		        4.0*(nx1*L2 + nx2*L1)*(*_states[3])[i] +
		        4.0*(nx2*L3 + nx3*L2)*(*_states[4])[i] +
		        4.0*(nx3*L1 + nx1*L3)*(*_states[5])[i])*coeffGrad;

    (*_gradients[i])[YY] = (ny1*(4.0*L1 - 1.0)*(*_states[0])[i] +
		        ny2*(4.0*L2 - 1.0)*(*_states[1])[i] +
		        ny3*(4.0*L3 - 1.0)*(*_states[2])[i] +
		        4.0*(ny1*L2 + ny2*L1)*(*_states[3])[i] +
		        4.0*(ny2*L3 + ny3*L2)*(*_states[4])[i] +
		        4.0*(ny3*L1 + ny1*L3)*(*_states[5])[i])*coeffGrad;
  }

  _normal[XX] = (node0[YY] - node2[YY]);
  _normal[YY] = (node2[XX] - node0[XX]);

  RealVector F1;
  F1.resize(nbEqs);
  F1 = _diffVar->getFlux(qdstates, _gradients,_normal , 0.0);




  for (CFuint i = 0 ; i < nbEqs; ++i){
     (*_phi_diff_bub)[i] = F1[i]*_cellVolume/
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

  for (CFuint i = 0; i< nbEqs; ++i){
     (*_gradients[i])[XX] = (nx1*(4.0*L1 - 1.0)*(*_states[0])[i] +
		        nx2*(4.0*L2 - 1.0)*(*_states[1])[i] +
		        nx3*(4.0*L3 - 1.0)*(*_states[2])[i] +
		        4.0*(nx1*L2 + nx2*L1)*(*_states[3])[i] +
		        4.0*(nx2*L3 + nx3*L2)*(*_states[4])[i] +
		        4.0*(nx3*L1 + nx1*L3)*(*_states[5])[i])*coeffGrad;

    (*_gradients[i])[YY] = (ny1*(4.0*L1 - 1.0)*(*_states[0])[i] +
		        ny2*(4.0*L2 - 1.0)*(*_states[1])[i] +
		        ny3*(4.0*L3 - 1.0)*(*_states[2])[i] +
		        4.0*(ny1*L2 + ny2*L1)*(*_states[3])[i] +
		        4.0*(ny2*L3 + ny3*L2)*(*_states[4])[i] +
		        4.0*(ny3*L1 + ny1*L3)*(*_states[5])[i])*coeffGrad;
  }

  _normal[XX] = (node1[YY] - node2[YY]);
  _normal[YY] = (node2[XX] - node1[XX]);

  F1 = _diffVar->getFlux(qdstates, _gradients,_normal , 0.0);

  for (CFuint i = 0 ; i < nbEqs; ++i){
     (*_phi_diff_bub)[i] += F1[i]*_cellVolume/
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

  for (CFuint i = 0; i< nbEqs; ++i){
     (*_gradients[i])[XX] = (nx1*(4.0*L1 - 1.0)*(*_states[0])[i] +
		        nx2*(4.0*L2 - 1.0)*(*_states[1])[i] +
		        nx3*(4.0*L3 - 1.0)*(*_states[2])[i] +
		        4.0*(nx1*L2 + nx2*L1)*(*_states[3])[i] +
		        4.0*(nx2*L3 + nx3*L2)*(*_states[4])[i] +
		        4.0*(nx3*L1 + nx1*L3)*(*_states[5])[i])*coeffGrad;

    (*_gradients[i])[YY] = (ny1*(4.0*L1 - 1.0)*(*_states[0])[i] +
		        ny2*(4.0*L2 - 1.0)*(*_states[1])[i] +
		        ny3*(4.0*L3 - 1.0)*(*_states[2])[i] +
		        4.0*(ny1*L2 + ny2*L1)*(*_states[3])[i] +
		        4.0*(ny2*L3 + ny3*L2)*(*_states[4])[i] +
		        4.0*(ny3*L1 + ny1*L3)*(*_states[5])[i])*coeffGrad;
  }

  _normal[XX] = (node0[YY] - node1[YY]);
  _normal[YY] = (node1[XX] - node0[XX]);

  F1 = _diffVar->getFlux(qdstates, _gradients, _normal , 0.0);

  for (CFuint i = 0 ; i < nbEqs; ++i){
      (*_phi_diff_bub)[i] += F1[i]*_cellVolume/
	                     (12.0*((node1[XX]*yg - xg*node1[YY]) - node0[XX]*(yg - node1[YY]) + node0[YY]*(xg - node1[XX])));

      }

}
//////////////////////////////////////////////////////////////////////////////

void ScalarDiffusionTermHO::setup()
{

  _values.resize(PhysicalModelStack::getActive()->getNbEq());
  _gradients.resize(PhysicalModelStack::getActive()->getNbEq());
  _avValues.resize(PhysicalModelStack::getActive()->getNbEq());
  _normal.resize(PhysicalModelStack::getActive()->getDim());
  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  const CFuint nbStatesInControlVolume =
    MeshDataStack::getActive()->Statistics().getMaxNbStatesInCell();
  const CFuint nbNodesInControlVolume =
    MeshDataStack::getActive()->Statistics().getMaxNbNodesInCell();

  _states.resize(nbStatesInControlVolume);

  for (CFuint i = 0; i< nbEqs; ++i) {
    _values[i] = new RealVector(nbStatesInControlVolume);
  }

  for (CFuint i = 0; i< nbEqs; ++i) {
    _gradients[i] = new RealVector(PhysicalModelStack::getActive()->getDim());
  }

  // Setup : should create a function to do all this setup
  const CFuint nbQdPts = 4;
  qd0.resize(nbQdPts); // quadrature points per face
  qd1.resize(nbQdPts); // quadrature points per face
  qd2.resize(nbQdPts);

  qd0[0] = 1.0/3.0;  qd1[0] = 1.0/3.0; qd2[0] = 1.0/3.0;
  qd0[1] = 0.6;  qd1[1] = 0.2; qd2[1] = 0.2;
  qd0[2] = 0.2;  qd1[2] = 0.6; qd2[2] = 0.2;
  qd0[3] = 0.2;  qd1[3] = 0.2; qd2[3] = 0.6;

  wqd.resize(nbQdPts); // 4 quadrature points per surface

  wqd[0] = -27.0/48.0;
  wqd[1] = 25.0/48.0;
  wqd[2] = 25.0/48.0;
  wqd[3] = 25.0/48.0;

  qdstates.resize(nbEqs);

  _phi_diff_bub = new RealVector(nbEqs);
  _phi_diff_bub_split.resize(nbNodesInControlVolume); // 3 residuals in each sub element

  _phi_diff_bub_split[0] = new RealVector(nbEqs);
  _phi_diff_bub_split[1] = new RealVector(nbEqs);
  _phi_diff_bub_split[2] = new RealVector(nbEqs);

  _phi_diff_gal.resize(nbStatesInControlVolume); // one galerkin function per state

  for (CFuint iState = 0; iState < nbStatesInControlVolume; ++ iState)
     _phi_diff_gal[iState] = new RealVector(nbEqs);

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

void ScalarDiffusionTermHO::configure ( Config::ConfigArgs& args )
{
  ComputeDiffusiveTerm::configure(args);
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
