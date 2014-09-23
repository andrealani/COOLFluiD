#include "NavierStokesTermHOlin.hh"
#include "FluctSplit/InwardNormalsData.hh"
#include "Environment/ObjectProvider.hh"
#include "Framework/GeometricEntity.hh"
#include "Framework/MeshData.hh"
#include "Framework/MethodStrategyProvider.hh"
#include "NavierStokes/NavierStokesVarSet.hh"
#include "NavierStokes/EulerVarSet.hh"
#include "FluctSplit/HONavierStokes/FluctSplitHONavierStokes.hh"
#include "FluctSplit/FluctuationSplitStrategy.hh"

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

MethodStrategyProvider<NavierStokesTermHOlin,
		       FluctuationSplitData,
		       ComputeDiffusiveTerm,
		       FluctSplitHONavierStokesModule>
navierStokesDiffusiveTermHOlinProvider("NavierStokesHOlin");

//////////////////////////////////////////////////////////////////////////////

void NavierStokesTermHOlin::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOption< bool >("isAxisymm","Flag telling if the case is axysimmetric.");
}

//////////////////////////////////////////////////////////////////////////////

NavierStokesTermHOlin::NavierStokesTermHOlin(const std::string& name) :
  ComputeDiffusiveTerm(name),
  _diffVar(CFNULL),
  _updateVar(CFNULL),
  _states(),
  _values(),
  _gradients(),
  _avValues(),
  _normal()
{
   addConfigOptionsTo(this);
  _isAxisymm = false;
   setParameter("isAxisymm",&_isAxisymm);


}

//////////////////////////////////////////////////////////////////////////////

NavierStokesTermHOlin::~NavierStokesTermHOlin()
{
  for (CFuint i = 0; i< _gradients.size(); ++i) {
    deletePtr(_gradients[i]);
  }
}

//////////////////////////////////////////////////////////////////////////////

void NavierStokesTermHOlin::setDiffusiveVarSet(SafePtr<DiffusiveVarSet> diffVar)
{
  _diffVar = diffVar.d_castTo<NavierStokesVarSet>();
}

//////////////////////////////////////////////////////////////////////////////

void NavierStokesTermHOlin::setUpdateVarSet(SafePtr<ConvectiveVarSet> updateVar)
{
  _updateVar = updateVar.d_castTo<EulerVarSet>();
}


//////////////////////////////////////////////////////////////////////////////

void NavierStokesTermHOlin::computeDiffusiveTerm
(GeometricEntity *const geo, vector<RealVector>& result, bool updateCoeffFlag)
{
 
  DataHandle< InwardNormalsData*> normals = socket_normals.getDataHandle();
  DataHandle< CFreal> updateCoeff = socket_updateCoeff.getDataHandle();
  CFreal onethird = 1.0/3.0;

  // Vector of the states of the cell
  m_cellStates = geo->getStates();
  const CFuint nbCellStates = m_cellStates->size();
  m_cellID = getMethodData().getDistributionData().cellID;
  
  // store the pointers to state in another array (of RealVector*)
  for (CFuint i = 0; i < nbCellStates; ++i) {
    _states[i] = (*m_cellStates)[i];
  }

  // compute vars that will be used to compute the gradients
  _diffVar->setGradientVars(_states, _values, geo->nbStates());

  ///@todo axisymetry case is not done
  // compute the radius (axysimmetric computations)
  //   CFreal radius = 0.0;

  //   if (_isAxisymm) {
  //     for (CFuint i = 0; i < nbCellStates; ++i) {
  //       const Node& node = *geo->getNode(i);
  //       radius += node[YY];
  //     }
  //     radius /= nbCellStates;
  //   }

  //Variables used to compute result
  _cellVolume = geo->computeVolume();
  const CFuint dim = PhysicalModelStack::getActive()->getDim();
  const CFreal dimCoeff = 1./dim;
  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();

  for (CFuint iEq = 0; iEq < nbEqs; ++iEq)
     for (CFuint iState = 0; iState < nbCellStates; ++iState)
        result[iState][iEq] = 0.0;


  // Compute coeficients necessary to update the updatecoefficient
  nbNodes = geo->nbNodes();
  const CFreal ovDimCoeff2 = 1./(dimCoeff*dimCoeff);
  NSTerm& model = _diffVar->getModel();
  const CFreal mu = model.getPhysicalData()[NSTerm::MU];
  vector<CFuint> index(nbNodes);
  
  // Triangle 0 : nodes 0-3-5

  CFuint i1 = 0;
  CFuint i2 = 3;
  CFuint i3 = 5;
  CFreal sign = 1.0;

  // Here we compute the two residual of the diffusive part
  // first one is the integrate of diffusion times linear basis function of the sub-element(corresponding to a galerkin)
  //the second is the integrate of diffusion times a bubble function.
  // This two integrals are done on each sub-elemement.
 
  fluctuation_diff_lingalerkin(i1, i2, i3, sign);

  fluctuation_diff_bubble(i1, i2, i3);


  // call beta of the sub-triangles, betas are always the one of the LDA schemes
  vector<RealMatrix> betasInTriag;
  betasInTriag.resize(nbNodes);
  for (CFuint i = 0 ; i< nbNodes; ++i)
    betasInTriag[i].resize(4,4);
  
  betasInTriag = getMethodData().getDistributionData().betaMats[0];
  
  // The "bubble" part should be distributed using the Beta of LDA, this for consistency
  for ( CFuint iNode = 0; iNode < nbNodes; ++iNode) {
     (*_phi_diff_bub_split[iNode]) = betasInTriag[iNode]*(*_phi_diff_bub);
  }

  // Diffusive resildual is distributed to each node of the element

  for (CFuint iEq = 0; iEq < nbEqs; ++iEq) {
      result[i1][iEq] += (*_phi_diff_gal[0])[iEq] + 3.0*(*_phi_diff_bub_split[0])[iEq] - onethird*(*_phi_diff_bub)[iEq];
      result[i2][iEq] += (*_phi_diff_gal[1])[iEq] + 3.0*(*_phi_diff_bub_split[1])[iEq] - onethird*(*_phi_diff_bub)[iEq];
      result[i3][iEq] += (*_phi_diff_gal[2])[iEq] + 3.0*(*_phi_diff_bub_split[2])[iEq] - onethird*(*_phi_diff_bub)[iEq];
  }
 

  // set the update coeff
  index[0] = i1;
  index[1] = i2;
  index[2] = i3;
  CFreal avRho = ((*_states[i1])[0] + (*_states[i2])[0] + (*_states[i3])[0])/3.0;   
  for (CFuint i = 0; i < nbNodes; ++i) {
    // this is not the unit normal !!
    for (CFuint iDim = 0; iDim < dim; ++iDim) {
      _normal[iDim] = normals[m_cellID]->getNodalNormComp(i,iDim);
      }
    if (updateCoeffFlag) {
      const CFreal faceArea = (normals[m_cellID]->getAreaNode(i))*0.5;
      const CFreal diffCoeff = mu / avRho;
      updateCoeff[geo->getState(index[i])->getLocalID()] += diffCoeff*faceArea*faceArea/(_cellVolume*ovDimCoeff2);
    }
  } 

  // Triangle 2 : nodes 3-1-4

  i1 = 3;
  i2 = 1;
  i3 = 4;

 
  fluctuation_diff_bubble(i1, i2, i3);
  fluctuation_diff_lingalerkin(i1, i2, i3, sign);
  betasInTriag = getMethodData().getDistributionData().betaMats[1];

  for ( CFuint iNode = 0; iNode < nbNodes; ++iNode) {
    (*_phi_diff_bub_split[iNode]) = betasInTriag[iNode]*(*_phi_diff_bub);
  }

  for (CFuint iEq = 0; iEq < nbEqs; ++iEq) {
    result[i1][iEq] += (*_phi_diff_gal[0])[iEq] + 3.0*(*_phi_diff_bub_split[0])[iEq] - onethird*(*_phi_diff_bub)[iEq];
    result[i2][iEq] += (*_phi_diff_gal[1])[iEq] + 3.0*(*_phi_diff_bub_split[1])[iEq] - onethird*(*_phi_diff_bub)[iEq];
    result[i3][iEq] += (*_phi_diff_gal[2])[iEq] + 3.0*(*_phi_diff_bub_split[2])[iEq] -  onethird*(*_phi_diff_bub)[iEq];
  }

  // set the update coeff
  index[0] = i1;  
  index[1] = i2;
  index[2] = i3;
  
  avRho = ((*_states[i1])[0] + (*_states[i2])[0] + (*_states[i3])[0])/3.0; 

  for (CFuint i = 0; i < nbNodes; ++i) {
    // this is not the unit normal !!
    for (CFuint iDim = 0; iDim < dim; ++iDim) {
      _normal[iDim] = normals[m_cellID]->getNodalNormComp(i,iDim);
      }
    if (updateCoeffFlag) {
      const CFreal faceArea = (normals[m_cellID]->getAreaNode(i))*0.5;
      const CFreal diffCoeff = mu / avRho;
      updateCoeff[geo->getState(index[i])->getLocalID()] += diffCoeff*faceArea*faceArea/(_cellVolume*ovDimCoeff2);
    }
  }

  // Triangle 3 : nodes 5-4-2
  i1 = 5;
  i2 = 4;
  i3 = 2;

 
  fluctuation_diff_bubble(i1, i2, i3);
  fluctuation_diff_lingalerkin(i1, i2, i3, sign);

  betasInTriag = getMethodData().getDistributionData().betaMats[2];

  for ( CFuint iNode = 0; iNode < nbNodes; ++iNode) {
    (*_phi_diff_bub_split[iNode]) = betasInTriag[iNode]*(*_phi_diff_bub);
  }


  for (CFuint iEq = 0; iEq < nbEqs; ++iEq) {
    result[i1][iEq] += (*_phi_diff_gal[0])[iEq] + 3.0*(*_phi_diff_bub_split[0])[iEq] - onethird*(*_phi_diff_bub)[iEq];
    result[i2][iEq] += (*_phi_diff_gal[1])[iEq] + 3.0*(*_phi_diff_bub_split[1])[iEq] - onethird*(*_phi_diff_bub)[iEq];
    result[i3][iEq] += (*_phi_diff_gal[2])[iEq] + 3.0*(*_phi_diff_bub_split[2])[iEq] -  onethird*(*_phi_diff_bub)[iEq];
  }

  // set the update coeff
  index[0] = i1;
  index[1] = i2;
  index[2] = i3;

  avRho = ((*_states[i1])[0] + (*_states[i2])[0] + (*_states[i3])[0])/3.0; 

  for (CFuint i = 0; i < nbNodes; ++i) {
    // this is not the unit normal !!
    for (CFuint iDim = 0; iDim < dim; ++iDim) {
      _normal[iDim] = normals[m_cellID]->getNodalNormComp(i,iDim);
      }
    if (updateCoeffFlag) {
      const CFreal faceArea = (normals[m_cellID]->getAreaNode(i))*0.5;
      const CFreal diffCoeff = mu / avRho;
      updateCoeff[geo->getState(index[i])->getLocalID()] += diffCoeff*faceArea*faceArea/(_cellVolume*ovDimCoeff2);
    }
  }
  
  // Triangle 4 : nodes 4-5-3
  i1 = 4;
  i2 = 5;
  i3 = 3;
  sign = -1.0;
 
  fluctuation_diff_bubble(i1, i2, i3);
  fluctuation_diff_lingalerkin(i1, i2, i3, sign);
  
  betasInTriag = getMethodData().getDistributionData().betaMats[3];

  for ( CFuint iNode = 0; iNode < nbNodes; ++iNode) {
     (*_phi_diff_bub_split[iNode]) = betasInTriag[iNode]*(*_phi_diff_bub);
  }

  for (CFuint iEq = 0; iEq < nbEqs; ++iEq) {
    result[i1][iEq] += (*_phi_diff_gal[0])[iEq] + 3.0*(*_phi_diff_bub_split[0])[iEq] - onethird*(*_phi_diff_bub)[iEq];
    result[i2][iEq] += (*_phi_diff_gal[1])[iEq] + 3.0*(*_phi_diff_bub_split[1])[iEq] - onethird*(*_phi_diff_bub)[iEq];
    result[i3][iEq] += (*_phi_diff_gal[2])[iEq] + 3.0*(*_phi_diff_bub_split[2])[iEq] -  onethird*(*_phi_diff_bub)[iEq];
  }

  // set the update coeff
  index[0] = i1;
  index[1] = i2;
  index[2] = i3;
  avRho = ((*_states[i1])[0] + (*_states[i2])[0] + (*_states[i3])[0])/3.0; 

  for (CFuint i = 0; i < nbNodes; ++i) {
    // this is not the unit normal !!
    for (CFuint iDim = 0; iDim < dim; ++iDim) {
      _normal[iDim] = normals[m_cellID]->getNodalNormComp(i,iDim);
    }
    if (updateCoeffFlag) {
      const CFreal faceArea = -(normals[m_cellID]->getAreaNode(i))*0.5;
      const CFreal diffCoeff = mu / avRho;
      updateCoeff[geo->getState(index[i])->getLocalID()] += diffCoeff*faceArea*faceArea/(_cellVolume*ovDimCoeff2);
    }
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

void NavierStokesTermHOlin::fluctuation_diff_lingalerkin(CFuint& i1, CFuint& i2, CFuint& i3, CFreal& sign){
  //This function compute the integrate on the sub-element( of vertex i1, i2, i3)
  //of the diffusion times the linear basis functions of the sub-element

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

  for (CFuint iStates = 0; iStates < nbNodes; ++ iStates ){
    for (CFuint iEq = 0; iEq < nbEqs; ++ iEq)
      (*_phi_diff_gal[iStates])[iEq] = 0.0;
  }
  for (CFuint iQd = 0; iQd < 4; ++iQd) {
     //point of quadrature
     const CFreal x = qd0[iQd] * node0[XX] + qd1[iQd] * node1[XX] + qd2[iQd] * node2[XX] ;
     const CFreal y = qd0[iQd] * node0[YY] + qd1[iQd] * node1[YY] + qd2[iQd] * node2[YY] ;

     // Linear basis function of the element
     CFreal L1 = 1.0 + 0.5*( ( x - x1_vert )*nx1 + ( y - y1_vert )*ny1 )*inv_volume ;
     CFreal L2 = 1.0 + 0.5*( ( x - x2_vert )*nx2 + ( y - y2_vert )*ny2 )*inv_volume ;
     CFreal L3 = 1.0 + 0.5*( ( x - x3_vert )*nx3 + ( y - y3_vert )*ny3 )*inv_volume ;

 
     RealVector grad_state_x;
     grad_state_x.resize(nbEqs);
     RealVector grad_state_y;
     grad_state_y.resize(nbEqs);

     for (CFuint i = 0; i< nbEqs; ++i){
         grad_state_x[i] = (nx1*(4.0*L1 - 1.0)*(*_states[0])[i] +
		            nx2*(4.0*L2 - 1.0)*(*_states[1])[i] +
		            nx3*(4.0*L3 - 1.0)*(*_states[2])[i] +
		            4.0*(nx1*L2 + nx2*L1)*(*_states[3])[i] +
		            4.0*(nx2*L3 + nx3*L2)*(*_states[4])[i] +
		            4.0*(nx3*L1 + nx1*L3)*(*_states[5])[i])*coeffGrad;

        grad_state_y[i] = (ny1*(4.0*L1 - 1.0)*(*_states[0])[i] +
		           ny2*(4.0*L2 - 1.0)*(*_states[1])[i] +
		           ny3*(4.0*L3 - 1.0)*(*_states[2])[i] +
		           4.0*(ny1*L2 + ny2*L1)*(*_states[3])[i] +
		           4.0*(ny2*L3 + ny3*L2)*(*_states[4])[i] +
		           4.0*(ny3*L1 + ny1*L3)*(*_states[5])[i])*coeffGrad;
        }


     CFreal grad_T_x;
     grad_T_x = (nx1*(4.0*L1 - 1.0)*_values(3,0)  +
	         nx2*(4.0*L2 - 1.0)*_values(3,1) +
	         nx3*(4.0*L3 - 1.0)*_values(3,2) +
	         4.0*(nx1*L2 + nx2*L1)*_values(3,3) +
	         4.0*(nx2*L3 + nx3*L2)*_values(3,4) +
	         4.0*(nx3*L1 + nx1*L3)*_values(3,5))*coeffGrad;
     
     CFreal grad_T_y;
     grad_T_y = (ny1*(4.0*L1 - 1.0)*_values(3,0)  +
		 ny2*(4.0*L2 - 1.0)*_values(3,1) +
		 ny3*(4.0*L3 - 1.0)*_values(3,2) +
		 4.0*(ny1*L2 + ny2*L1)*_values(3,3) +
		 4.0*(ny2*L3 + ny3*L2)*_values(3,4) +
		 4.0*(ny3*L1 + ny1*L3)*_values(3,5))*coeffGrad;
     

     (qdstates) = (L1*( 2.0*L1 - 1.0 ) * (*_states[0])) +
       (L2*( 2.0*L2 - 1.0 ) * (*_states[1])) +
       (L3*( 2.0*L3 - 1.0 ) * (*_states[2])) +
       (4.0*L1*L2           * (*_states[3])) +
       (4.0*L3*L2           * (*_states[4])) +
       (4.0*L1*L3           * (*_states[5]));
     
     (*_gradients[0])[XX] = grad_state_x[0];
     (*_gradients[0])[YY] = grad_state_y[0];
     (*_gradients[1])[XX] = -(qdstates[1]/(qdstates[0]*qdstates[0]))*grad_state_x[0] + (1.0/qdstates[0])*grad_state_x[1];
     (*_gradients[1])[YY] = -(qdstates[1]/(qdstates[0]*qdstates[0]))*grad_state_y[0] + (1.0/qdstates[0])*grad_state_y[1];
     (*_gradients[2])[XX] = -(qdstates[2]/(qdstates[0]*qdstates[0]))*grad_state_x[0] + (1.0/qdstates[0])*grad_state_x[2];
     (*_gradients[2])[YY] = -(qdstates[2]/(qdstates[0]*qdstates[0]))*grad_state_y[0] + (1.0/qdstates[0])*grad_state_y[2];
     (*_gradients[3])[XX] = grad_T_x;
     (*_gradients[3])[YY] = grad_T_y;

     _normal[XX] = sign*nx1*0.5;
     _normal[YY] = sign*ny1*0.5;

     RealVector F1;
     F1.resize(nbEqs);
     F1 = _diffVar->getFlux(qdstates, _gradients,_normal , 0.0);

     _normal[XX] = sign*nx2*0.5;
     _normal[YY] = sign*ny2*0.5;

     RealVector F2;
     F2.resize(nbEqs);
     F2 = _diffVar->getFlux(qdstates, _gradients,_normal , 0.0);

     _normal[XX] = sign*nx3*0.5;
     _normal[YY] = sign*ny3*0.5;

     RealVector F3;
     F3.resize(nbEqs);
     F3 = _diffVar->getFlux(qdstates, _gradients,_normal , 0.0);

     (*_phi_diff_gal[0])[0] = 0.0;
     (*_phi_diff_gal[1])[0] = 0.0;
     (*_phi_diff_gal[2])[0] = 0.0;

      for (CFuint i = 1 ; i < nbEqs; ++i)
	{
          (*_phi_diff_gal[0])[i] += F1[i]*0.5*wqd[iQd];
          (*_phi_diff_gal[1])[i] += F2[i]*0.5*wqd[iQd];
          (*_phi_diff_gal[2])[i] += F3[i]*0.5*wqd[iQd];
       }

  }

}

void NavierStokesTermHOlin::fluctuation_diff_bubble(CFuint& i1, CFuint& i2, CFuint& i3)
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
  // Then, we just need one node per constant part. I choose to take one point at the gravity center of the 3 triangkes
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

  for (CFuint i = 0; i<nbEqs; ++i){
     grad_state_x[i] = (nx1*(4.0*L1 - 1.0)*(*_states[0])[i] +
		        nx2*(4.0*L2 - 1.0)*(*_states[1])[i] +
		        nx3*(4.0*L3 - 1.0)*(*_states[2])[i] +
		        4.0*(nx1*L2 + nx2*L1)*(*_states[3])[i] +
		        4.0*(nx2*L3 + nx3*L2)*(*_states[4])[i] +
		        4.0*(nx3*L1 + nx1*L3)*(*_states[5])[i])*coeffGrad;

     grad_state_y[i] = (ny1*(4.0*L1 - 1.0)*(*_states[0])[i] +
		        ny2*(4.0*L2 - 1.0)*(*_states[1])[i] +
		        ny3*(4.0*L3 - 1.0)*(*_states[2])[i] +
		        4.0*(ny1*L2 + ny2*L1)*(*_states[3])[i] +
		        4.0*(ny2*L3 + ny3*L2)*(*_states[4])[i] +
		        4.0*(ny3*L1 + ny1*L3)*(*_states[5])[i])*coeffGrad;
  }



  CFreal grad_T_x;
  CFreal grad_T_y;

  grad_T_x = (nx1*(4.0*L1 - 1.0)*_values(3,0) +
	      nx2*(4.0*L2 - 1.0)*_values(3,1) +
	      nx3*(4.0*L3 - 1.0)*_values(3,2) +
	      4.0*(nx1*L2 + nx2*L1)*_values(3,3) +
	      4.0*(nx2*L3 + nx3*L2)*_values(3,4) +
	      4.0*(nx3*L1 + nx1*L3)*_values(3,5))*coeffGrad;

  grad_T_y = (ny1*(4.0*L1 - 1.0)*_values(3,0) +
	      ny2*(4.0*L2 - 1.0)*_values(3,1) +
	      ny3*(4.0*L3 - 1.0)*_values(3,2) +
	      4.0*(ny1*L2 + ny2*L1)*_values(3,3) +
	      4.0*(ny2*L3 + ny3*L2)*_values(3,4) +
	      4.0*(ny3*L1 + ny1*L3)*_values(3,5))*coeffGrad;


  qdstates = (L1*( 2.0*L1 - 1.0 ) * (*_states[0])) +
             (L2*( 2.0*L2 - 1.0 ) * (*_states[1])) +
             (L3*( 2.0*L3 - 1.0 ) * (*_states[2])) +
             (4.0*L1*L2           * (*_states[3])) +
             (4.0*L3*L2           * (*_states[4])) +
             (4.0*L1*L3           * (*_states[5]));

  (*_gradients[0])[XX] = grad_state_x[0];
  (*_gradients[0])[YY] = grad_state_y[0];
  (*_gradients[1])[XX] = -(qdstates[1]/(qdstates[0]*qdstates[0]))*grad_state_x[0] + (1.0/qdstates[0])*grad_state_x[1];
  (*_gradients[1])[YY] = -(qdstates[1]/(qdstates[0]*qdstates[0]))*grad_state_y[0] + (1.0/qdstates[0])*grad_state_y[1];
  (*_gradients[2])[XX] = -(qdstates[2]/(qdstates[0]*qdstates[0]))*grad_state_x[0] + (1.0/qdstates[0])*grad_state_x[2];
  (*_gradients[2])[YY] = -(qdstates[2]/(qdstates[0]*qdstates[0]))*grad_state_y[0] + (1.0/qdstates[0])*grad_state_y[2];
  (*_gradients[3])[XX] = grad_T_x;
  (*_gradients[3])[YY] = grad_T_y;

  _normal[XX] = (node0[YY] - node2[YY]);
  _normal[YY] = (node2[XX] - node0[XX]);

  RealVector F1;
  F1.resize(nbEqs);
  F1 = _diffVar->getFlux(qdstates, _gradients,_normal , 0.0);




  for (CFuint i = 1 ; i < nbEqs; ++i){
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

  grad_state_x = (nx1*(4.0*L1 - 1.0)*(*_states[0]) +
		  nx2*(4.0*L2 - 1.0)*(*_states[1]) +
		  nx3*(4.0*L3 - 1.0)*(*_states[2]) +
		  4.0*(nx1*L2 + nx2*L1)*(*_states[3]) +
		  4.0*(nx2*L3 + nx3*L2)*(*_states[4]) +
		  4.0*(nx3*L1 + nx1*L3)*(*_states[5]))*coeffGrad;

  grad_state_y = (ny1*(4.0*L1 - 1.0)*(*_states[0]) +
		  ny2*(4.0*L2 - 1.0)*(*_states[1]) +
		  ny3*(4.0*L3 - 1.0)*(*_states[2]) +
		  4.0*(ny1*L2 + ny2*L1)*(*_states[3]) +
		  4.0*(ny2*L3 + ny3*L2)*(*_states[4]) +
		  4.0*(ny3*L1 + ny1*L3)*(*_states[5]))*coeffGrad;

  grad_T_x = (nx1*(4.0*L1 - 1.0)*_values(3,0) +
	      nx2*(4.0*L2 - 1.0)*_values(3,1) +
	      nx3*(4.0*L3 - 1.0)*_values(3,2) +
	      4.0*(nx1*L2 + nx2*L1)*_values(3,3) +
	      4.0*(nx2*L3 + nx3*L2)*_values(3,4) +
	      4.0*(nx3*L1 + nx1*L3)*_values(3,5))*coeffGrad;

  grad_T_y = (ny1*(4.0*L1 - 1.0)*_values(3,0) +
	      ny2*(4.0*L2 - 1.0)*_values(3,1) +
	      ny3*(4.0*L3 - 1.0)*_values(3,2) +
	      4.0*(ny1*L2 + ny2*L1)*_values(3,3) +
	      4.0*(ny2*L3 + ny3*L2)*_values(3,4) +
	      4.0*(ny3*L1 + ny1*L3)*_values(3,5))*coeffGrad;


  qdstates = (L1*( 2.0*L1 - 1.0 ) * (*_states[0])) +
             (L2*( 2.0*L2 - 1.0 ) * (*_states[1])) +
             (L3*( 2.0*L3 - 1.0 ) * (*_states[2])) +
             (4.0*L1*L2           * (*_states[3])) +
             (4.0*L3*L2           * (*_states[4])) +
             (4.0*L1*L3           * (*_states[5]));


  (*_gradients[0])[XX] = grad_state_x[0];
  (*_gradients[0])[YY] = grad_state_y[0];
  (*_gradients[1])[XX] = -(qdstates[1]/(qdstates[0]*qdstates[0]))*grad_state_x[0] + (1.0/qdstates[0])*grad_state_x[1];
  (*_gradients[1])[YY] = -(qdstates[1]/(qdstates[0]*qdstates[0]))*grad_state_y[0] + (1.0/qdstates[0])*grad_state_y[1];
  (*_gradients[2])[XX] = -(qdstates[2]/(qdstates[0]*qdstates[0]))*grad_state_x[0] + (1.0/qdstates[0])*grad_state_x[2];
  (*_gradients[2])[YY] = -(qdstates[2]/(qdstates[0]*qdstates[0]))*grad_state_y[0] + (1.0/qdstates[0])*grad_state_y[2];
  (*_gradients[3])[XX] = grad_T_x;
  (*_gradients[3])[YY] = grad_T_y;

  _normal[XX] = (node1[YY] - node2[YY]);
  _normal[YY] = (node2[XX] - node1[XX]);

  F1 = _diffVar->getFlux(qdstates, _gradients,_normal , 0.0);

  for (CFuint i = 1 ; i < nbEqs; ++i){
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

  grad_state_x = (nx1*(4.0*L1 - 1.0)*(*_states[0]) +
		  nx2*(4.0*L2 - 1.0)*(*_states[1]) +
		  nx3*(4.0*L3 - 1.0)*(*_states[2]) +
		  4.0*(nx1*L2 + nx2*L1)*(*_states[3]) +
		  4.0*(nx2*L3 + nx3*L2)*(*_states[4]) +
		  4.0*(nx3*L1 + nx1*L3)*(*_states[5]))*coeffGrad;

  grad_state_y = (ny1*(4.0*L1 - 1.0)*(*_states[0]) +
		  ny2*(4.0*L2 - 1.0)*(*_states[1]) +
		  ny3*(4.0*L3 - 1.0)*(*_states[2]) +
		  4.0*(ny1*L2 + ny2*L1)*(*_states[3]) +
		  4.0*(ny2*L3 + ny3*L2)*(*_states[4]) +
		  4.0*(ny3*L1 + ny1*L3)*(*_states[5]))*coeffGrad;

  grad_T_x = (nx1*(4.0*L1 - 1.0)*_values(3,0) +
	      nx2*(4.0*L2 - 1.0)*_values(3,1) +
	      nx3*(4.0*L3 - 1.0)*_values(3,2) +
	      4.0*(nx1*L2 + nx2*L1)*_values(3,3) +
	      4.0*(nx2*L3 + nx3*L2)*_values(3,4) +
	      4.0*(nx3*L1 + nx1*L3)*_values(3,5))*coeffGrad;

  grad_T_y = (ny1*(4.0*L1 - 1.0)*_values(3,0) +
	      ny2*(4.0*L2 - 1.0)*_values(3,1) +
	      ny3*(4.0*L3 - 1.0)*_values(3,2) +
	      4.0*(ny1*L2 + ny2*L1)*_values(3,3) +
	      4.0*(ny2*L3 + ny3*L2)*_values(3,4) +
	      4.0*(ny3*L1 + ny1*L3)*_values(3,5))*coeffGrad;


  qdstates = (L1*( 2.0*L1 - 1.0 ) * (*_states[0])) +
             (L2*( 2.0*L2 - 1.0 ) * (*_states[1])) +
             (L3*( 2.0*L3 - 1.0 ) * (*_states[2])) +
             (4.0*L1*L2           * (*_states[3])) +
             (4.0*L3*L2           * (*_states[4])) +
             (4.0*L1*L3           * (*_states[5]));


  (*_gradients[0])[XX] = grad_state_x[0];
  (*_gradients[0])[YY] = grad_state_y[0];
  (*_gradients[1])[XX] = -(qdstates[1]/(qdstates[0]*qdstates[0]))*grad_state_x[0] + (1.0/qdstates[0])*grad_state_x[1];
  (*_gradients[1])[YY] = -(qdstates[1]/(qdstates[0]*qdstates[0]))*grad_state_y[0] + (1.0/qdstates[0])*grad_state_y[1];
  (*_gradients[2])[XX] = -(qdstates[2]/(qdstates[0]*qdstates[0]))*grad_state_x[0] + (1.0/qdstates[0])*grad_state_x[2];
  (*_gradients[2])[YY] = -(qdstates[2]/(qdstates[0]*qdstates[0]))*grad_state_y[0] + (1.0/qdstates[0])*grad_state_y[2];
  (*_gradients[3])[XX] = grad_T_x;
  (*_gradients[3])[YY] = grad_T_y;

  _normal[XX] = (node0[YY] - node1[YY]);
  _normal[YY] = (node1[XX] - node0[XX]);

  F1 = _diffVar->getFlux(qdstates, _gradients, _normal , 0.0);

  for (CFuint i = 1 ; i < nbEqs; ++i){
      (*_phi_diff_bub)[i] += F1[i]*_cellVolume/
	                     (12.0*((node1[XX]*yg - xg*node1[YY]) - node0[XX]*(yg - node1[YY]) + node0[YY]*(xg - node1[XX])));

      }

}
//////////////////////////////////////////////////////////////////////////////

void NavierStokesTermHOlin::setup()
{
 const CFuint nbStatesInControlVolume =
    MeshDataStack::getActive()->Statistics().getMaxNbStatesInCell();
  const CFuint nbNodesInControlVolume =
    MeshDataStack::getActive()->Statistics().getMaxNbNodesInCell();
  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  
  _states.resize(nbStatesInControlVolume);
  _values.resize(nbEqs, nbStatesInControlVolume);
  
  _gradients.resize(nbEqs);
  for (CFuint i = 0; i< nbEqs; ++i) {
    _gradients[i] = new RealVector(PhysicalModelStack::getActive()->getDim());
  }
  
  _avValues.resize(nbEqs);
  _normal.resize(PhysicalModelStack::getActive()->getDim());
  
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

  _phi_diff_gal.resize(nbNodesInControlVolume); // one galerkin function per state

  for (CFuint iState = 0; iState < 3; ++ iState)
     _phi_diff_gal[iState] = new RealVector(nbEqs);

 

}

//////////////////////////////////////////////////////////////////////////////

void NavierStokesTermHOlin::configure ( Config::ConfigArgs& args )
{
  ComputeDiffusiveTerm::configure(args);
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
