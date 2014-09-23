#include "Environment/ObjectProvider.hh"

#include "Framework/GeometricEntity.hh"
#include "Framework/MeshData.hh"
#include "Framework/MethodStrategyProvider.hh"

#include "NavierStokes/NavierStokesVarSet.hh"
#include "NavierStokes/EulerVarSet.hh"

#include "FluctSplit/HONavierStokes/FluctSplitHONavierStokes.hh"
#include "FluctSplit/FluctuationSplitStrategy.hh"
#include "FluctSplit/HONavierStokes/NavierStokesTermHO3D.hh"
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

MethodStrategyProvider<NavierStokesTermHO3D,
		       FluctuationSplitData,
		       ComputeDiffusiveTerm,
		       FluctSplitHONavierStokesModule>
navierStokesDiffusiveTermHO3DProvider("NavierStokesHO3D");

//////////////////////////////////////////////////////////////////////////////

NavierStokesTermHO3D::NavierStokesTermHO3D(const std::string& name) :
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

NavierStokesTermHO3D::~NavierStokesTermHO3D()
{
  for (CFuint i = 0; i< _gradients.size(); ++i) {
    deletePtr(_gradients[i]);
  }
}

//////////////////////////////////////////////////////////////////////////////

void NavierStokesTermHO3D::setDiffusiveVarSet(SafePtr<DiffusiveVarSet> diffVar)
{
  _diffVar = diffVar.d_castTo<NavierStokesVarSet>();
}

//////////////////////////////////////////////////////////////////////////////

void NavierStokesTermHO3D::setUpdateVarSet(SafePtr<ConvectiveVarSet> updateVar)
{
  _updateVar = updateVar.d_castTo<EulerVarSet>();
}


//////////////////////////////////////////////////////////////////////////////

void NavierStokesTermHO3D::computeDiffusiveTerm
(GeometricEntity *const geo, vector<RealVector>& result, bool updateCoeffFlag)
{
  DataHandle< InwardNormalsData*> normals = socket_normals.getDataHandle();
  DataHandle< CFreal> updateCoeff = socket_updateCoeff.getDataHandle();
  DistributionData& distdata = getMethodData().getDistributionData();
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

  NSTerm& model = _diffVar->getModel();
const CFreal mu = model.getPhysicalData()[NSTerm::MU];
// unused // const CFreal ovDimCoeff = 1./(dimCoeff);

 CFreal betacoef = (105.0/32.0);

 FacesMat = distdata.FacesMat_tet;
  // Tetra 0 : nodes 1-5-4-7
// CF_DEBUG_POINT;
  CFuint i1 = 1;
  CFuint i2 = 5;
  CFuint i3 = 4;
  CFuint i4 = 7;


  matrix_node_norms (0,XX) = FacesMat(16,XX);
  matrix_node_norms (0,YY) = FacesMat(16,YY);
  matrix_node_norms (0,ZZ) = FacesMat(16,ZZ);
  
  matrix_node_norms (1,XX) = -FacesMat(1,XX);
  matrix_node_norms (1,YY) = -FacesMat(1,YY);
  matrix_node_norms (1,ZZ) = -FacesMat(1,ZZ);

  matrix_node_norms (2,XX) = -FacesMat(14,XX);
  matrix_node_norms (2,YY) = -FacesMat(14,YY);
  matrix_node_norms (2,ZZ) = -FacesMat(14,ZZ);

  matrix_node_norms (3,XX) = -FacesMat(10,XX);
  matrix_node_norms (3,YY) = -FacesMat(10,YY);
  matrix_node_norms (3,ZZ) = -FacesMat(10,ZZ);
    

// CF_DEBUG_POINT;
  // Here we compute the two residual of the diffusive part
  // first one is the integrate of diffusion times basis function (corresponding to a galerkin)
  //the second is the integrate of diffusion times a bubble function.
  // This two integrals are done on each sub-elemement.
  // @todo : possible to compute fluct_diff_galerkin only one time for the whole triangle

  fluctuation_diff_galerkin(i1, i2, i3,i4, false);
  fluctuation_diff_bubble(i1, i2, i3, i4, false);


  // call beta of the sub-triangles, betas are always the one of the LDA schemes
  vector<RealMatrix> betasInTriag;
  betasInTriag.resize(4);
  for (CFuint i = 0 ; i< 4; ++i)
    betasInTriag[i].resize(5,5);
  betasInTriag = getMethodData().getDistributionData().betaMats[0];
  
  // The "bubble" part should be distributed using the Beta of LDA, this for consistency
  for ( CFuint iNode = 0; iNode < 4; ++iNode) {
     (*_phi_diff_bub_split[iNode]) = betasInTriag[iNode]*(*_phi_diff_bub);
  }

  // Diffusive resildual is distributed to each node of the element

  for (CFuint iEq = 0; iEq < nbEqs; ++iEq) {
      result[0][iEq] += (*_phi_diff_gal[0])[iEq] + kappa(0,0)*(*_phi_diff_bub)[iEq];
      result[1][iEq] += (*_phi_diff_gal[1])[iEq] + betacoef*(*_phi_diff_bub_split[0])[iEq] + kappa(0,1)*(*_phi_diff_bub)[iEq];
      result[2][iEq] += (*_phi_diff_gal[2])[iEq] + kappa(0,2)*(*_phi_diff_bub)[iEq];
      result[3][iEq] += (*_phi_diff_gal[3])[iEq] + kappa(0,3)*(*_phi_diff_bub)[iEq];
      result[4][iEq] += (*_phi_diff_gal[4])[iEq] + betacoef*(*_phi_diff_bub_split[2])[iEq] + kappa(0,4)*(*_phi_diff_bub)[iEq];
      result[5][iEq] += (*_phi_diff_gal[5])[iEq] + betacoef*(*_phi_diff_bub_split[1])[iEq] + kappa(0,5)*(*_phi_diff_bub)[iEq];
      result[6][iEq] += (*_phi_diff_gal[6])[iEq] + kappa(0,6)*(*_phi_diff_bub)[iEq];
      result[7][iEq] += (*_phi_diff_gal[7])[iEq] + betacoef*(*_phi_diff_bub_split[3])[iEq] + kappa(0,7)*(*_phi_diff_bub)[iEq];
      result[8][iEq] += (*_phi_diff_gal[8])[iEq] + kappa(0,8)*(*_phi_diff_bub)[iEq];
      result[9][iEq] += (*_phi_diff_gal[9])[iEq] + kappa(0,9)*(*_phi_diff_bub)[iEq];

  }


  // Tetra 1 : nodes 5-2-6-8

  i1 = 5;
  i2 = 2;
  i3 = 6;
  i4 = 8;

  matrix_node_norms (0,XX) = -FacesMat(6,XX);
  matrix_node_norms (0,YY) = -FacesMat(6,YY);
  matrix_node_norms (0,ZZ) = -FacesMat(6,ZZ);
  
  matrix_node_norms (1,XX) = FacesMat(18,XX);
  matrix_node_norms (1,YY) = FacesMat(18,YY);
  matrix_node_norms (1,ZZ) = FacesMat(18,ZZ);
  
  matrix_node_norms (2,XX) = -FacesMat(13,XX);
  matrix_node_norms (2,YY) = -FacesMat(13,YY);
  matrix_node_norms (2,ZZ) = -FacesMat(13,ZZ);

  matrix_node_norms (3,XX) = -FacesMat(9,XX);
  matrix_node_norms (3,YY) = -FacesMat(9,YY);
  matrix_node_norms (3,ZZ) = -FacesMat(9,ZZ);


  fluctuation_diff_galerkin(i1, i2, i3,i4, false);
  fluctuation_diff_bubble(i1, i2, i3,i4, false);

  betasInTriag = getMethodData().getDistributionData().betaMats[1];

  for ( CFuint iNode = 0; iNode < 4; ++iNode) {
     (*_phi_diff_bub_split[iNode]) = betasInTriag[iNode]*(*_phi_diff_bub);
  }

  for (CFuint iEq = 0; iEq < nbEqs; ++iEq) {
     result[0][iEq] += (*_phi_diff_gal[0])[iEq] + kappa(1,0)*(*_phi_diff_bub)[iEq];
     result[1][iEq] += (*_phi_diff_gal[1])[iEq] + kappa(1,1)*(*_phi_diff_bub)[iEq];
     result[2][iEq] += (*_phi_diff_gal[2])[iEq] + betacoef*(*_phi_diff_bub_split[1])[iEq] + kappa(1,2)*(*_phi_diff_bub)[iEq];
     result[3][iEq] += (*_phi_diff_gal[3])[iEq] + kappa(1,3)*(*_phi_diff_bub)[iEq];
     result[4][iEq] += (*_phi_diff_gal[4])[iEq] + kappa(1,4)*(*_phi_diff_bub)[iEq];
     result[5][iEq] += (*_phi_diff_gal[5])[iEq] + betacoef*(*_phi_diff_bub_split[0])[iEq] + kappa(1,5)*(*_phi_diff_bub)[iEq];
     result[6][iEq] += (*_phi_diff_gal[6])[iEq] + betacoef*(*_phi_diff_bub_split[2])[iEq] + kappa(1,6)*(*_phi_diff_bub)[iEq];
     result[7][iEq] += (*_phi_diff_gal[7])[iEq] + kappa(1,7)*(*_phi_diff_bub)[iEq];
     result[8][iEq] += (*_phi_diff_gal[8])[iEq] + betacoef*(*_phi_diff_bub_split[3])[iEq] + kappa(1,8)*(*_phi_diff_bub)[iEq];
     result[9][iEq] += (*_phi_diff_gal[9])[iEq] + kappa(1,9)*(*_phi_diff_bub)[iEq];
  }


  // Tetra 2 : nodes 4-5-0-7
  i1 = 4;
  i2 = 5;
  i3 = 0;
  i4 = 7;


  matrix_node_norms (0,XX) = -FacesMat(19,XX);
  matrix_node_norms (0,YY) = -FacesMat(19,YY);
  matrix_node_norms (0,ZZ) = -FacesMat(19,ZZ);

  matrix_node_norms (1,XX) = -FacesMat(0,XX);
  matrix_node_norms (1,YY) = -FacesMat(0,YY);
  matrix_node_norms (1,ZZ) = -FacesMat(0,ZZ);

  matrix_node_norms (2,XX) = -FacesMat(16,XX);
  matrix_node_norms (2,YY) = -FacesMat(16,YY);
  matrix_node_norms (2,ZZ) = -FacesMat(16,ZZ);

  matrix_node_norms (3,XX) = -FacesMat(11,XX);
  matrix_node_norms (3,YY) = -FacesMat(11,YY);
  matrix_node_norms (3,ZZ) = -FacesMat(11,ZZ);


  fluctuation_diff_galerkin(i1, i2, i3, i4, false);
  fluctuation_diff_bubble(i1, i2, i3, i4, false);


  betasInTriag = getMethodData().getDistributionData().betaMats[2];

  for ( CFuint iNode = 0; iNode < 4; ++iNode) {
     (*_phi_diff_bub_split[iNode]) = betasInTriag[iNode]*(*_phi_diff_bub);
  }


  for (CFuint iEq = 0; iEq < nbEqs; ++iEq) {
    result[0][iEq] += (*_phi_diff_gal[0])[iEq] + betacoef*(*_phi_diff_bub_split[2])[iEq] + kappa(2,0)*(*_phi_diff_bub)[iEq];
    result[1][iEq] += (*_phi_diff_gal[1])[iEq] + kappa(2,1)*(*_phi_diff_bub)[iEq];
    result[2][iEq] += (*_phi_diff_gal[2])[iEq] + kappa(2,2)*(*_phi_diff_bub)[iEq];
    result[3][iEq] += (*_phi_diff_gal[3])[iEq] + kappa(2,3)*(*_phi_diff_bub)[iEq];
    result[4][iEq] += (*_phi_diff_gal[4])[iEq] + betacoef*(*_phi_diff_bub_split[0])[iEq] + kappa(2,4)*(*_phi_diff_bub)[iEq];
    result[5][iEq] += (*_phi_diff_gal[5])[iEq] + betacoef*(*_phi_diff_bub_split[1])[iEq] + kappa(2,5)*(*_phi_diff_bub)[iEq];
    result[6][iEq] += (*_phi_diff_gal[6])[iEq] + kappa(2,6)*(*_phi_diff_bub)[iEq];
    result[7][iEq] += (*_phi_diff_gal[7])[iEq] + betacoef*(*_phi_diff_bub_split[3])[iEq] + kappa(2,7)*(*_phi_diff_bub)[iEq];
    result[8][iEq] += (*_phi_diff_gal[8])[iEq] + kappa(2,8)*(*_phi_diff_bub)[iEq];
    result[9][iEq] += (*_phi_diff_gal[9])[iEq] + kappa(2,9)*(*_phi_diff_bub)[iEq];
  }

  // Tetra 3 : nodes 0-5-6-8
  i1 = 0;
  i2 = 5;
  i3 = 6;
  i4 = 8;

  matrix_node_norms (0,XX) = -FacesMat(18,XX);
  matrix_node_norms (0,YY) = -FacesMat(18,YY);
  matrix_node_norms (0,ZZ) = -FacesMat(18,ZZ);

  matrix_node_norms (1,XX) = -FacesMat(7,XX);
  matrix_node_norms (1,YY) = -FacesMat(7,YY);
  matrix_node_norms (1,ZZ) = -FacesMat(7,ZZ);

  matrix_node_norms (2,XX) = -FacesMat(21,XX);
  matrix_node_norms (2,YY) = -FacesMat(21,YY);
  matrix_node_norms (2,ZZ) = -FacesMat(21,ZZ);

  matrix_node_norms (3,XX) = -FacesMat(8,XX);
  matrix_node_norms (3,YY) = -FacesMat(8,YY);
  matrix_node_norms (3,ZZ) = -FacesMat(8,ZZ);


  fluctuation_diff_galerkin(i1, i2, i3, i4, false);
  fluctuation_diff_bubble(i1, i2, i3, i4, false);

  betasInTriag = getMethodData().getDistributionData().betaMats[3];

  for ( CFuint iNode = 0; iNode < 4; ++iNode) {
     (*_phi_diff_bub_split[iNode]) = betasInTriag[iNode]*(*_phi_diff_bub);
  }

  for (CFuint iEq = 0; iEq < nbEqs; ++iEq) {
    result[0][iEq] += (*_phi_diff_gal[0])[iEq] + betacoef*(*_phi_diff_bub_split[0])[iEq] + kappa(3,0)*(*_phi_diff_bub)[iEq];
    result[1][iEq] += (*_phi_diff_gal[1])[iEq] + kappa(3,1)*(*_phi_diff_bub)[iEq];
    result[2][iEq] += (*_phi_diff_gal[2])[iEq] + kappa(3,2)*(*_phi_diff_bub)[iEq];
    result[3][iEq] += (*_phi_diff_gal[3])[iEq] + kappa(3,3)*(*_phi_diff_bub)[iEq];
    result[4][iEq] += (*_phi_diff_gal[4])[iEq] + kappa(3,4)*(*_phi_diff_bub)[iEq];
    result[5][iEq] += (*_phi_diff_gal[5])[iEq] + betacoef*(*_phi_diff_bub_split[1])[iEq] + kappa(3,5)*(*_phi_diff_bub)[iEq];
    result[6][iEq] += (*_phi_diff_gal[6])[iEq] + betacoef*(*_phi_diff_bub_split[2])[iEq] + kappa(3,6)*(*_phi_diff_bub)[iEq];
    result[7][iEq] += (*_phi_diff_gal[7])[iEq] + kappa(3,7)*(*_phi_diff_bub)[iEq];
    result[8][iEq] += (*_phi_diff_gal[8])[iEq] + betacoef*(*_phi_diff_bub_split[3])[iEq] + kappa(3,8)*(*_phi_diff_bub)[iEq];
    result[9][iEq] += (*_phi_diff_gal[9])[iEq] + kappa(3,9)*(*_phi_diff_bub)[iEq];
  }


// Tetra 4 : nodes 5-7-8-0
  i1 = 5;
  i2 = 7;
  i3 = 8;
  i4 = 0;


  matrix_node_norms (0,XX) =  FacesMat(20,XX);
  matrix_node_norms (0,YY) =  FacesMat(20,YY);
  matrix_node_norms (0,ZZ) =  FacesMat(20,ZZ);

  matrix_node_norms (1,XX) =  FacesMat(21,XX);
  matrix_node_norms (1,YY) =  FacesMat(21,YY);
  matrix_node_norms (1,ZZ) =  FacesMat(21,ZZ);

  matrix_node_norms (2,XX) =  FacesMat(19,XX);
  matrix_node_norms (2,YY) =  FacesMat(19,YY);
  matrix_node_norms (2,ZZ) =  FacesMat(19,ZZ);

  matrix_node_norms (3,XX) = -FacesMat(15,XX);
  matrix_node_norms (3,YY) = -FacesMat(15,YY);
  matrix_node_norms (3,ZZ) = -FacesMat(15,ZZ);

  fluctuation_diff_galerkin(i1, i2, i3, i4, true);
  fluctuation_diff_bubble(i1, i2, i3, i4,true);

  betasInTriag = getMethodData().getDistributionData().betaMats[4];

  for ( CFuint iNode = 0; iNode < 4; ++iNode) {
     (*_phi_diff_bub_split[iNode]) = betasInTriag[iNode]*(*_phi_diff_bub);
  }

  for (CFuint iEq = 0; iEq < nbEqs; ++iEq) {
    result[0][iEq] += (*_phi_diff_gal[0])[iEq] + betacoef*(*_phi_diff_bub_split[3])[iEq] + kappa(4,0)*(*_phi_diff_bub)[iEq];
    result[1][iEq] += (*_phi_diff_gal[1])[iEq] + kappa(4,1)*(*_phi_diff_bub)[iEq];
    result[2][iEq] += (*_phi_diff_gal[2])[iEq] + kappa(4,2)*(*_phi_diff_bub)[iEq];
    result[3][iEq] += (*_phi_diff_gal[3])[iEq] + kappa(4,3)*(*_phi_diff_bub)[iEq];
    result[4][iEq] += (*_phi_diff_gal[4])[iEq] + kappa(4,4)*(*_phi_diff_bub)[iEq];
    result[5][iEq] += (*_phi_diff_gal[5])[iEq] + betacoef*(*_phi_diff_bub_split[0])[iEq] + kappa(4,5)*(*_phi_diff_bub)[iEq];
    result[6][iEq] += (*_phi_diff_gal[6])[iEq] +  kappa(4,6)*(*_phi_diff_bub)[iEq];
    result[7][iEq] += (*_phi_diff_gal[7])[iEq] + betacoef*(*_phi_diff_bub_split[1])[iEq] + kappa(4,7)*(*_phi_diff_bub)[iEq];
    result[8][iEq] += (*_phi_diff_gal[8])[iEq] + betacoef*(*_phi_diff_bub_split[2])[iEq] + kappa(4,8)*(*_phi_diff_bub)[iEq];
    result[9][iEq] += (*_phi_diff_gal[9])[iEq] + kappa(4,9)*(*_phi_diff_bub)[iEq];
  }





 // Tetra 5 : nodes 7-8-9-3
  i1 = 7;
  i2 = 8;
  i3 = 9;
  i4 = 3;

  matrix_node_norms (0,XX) = -FacesMat(5,XX);
  matrix_node_norms (0,YY) = -FacesMat(5,YY);
  matrix_node_norms (0,ZZ) = -FacesMat(5,ZZ);
  
  matrix_node_norms (1,XX) = -FacesMat(2,XX);
  matrix_node_norms (1,YY) = -FacesMat(2,YY);
  matrix_node_norms (1,ZZ) = -FacesMat(2,ZZ);
  
  matrix_node_norms (2,XX) = -FacesMat(12,XX);
  matrix_node_norms (2,YY) = -FacesMat(12,YY);
  matrix_node_norms (2,ZZ) = -FacesMat(12,ZZ);
  
  matrix_node_norms (3,XX) = FacesMat(17,XX);
  matrix_node_norms (3,YY) = FacesMat(17,YY);
  matrix_node_norms (3,ZZ) = FacesMat(17,ZZ);


  fluctuation_diff_galerkin(i1, i2, i3, i4, false);
  fluctuation_diff_bubble(i1, i2, i3, i4,false);

  betasInTriag = getMethodData().getDistributionData().betaMats[5];

  for ( CFuint iNode = 0; iNode < 4; ++iNode) {
     (*_phi_diff_bub_split[iNode]) = betasInTriag[iNode]*(*_phi_diff_bub);
  }

  for (CFuint iEq = 0; iEq < nbEqs; ++iEq) {
    result[0][iEq] += (*_phi_diff_gal[0])[iEq] + kappa(5,0)*(*_phi_diff_bub)[iEq];
    result[1][iEq] += (*_phi_diff_gal[1])[iEq] + kappa(5,1)*(*_phi_diff_bub)[iEq];
    result[2][iEq] += (*_phi_diff_gal[2])[iEq] + kappa(5,2)*(*_phi_diff_bub)[iEq];
    result[3][iEq] += (*_phi_diff_gal[3])[iEq] + betacoef*(*_phi_diff_bub_split[3])[iEq] + kappa(5,3)*(*_phi_diff_bub)[iEq];
    result[4][iEq] += (*_phi_diff_gal[4])[iEq] + kappa(5,4)*(*_phi_diff_bub)[iEq];
    result[5][iEq] += (*_phi_diff_gal[5])[iEq] + kappa(5,5)*(*_phi_diff_bub)[iEq];
    result[6][iEq] += (*_phi_diff_gal[6])[iEq] + kappa(5,6)*(*_phi_diff_bub)[iEq];
    result[7][iEq] += (*_phi_diff_gal[7])[iEq] + betacoef*(*_phi_diff_bub_split[0])[iEq] + kappa(5,7)*(*_phi_diff_bub)[iEq];
    result[8][iEq] += (*_phi_diff_gal[8])[iEq] + betacoef*(*_phi_diff_bub_split[1])[iEq] + kappa(5,8)*(*_phi_diff_bub)[iEq];
    result[9][iEq] += (*_phi_diff_gal[9])[iEq] + betacoef*(*_phi_diff_bub_split[2])[iEq] + kappa(5,9)*(*_phi_diff_bub)[iEq];
  }

 // Tetra 6 : nodes 0-8-9-7
  i1 = 0;
  i2 = 8;
  i3 = 9;
  i4 = 7;


  matrix_node_norms (0,XX) = -FacesMat(17,XX);
  matrix_node_norms (0,YY) = -FacesMat(17,YY);
  matrix_node_norms (0,ZZ) = -FacesMat(17,ZZ);

  matrix_node_norms (1,XX) = -FacesMat(3,XX);
  matrix_node_norms (1,YY) = -FacesMat(3,YY);
  matrix_node_norms (1,ZZ) = -FacesMat(3,ZZ);

  matrix_node_norms (2,XX) = -FacesMat(20,XX);
  matrix_node_norms (2,YY) = -FacesMat(20,YY);
  matrix_node_norms (2,ZZ) = -FacesMat(20,ZZ);

  matrix_node_norms (3,XX) = -FacesMat(4,XX);
  matrix_node_norms (3,YY) = -FacesMat(4,YY);
  matrix_node_norms (3,ZZ) = -FacesMat(4,ZZ);

  fluctuation_diff_galerkin(i1, i2, i3, i4, false);
  fluctuation_diff_bubble(i1, i2, i3, i4, false);

  betasInTriag = getMethodData().getDistributionData().betaMats[6];

  for ( CFuint iNode = 0; iNode < 4; ++iNode) {
     (*_phi_diff_bub_split[iNode]) = betasInTriag[iNode]*(*_phi_diff_bub);
  }

  for (CFuint iEq = 0; iEq < nbEqs; ++iEq) {
    result[0][iEq] += (*_phi_diff_gal[0])[iEq] + betacoef*(*_phi_diff_bub_split[0])[iEq] + kappa(6,0)*(*_phi_diff_bub)[iEq];
    result[1][iEq] += (*_phi_diff_gal[1])[iEq] + kappa(6,1)*(*_phi_diff_bub)[iEq];
    result[2][iEq] += (*_phi_diff_gal[2])[iEq] + kappa(6,2)*(*_phi_diff_bub)[iEq];
    result[3][iEq] += (*_phi_diff_gal[3])[iEq] + kappa(6,3)*(*_phi_diff_bub)[iEq];
    result[4][iEq] += (*_phi_diff_gal[4])[iEq] + kappa(6,4)*(*_phi_diff_bub)[iEq];
    result[5][iEq] += (*_phi_diff_gal[5])[iEq] + kappa(6,5)*(*_phi_diff_bub)[iEq];
    result[6][iEq] += (*_phi_diff_gal[6])[iEq] + kappa(6,6)*(*_phi_diff_bub)[iEq];
    result[7][iEq] += (*_phi_diff_gal[7])[iEq] + betacoef*(*_phi_diff_bub_split[3])[iEq] + kappa(6,7)*(*_phi_diff_bub)[iEq];
    result[8][iEq] += (*_phi_diff_gal[8])[iEq] + betacoef*(*_phi_diff_bub_split[1])[iEq] + kappa(6,8)*(*_phi_diff_bub)[iEq];
    result[9][iEq] += (*_phi_diff_gal[9])[iEq] + betacoef*(*_phi_diff_bub_split[2])[iEq] + kappa(6,9)*(*_phi_diff_bub)[iEq];
  }


    //update coefficient
    CFreal avRho = ((*_states[0])[0] + (*_states[1])[0] + (*_states[2])[0])/3.0;   
    const CFreal diffCoeff = mu / avRho;

    CFreal coeff = diffCoeff/(_cellVolume*36.0);
    if (updateCoeffFlag) {
      for (CFuint i = 0 ; i < 4; ++i) { 
        const CFreal faceArea = (normals[m_cellID]->getAreaNode(i));

        updateCoeff[geo->getState(i)->getLocalID()] += coeff*faceArea*faceArea;
      } 

      // coeff = 2.0*diffCoeff/(3.0*_cellVolume);
      // CFreal dot_prod = (normals[m_cellID]->getNodalNormComp(0,XX))*(normals[m_cellID]->getNodalNormComp(1,XX)) +
      //                   (normals[m_cellID]->getNodalNormComp(0,YY))*(normals[m_cellID]->getNodalNormComp(1,YY));
  
      // CFreal faceArea1 = (normals[m_cellID]->getAreaNode(0));
      // CFreal faceArea2 = (normals[m_cellID]->getAreaNode(1));
      // updateCoeff[geo->getState(3)->getLocalID()] += coeff*(faceArea1*faceArea1 + faceArea2*faceArea2 + dot_prod);
      //                                                ;

      // faceArea1 = (normals[m_cellID]->getAreaNode(1));
      // faceArea2 = (normals[m_cellID]->getAreaNode(2));
      // dot_prod = (normals[m_cellID]->getNodalNormComp(1,XX))*(normals[m_cellID]->getNodalNormComp(2,XX)) +
      //            (normals[m_cellID]->getNodalNormComp(1,YY))*(normals[m_cellID]->getNodalNormComp(2,YY));
      // updateCoeff[geo->getState(4)->getLocalID()] += coeff*(faceArea1*faceArea1 + faceArea2*faceArea2 + dot_prod);

      // faceArea1 = (normals[m_cellID]->getAreaNode(2));
      // faceArea2 = (normals[m_cellID]->getAreaNode(0));
      // dot_prod = (normals[m_cellID]->getNodalNormComp(0,XX))*(normals[m_cellID]->getNodalNormComp(2,XX)) +
      //            (normals[m_cellID]->getNodalNormComp(0,YY))*(normals[m_cellID]->getNodalNormComp(2,YY));
      // updateCoeff[geo->getState(5)->getLocalID()] += coeff*(faceArea1*faceArea1 + faceArea2*faceArea2 + dot_prod);
}

  for (CFuint iEq = 0; iEq < nbEqs; ++iEq) {
     result[0][iEq] = -result[0][iEq];
     result[1][iEq] = -result[1][iEq];
     result[2][iEq] = -result[2][iEq];
     result[3][iEq] = -result[3][iEq];
     result[4][iEq] = -result[4][iEq];
     result[5][iEq] = -result[5][iEq];
     result[6][iEq] = -result[6][iEq];
     result[7][iEq] = -result[7][iEq];
     result[8][iEq] = -result[8][iEq];
     result[9][iEq] = -result[9][iEq];
	




  }


}

//////////////////////////////////////////////////////////////////////////////
void NavierStokesTermHO3D::fluctuation_diff_galerkin(CFuint& i1, CFuint& i2 , CFuint& i3, CFuint & i4, bool triag){


  // Coordinate of the vertex of the element
  const Node& node0_vert = (*m_cellStates)[0]->getCoordinates();
  const Node& node1_vert = (*m_cellStates)[1]->getCoordinates();
  const Node& node2_vert = (*m_cellStates)[2]->getCoordinates();
  const Node& node3_vert = (*m_cellStates)[3]->getCoordinates();


  const CFreal x0_vert = node0_vert[XX];
  const CFreal x1_vert = node1_vert[XX];
  const CFreal x2_vert = node2_vert[XX];
  const CFreal x3_vert = node3_vert[XX];


  const CFreal y0_vert = node0_vert[YY];
  const CFreal y1_vert = node1_vert[YY];
  const CFreal y2_vert = node2_vert[YY];
  const CFreal y3_vert = node3_vert[YY];

  const CFreal z0_vert = node0_vert[ZZ];
  const CFreal z1_vert = node1_vert[ZZ];
  const CFreal z2_vert = node2_vert[ZZ];
  const CFreal z3_vert = node3_vert[ZZ];



  //Coordinates of the node defining the surface where the integral is computed
  //  This correspond to vertex of the sub-element on which we integrate
  Node& node0 = (*m_cellStates)[i1]->getCoordinates();
  Node& node1 = (*m_cellStates)[i2]->getCoordinates();
  Node& node2 = (*m_cellStates)[i3]->getCoordinates();
  Node& node3 = (*m_cellStates)[i4]->getCoordinates();

  const CFuint nbStates = _states.size();
  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();

  InwardNormalsData& cellnormals = (*socket_normals.getDataHandle()[m_cellID]);

  const CFreal nx0 = cellnormals.getNodalNormComp(0,XX);
  const CFreal nx1 = cellnormals.getNodalNormComp(1,XX);
  const CFreal nx2 = cellnormals.getNodalNormComp(2,XX);
  const CFreal nx3 = cellnormals.getNodalNormComp(3,XX);



  const CFreal ny0 = cellnormals.getNodalNormComp(0,YY);
  const CFreal ny1 = cellnormals.getNodalNormComp(1,YY);
  const CFreal ny2 = cellnormals.getNodalNormComp(2,YY);
  const CFreal ny3 = cellnormals.getNodalNormComp(3,YY);


  const CFreal nz0 = cellnormals.getNodalNormComp(0,ZZ);
  const CFreal nz1 = cellnormals.getNodalNormComp(1,ZZ);
  const CFreal nz2 = cellnormals.getNodalNormComp(2,ZZ);
  const CFreal nz3 = cellnormals.getNodalNormComp(3,ZZ);


  const CFuint dim = PhysicalModelStack::getActive()->getDim();
  const CFreal inv_sixvolume = 1.0/(6.0*_cellVolume);
 // const CFreal one_eighth = 1.0/8.0;


  for (CFuint iStates = 0; iStates < nbStates; ++ iStates ){
     for (CFuint iEq = 0; iEq < nbEqs; ++ iEq)
           (*_phi_diff_gal[iStates])[iEq] = 0.0;
  }

  for (CFuint iQd = 0; iQd < 4; ++iQd) {
     //point od quadrature
     const CFreal x = qd0[iQd] * node0[XX] + qd1[iQd] * node1[XX] + qd2[iQd] * node2[XX] + qd3[iQd]*node3[XX];
     const CFreal y = qd0[iQd] * node0[YY] + qd1[iQd] * node1[YY] + qd2[iQd] * node2[YY] + qd3[iQd]*node3[YY];
     const CFreal z = qd0[iQd] * node0[ZZ] + qd1[iQd] * node1[ZZ] + qd2[iQd] * node2[ZZ] + qd3[iQd]*node3[ZZ];



     // Linear basis function
     CFreal L0 = 1.0 + ( ( x - x0_vert )*nx0 + ( y - y0_vert )*ny0 + (z - z0_vert)*nz0 )*inv_sixvolume ;
     CFreal L1 = 1.0 + ( ( x - x1_vert )*nx1 + ( y - y1_vert )*ny1 + (z - z1_vert)*nz1 )*inv_sixvolume ;
     CFreal L2 = 1.0 + ( ( x - x2_vert )*nx2 + ( y - y2_vert )*ny2 + (z - z2_vert)*nz2 )*inv_sixvolume ;
     CFreal L3 = 1.0 + ( ( x - x3_vert )*nx3 + ( y - y3_vert )*ny3 + (z - z3_vert)*nz3 )*inv_sixvolume ;




     RealVector grad_state_x;
     grad_state_x.resize(5);
     RealVector grad_state_y;
     grad_state_y.resize(5);
     RealVector grad_state_z;
     grad_state_z.resize(5);



     for (CFuint i = 0; i<5; ++i){
         grad_state_x[i] = (nx0*(4.0*L0 - 1.0)*(*_states[0])[i] +
		            nx1*(4.0*L1 - 1.0)*(*_states[1])[i] +
		            nx2*(4.0*L2 - 1.0)*(*_states[2])[i] +
			    nx3*(4.0*L3 - 1.0)*(*_states[3])[i] +
		            4.0*(nx1*L0 + nx0*L1)*(*_states[4])[i] +
		            4.0*(nx2*L1 + nx1*L2)*(*_states[5])[i] +
		            4.0*(nx0*L2 + nx2*L0)*(*_states[6])[i] +
                            4.0*(nx1*L3 + nx3*L1)*(*_states[7])[i] +
                            4.0*(nx2*L3 + nx3*L2)*(*_states[8])[i] +
                            4.0*(nx0*L3 + nx3*L0)*(*_states[9])[i])*inv_sixvolume;

        grad_state_y[i] = (ny0*(4.0*L0 - 1.0)*(*_states[0])[i] +
                            ny1*(4.0*L1 - 1.0)*(*_states[1])[i] +
                            ny2*(4.0*L2 - 1.0)*(*_states[2])[i] +
                            ny3*(4.0*L3 - 1.0)*(*_states[3])[i] +
                            4.0*(ny1*L0 + ny0*L1)*(*_states[4])[i] +
                            4.0*(ny2*L1 + ny1*L2)*(*_states[5])[i] +
                            4.0*(ny0*L2 + ny2*L0)*(*_states[6])[i] +
                            4.0*(ny1*L3 + ny3*L1)*(*_states[7])[i] +
                            4.0*(ny2*L3 + ny3*L2)*(*_states[8])[i] +
                            4.0*(ny0*L3 + ny3*L0)*(*_states[9])[i])*inv_sixvolume;

         grad_state_z[i] = (nz0*(4.0*L0 - 1.0)*(*_states[0])[i] +
                            nz1*(4.0*L1 - 1.0)*(*_states[1])[i] +
                            nz2*(4.0*L2 - 1.0)*(*_states[2])[i] +
                            nz3*(4.0*L3 - 1.0)*(*_states[3])[i] +
                            4.0*(nz1*L0 + nz0*L1)*(*_states[4])[i] +
                            4.0*(nz2*L1 + nz1*L2)*(*_states[5])[i] +
                            4.0*(nz0*L2 + nz2*L0)*(*_states[6])[i] +
                            4.0*(nz1*L3 + nz3*L1)*(*_states[7])[i] +
                            4.0*(nz2*L3 + nz3*L2)*(*_states[8])[i] +
                            4.0*(nz0*L3 + nz3*L0)*(*_states[9])[i])*inv_sixvolume;




        }



     CFreal grad_T_x;
     grad_T_x = (nx0*(4.0*L0 - 1.0)*_values(4,0) +
		 nx1*(4.0*L1 - 1.0)*_values(4,1) +
                 nx2*(4.0*L2 - 1.0)*_values(4,2) +
                 nx3*(4.0*L3 - 1.0)*_values(4,3) + 
                 4.0*(nx1*L0 + nx0*L1)*_values(4,4)  +
                 4.0*(nx2*L1 + nx1*L2)*_values(4,5)  +
                 4.0*(nx0*L2 + nx2*L0)*_values(4,6)  +
                 4.0*(nx1*L3 + nx3*L1)*_values(4,7)  +
                 4.0*(nx2*L3 + nx3*L2)*_values(4,8)  +
                 4.0*(nx0*L3 + nx3*L0)*_values(4,9) )*inv_sixvolume;

     CFreal grad_T_y;
     grad_T_y = (ny0*(4.0*L0 - 1.0)*_values(4,0) +
                 ny1*(4.0*L1 - 1.0)*_values(4,1) +
                 ny2*(4.0*L2 - 1.0)*_values(4,2) +
                 ny3*(4.0*L3 - 1.0)*_values(4,3) +
                 4.0*(ny1*L0 + ny0*L1)*_values(4,4)  +
                 4.0*(ny2*L1 + ny1*L2)*_values(4,5)  +
                 4.0*(ny0*L2 + ny2*L0)*_values(4,6)  +
                 4.0*(ny1*L3 + ny3*L1)*_values(4,7)  +
                 4.0*(ny2*L3 + ny3*L2)*_values(4,8)  +
                 4.0*(ny0*L3 + ny3*L0)*_values(4,9) )*inv_sixvolume;

     CFreal grad_T_z;
     grad_T_z = (nz0*(4.0*L0 - 1.0)*_values(4,0) +
                 nz1*(4.0*L1 - 1.0)*_values(4,1) +
                 nz2*(4.0*L2 - 1.0)*_values(4,2) +
                 nz3*(4.0*L3 - 1.0)*_values(4,3) +
                 4.0*(nz1*L0 + nz0*L1)*_values(4,4)  +
                 4.0*(nz2*L1 + nz1*L2)*_values(4,5)  +
                 4.0*(nz0*L2 + nz2*L0)*_values(4,6)  +
                 4.0*(nz1*L3 + nz3*L1)*_values(4,7)  +
                 4.0*(nz2*L3 + nz3*L2)*_values(4,8)  +
                 4.0*(nz0*L3 + nz3*L0)*_values(4,9) )*inv_sixvolume;


     (qdstates) = (L0*( 2.0*L0 - 1.0 ) * (*_states[0])) +
                  (L1*( 2.0*L1 - 1.0 ) * (*_states[1])) +
                  (L2*( 2.0*L2 - 1.0 ) * (*_states[2])) +
		  (L3*( 2.0*L3 - 1.0 ) * (*_states[3])) +
                  (4.0*L0*L1           * (*_states[4])) +
                  (4.0*L1*L2           * (*_states[5])) +
                  (4.0*L0*L2           * (*_states[6])) +
	          (4.0*L1*L3           * (*_states[7])) +
                  (4.0*L3*L2           * (*_states[8])) +
                  (4.0*L0*L3           * (*_states[9]));



     (*_gradients[0])[XX] = grad_state_x[0];
     (*_gradients[0])[YY] = grad_state_y[0];
     (*_gradients[0])[ZZ] = grad_state_z[0];
   
     (*_gradients[1])[XX] = -(qdstates[1]/(qdstates[0]*qdstates[0]))*grad_state_x[0] + (1.0/qdstates[0])*grad_state_x[1];
     (*_gradients[1])[YY] = -(qdstates[1]/(qdstates[0]*qdstates[0]))*grad_state_y[0] + (1.0/qdstates[0])*grad_state_y[1];
     (*_gradients[1])[ZZ] = -(qdstates[1]/(qdstates[0]*qdstates[0]))*grad_state_z[0] + (1.0/qdstates[0])*grad_state_z[1];
   
     (*_gradients[2])[XX] = -(qdstates[2]/(qdstates[0]*qdstates[0]))*grad_state_x[0] + (1.0/qdstates[0])*grad_state_x[2];
     (*_gradients[2])[YY] = -(qdstates[2]/(qdstates[0]*qdstates[0]))*grad_state_y[0] + (1.0/qdstates[0])*grad_state_y[2];
     (*_gradients[2])[ZZ] = -(qdstates[2]/(qdstates[0]*qdstates[0]))*grad_state_z[0] + (1.0/qdstates[0])*grad_state_z[2];     
   
     (*_gradients[3])[XX] = -(qdstates[3]/(qdstates[0]*qdstates[0]))*grad_state_x[0] + (1.0/qdstates[0])*grad_state_x[3];
     (*_gradients[3])[YY] = -(qdstates[3]/(qdstates[0]*qdstates[0]))*grad_state_y[0] + (1.0/qdstates[0])*grad_state_y[3];
     (*_gradients[3])[ZZ] = -(qdstates[3]/(qdstates[0]*qdstates[0]))*grad_state_z[0] + (1.0/qdstates[0])*grad_state_z[3];
   
     (*_gradients[4])[XX] = grad_T_x;
     (*_gradients[4])[YY] = grad_T_y;
     (*_gradients[4])[ZZ] = grad_T_z;


     _normal[XX] = nx1;
     _normal[YY] = ny1;
     _normal[ZZ] = nz1;

     RealVector F1;
     F1.resize(5);
     F1 = _diffVar->getFlux(qdstates, _gradients,_normal , 0.0);

     _normal[XX] = nx2;
     _normal[YY] = ny2;
     _normal[ZZ] = nz2;	

     RealVector F2;
     F2.resize(5);
     F2 = _diffVar->getFlux(qdstates, _gradients,_normal , 0.0);

     _normal[XX] = nx3;
     _normal[YY] = ny3;
     _normal[ZZ] = nz3;	

     RealVector F3;
     F3.resize(5);
     F3 = _diffVar->getFlux(qdstates, _gradients,_normal , 0.0);

     _normal[XX] = nx0;
     _normal[YY] = ny0;
     _normal[ZZ] = nz0;

     RealVector F0;
     F0.resize(5);
     F0 = _diffVar->getFlux(qdstates, _gradients,_normal , 0.0);
	


     (*_phi_diff_gal[0])[0] = 0.0;
     (*_phi_diff_gal[1])[0] = 0.0;
     (*_phi_diff_gal[2])[0] = 0.0;
     (*_phi_diff_gal[3])[0] = 0.0;
     (*_phi_diff_gal[4])[0] = 0.0;
     (*_phi_diff_gal[5])[0] = 0.0;

	CFreal coeff;
	if (!triag) coeff = 1.0/48.8;
        else coeff = 1.0/24.0;


      for (CFuint i = 1 ; i < 5; ++i)
	{
          (*_phi_diff_gal[0])[i] += F0[i]*(4.0*L0 - 1.0)*coeff*wqd[iQd];
          (*_phi_diff_gal[1])[i] += F1[i]*(4.0*L1 - 1.0)*coeff*wqd[iQd];
          (*_phi_diff_gal[2])[i] += F2[i]*(4.0*L2 - 1.0)*coeff*wqd[iQd];
	  (*_phi_diff_gal[3])[i] += F3[i]*(4.0*L3 - 1.0)*coeff*wqd[iQd];
          (*_phi_diff_gal[4])[i] += 4.0*(F0[i]*L1 + F1[i]*L0)*coeff*wqd[iQd];
          (*_phi_diff_gal[5])[i] += 4.0*(F2[i]*L1 + F1[i]*L2)*coeff*wqd[iQd];
          (*_phi_diff_gal[6])[i] += 4.0*(F0[i]*L2 + F2[i]*L0)*coeff*wqd[iQd];
          (*_phi_diff_gal[7])[i] += 4.0*(F1[i]*L3 + F3[i]*L1)*coeff*wqd[iQd];
          (*_phi_diff_gal[8])[i] += 4.0*(F2[i]*L3 + F3[i]*L2)*coeff*wqd[iQd];
          (*_phi_diff_gal[9])[i] += 4.0*(F0[i]*L3 + F3[i]*L0)*coeff*wqd[iQd];

       }

  }

}

//////////////////////////////////////////////////////////////////////////////
void NavierStokesTermHO3D::fluctuation_diff_bubble(CFuint& i1, CFuint& i2, CFuint& i3, CFuint& i4, bool triag)
{

  Node& node0 = (*m_cellStates)[i1]->getCoordinates();
  Node& node1 = (*m_cellStates)[i2]->getCoordinates();
  Node& node2 = (*m_cellStates)[i3]->getCoordinates();
  Node& node3 = (*m_cellStates)[i4]->getCoordinates();

  CFreal x0 = node0[XX];
  CFreal y0 = node0[YY];
  CFreal z0 = node0[ZZ];
  
  CFreal x1 = node1[XX];
  CFreal y1 = node1[YY];
  CFreal z1 = node1[ZZ];

  CFreal x2 = node2[XX];
  CFreal y2 = node2[YY];
  CFreal z2 = node2[ZZ];

  CFreal x3 = node3[XX];
  CFreal y3 = node3[YY];
  CFreal z3 = node3[ZZ];

  Node& node0_vert = (*m_cellStates)[0]->getCoordinates();
  Node& node1_vert = (*m_cellStates)[1]->getCoordinates();
  Node& node2_vert = (*m_cellStates)[2]->getCoordinates();
  Node& node3_vert = (*m_cellStates)[3]->getCoordinates();

  CFreal x0_vert = node0_vert[XX];
  CFreal y0_vert = node0_vert[YY];
  CFreal z0_vert = node0_vert[ZZ];
  

  CFreal x1_vert = node1_vert[XX];
  CFreal y1_vert = node1_vert[YY];
  CFreal z1_vert = node1_vert[ZZ];

  CFreal x2_vert = node2_vert[XX];
  CFreal y2_vert = node2_vert[YY];
  CFreal z2_vert = node2_vert[ZZ];

  CFreal x3_vert = node3_vert[XX];
  CFreal y3_vert = node3_vert[YY];
  CFreal z3_vert = node3_vert[ZZ];

  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();


  const CFreal nx0_bub = matrix_node_norms(0,XX);
  const CFreal nx1_bub = matrix_node_norms(1,XX);
  const CFreal nx2_bub = matrix_node_norms(2,XX);
  const CFreal nx3_bub = matrix_node_norms(3,XX);

  const CFreal ny0_bub = matrix_node_norms(0,YY);
  const CFreal ny1_bub = matrix_node_norms(1,YY);
  const CFreal ny2_bub = matrix_node_norms(2,YY);
  const CFreal ny3_bub = matrix_node_norms(3,YY);


  const CFreal nz0_bub = matrix_node_norms(0,ZZ);
  const CFreal nz1_bub = matrix_node_norms(1,ZZ);
  const CFreal nz2_bub = matrix_node_norms(2,ZZ);
  const CFreal nz3_bub = matrix_node_norms(3,ZZ);

  
  InwardNormalsData& cellnormals = (*socket_normals.getDataHandle()[m_cellID]);
// CF_DEBUG_POINT;
  const CFreal nx0 = cellnormals.getNodalNormComp(0,XX);
  const CFreal nx1 = cellnormals.getNodalNormComp(1,XX);
  const CFreal nx2 = cellnormals.getNodalNormComp(2,XX);
  const CFreal nx3 = cellnormals.getNodalNormComp(3,XX);



  const CFreal ny0 = cellnormals.getNodalNormComp(0,YY);
  const CFreal ny1 = cellnormals.getNodalNormComp(1,YY);
  const CFreal ny2 = cellnormals.getNodalNormComp(2,YY);
  const CFreal ny3 = cellnormals.getNodalNormComp(3,YY);


  const CFreal nz0 = cellnormals.getNodalNormComp(0,ZZ);
  const CFreal nz1 = cellnormals.getNodalNormComp(1,ZZ);
  const CFreal nz2 = cellnormals.getNodalNormComp(2,ZZ);
  const CFreal nz3 = cellnormals.getNodalNormComp(3,ZZ);

  CFreal inv_sixvolume = 1.0/(6.0*_cellVolume);
  CFreal inv_sixvolumebub;
  if (!triag) inv_sixvolumebub = 8.0/(6.0*_cellVolume);
  else inv_sixvolumebub = 4.0/(6.0*_cellVolume);

  // The buble function is linear peacewise linear. Which means that its gradient is not continuous
  // Then, we just need one node per constant part. i choose to take one point at the gravity center of the 3 triangkes
  // composed by nodes and gravity center of the sub-element.
  for (CFuint iEq = 0; iEq < nbEqs; ++ iEq)
     (*_phi_diff_bub)[iEq] = 0.0;

  for (CFuint iQd = 0; iQd < 4; ++iQd) {
    //point od quadrature
    const CFreal x = qd0[iQd] * node0[XX] + qd1[iQd] * node1[XX] + qd2[iQd] * node2[XX] + qd3[iQd]*node3[XX];
    const CFreal y = qd0[iQd] * node0[YY] + qd1[iQd] * node1[YY] + qd2[iQd] * node2[YY] + qd3[iQd]*node3[YY];
    const CFreal z = qd0[iQd] * node0[ZZ] + qd1[iQd] * node1[ZZ] + qd2[iQd] * node2[ZZ] + qd3[iQd]*node3[ZZ];
    
     
     
    // Linear basis function
    CFreal L0bub = 1.0 + ( ( x - x0 )*nx0_bub + ( y - y0 )*ny0_bub + (z - z0)*nz0_bub )*inv_sixvolumebub ;
    CFreal L1bub = 1.0 + ( ( x - x1 )*nx1_bub + ( y - y1 )*ny1_bub + (z - z1)*nz1_bub )*inv_sixvolumebub ;
    CFreal L2bub = 1.0 + ( ( x - x2 )*nx2_bub + ( y - y2 )*ny2_bub + (z - z2)*nz2_bub )*inv_sixvolumebub ;
    CFreal L3bub = 1.0 + ( ( x - x3 )*nx3_bub + ( y - y3 )*ny3_bub + (z - z3)*nz3_bub )*inv_sixvolumebub ;
          



    CFreal L0 = 1.0 + ( ( x - x0_vert )*nx0 + ( y - y0_vert )*ny0 + (z - z0_vert)*nz0 )*inv_sixvolume ;
    CFreal L1 = 1.0 + ( ( x - x1_vert )*nx1 + ( y - y1_vert )*ny1 + (z - z1_vert)*nz1 )*inv_sixvolume ;
    CFreal L2 = 1.0 + ( ( x - x2_vert )*nx2 + ( y - y2_vert )*ny2 + (z - z2_vert)*nz2 )*inv_sixvolume ;
    CFreal L3 = 1.0 + ( ( x - x3_vert )*nx3 + ( y - y3_vert )*ny3 + (z - z3_vert)*nz3 )*inv_sixvolume ;

    
    RealVector grad_state_x;
    grad_state_x.resize(5);
    RealVector grad_state_y;
    grad_state_y.resize(5);
    RealVector grad_state_z;
    grad_state_z.resize(5);
    
 for (CFuint i = 0; i<5; ++i){
         grad_state_x[i] = (nx0*(4.0*L0 - 1.0)*(*_states[0])[i] +
                            nx1*(4.0*L1 - 1.0)*(*_states[1])[i] +
                            nx2*(4.0*L2 - 1.0)*(*_states[2])[i] +
                            nx3*(4.0*L3 - 1.0)*(*_states[3])[i] +
                            4.0*(nx1*L0 + nx0*L1)*(*_states[4])[i] +
                            4.0*(nx2*L1 + nx1*L2)*(*_states[5])[i] +
                            4.0*(nx0*L2 + nx2*L0)*(*_states[6])[i] +
                            4.0*(nx1*L3 + nx3*L1)*(*_states[7])[i] +
                            4.0*(nx2*L3 + nx3*L2)*(*_states[8])[i] +
                            4.0*(nx0*L3 + nx3*L0)*(*_states[9])[i])*inv_sixvolume;

        grad_state_y[i] = (ny0*(4.0*L0 - 1.0)*(*_states[0])[i] +
                            ny1*(4.0*L1 - 1.0)*(*_states[1])[i] +
                            ny2*(4.0*L2 - 1.0)*(*_states[2])[i] +
                            ny3*(4.0*L3 - 1.0)*(*_states[3])[i] +
                            4.0*(ny1*L0 + ny0*L1)*(*_states[4])[i] +
                            4.0*(ny2*L1 + ny1*L2)*(*_states[5])[i] +
                            4.0*(ny0*L2 + ny2*L0)*(*_states[6])[i] +
                            4.0*(ny1*L3 + ny3*L1)*(*_states[7])[i] +
                            4.0*(ny2*L3 + ny3*L2)*(*_states[8])[i] +
                            4.0*(ny0*L3 + ny3*L0)*(*_states[9])[i])*inv_sixvolume;

         grad_state_z[i] = (nz0*(4.0*L0 - 1.0)*(*_states[0])[i] +
                            nz1*(4.0*L1 - 1.0)*(*_states[1])[i] +
                            nz2*(4.0*L2 - 1.0)*(*_states[2])[i] +
                            nz3*(4.0*L3 - 1.0)*(*_states[3])[i] +
                            4.0*(nz1*L0 + nz0*L1)*(*_states[4])[i] +
                            4.0*(nz2*L1 + nz1*L2)*(*_states[5])[i] +
                            4.0*(nz0*L2 + nz2*L0)*(*_states[6])[i] +
                            4.0*(nz1*L3 + nz3*L1)*(*_states[7])[i] +
                            4.0*(nz2*L3 + nz3*L2)*(*_states[8])[i] +
                            4.0*(nz0*L3 + nz3*L0)*(*_states[9])[i])*inv_sixvolume;




        }
     CFreal grad_T_x;
     grad_T_x = (nx0*(4.0*L0 - 1.0)*_values(4,0) +
                 nx1*(4.0*L1 - 1.0)*_values(4,1) +
                 nx2*(4.0*L2 - 1.0)*_values(4,2) +
                 nx3*(4.0*L3 - 1.0)*_values(4,3) +
                 4.0*(nx1*L0 + nx0*L1)*_values(4,4)  +
                 4.0*(nx2*L1 + nx1*L2)*_values(4,5)  +
                 4.0*(nx0*L2 + nx2*L0)*_values(4,6)  +
                 4.0*(nx1*L3 + nx3*L1)*_values(4,7)  +
                 4.0*(nx2*L3 + nx3*L2)*_values(4,8)  +
                 4.0*(nx0*L3 + nx3*L0)*_values(4,9) )*inv_sixvolume;

     CFreal grad_T_y;
     grad_T_y = (ny0*(4.0*L0 - 1.0)*_values(4,0) +
                 ny1*(4.0*L1 - 1.0)*_values(4,1) +
                 ny2*(4.0*L2 - 1.0)*_values(4,2) +
                 ny3*(4.0*L3 - 1.0)*_values(4,3) +
                 4.0*(ny1*L0 + ny0*L1)*_values(4,4)  +
                 4.0*(ny2*L1 + ny1*L2)*_values(4,5)  +
                 4.0*(ny0*L2 + ny2*L0)*_values(4,6)  +
                 4.0*(ny1*L3 + ny3*L1)*_values(4,7)  +
                 4.0*(ny2*L3 + ny3*L2)*_values(4,8)  +
                 4.0*(ny0*L3 + ny3*L0)*_values(4,9) )*inv_sixvolume;

     CFreal grad_T_z;
     grad_T_z = (nz0*(4.0*L0 - 1.0)*_values(4,0) +
                 nz1*(4.0*L1 - 1.0)*_values(4,1) +
                 nz2*(4.0*L2 - 1.0)*_values(4,2) +
                 nz3*(4.0*L3 - 1.0)*_values(4,3) +
                 4.0*(nz1*L0 + nz0*L1)*_values(4,4)  +
                 4.0*(nz2*L1 + nz1*L2)*_values(4,5)  +
                 4.0*(nz0*L2 + nz2*L0)*_values(4,6)  +
                 4.0*(nz1*L3 + nz3*L1)*_values(4,7)  +
                 4.0*(nz2*L3 + nz3*L2)*_values(4,8)  +
                 4.0*(nz0*L3 + nz3*L0)*_values(4,9) )*inv_sixvolume;

  (qdstates) = (L0*( 2.0*L0 - 1.0 ) * (*_states[0])) +
                  (L1*( 2.0*L1 - 1.0 ) * (*_states[1])) +
                  (L2*( 2.0*L2 - 1.0 ) * (*_states[2])) +
                  (L3*( 2.0*L3 - 1.0 ) * (*_states[3])) +
                  (4.0*L0*L1           * (*_states[4])) +
                  (4.0*L1*L2           * (*_states[5])) +
                  (4.0*L0*L2           * (*_states[6])) +
                  (4.0*L1*L3           * (*_states[7])) +
                  (4.0*L3*L2           * (*_states[8])) +
                  (4.0*L0*L3           * (*_states[9]));






     (*_gradients[0])[XX] = grad_state_x[0];
     (*_gradients[0])[YY] = grad_state_y[0];
     (*_gradients[0])[ZZ] = grad_state_z[0];
     (*_gradients[1])[XX] = -(qdstates[1]/(qdstates[0]*qdstates[0]))*grad_state_x[0] + (1.0/qdstates[0])*grad_state_x[1];
     (*_gradients[1])[YY] = -(qdstates[1]/(qdstates[0]*qdstates[0]))*grad_state_y[0] + (1.0/qdstates[0])*grad_state_y[1];
     (*_gradients[1])[ZZ] = -(qdstates[1]/(qdstates[0]*qdstates[0]))*grad_state_z[0] + (1.0/qdstates[0])*grad_state_z[1];
     (*_gradients[2])[XX] = -(qdstates[2]/(qdstates[0]*qdstates[0]))*grad_state_x[0] + (1.0/qdstates[0])*grad_state_x[2];
     (*_gradients[2])[YY] = -(qdstates[2]/(qdstates[0]*qdstates[0]))*grad_state_y[0] + (1.0/qdstates[0])*grad_state_y[2];
     (*_gradients[2])[ZZ] = -(qdstates[2]/(qdstates[0]*qdstates[0]))*grad_state_z[0] + (1.0/qdstates[0])*grad_state_z[2];
     (*_gradients[3])[XX] = -(qdstates[3]/(qdstates[0]*qdstates[0]))*grad_state_x[0] + (1.0/qdstates[0])*grad_state_x[3];
     (*_gradients[3])[YY] = -(qdstates[3]/(qdstates[0]*qdstates[0]))*grad_state_y[0] + (1.0/qdstates[0])*grad_state_y[3];
     (*_gradients[3])[ZZ] = -(qdstates[3]/(qdstates[0]*qdstates[0]))*grad_state_z[0] + (1.0/qdstates[0])*grad_state_z[3];
     (*_gradients[4])[XX] = grad_T_x;
     (*_gradients[4])[YY] = grad_T_y;
     (*_gradients[4])[ZZ] = grad_T_z;




  _normal[XX] = nx0_bub;;
  _normal[YY] = ny0_bub;
  _normal[ZZ] = nz0_bub;

  RealVector F0;
  F0.resize(nbEqs);
  F0 = _diffVar->getFlux(qdstates, _gradients,_normal , 0.0);

  _normal[XX] = nx1_bub;
  _normal[YY] = ny1_bub;
  _normal[ZZ] = nz1_bub;

  RealVector F1;
  F1.resize(nbEqs);
  F1 = _diffVar->getFlux(qdstates, _gradients,_normal , 0.0);



  _normal[XX] = nx2_bub;
  _normal[YY] = ny2_bub;
  _normal[ZZ] = nz2_bub;

  RealVector F2;
  F2.resize(nbEqs);
  F2 = _diffVar->getFlux(qdstates, _gradients,_normal , 0.0);

  _normal[XX] = nx3_bub;
  _normal[YY] = ny3_bub;
  _normal[ZZ] = nz3_bub;

  RealVector F3;
  F3.resize(nbEqs);
  F3 = _diffVar->getFlux(qdstates, _gradients,_normal , 0.0);

     CFreal coeff = 1.0/6.0;

      for (CFuint i = 1 ; i < 5; ++i)
        {
          (*_phi_diff_bub)[i] += 256.0*(F0[i]*L1bub*L2bub*L3bub + F1[i]*L0bub*L2bub*L3bub + F2[i]*L0bub*L1bub*L3bub + F3[i]*L0bub*L1bub*L2bub)*coeff*wqd[iQd];

}
}


}
//////////////////////////////////////////////////////////////////////////////

void NavierStokesTermHO3D::setup()
{
  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  const CFuint nbStatesInControlVolume =
    MeshDataStack::getActive()->Statistics().getMaxNbStatesInCell();
  const CFuint nbNodesInControlVolume =
    MeshDataStack::getActive()->Statistics().getMaxNbNodesInCell();
  
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
  qd3.resize(nbQdPts);


CFreal alpha = 0.58541020;
CFreal beta = 0.13819660;
  qd0[0] = alpha;  qd1[0] = beta;  qd2[0] = beta;   qd3[0] = beta;
  qd0[1] = beta;   qd1[1] = alpha; qd2[1] = beta;   qd3[1] = beta; 
  qd0[2] = beta;   qd1[2] = beta;  qd2[2] = alpha;  qd3[2] = beta;
  qd0[3] = beta;   qd1[3] = beta;  qd2[3] = beta;   qd3[3] = alpha;

  wqd.resize(nbQdPts); // 4 quadrature points per surface

  wqd[0] = 0.25;
  wqd[1] = 0.25;
  wqd[2] = 0.25;
  wqd[3] = 0.25;

  qdstates.resize(nbEqs);

  _phi_diff_bub = new RealVector(nbEqs);
  _phi_diff_bub_split.resize(nbNodesInControlVolume); //4  residuals in each sub element

  _phi_diff_bub_split[0] = new RealVector(nbEqs);
  _phi_diff_bub_split[1] = new RealVector(nbEqs);
  _phi_diff_bub_split[2] = new RealVector(nbEqs);
  _phi_diff_bub_split[3] = new RealVector(nbEqs);	

  _phi_diff_gal.resize(nbStatesInControlVolume); // one galerkin function per state

  for (CFuint iState = 0; iState < nbStatesInControlVolume; ++ iState)
     _phi_diff_gal[iState] = new RealVector(nbEqs);

  kappa.resize(7,10); // 7 sub-elements 10 nodes

  CFreal one_fourth = 1.0/4.0;
  alpha = 63.0/256.0;
  beta = -147.0/256.0;
  CFreal gamma = -63.0/64.0;
  CFreal delta = -21.0/128.0;
  CFreal lambda = 21.0/256.0;
  CFreal eta = -105.0/64.0;
  CFreal theta = -21.0/32.0;
  CFreal epsilon = -63.0/128.0;
  CFreal phi = 21.0/64.0;
  CFreal psi = -105.0/128.0;


  // kappa is used to "distribute" the part of the fluctuation with the bubble function
  kappa(0,0) = alpha;
  kappa(0,1) = beta;
  kappa(0,2) = alpha;
  kappa(0,3) = alpha;
  kappa(0,4) = gamma;
  kappa(0,5) = gamma;
  kappa(0,6) = delta;
  kappa(0,7) = gamma;
  kappa(0,8) = delta;
  kappa(0,9) = delta;

  //OK


  kappa(1,0) = alpha;
  kappa(1,1) = alpha;
  kappa(1,2) = beta;
  kappa(1,3) = alpha;
  kappa(1,4) = delta;
  kappa(1,5) = gamma;
  kappa(1,6) = gamma;
  kappa(1,7) = delta;
  kappa(1,8) = gamma;
  kappa(1,9) = delta;

  //OK


  kappa(2,0) = lambda;
  kappa(2,1) = alpha;
  kappa(2,2) = alpha;
  kappa(2,3) = alpha;
  kappa(2,4) = eta;
  kappa(2,5) = theta;
  kappa(2,6) = epsilon;
  kappa(2,7) = theta;
  kappa(2,8) = delta;
  kappa(2,9) = epsilon;

  //OK

  kappa(3,0) = lambda;
  kappa(3,1) = alpha;
  kappa(3,2) = alpha;
  kappa(3,3) = alpha;
  kappa(3,4) = epsilon;
  kappa(3,5) = theta;
  kappa(3,6) = eta;
  kappa(3,7) = delta;
  kappa(3,8) = theta;
  kappa(3,9) = epsilon;
//OK

  kappa(4,0) = -delta;
  kappa(4,1) = phi;
  kappa(4,2) = phi;
  kappa(4,3) = phi;
  kappa(4,4) = theta;
  kappa(4,5) = psi;
  kappa(4,6) = theta;
  kappa(4,7) = psi;
  kappa(4,8) = psi;
  kappa(4,9) = theta;


   kappa(5,0) = alpha;
   kappa(5,1) = alpha;
   kappa(5,2) = alpha;
   kappa(5,3) = beta;
   kappa(5,4) = delta;
   kappa(5,5) = delta;
   kappa(5,6) = delta;
   kappa(5,7) = gamma;
   kappa(5,8) = gamma;
   kappa(5,9) = gamma;
//OK

  kappa(6,0) = lambda;
  kappa(6,1) = alpha;
  kappa(6,2) = alpha;
  kappa(6,3) = alpha;
  kappa(6,4) = epsilon;
  kappa(6,5) = delta;
  kappa(6,6) = epsilon;
  kappa(6,7) = theta;
  kappa(6,8) = theta;
  kappa(6,9) = eta;
//OK

  getMethodData().getDistributionData().computeBetas = true;

  FacesMat.resize(22,3);
  matrix_node_norms.resize(4,3); 
}

//////////////////////////////////////////////////////////////////////////////

void NavierStokesTermHO3D::configure ( Config::ConfigArgs& args )
{
  ComputeDiffusiveTerm::configure(args);
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
