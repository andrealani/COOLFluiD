#include "Framework/CFL.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/LSSMatrix.hh"
#include "Framework/BlockAccumulator.hh"
#include "Framework/MeshData.hh"

#include "FiniteVolume/FiniteVolume.hh"
#include "FiniteVolume/FVMCC_CrankNichLimComputeRhs.hh"
#include "FiniteVolume/FVMCC_BC.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<FVMCC_CrankNichLimComputeRhs,
		      CellCenterFVMData,
		      FiniteVolumeModule>
FVMCC_CrankNichLimComputeRhs("CrankNichLim");

//////////////////////////////////////////////////////////////////////////////

FVMCC_CrankNichLimComputeRhs::FVMCC_CrankNichLimComputeRhs(const std::string& name) :
  FVMCC_ComputeRhsJacob(name),
  socket_timeLimiter("timeLimiter"),
  _upFactorVec(2),
  _upStFactorVec(2)  
{
}

//////////////////////////////////////////////////////////////////////////////

FVMCC_CrankNichLimComputeRhs::~FVMCC_CrankNichLimComputeRhs()
{
}

//////////////////////////////////////////////////////////////////////////////
      
vector<SafePtr<BaseDataSocketSink> > FVMCC_CrankNichLimComputeRhs::needsSockets()
{
  vector<SafePtr<BaseDataSocketSink> > result = FVMCC_ComputeRhsJacob::needsSockets();
  
  result.push_back(&socket_timeLimiter);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

void FVMCC_CrankNichLimComputeRhs::setup()
{
  FVMCC_ComputeRhsJacob::setup();
  
  for (CFuint i = 0; i < 2; ++i) {
    _upFactorVec[i].resize(PhysicalModelStack::getActive()->getNbEq());
  }
  
  for (CFuint i = 0; i < 2; ++i) {
    _upStFactorVec[i].resize(PhysicalModelStack::getActive()->getNbEq());
  }
}

//////////////////////////////////////////////////////////////////////////////

void FVMCC_CrankNichLimComputeRhs::unsetup()
{
  
  for (CFuint i = 0; i < 2; ++i) {
    _upFactorVec[i].resize(0);
  }
  _upFactorVec.resize(0);
  
  for (CFuint i = 0; i < 2; ++i) {
    _upStFactorVec[i].resize(0);
  }
  _upStFactorVec.resize(0);
  
  FVMCC_ComputeRhsJacob::unsetup();
}

//////////////////////////////////////////////////////////////////////////////

void FVMCC_CrankNichLimComputeRhs::updateRHS()
{
  if (getMethodData().isAxisymmetric()) {
    cf_assert(_currFace->nbNodes() == 2);
    const Node *const  node0 = _currFace->getNode(0);
    const Node *const  node1 = _currFace->getNode(1);
    
    // distance of the face mid point to the axis
    // _rMid = 0 on a centerline boundary face
    // _rMid == average y between the two nodes the face
    _rMid = 0.5*std::abs((*node0)[YY] + (*node1)[YY]);
    _invr[0] = 1./_currFace->getState(0)->getCoordinates()[YY];
    _invr[1] = 1./_currFace->getState(1)->getCoordinates()[YY];
  }
  
  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  for (CFuint iEq = 0; iEq < nbEqs; ++iEq) {
    _rFlux = (getResFactorCN(iEq)*_rMid)*_flux[iEq];
  }
  
  // distribute the computed flux to the two neighbor states
  // of the corresponding face, subtracting the residual to the
  // current state and summing it to the neighbor state
  State *const firstState = _currFace->getState(0);
  const CFuint firstStateID = firstState->getLocalID();
  DataHandle<CFreal> rhs = socket_rhs.getDataHandle();
  for (CFuint iEq = 0; iEq < nbEqs; ++iEq) {
    rhs(firstStateID, iEq, nbEqs) -= _rFlux[iEq]*_invr[0];
    CFLogDebugMed("rhs " << rhs(firstStateID, iEq, nbEqs) << "\n");
  }
  
  CFLogDebugMed("updateCoeff[firstStateID] = "
                << socket_updateCoeff.getDataHandle()[firstStateID] << "\n");
  
  // cout.precision(18); cout << firstStateID << ", "   << socket_updateCoeff.getDataHandle()[firstStateID] << endl;
  
  State *const lastState = _currFace->getState(1);
  if (!lastState->isGhost()) {
    const CFuint lastStateID = lastState->getLocalID();
    cf_assert(firstStateID < lastStateID);
    for (CFuint iEq = 0; iEq < nbEqs; ++iEq) {
      rhs(lastStateID, iEq, nbEqs) += _rFlux[iEq]*_invr[1];
      CFLogDebugMed("rhs " << rhs(lastStateID, iEq, nbEqs) << "\n");
    }
    
    CFLogDebugMed("updateCoeff[lastStateID] = "
	       << socket_updateCoeff.getDataHandle()[lastStateID]  << "\n");
    
    // cout.precision(18); cout << lastStateID << ", "   << socket_updateCoeff.getDataHandle()[lastStateID] << endl;
  }
  
  // set the cell flags to true if the cells have already been used
  DataHandle<bool> cellFlag = socket_cellFlag.getDataHandle();
  for (CFuint i = 0; i < 2; ++i) { 
    if (!_currFace->getState(i)->isGhost()) {
      cellFlag[_currFace->getState(i)->getLocalID()] = true;
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void FVMCC_CrankNichLimComputeRhs::computeBothJacobTerms()
{
  // set on the perturbation flag
  getMethodData().setIsPerturb(true);
  
  // set the index of the block corresponding to the current
  // state in the jacobian matrix
  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  
  State& state0 = *_currFace->getState(0);
  State& state1 = *_currFace->getState(1);
    
  _acc->setRowColIndex(0, state0.getLocalID());
  _acc->setRowColIndex(1, state1.getLocalID());
  
  // first node contribution
  for (CFuint iVar = 0; iVar < nbEqs; ++iVar) {
    CFLogDebugMax( "Perturbing iVar = " << iVar << "\n");
    
    // perturb the given component of the state vector
    _numericalJacob->perturb(iVar, state0[iVar]); 
    (!getMethodData().isAxisymmetric()) ? computeNoAxiUpFactorsVec() : computeAxiUpFactorsVec();
    computeConvDiffFluxes(iVar, LEFT);  
    addBothJacobTermsCN(iVar, LEFT);
    
    // restore the unperturbed value
    _numericalJacob->restore(state0[iVar]);
    _polyRec->restoreValues(iVar, LEFT);
    
    // perturb the given component of the state vector
    _numericalJacob->perturb(iVar, state1[iVar]);
    (!getMethodData().isAxisymmetric()) ? computeNoAxiUpFactorsVec() : computeAxiUpFactorsVec();
    computeConvDiffFluxes(iVar, RIGHT);
    addBothJacobTermsCN(iVar, RIGHT);
    
    // restore the unperturbed value
    _numericalJacob->restore(state1[iVar]);
    _polyRec->restoreValues(iVar, RIGHT);
  }
  
  // compute analytical jacobian for source term 
  if (computeSourceTermJacob(LEFT,_stAnJacobIDs)) {
    addAnalyticSourceTermJacobCN(LEFT, _acc.get());
  }
  
  if (computeSourceTermJacob(RIGHT,_stAnJacobIDs)) {
    addAnalyticSourceTermJacobCN(RIGHT, _acc.get());
  }
  
  // add the values in the jacobian matrix
  _lss->getMatrix()->addValues(*_acc);
  
  // reset to zero the entries in the block accumulator
  _acc->reset(); 
  _sourceJacobOnCell[LEFT] = _sourceJacobOnCell[RIGHT] = false;
}

//////////////////////////////////////////////////////////////////////////////

void FVMCC_CrankNichLimComputeRhs::addBothJacobTermsCN(CFuint iVar, CFuint iCell)
{
  //multiply by the residual factor
  _fluxDiff *= _upFactorVec[0]; 
  
  _acc->addValues(0, iCell, iVar, &_fluxDiff[0]);
  
  // flux is opposite in sign for the other state
  _fluxDiff *= _upFactorVec[1];
  
  _acc->addValues(1, iCell, iVar, &_fluxDiff[0]);
  
  if (computeSourceTermJacob(iCell,_stNumJacobIDs)) {
    _sourceDiffSum = 0.0;
    addSourceTermNumJacob(_currFace->getNeighborGeo(iCell), iCell);
    _sourceDiffSum *= _upStFactorVec[iCell];
    _acc->addValues(iCell, iCell, iVar, &_sourceDiffSum[0]);
  }
}

//////////////////////////////////////////////////////////////////////////////

void FVMCC_CrankNichLimComputeRhs::computeJacobTerm(const CFuint idx)
{
  cf_assert(_currFace->getState(idx)->isParUpdatable());

  // set on the perturbation flag
  getMethodData().setIsPerturb(true);
  
  // set the index of the block corresponding to the current
  // state in the jacobian matrix
  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  const bool isAxi = getMethodData().isAxisymmetric();
  
  _acc->setRowColIndex(0, _currFace->getState(0)->getLocalID());
  _acc->setRowColIndex(1, _currFace->getState(1)->getLocalID());
  
  // first node contribution
  for (CFuint iCell =0 ; iCell < 2; ++iCell) {
    for(CFuint iEq=0; iEq<nbEqs; ++iEq){
      _upStFactorVec[iCell][iEq] = (!isAxi) ? -getResFactorCN(iEq) :-getResFactorCN(iEq)*_invr[iCell];
    }
    
    State& state = *_currFace->getState(iCell);
    
    for (CFuint iVar = 0; iVar < nbEqs; ++iVar) {
      CFLogDebugMax( "Perturbing iVar = " << iVar << "\n");
      
      _numericalJacob->perturb(iVar, state[iVar]);
      
      for(CFuint iEq=0; iEq<nbEqs; ++iEq){
	const CFreal factor = pow(-1.,static_cast<CFreal>(idx))*getResFactorCN(iEq);
	_upFactorVec[idx][iEq] = (!isAxi) ? factor : factor*_rMid*_invr[idx];
      }
      
      computeConvDiffFluxes(iVar, iCell);
      addJacobTermCN(idx, iVar, iCell, _acc.get());
      
      // restore the unperturbed value
      _numericalJacob->restore(state[iVar]);
      _polyRec->restoreValues(iVar, iCell);
    }
  }
  
  // compute analytical jacobian for source term 
  if (computeSourceTermJacob(idx,_stAnJacobIDs)) {
    addAnalyticSourceTermJacobCN(idx, _acc.get());
  }
  
  // add the values in the jacobian matrix
  _lss->getMatrix()->addValues(*_acc);
  
  // reset to zero the entries in the block accumulator
  _acc->reset();
  _sourceJacobOnCell[idx] = false;
}      

//////////////////////////////////////////////////////////////////////////////

void FVMCC_CrankNichLimComputeRhs::addJacobTermCN(CFuint idx, CFuint iVar, CFuint iCell,
						  BlockAccumulator *const acc)
{ 
  //multiply by the residual factor
  _fluxDiff *= _upFactorVec[idx]; 
  
  acc->addValues(idx, iCell, iVar, &_fluxDiff[0]);
  
  if (computeSourceTermJacob(idx, iCell,_stNumJacobIDs)) {
    _sourceDiffSum = 0.0;
    addSourceTermNumJacob(_currFace->getNeighborGeo(idx), idx);
    _sourceDiffSum *= _upStFactorVec[idx];
    acc->addValues(iCell, iCell, iVar, &_sourceDiffSum[0]);
  }  
}
      
//////////////////////////////////////////////////////////////////////////////

void FVMCC_CrankNichLimComputeRhs::computeBoundaryJacobianTerm()
{
  // set the index of the block corresponding to the current
  // state in the jacobian matrix
  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  State& currState = *_currFace->getState(0);
  State& ghostState = *_currFace->getState(1);
  cf_assert(ghostState.isGhost());
  
  if (currState.isParUpdatable()) {
    const bool isAxi = getMethodData().isAxisymmetric();
    // copy the original value of the ghost state
    _origState = ghostState;

    _bAcc->setRowColIndex(0, currState.getLocalID());

    for (CFuint iVar = 0; iVar < nbEqs; ++iVar) {
      CFLogDebugMax( "Perturbing iVar = " << iVar << "\n");

      // perturb the given component of the state vector
      _numericalJacob->perturb(iVar, currState[iVar]);

      // compute the ghost state in the perturbed inner state
      _currBC->setGhostState(_currFace);

      // extrapolate (and LIMIT, if the reconstruction is linear or more)
      // the solution in the quadrature points
      _polyRec->extrapolate(_currFace);

      // compute the physical data for each left and right reconstructed
      // state and in the left and right cell centers
      computeStatesData();

      // set the perturbed variable
      getMethodData().setIPerturbVar(iVar);
      
      _currBC->computeFlux(_pertFlux);
      
      if (_hasDiffusiveTerm) {
	//_nodalExtrapolator->extrapolateInNodes(*_currFace->getNodes());
        _diffusiveFlux->computeFlux(_dFlux);
        _pertFlux -= _dFlux;
      }

      // compute the finite difference derivative of the flux
      _numericalJacob->computeDerivative(getJacobianFlux(),_pertFlux,_fluxDiff);
    
      for(CFuint iEq=0; iEq<nbEqs; ++iEq){
	    const CFreal factor = getResFactorCN(iEq);
	    _upFactorVec[LEFT][iEq] = (!isAxi) ? factor : factor*(_rMid*_invr[0]);
        _upStFactorVec[LEFT][iEq] = (!isAxi) ? -factor :-factor*_invr[0];
      }

      addJacobTermCN(0, iVar, 0, _bAcc.get());

      // restore the unperturbed value
      _numericalJacob->restore(currState[iVar]);

      // restore the original ghost state
      ghostState = _origState;
    }
    
    // compute analytical jacobian for source term 
     if (computeSourceTermJacob(LEFT,_stAnJacobIDs)) {
       addAnalyticSourceTermJacob(LEFT, _bAcc.get());
     }
    
    // add the values in the jacobian matrix
    _lss->getMatrix()->addValues(*_bAcc);

    // reset to zero the entries in the block accumulator
    _bAcc->reset();
    _sourceJacobOnCell[LEFT] = false;
  }
}

//////////////////////////////////////////////////////////////////////////////

} // namespace FiniteVolume

} // namespace Numerics

} // namespace COOLFluiD
