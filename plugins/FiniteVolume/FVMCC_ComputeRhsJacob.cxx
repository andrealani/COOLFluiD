#include "Framework/CFL.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/LSSMatrix.hh"
#include "Framework/BlockAccumulator.hh"
#include "Framework/MeshData.hh"

#include "FiniteVolume/FiniteVolume.hh"
#include "FiniteVolume/FVMCC_ComputeRhsJacob.hh"
#include "FiniteVolume/FVMCC_BC.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<FVMCC_ComputeRhsJacob,
		      CellCenterFVMData,
		      FiniteVolumeModule>
fvmcc_computeRhsJacob("NumJacob");

//////////////////////////////////////////////////////////////////////////////

FVMCC_ComputeRhsJacob::FVMCC_ComputeRhsJacob(const std::string& name) :
  FVMCC_ComputeRHS(name),
  _lss(CFNULL),
  _numericalJacob(CFNULL),
  _pertFlux(),
  _fluxDiff(),
  _origState(),
  _acc(CFNULL),
  _bAcc(CFNULL),
  _pertSource(),
  _sourceDiff(),
  _sourceDiffSum(),
  _dummyJacob()
{
  addConfigOptionsTo(this);
}

//////////////////////////////////////////////////////////////////////////////

FVMCC_ComputeRhsJacob::~FVMCC_ComputeRhsJacob()
{
}

//////////////////////////////////////////////////////////////////////////////

void FVMCC_ComputeRhsJacob::defineConfigOptions(Config::OptionList& options)
{
}

//////////////////////////////////////////////////////////////////////////////

void FVMCC_ComputeRhsJacob::configure ( Config::ConfigArgs& args )
{
  FVMCC_ComputeRHS::configure(args);
}

//////////////////////////////////////////////////////////////////////////////

void FVMCC_ComputeRhsJacob::setup()
{
  FVMCC_ComputeRHS::setup();

  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  _pertFlux.resize(nbEqs);
  _fluxDiff.resize(nbEqs);
  _origState.resize(nbEqs);
  
  const CFuint nbSourceTerms = _stComputers->size();
  _pertSource.resize(nbSourceTerms);
  for (CFuint i = 0; i < nbSourceTerms; ++i) {
    _pertSource[i].resize(nbEqs);
    _pertSource[i] = 0.0;
  }
  
  _sourceDiff.resize(nbEqs);
  _sourceDiffSum.resize(nbEqs);
  _dummyJacob.resize(nbEqs,nbEqs);
  
  // linear system solver
  _lss = getMethodData().getLinearSystemSolver()[0];

  // numerical jacobian
  _numericalJacob = &getMethodData().getNumericalJacobian();

  _acc.reset(_lss->createBlockAccumulator(2, 2, nbEqs));
  _bAcc.reset(_lss->createBlockAccumulator(1, 1, nbEqs));
}

//////////////////////////////////////////////////////////////////////////////

void FVMCC_ComputeRhsJacob::computeJacobianTerm()
{ 
  CFLog(DEBUG_MIN, "FVMCC_ComputeRhsJacob::computeJacobianTerm() START\n"); 
  
  const bool upLeft = _currFace->getState(0)->isParUpdatable();
  const bool upRight = _currFace->getState(1)->isParUpdatable();
  
  if (upLeft && upRight) {
    computeBothJacobTerms();
  }
  else if (upLeft) {
    computeJacobTerm(0);
  }
  else if (upRight) {
    computeJacobTerm(1);
  }
  
  CFLog(DEBUG_MIN, "FVMCC_ComputeRhsJacob::computeJacobianTerm() END\n");
}
      
//////////////////////////////////////////////////////////////////////////////

void FVMCC_ComputeRhsJacob::computeBothJacobTerms()
{
  //cout.precision(14); cout.setf(ios::scientific,ios::floatfield); cout << "_flux = " << _flux << endl;
  
  // set on the perturbation flag
  getMethodData().setIsPerturb(true);

  // set the index of the block corresponding to the current
  // state in the jacobian matrix
  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  (!getMethodData().isAxisymmetric()) ? computeNoAxiUpFactors() : computeAxiUpFactors();
  
  _acc->setRowColIndex(0, _currFace->getState(0)->getLocalID());
  _acc->setRowColIndex(1, _currFace->getState(1)->getLocalID());
  
  // first node contribution
  for (CFuint iCell = 0; iCell < 2; ++iCell) {
    State& state = *_currFace->getState(iCell); 
    cf_assert(!state.isGhost());
    
    for (CFuint iVar = 0; iVar < nbEqs; ++iVar) {
      CFLogDebugMax("Perturbing iVar = " << iVar << "\n");
    
      // set the perturbed variable
      getMethodData().setIPerturbVar(iVar);
      
      // perturb the given component of the state vector
      _numericalJacob->perturb(iVar, state[iVar]); 
      
      computeConvDiffFluxes(iVar, iCell);  
      addBothJacobTerms(iVar, iCell);
      
      // restore the unperturbed value
      _numericalJacob->restore(state[iVar]);
      _polyRec->restoreValues(iVar, iCell);
    }
    
    // restoring right states is NOT needed because once you get out of here a new face is processed
    if (iCell == 0) restoreState(0);
    
    // compute analytical jacobian for source term 
    if (computeSourceTermJacob(iCell,_stAnJacobIDs)) {
      addAnalyticSourceTermJacob(iCell, _acc.get());
    }
  }
    
  // add the values in the jacobian matrix
  _lss->getMatrix()->addValues(*_acc);
  
  // _acc->print(); EXIT_AT(1);
  
  // reset to zero the entries in the block accumulator
  _acc->reset(); 
  _sourceJacobOnCell[LEFT] = _sourceJacobOnCell[RIGHT] = false;
}

//////////////////////////////////////////////////////////////////////////////

void FVMCC_ComputeRhsJacob::computeJacobTerm(const CFuint idx)
{
  cf_assert(_currFace->getState(idx)->isParUpdatable());

  // set on the perturbation flag
  getMethodData().setIsPerturb(true);
  
  // set the index of the block corresponding to the current
  // state in the jacobian matrix
  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  const bool isAxi = getMethodData().isAxisymmetric();
  const CFreal factor = pow(-1.,static_cast<CFreal>(idx))*getResFactor();
  _upFactor[idx] = (!isAxi) ? factor : factor*_rMid*_invr[idx];
  
  _acc->setRowColIndex(0, _currFace->getState(0)->getLocalID());
  _acc->setRowColIndex(1, _currFace->getState(1)->getLocalID());
  
  // first node contribution
  for (CFuint iCell = 0; iCell < 2; ++iCell) {
    _upStFactor[iCell] = (!isAxi) ? -getResFactor() :-getResFactor()*_invr[iCell];
    State& state = *_currFace->getState(iCell);
    
    for (CFuint iVar = 0; iVar < nbEqs; ++iVar) {
      CFLogDebugMax( "Perturbing iVar = " << iVar << "\n");
      
      // set the perturbed variable
      getMethodData().setIPerturbVar(iVar);
      
      _numericalJacob->perturb(iVar, state[iVar]);
      computeConvDiffFluxes(iVar, iCell);
      addJacobTerm(idx, iVar, iCell, _acc.get());
      
      // restore the unperturbed value
      _numericalJacob->restore(state[iVar]);
      _polyRec->restoreValues(iVar, iCell);
    }
    
    // restoring right states is NOT needed because once you get out of here a new face is processed
    if (iCell == 0) restoreState(0);
  }
  
  // compute analytical jacobian for source term 
  if (computeSourceTermJacob(idx,_stAnJacobIDs)) {
    addAnalyticSourceTermJacob(idx, _acc.get());
  }
  
  // add the values in the jacobian matrix
  _lss->getMatrix()->addValues(*_acc);
  
  // reset to zero the entries in the block accumulator
  _acc->reset();
  _sourceJacobOnCell[idx] = false;
}      

//////////////////////////////////////////////////////////////////////////////

void FVMCC_ComputeRhsJacob::computeBoundaryJacobianTerm()
{
  getMethodData().setIsPerturb(true);

  // set the index of the block corresponding to the current
  // state in the jacobian matrix
  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  State& currState = *_currFace->getState(0);
  State& ghostState = *_currFace->getState(1);
  cf_assert(ghostState.isGhost());

  if (currState.isParUpdatable()) {
    const bool isAxi = getMethodData().isAxisymmetric();
    _upFactor[LEFT] = (!isAxi) ? getResFactor() : getResFactor()*(_rMid*_invr[0]);
    _upStFactor[LEFT] = (!isAxi) ? -getResFactor() :-getResFactor()*_invr[0];
    
    // copy the original value of the ghost state
    _origState = ghostState;
    
    _bAcc->setRowColIndex(0, currState.getLocalID());
    
    for (CFuint iVar = 0; iVar < nbEqs; ++iVar) {
      CFLogDebugMax( "Perturbing iVar = " << iVar << "\n");
   
      // set the perturbed variable
      getMethodData().setIPerturbVar(iVar);
      
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

      _pertFlux = 0.;
      _currBC->computeFlux(_pertFlux);

      if (_hasDiffusiveTerm && _isDiffusionActive) {
	//_nodalExtrapolator->extrapolateInNodes(*_currFace->getNodes());
        _diffusiveFlux->computeFlux(_dFlux);
        _pertFlux -= _dFlux;
      }

      // compute the finite difference derivative of the flux
      _numericalJacob->computeDerivative(getJacobianFlux(),_pertFlux,_fluxDiff);

      addJacobTerm(0, iVar, 0, _bAcc.get());

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
    // cout << "BAC" << endl;_bAcc->print();

    // reset to zero the entries in the block accumulator
    _bAcc->reset();
    _sourceJacobOnCell[LEFT] = false;
  } 
  
  // cout << "BAC" << endl;_bAcc->print();
}

//////////////////////////////////////////////////////////////////////////////

void FVMCC_ComputeRhsJacob::computeConvDiffFluxes(CFuint iVar, CFuint iCell)
{  
  // extrapolate (and LIMIT, if the reconstruction is linear or more)
  // the solution in the quadrature points
  _polyRec->extrapolate(_currFace, iVar, iCell);
  
  // compute the physical data for each left and right reconstructed
  // state and in the left and right cell centers
  computeStatesData();
  
  // linearization will be done in the flux splitter if needed
  _pertFlux = 0.;
  _fluxSplitter->computeFlux(_pertFlux);
  
  if (_hasDiffusiveTerm && _isDiffusionActive) {
    //_nodalExtrapolator->extrapolateInNodes(*_currFace->getNodes());
    _diffusiveFlux->computeFlux(_dFlux);
    _pertFlux -= _dFlux;
  }
  
  // compute the finite difference derivative of the flux
  _numericalJacob->computeDerivative(getJacobianFlux(),_pertFlux,_fluxDiff);
}

//////////////////////////////////////////////////////////////////////////////

void FVMCC_ComputeRhsJacob::addBothJacobTerms(CFuint iVar, CFuint iCell)
{
  //multiply by the residual factor
  _fluxDiff *= _upFactor[0]; 
  
  _acc->addValues(0, iCell, iVar, &_fluxDiff[0]);
  
  // flux is opposite in sign for the other state
  _fluxDiff *= _upFactor[1];
  
  _acc->addValues(1, iCell, iVar, &_fluxDiff[0]);

  if (computeSourceTermJacob(iCell,_stNumJacobIDs)) {
    _sourceDiffSum = 0.0;
    addSourceTermNumJacob(_currFace->getNeighborGeo(iCell), iCell);
    _sourceDiffSum *= _upStFactor[iCell];
    _acc->addValues(iCell, iCell, iVar, &_sourceDiffSum[0]);
  }
  

  
  // if (computeSourceTermJacob(iCell,_stNumJacobIDs) && _currFace->getNeighborGeo(iCell)->getState(0)->getLocalID() == 100) {
//     // // output source term jacobian on screen
//     static CFuint nbEqs = _fluxDiff.size();
//     static CFuint count = 0;
//     static BlockAccumulator* b = _lss->createBlockAccumulator(1, 1, nbEqs);
    
//     //   //   // debugging
//     if (count < nbEqs && iCell == 1) { 
//       if (computeSourceTermJacob(1,_stNumJacobIDs)) {
// 	b->addValues(0, 0, iVar, &_sourceDiffSum[0]);
// 	count++;
//       }
//     }
    
//     if (count == nbEqs && iCell == 1) {
//       b->print();
//       count = 0;
//     }
//   }  
 
}
      
//////////////////////////////////////////////////////////////////////////////

void FVMCC_ComputeRhsJacob::addJacobTerm(CFuint idx, CFuint iVar, CFuint iCell,
					 BlockAccumulator *const acc)
{ 
  //multiply by the residual factor
  _fluxDiff *= _upFactor[idx]; 
  
  acc->addValues(idx, iCell, iVar, &_fluxDiff[0]);
  
  if (computeSourceTermJacob(idx, iCell,_stNumJacobIDs)) {
    _sourceDiffSum = 0.0;
    addSourceTermNumJacob(_currFace->getNeighborGeo(idx), idx);
    _sourceDiffSum *= _upStFactor[idx];
    acc->addValues(iCell, iCell, iVar, &_sourceDiffSum[0]);
  }  
}

//////////////////////////////////////////////////////////////////////////////

void FVMCC_ComputeRhsJacob::initializeComputationRHS()
{
  // reset rhs to 0
  socket_rhs.getDataHandle() = 0.0;
  socket_cellFlag.getDataHandle() = false;
  
  // reset to zero all non zero entries in the jacobian
  if (getMethodData().doComputeJacobian()) {
    _lss->getMatrix()->resetToZeroEntries();
  }
  
  // try to see if this fixes the LS with Roe
  _polyRec->updateWeights();
  
  // gradients and limiters are computed on all the variables at once
  _polyRec->computeGradients();
  
  // extrapolate the solution from cell centers to all vertices
  _nodalExtrapolator->extrapolateInAllNodes();
}
      
//////////////////////////////////////////////////////////////////////////////
 
void FVMCC_ComputeRhsJacob::computeRHSJacobian()
{
  if (getMethodData().doComputeJacobian()) {
    _diffVar->setFreezeCoeff(_freezeDiffCoeff);
    const bool isBFace = _currFace->getState(1)->isGhost();
    (!isBFace) ? computeJacobianTerm() : computeBoundaryJacobianTerm();
  }
  
  // set off the perturbation flag 
  getMethodData().setIsPerturb(false);
}
      
//////////////////////////////////////////////////////////////////////////////

void FVMCC_ComputeRhsJacob::finalizeComputationRHS()
{
  // reset the flag for freezing the transport properties
  _diffVar->setFreezeCoeff(false);
}

//////////////////////////////////////////////////////////////////////////////
 
} // namespace FiniteVolume

} // namespace Numerics

} // namespace COOLFluiD
