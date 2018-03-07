#include "FiniteVolume/FiniteVolume.hh"
#include "FVMCC_ComputeRhsJacobCoupling.hh"
#include "FiniteVolume/FVMCC_BC.hh"
#include "Framework/CFL.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/LSSMatrix.hh"
#include "Framework/MeshData.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<FVMCC_ComputeRhsJacobCoupling, CellCenterFVMData, FiniteVolumeModule>
FVMCC_computeRhsJacobCoupling("NumJacobCoupling");

//////////////////////////////////////////////////////////////////////////////

FVMCC_ComputeRhsJacobCoupling::FVMCC_ComputeRhsJacobCoupling
(const std::string& name) :
  FVMCC_ComputeRHS(name),
  _numericalJacob(CFNULL),
  _nbEqs(0),
  _start(0),
  _iLSS(0),
  _pertFlux(),
  _fluxDiff(),
  _origState(),
  _lss(),
  _acc(),
  _bAcc(),
  _equations(),
  _pertSource(),
  _sourceDiff(),
  _sourceDiffSum(),
  _dummyJacob(),
  _tmpAnalytMatrix()
{
}

//////////////////////////////////////////////////////////////////////////////

FVMCC_ComputeRhsJacobCoupling::~FVMCC_ComputeRhsJacobCoupling()
{
  for (CFuint i = 0; i < _acc.size(); ++i) {
    deletePtr(_acc[i]);
  }

  for (CFuint i = 0; i < _bAcc.size(); ++i) {
    deletePtr(_bAcc[i]);
  }
}

//////////////////////////////////////////////////////////////////////////////

void FVMCC_ComputeRhsJacobCoupling::setup()
{
  FVMCC_ComputeRHS::setup();

  _pertFlux.resize(PhysicalModelStack::getActive()->getNbEq());
  _fluxDiff.resize(PhysicalModelStack::getActive()->getNbEq());
  _origState.resize(PhysicalModelStack::getActive()->getNbEq());

  // numerical jacobian
  _numericalJacob = &getMethodData().getNumericalJacobian();

  const CFuint nbLSS = getMethodData().getLinearSystemSolver().size();

  // set the total number of equation subsystems
  PhysicalModelStack::getActive()->setTotalNbEqSS(nbLSS);

  // acquaintance of the linear systems
  _lss.resize(nbLSS);
  for (CFuint i = 0; i < nbLSS; ++i) {
    _lss[i] = getMethodData().getLinearSystemSolver()[i];
  }

  _acc.resize(nbLSS);
  for (CFuint i = 0; i < nbLSS; ++i) {
    const CFuint nbSysEq = _lss[i]->getNbSysEqs();
    _acc[i] = _lss[i]->createBlockAccumulator(2, 2, nbSysEq);
  }
  
  _bAcc.resize(nbLSS);
  for (CFuint i = 0; i < nbLSS; ++i) {
    const CFuint nbSysEq = _lss[i]->getNbSysEqs();
    _bAcc[i] = _lss[i]->createBlockAccumulator(1, 1, nbSysEq);
  }

  _equations.resize(nbLSS);
  for (CFuint i = 0; i < nbLSS; ++i) {
    // acquaintance of the equation IDs to solve in each LSS
    _equations[i] = _lss[i]->getEquationIDs();
    cf_assert(_equations[i]->size() == _lss[i]->getNbSysEqs());
  }

  vector<vector<CFuint> > eqVarPatterns(nbLSS);
  for (CFuint i = 0; i < nbLSS; ++i) {
    eqVarPatterns[i] = *_equations[i];
    cf_assert(eqVarPatterns[i].size() == _equations[i]->size());
  }

  PhysicalModelStack::getActive()->setEquationVarPatterns(eqVarPatterns);
  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  
  const CFuint nbSourceTerms = _stComputers->size();
  _pertSource.resize(nbSourceTerms);
  for (CFuint i = 0; i < nbSourceTerms; ++i) {
    _pertSource[i].resize(nbEqs);
    _pertSource[i] = 0.0;
  }
  
  _sourceDiff.resize(nbEqs);
  _sourceDiffSum.resize(nbEqs);
  _dummyJacob.resize(nbEqs,nbEqs);
  
  // temporary storage for analytical matrices
  _tmpAnalytMatrix.resize(nbLSS);
  for (CFuint i = 0; i < nbLSS; ++i) {
    const CFuint nbSysEq = _lss[i]->getNbSysEqs();
    _tmpAnalytMatrix[i].resize(nbSysEq, nbSysEq);
    _tmpAnalytMatrix[i] = 0.0;
  } 
  
  // initialization
  _nbEqs = nbEqs;
  _start = 0;
  _iLSS = 0;
  PhysicalModelStack::getActive()->setEquationSubSysDescriptor(_start, _nbEqs, 0);
}
      
//////////////////////////////////////////////////////////////////////////////

void FVMCC_ComputeRhsJacobCoupling::computeJacobianTerm()
{
  CFLog(DEBUG_MED, "FVMCC_ComputeRhsJacobCoupling::computeJacobianTerm() START\n");
  
  const CFuint nbLSS = _lss.size();
  const bool upLeft  = _currFace->getState(0)->isParUpdatable();
  const bool upRight = _currFace->getState(1)->isParUpdatable();
  
  if (upLeft && upRight) {
    // the current index _iLSS of the LSS is set here not to have to carry it around
    for (_iLSS = 0; _iLSS < nbLSS; ++_iLSS) {
      setCurrentSubsystemData(*_equations[_iLSS]);
      computeBothJacobTerms();
    }
  }
  else if (upLeft) {
    // the current index _iLSS of the LSS is set here not to have to carry it around 
    for (_iLSS = 0; _iLSS < nbLSS; ++_iLSS) {
      setCurrentSubsystemData(*_equations[_iLSS]);
      computeJacobTerm(0);
    }
  }
  else if (upRight) {
    // the current index _iLSS of the LSS is set here not to have to carry it around 
    for (_iLSS = 0; _iLSS < nbLSS; ++_iLSS) {
      setCurrentSubsystemData(*_equations[_iLSS]);
      computeJacobTerm(1);
    }
  } 
  
  CFLog(DEBUG_MED, "FVMCC_ComputeRhsJacobCoupling::computeJacobianTerm() END\n");
}

//////////////////////////////////////////////////////////////////////////////

void FVMCC_ComputeRhsJacobCoupling::computeBothJacobTerms()
{
  CFLog(DEBUG_MED, "FVMCC_ComputeRhsJacobCoupling::computeBothJacobTerms() START\n");
  //cout.precision(14); cout.setf(ios::scientific,ios::floatfield); cout << "_flux = " << _flux << endl;
  
  // set on the perturbation flag
  getMethodData().setIsPerturb(true);

  // set the index of the block corresponding to the current
  // state in the jacobian matrix
  (!getMethodData().isAxisymmetric()) ? computeNoAxiUpFactors() : computeAxiUpFactors();
  
  // set the index of the block corresponding to the current
  // state in the jacobian matrix  
  BlockAccumulator& acc = *_acc[_iLSS];
  SafePtr<vector<CFuint> > equations = _equations[_iLSS];
  
  // contribution to the jacobian of the first LSS
  // set row and column index for the first block accumulator (LSS 0)
  acc.setRowColIndex(0, _currFace->getState(0)->getLocalID());
  acc.setRowColIndex(1, _currFace->getState(1)->getLocalID());
  
  // first node contribution
  for (CFuint iCell = 0; iCell < 2; ++iCell) {
    State& state = *_currFace->getState(iCell); 
    cf_assert(!state.isGhost());
    
    for (CFuint iEq = 0; iEq < _nbEqs; ++iEq) {
      // variable to perturb
      const CFuint iVar = (*equations)[iEq];
      cf_assert(iEq <= iVar);
      
      CFLogDebugMax("Perturbing iVar = " << iVar << "\n");
    
      // set the perturbed variable
      getMethodData().setIPerturbVar(iVar);
      
      // perturb the given component of the state vector
      _numericalJacob->perturb(iVar, state[iVar]); 
      
      computeConvDiffFluxes(iVar, iCell);  
      // iEq is the equation ID within the current equation subsystem
      addBothJacobTerms(iEq, iCell);
      
      // restore the unperturbed value
      _numericalJacob->restore(state[iVar]);
      _polyRec->restoreValues(iVar, iCell);
    }
    
    // restoring right states is NOT needed because once you get out of here a new face is processed
    if (iCell == 0) restoreState(0);
    
    // compute analytical jacobian for source term 
    if (computeSourceTermJacob(iCell,_stAnJacobIDs)) {
      addAnalyticSourceTermJacob(iCell, &acc);
    }
  }
  
  // add the values in the jacobian matrix
  _lss[_iLSS]->getMatrix()->addValues(acc);
 
  // acc.print(); EXIT_AT(1);
  
  // reset to zero the entries in the block accumulator
  acc.reset();
  
  if (_iLSS == _lss.size() - 1) {
    _sourceJacobOnCell[LEFT] = _sourceJacobOnCell[RIGHT] = false;
  }
  
  CFLog(DEBUG_MED, "FVMCC_ComputeRhsJacobCoupling::computeBothJacobTerms() END\n");
}

//////////////////////////////////////////////////////////////////////////////

void FVMCC_ComputeRhsJacobCoupling::computeJacobTerm(CFuint idx)
{
  CFLog(DEBUG_MED, "FVMCC_ComputeRhsJacobCoupling::computeJacobTerm() START\n");
  
  cf_assert(_currFace->getState(idx)->isParUpdatable());
  
  // set on the perturbation flag
  getMethodData().setIsPerturb(true);
  
  BlockAccumulator& acc = *_acc[_iLSS];
  SafePtr<vector<CFuint> > equations = _equations[_iLSS];
  
  // set the index of the block corresponding to the current
  // state in the jacobian matrix
  const bool isAxi = getMethodData().isAxisymmetric();
  const CFreal factor = pow(-1.,static_cast<CFreal>(idx))*getResFactor();
  _upFactor[idx] = (!isAxi) ? factor : factor*_rMid*_invr[idx];
  
  acc.setRowColIndex(0, _currFace->getState(0)->getLocalID());
  acc.setRowColIndex(1, _currFace->getState(1)->getLocalID());
  
  // first node contribution
  for (CFuint iCell = 0; iCell < 2; ++iCell) {
    _upStFactor[iCell] = (!isAxi) ? -getResFactor() :-getResFactor()*_invr[iCell];
    State& state = *_currFace->getState(iCell);
    
    for (CFuint iEq = 0; iEq < _nbEqs; ++iEq) {
      // variable to perturb
      const CFuint iVar = (*equations)[iEq];
      cf_assert(iEq <= iVar);
      CFLogDebugMax( "Perturbing iVar = " << iVar << "\n");
      
      // set the perturbed variable
      getMethodData().setIPerturbVar(iVar);

      // perturb the given component of the state vector
      _numericalJacob->perturb(iVar, state[iVar]);
      computeConvDiffFluxes(iVar, iCell);
      // iEq is the equation ID within the current equation subsystem
      addJacobTerm(idx, iEq, iCell, &acc);
      
      // restore the unperturbed value
      _numericalJacob->restore(state[iVar]);
      _polyRec->restoreValues(iVar, iCell);
    }
  }
  
  // compute analytical jacobian for source term 
  if (computeSourceTermJacob(idx,_stAnJacobIDs)) {
    addAnalyticSourceTermJacob(idx, &acc);
  }
  
  // add the values in the jacobian matrix
  _lss[_iLSS]->getMatrix()->addValues(acc);
  
  // reset to zero the entries in the block accumulator
  acc.reset();
  
  if (_iLSS == _lss.size() - 1) {
    _sourceJacobOnCell[idx] = false;
  }
  
  CFLog(DEBUG_MED, "FVMCC_ComputeRhsJacobCoupling::computeJacobTerm() END\n");
}

//////////////////////////////////////////////////////////////////////////////

void FVMCC_ComputeRhsJacobCoupling::computeBoundaryJacobianTerm()
{
 CFLog(DEBUG_MED, "FVMCC_ComputeRhsJacobCoupling::computeBoundaryJacobianTerm() START\n");
 
 getMethodData().setIsPerturb(true);
 
 const CFuint nbLSS = _lss.size();
 
 // the current index of the LSS is set here not to have to carry it around
 for (_iLSS = 0; _iLSS < nbLSS; ++_iLSS) {
   setCurrentSubsystemData(*_equations[_iLSS]);    
   
   // set the index of the block corresponding to the current
   // state in the jacobian matrix
   State& currState = *_currFace->getState(0);
   State& ghostState = *_currFace->getState(1);
   cf_assert(ghostState.isGhost());
   
   BlockAccumulator& bAcc = *_bAcc[_iLSS];
   SafePtr<LinearSystemSolver> lss = _lss[_iLSS];
   // we work with a slice [_start, _start+_nbEqs)
   SafePtr<vector<CFuint> > equations = _equations[_iLSS];
   
   if (currState.isParUpdatable()) {
     const bool isAxi = getMethodData().isAxisymmetric();
     _upFactor[LEFT] = (!isAxi) ? getResFactor() : getResFactor()*(_rMid*_invr[0]);
     _upStFactor[LEFT] = (!isAxi) ? -getResFactor() :-getResFactor()*_invr[0];

     // copy the original value of the ghost state
     _origState.slice(_start, _nbEqs) = ghostState.slice(_start, _nbEqs);
     
     bAcc.setRowColIndex(0, currState.getLocalID());
     
     for (CFuint iEq = 0; iEq < _nbEqs; ++iEq) {
       // variable to perturb
       const CFuint iVar = (*equations)[iEq];
       cf_assert(iEq <= iVar);
       
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
       
       // AL: this initialization is fundamental, especially for cases with coupling
       // where some equation subsystems don't have convective terms
        
       // at this level, we are working only inside each subsystem
       // we first reset to 0 a flux slice and then compute convective + diffusive term for that slice
       // it can happen that the last slices of the flux array are != 0 even if the convective part is 0:
       // this is due to the previous diffusive contribution to those slices     
       _pertFlux.slice(_start, _nbEqs) = 0.;
       _currBC->computeFlux(_pertFlux);
       
       if (_hasDiffusiveTerm) {
	 if (_extrapolateInNodes) {
	   _nodalExtrapolator->extrapolateInNodes(*_currFace->getNodes());
	 }
	 
	 // this initialization is fundamental, especially for cases with coupling
	 // where some equation subsystems don't have diffusive terms
	 _dFlux.slice(_start, _nbEqs) = 0.;
	 _diffusiveFlux->computeFlux(_dFlux);
	 _pertFlux.slice(_start, _nbEqs) -= _dFlux.slice(_start, _nbEqs);
       }
       
       // compute the finite difference derivative of the flux
       _numericalJacob->computeDerivative(getJacobianFlux().slice(_start, _nbEqs),
					  _pertFlux.slice(_start, _nbEqs),
					  _fluxDiff.slice(_start, _nbEqs));
       
       addJacobTerm(0, iEq, 0, &bAcc);
       
       // restore the unperturbed value
       _numericalJacob->restore(currState[iVar]);
       
       // restore the original ghost state
       ghostState.slice(_start, _nbEqs) = _origState.slice(_start, _nbEqs);
     }
     
     // compute analytical jacobian for source term 
     if (computeSourceTermJacob(LEFT,_stAnJacobIDs)) {
       addAnalyticSourceTermJacob(LEFT, &bAcc);
     }
     
     // add the values in the jacobian matrix
     lss->getMatrix()->addValues(bAcc);
   
     // bAcc.print();
     
     // static int count = 0; if (++count == 100) abort();
     
     // reset to zero the entries in the block accumulator
     bAcc.reset();
     
     if (_iLSS == nbLSS - 1) {
       _sourceJacobOnCell[LEFT] = false;
     }
   }
 }
 
 CFLog(DEBUG_MED, "FVMCC_ComputeRhsJacobCoupling::computeBoundaryJacobianTerm() END\n");
}
      
//////////////////////////////////////////////////////////////////////////////

void FVMCC_ComputeRhsJacobCoupling::computeConvDiffFluxes(CFuint iVar, 
							  CFuint iCell)
{  
  // extrapolate (and LIMIT, if the reconstruction is linear or more)
  // the solution in the quadrature points
  _polyRec->extrapolate(_currFace, iVar, iCell);
  
  // compute the physical data for each left and right reconstructed
  // state and in the left and right cell centers
  // here not all the quantities should be computed ...
  computeStatesData();
  
  // compute the flux corresponding to the first _nbEqs components
  _pertFlux.slice(_start, _nbEqs) = 0.;
  _fluxSplitter->computeFlux(_pertFlux);

  if (_hasDiffusiveTerm) {
    if (_extrapolateInNodes) {
      _nodalExtrapolator->extrapolateInNodes(*_currFace->getNodes());
    }
    
    _dFlux.slice(_start, _nbEqs) = 0.;
    _diffusiveFlux->computeFlux(_dFlux); 
    
    // modify only the first _nbEqs components
    _pertFlux.slice(_start, _nbEqs) -= _dFlux.slice(_start, _nbEqs);
  }
  
  // compute the finite difference derivative of the flux
  _numericalJacob->computeDerivative(getJacobianFlux().slice(_start, _nbEqs),
				     _pertFlux.slice(_start, _nbEqs),
				     _fluxDiff.slice(_start, _nbEqs));
}

//////////////////////////////////////////////////////////////////////////////

void FVMCC_ComputeRhsJacobCoupling::addBothJacobTerms(CFuint iVar, CFuint iCell)
{
  cf_assert(iVar < _nbEqs);
  //multiply by the residual factor
  _fluxDiff.slice(_start, _nbEqs) *= _upFactor[0];

  BlockAccumulator& acc = *_acc[_iLSS];
  acc.addValues(0, iCell, iVar, &_fluxDiff[_start]);
  
  // flux is opposite in sign for the other state
  _fluxDiff.slice(_start, _nbEqs) *= _upFactor[1];
  
  acc.addValues(1, iCell, iVar, &_fluxDiff[_start]);

  if (computeSourceTermJacob(iCell,_stNumJacobIDs)) {
    _sourceDiffSum.slice(_start, _nbEqs) = 0.0;
    addSourceTermNumJacob(_currFace->getNeighborGeo(iCell), iCell);
    _sourceDiffSum.slice(_start, _nbEqs) *= _upStFactor[iCell];
    acc.addValues(iCell, iCell, iVar, &_sourceDiffSum[_start]);
  }
}
      
//////////////////////////////////////////////////////////////////////////////

void FVMCC_ComputeRhsJacobCoupling::addJacobTerm(CFuint idx, CFuint iVar, CFuint iCell,
						 BlockAccumulator *const acc)
{ 
  cf_assert(iVar < _nbEqs);
  
  //multiply by the residual factor
  _fluxDiff.slice(_start, _nbEqs) *= _upFactor[idx]; 
  
  acc->addValues(idx, iCell, iVar, &_fluxDiff[_start]);
  
  if (computeSourceTermJacob(idx, iCell,_stNumJacobIDs)) {
    _sourceDiffSum.slice(_start, _nbEqs) = 0.0;
    addSourceTermNumJacob(_currFace->getNeighborGeo(idx), idx);
    _sourceDiffSum.slice(_start, _nbEqs) *= _upStFactor[idx];
    acc->addValues(iCell, iCell, iVar, &_sourceDiffSum[_start]);
  }  
}

//////////////////////////////////////////////////////////////////////////////

void FVMCC_ComputeRhsJacobCoupling::initializeComputationRHS()
{
  // reset rhs to 0
  socket_rhs.getDataHandle() = 0.0;
  socket_cellFlag.getDataHandle() = false;
  
  // reset to zero all non zero entries in the jacobian
  if (getMethodData().doComputeJacobian()) {
    const CFuint nbLSS = _lss.size();
    for (CFuint i = 0; i < nbLSS; ++i) {
      CFLog(VERBOSE, "FVMCC_ComputeRhsJacobCoupling::initializeComputationRHS() => before zero\n");
      _lss[i]->getMatrix()->resetToZeroEntries();
      CFLog(VERBOSE, "FVMCC_ComputeRhsJacobCoupling::initializeComputationRHS() => after zero\n");
    }
  }
  
  // try to see if this fixes the LS with Roe
  _polyRec->updateWeights();
  
  // gradients and limiters are computed on all the variables at once
  _polyRec->computeGradients();
  
  // extrapolate the solution from cell centers to all vertices
  _nodalExtrapolator->extrapolateInAllNodes();
}
      
//////////////////////////////////////////////////////////////////////////////
 
void FVMCC_ComputeRhsJacobCoupling::computeRHSJacobian()
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

void FVMCC_ComputeRhsJacobCoupling::finalizeComputationRHS()
{
  // reset the flag for freezing the transport properties
  _diffVar->setFreezeCoeff(false);
  
  // reset the equation subsystem descriptor
  PhysicalModelStack::getActive()->resetEquationSubSysDescriptor();
}  

//////////////////////////////////////////////////////////////////////////////
 
} // namespace FiniteVolume

} // namespace Numerics

} // namespace COOLFluiD
