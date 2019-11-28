#include "Framework/MethodCommandProvider.hh"
#include "FiniteVolume/FiniteVolume.hh"
#include "FiniteVolume/FVMCC_ComputeRhsBlockJacob.hh"
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

MethodCommandProvider<FVMCC_ComputeRhsBlockJacob,
		      CellCenterFVMData,
		      FiniteVolumeModule>
fvmcc_computeRhsBlockJacob("NumBlockJacob");

//////////////////////////////////////////////////////////////////////////////

FVMCC_ComputeRhsBlockJacob::FVMCC_ComputeRhsBlockJacob(const std::string& name) :
  FVMCC_ComputeRHS(name),
  socket_diagMatrix("diagMatrix"),
  _acc(CFNULL),
  _numericalJacob(CFNULL),
  _pertFlux(),
  _fluxDiff(),
  _origState(),
  _pertSource(),
  _sourceDiff(),
  _sourceDiffSum(),
  _dummyJacob()
{
  addConfigOptionsTo(this);
}

//////////////////////////////////////////////////////////////////////////////

FVMCC_ComputeRhsBlockJacob::~FVMCC_ComputeRhsBlockJacob()
{
}

//////////////////////////////////////////////////////////////////////////////

void FVMCC_ComputeRhsBlockJacob::defineConfigOptions(Config::OptionList& options)
{
}

//////////////////////////////////////////////////////////////////////////////

void FVMCC_ComputeRhsBlockJacob::configure ( Config::ConfigArgs& args )
{
  FVMCC_ComputeRHS::configure(args);
}

//////////////////////////////////////////////////////////////////////////////

void FVMCC_ComputeRhsBlockJacob::setup()
{
  FVMCC_ComputeRHS::setup();
  
  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  
  _acc.reset(new BlockAccumulatorBase(1, 1, nbEqs, CFNULL));
  
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
  
  // numerical jacobian
  _numericalJacob = &getMethodData().getNumericalJacobian();

  // allocate the diagonal block matrix storage
  const CFuint nbStates = socket_states.getDataHandle().size();
  socket_diagMatrix.getDataHandle().resize(nbStates*nbEqs*nbEqs);
  socket_diagMatrix.getDataHandle() = 0.;
}
      
//////////////////////////////////////////////////////////////////////////////

void FVMCC_ComputeRhsBlockJacob::computeJacobianTerm()
{ 
  CFLog(DEBUG_MIN, "FVMCC_ComputeRhsBlockJacob::computeJacobianTerm() START\n"); 
  
  // set the index of the block corresponding to the current
  // state in the jacobian matrix
  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  const CFuint nbEqsSq = nbEqs*nbEqs;
  (!getMethodData().isAxisymmetric()) ? computeNoAxiUpFactors() : computeAxiUpFactors();
  
  DataHandle<CFreal> diagMat = socket_diagMatrix.getDataHandle();
  // first node contribution
  for (CFuint iCell = 0; iCell < 2; ++iCell) {
    State& state = *_currFace->getState(iCell);
    cf_assert(!state.isGhost());
    const CFuint stateID = state.getLocalID();
    _acc->reset(1, 1, nbEqs, &diagMat[stateID*nbEqsSq]);
    
    for (CFuint iVar = 0; iVar < nbEqs; ++iVar) {
      CFLogDebugMax("Perturbing iVar = " << iVar << "\n");
    
      // set the perturbed variable
      getMethodData().setIPerturbVar(iVar);
      
      // perturb the given component of the state vector
      _numericalJacob->perturb(iVar, state[iVar]); 
      
      computeConvDiffFluxes(iVar, iCell);  
      addJacobTerm(0, iCell, iVar, 0, _acc.get());
      
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
  
  // reset to zero the entries in the block accumulator
  _sourceJacobOnCell[LEFT] = _sourceJacobOnCell[RIGHT] = false;
  
  CFLog(DEBUG_MIN, "FVMCC_ComputeRhsBlockJacob::computeJacobianTerm() END\n"); 
}
      
//////////////////////////////////////////////////////////////////////////////

void FVMCC_ComputeRhsBlockJacob::computeBoundaryJacobianTerm()
{
  getMethodData().setIsPerturb(true);
  
  // set the index of the block corresponding to the current
  // state in the jacobian matrix
  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  const CFuint nbEqsSq = nbEqs*nbEqs;
  State& currState = *_currFace->getState(0);
  State& ghostState = *_currFace->getState(1);
  cf_assert(ghostState.isGhost());
  
  if (currState.isParUpdatable()) {
    const bool isAxi = getMethodData().isAxisymmetric();
    _upFactor[LEFT] = (!isAxi) ? getResFactor() : getResFactor()*(_rMid*_invr[0]);
    _upStFactor[LEFT] = (!isAxi) ? -getResFactor() :-getResFactor()*_invr[0];
    
    // copy the original value of the ghost state
    _origState = ghostState;

    const CFuint stateID = currState.getLocalID();
    DataHandle<CFreal> diagMat = socket_diagMatrix.getDataHandle();
    _acc->reset(1, 1, nbEqs, &diagMat[stateID*nbEqsSq]);
    
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
      
      addJacobTerm(0, 0, iVar, 0, _acc.get());
      
      // restore the unperturbed value
      _numericalJacob->restore(currState[iVar]);

      // restore the original ghost state
      ghostState = _origState;
    }

    // compute analytical jacobian for source term 
    if (computeSourceTermJacob(LEFT,_stAnJacobIDs)) {
      addAnalyticSourceTermJacob(LEFT, _acc.get());
    }
    _sourceJacobOnCell[LEFT] = false;
  } 
}

//////////////////////////////////////////////////////////////////////////////

void FVMCC_ComputeRhsBlockJacob::computeConvDiffFluxes(CFuint iVar, CFuint iCell)
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

void FVMCC_ComputeRhsBlockJacob::addJacobTerm(CFuint idx, CFuint idxUpFactor,
					      CFuint iVar, CFuint iCell,
					      BlockAccumulatorBase *const acc)
{ 
  //multiply by the residual factor
  _fluxDiff *= _upFactor[idxUpFactor]; 
  
  acc->addValues(idx, iCell, iVar, &_fluxDiff[0]);
  
  if (computeSourceTermJacob(idx, iCell,_stNumJacobIDs)) {
    _sourceDiffSum = 0.0;
    addSourceTermNumJacob(_currFace->getNeighborGeo(idx), idx);
    _sourceDiffSum *= _upStFactor[idx];
    acc->addValues(iCell, iCell, iVar, &_sourceDiffSum[0]);
  }  
}

//////////////////////////////////////////////////////////////////////////////
      
void FVMCC_ComputeRhsBlockJacob::initializeComputationRHS()
{
  // reset rhs to 0
  socket_rhs.getDataHandle() = 0.0;
  socket_cellFlag.getDataHandle() = false;
  
  // reset to zero all non zero entries in the jacobian
  if (getMethodData().doComputeJacobian()) {
    socket_diagMatrix.getDataHandle() = 0.;
  }
  
  // try to see if this fixes the LS with Roe
  _polyRec->updateWeights();
  
  // gradients and limiters are computed on all the variables at once
  _polyRec->computeGradients();
  
  // extrapolate the solution from cell centers to all vertices
  _nodalExtrapolator->extrapolateInAllNodes();
}
      
//////////////////////////////////////////////////////////////////////////////
 
void FVMCC_ComputeRhsBlockJacob::computeRHSJacobian()
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

void FVMCC_ComputeRhsBlockJacob::finalizeComputationRHS()
{
  // reset the flag for freezing the transport properties
  _diffVar->setFreezeCoeff(false);
}

//////////////////////////////////////////////////////////////////////////////

vector<SafePtr<BaseDataSocketSource> > FVMCC_ComputeRhsBlockJacob::providesSockets()
{  
  CFAUTOTRACE;
  
  vector<SafePtr<BaseDataSocketSource> > result = FVMCC_ComputeRHS::providesSockets();
  result.push_back(&socket_diagMatrix);
  return result;
}

//////////////////////////////////////////////////////////////////////////////

} // namespace FiniteVolume

} // namespace Numerics

} // namespace COOLFluiD
