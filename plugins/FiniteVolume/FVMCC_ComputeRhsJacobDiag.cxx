#include "FiniteVolume/FiniteVolume.hh"
#include "FVMCC_ComputeRhsJacobDiag.hh"
#include "FiniteVolume/FVMCC_BC.hh"
#include "Framework/CFL.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/MeshData.hh"

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

MethodCommandProvider<FVMCC_ComputeRhsJacobDiag,
		      CellCenterFVMData,
		      FiniteVolumeModule>
fvmcc_computeRhsJacobDiag("NumJacobDiag");

//////////////////////////////////////////////////////////////////////////////

FVMCC_ComputeRhsJacobDiag::FVMCC_ComputeRhsJacobDiag
(const std::string& name) :
  FVMCC_ComputeRhsJacob(name),
  socket_diagMatrices("diagMatrices"),
  socket_upLocalIDsAll("upLocalIDsAll"),
  _matIter0(),
  _matIter1()
{
}
    
//////////////////////////////////////////////////////////////////////////////

FVMCC_ComputeRhsJacobDiag::~FVMCC_ComputeRhsJacobDiag()
{
}

//////////////////////////////////////////////////////////////////////////////

void FVMCC_ComputeRhsJacobDiag::defineConfigOptions(Config::OptionList& options)
{
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSource> >
FVMCC_ComputeRhsJacobDiag::providesSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSource> > result =  
    FVMCC_ComputeRhsJacob::providesSockets();
  
  result.push_back(&socket_diagMatrices);
  result.push_back(&socket_upLocalIDsAll);
  
  return result;
}

//////////////////////////////////////////////////////////////////////////////

void FVMCC_ComputeRhsJacobDiag::setup()
{
  FVMCC_ComputeRhsJacob::setup();
  
  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
 
  DataHandle<State*, GLOBAL> states = socket_states.getDataHandle();
  DataHandle<CFint> upLocalIDsAll = socket_upLocalIDsAll.getDataHandle();
  upLocalIDsAll.resize(states.size());
  
  CFuint countUpdatables = 0;
  for (CFuint i = 0; i < states.size(); ++i) {
    if (states[i]->isParUpdatable()) {
      upLocalIDsAll[i] = countUpdatables;
      countUpdatables++;
    }
    else {
      upLocalIDsAll[i] = -1;
    }
  }  
  
  // only block matrices corresponding to updatable states are stored
  // in order to save memory in parallel
  DataHandle<CFreal> diagMatrices = socket_diagMatrices.getDataHandle();
  diagMatrices.resize(countUpdatables*nbEqs*nbEqs);
  
  _matIter0.resize(nbEqs,nbEqs,false);
  _matIter1.resize(nbEqs,nbEqs,false);
}
      
//////////////////////////////////////////////////////////////////////////////

void FVMCC_ComputeRhsJacobDiag::computeBothJacobTerms()
{
  // set the index of the block corresponding to the current
  // state in the jacobian matrix
  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq(); 
  const CFuint nbEqs2 = nbEqs*nbEqs;
  
  State& state0 = *_currFace->getState(0);
  State& state1 = *_currFace->getState(1);
  
  DataHandle<CFreal> diagMatrices = socket_diagMatrices.getDataHandle();
  DataHandle<CFint> upLocalIDsAll = socket_upLocalIDsAll.getDataHandle();
  
  const CFint stateID0 = upLocalIDsAll[state0.getLocalID()];
  const CFint stateID1 = upLocalIDsAll[state1.getLocalID()];
  cf_assert(stateID0 >= 0);
  cf_assert(stateID1 >= 0);
  
  _matIter0.resetPtr(&diagMatrices[stateID0*nbEqs2]);
  _matIter1.resetPtr(&diagMatrices[stateID1*nbEqs2]);
  
  // first node contribution
  for (CFuint iVar = 0; iVar < nbEqs; ++iVar) {
    CFLogDebugMax( "Perturbing iVar = " << iVar << "\n");

    // perturb the given component of the state vector
    _numericalJacob->perturb(iVar, state0[iVar]);

    // extrapolate (and LIMIT, if the reconstruction is linear or more)
    // the solution in the quadrature points
    _polyRec->extrapolate(_currFace, iVar, LEFT);

    // compute the physical data for each left and right reconstructed
    // state and in the left and right cell centers
    computeStatesData();

    // set the perturbed variable
    _fluxData->iPerturbVar = iVar;

    // linearization will be done in the flux splitter if needed
    _fluxSplitter->computeFlux(_pertFlux);

    if (_fluxData->hasDiffusiveTerm) {
      // _nodalExtrapolator->extrapolateInNodes(*_currFace->getNodes());
      _diffusiveFlux->computeFlux(_dFlux);
      _pertFlux -= _dFlux;
    }

    // compute the finite difference derivative of the flux
    _numericalJacob->computeDerivative(getJacobianFlux(),_pertFlux,_fluxDiff);

    _fluxDiff *= getResFactor();
    
    if (!getMethodData().isAxisymmetric()) {
      _matIter0.addColumn(_fluxDiff, iVar);
    }
    else {
      const CFreal invR0 = _rMid/(state0.getCoordinates())[YY];
      _axiFlux = _fluxDiff*invR0;
      _matIter0.addColumn(_axiFlux, iVar);
    }

    // restore the unperturbed value
    _numericalJacob->restore(state0[iVar]);
    _polyRec->restoreValues(iVar, LEFT);

    // perturb the given component of the second state vector
    _numericalJacob->perturb(iVar, state1[iVar]);

    // extrapolate (and LIMIT, if the reconstruction is linear or more)
    // the solution in the quadrature points
    _polyRec->extrapolate(_currFace, iVar, RIGHT);

    // compute the physical data for each left and right reconstructed
    // state and in the left and right cell centers
    computeStatesData();

    // set the perturbed variable
    _fluxData->iPerturbVar = iVar;

    // linearization will be done in the flux splitter if needed
    _fluxSplitter->computeFlux(_pertFlux);

    if (_fluxData->hasDiffusiveTerm) {
      // _nodalExtrapolator->extrapolateInNodes(*_currFace->getNodes());
      _diffusiveFlux->computeFlux(_dFlux);
      _pertFlux -= _dFlux;
    }

    // compute the finite difference derivative of the flux
    _numericalJacob->computeDerivative(getJacobianFlux(),
                                       _pertFlux,
                                       _fluxDiff);

    // multiply for the residual factor
    _fluxDiff *= getResFactor();

    if (!getMethodData().isAxisymmetric()) {
      // flux is opposite in sign for the other state
      _fluxDiff *= -1.;
      _matIter1.addColumn(_fluxDiff, iVar);
    }
    else {
      const CFreal invR1 = _rMid/(state1.getCoordinates())[YY];
      _axiFlux = _fluxDiff*(-invR1);
      _matIter1.addColumn(_axiFlux, iVar);
    }

          
    // restore the unperturbed value
    _numericalJacob->restore(state1[iVar]);
    _polyRec->restoreValues(iVar, RIGHT);
  }
}
      
//////////////////////////////////////////////////////////////////////////////

void FVMCC_ComputeRhsJacobDiag::computeJacobTerm(const CFuint idx)
{
  cf_assert(_currFace->getState(idx)->isParUpdatable());
  DataHandle<CFreal> diagMatrices = socket_diagMatrices.getDataHandle();
  DataHandle<CFint> upLocalIDsAll = socket_upLocalIDsAll.getDataHandle();
  
  // set the index of the block corresponding to the current
  // state in the jacobian matrix
  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  const CFuint nbEqs2 = nbEqs*nbEqs;
  const CFreal factor = pow(-1.,static_cast<CFreal>(idx))*getResFactor();
  
  if (idx == 0) {
    State& state0 = *_currFace->getState(0);
    const CFint stateID0 = upLocalIDsAll[state0.getLocalID()];
    cf_assert(stateID0 >= 0);
    
    _matIter0.resetPtr(&diagMatrices[stateID0*nbEqs2]);
    
    // first node contribution
    for (CFuint iVar = 0; iVar < nbEqs; ++iVar) {
      CFLogDebugMax( "Perturbing iVar = " << iVar << "\n");
      
      // perturb the given component of the state vector
      _numericalJacob->perturb(iVar, state0[iVar]);
      
      // extrapolate (and LIMIT, if the reconstruction is linear or more)
      // the solution in the quadrature points
      _polyRec->extrapolate(_currFace, iVar, LEFT);
      
      // compute the physical data for each left and right reconstructed
      // state and in the left and right cell centers
      computeStatesData();
      
      // set the perturbed variable
      _fluxData->iPerturbVar = iVar;
      
      // linearization will be done in the flux splitter if needed
      _fluxSplitter->computeFlux(_pertFlux);
      
      if (_fluxData->hasDiffusiveTerm) {
	// _nodalExtrapolator->extrapolateInNodes(*_currFace->getNodes());
	_diffusiveFlux->computeFlux(_dFlux);
	_pertFlux -= _dFlux;
      }
      
      // compute the finite difference derivative of the flux
      _numericalJacob->computeDerivative(getJacobianFlux(),_pertFlux,_fluxDiff);
      
      //multiply for the residual factor
      _fluxDiff *= factor;
      
      if (!getMethodData().isAxisymmetric()) {
	_matIter0.addColumn(_fluxDiff, iVar);
      }
      else {
	const CFreal invR  = _rMid/(_currFace->getState(idx)->getCoordinates())[YY];
	_axiFlux = _fluxDiff*invR;
	_matIter0.addColumn(_axiFlux, iVar);
      }
      
      // restore the unperturbed value
      _numericalJacob->restore(state0[iVar]);
      _polyRec->restoreValues(iVar, LEFT);
    }
  }
  
  if (idx == 1) {
    State& state1 = *_currFace->getState(1);
    const CFint stateID1 = upLocalIDsAll[state1.getLocalID()];
    cf_assert(stateID1 >= 0);
    
    _matIter1.resetPtr(&diagMatrices[stateID1*nbEqs2]);
    
    // second node contribution
    for (CFuint iVar = 0; iVar < nbEqs; ++iVar) {
      CFLogDebugMax( "Perturbing iVar = " << iVar << "\n");
      
      // perturb the given component of the state vector
      _numericalJacob->perturb(iVar, state1[iVar]);
      // extrapolate (and LIMIT, if the reconstruction is linear or more)
      // the solution in the quadrature points
      _polyRec->extrapolate(_currFace, iVar, RIGHT);
      
      // compute the physical data for each left and right reconstructed
      // state and in the left and right cell centers
      computeStatesData();
      
      // set the perturbed variable
      _fluxData->iPerturbVar = iVar;
      
      // linearization will be done in the flux splitter if needed
      _fluxSplitter->computeFlux(_pertFlux);
      
      if (_fluxData->hasDiffusiveTerm) {
	// _nodalExtrapolator->extrapolateInNodes(*_currFace->getNodes());
	_diffusiveFlux->computeFlux(_dFlux);
	_pertFlux -= _dFlux;
      }
      
      // compute the finite difference derivative of the flux
      _numericalJacob->computeDerivative(getJacobianFlux(),
					 _pertFlux,
					 _fluxDiff);
      
      // multiply for the residual factor
      _fluxDiff *= factor;
      
      if (!getMethodData().isAxisymmetric()) {
	_matIter1.addColumn(_fluxDiff, iVar);
      }
      else {  
	const CFreal invR   = _rMid/(_currFace->getState(idx)->getCoordinates())[YY];
	_axiFlux = _fluxDiff*invR;
	_matIter1.addColumn(_axiFlux, iVar);
      }
      
      // restore the unperturbed value
      _numericalJacob->restore(state1[iVar]);
      _polyRec->restoreValues(iVar, RIGHT);
    }
  }
}


//////////////////////////////////////////////////////////////////////////////

void FVMCC_ComputeRhsJacobDiag::computeBoundaryJacobianTerm()
{
  // set the index of the block corresponding to the current
  // state in the jacobian matrix
  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  const CFuint nbEqs2 = nbEqs*nbEqs;
  
  State& currState = *_currFace->getState(0);
  State& ghostState = *_currFace->getState(1);
  cf_assert(ghostState.isGhost());
  
  DataHandle<CFreal> diagMatrices = socket_diagMatrices.getDataHandle();
  DataHandle<CFint> upLocalIDsAll = socket_upLocalIDsAll.getDataHandle();
  
  if (currState.isParUpdatable()) {
    const CFint stateID0 = upLocalIDsAll[currState.getLocalID()];
    cf_assert(stateID0 >= 0);
    
    _matIter0.resetPtr(&diagMatrices[stateID0*nbEqs2]); 
   
    // copy the original value of the ghost state
    _origState = ghostState;
    
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
      _fluxData->iPerturbVar = iVar;
      
      _currBC->computeFlux(_pertFlux);

      if (_fluxData->hasDiffusiveTerm) {
	// _nodalExtrapolator->extrapolateInNodes(*_currFace->getNodes());
        _diffusiveFlux->computeFlux(_dFlux);
        _pertFlux -= _dFlux;
      }
      
      // compute the finite difference derivative of the flux
      _numericalJacob->computeDerivative(getJacobianFlux(),
                                         _pertFlux,
                                         _fluxDiff);

      // multiply for the residual factor
      _fluxDiff *= getResFactor();

      if (!getMethodData().isAxisymmetric()) {
        _matIter0.addColumn(_fluxDiff, iVar);
      }
      else {
	const CFreal invR0 = _rMid/(currState.getCoordinates())[YY];
	_axiFlux = _fluxDiff*invR0;
	_matIter0.addColumn(_axiFlux, iVar);
      }

      // restore the unperturbed value
      _numericalJacob->restore(currState[iVar]);
      
      // restore the original ghost state
      ghostState = _origState;
    }
  }
}
      
//////////////////////////////////////////////////////////////////////////////

void FVMCC_ComputeRhsJacobDiag::computeSourceTermAnalytJacob(CFuint ist)
{
  SafePtr<LSSMatrix> jacobMatrix = _lss->getMatrix();

  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  const CFuint nbEqs2 = nbEqs*nbEqs;
  
  const bool isAxisymm = getMethodData().isAxisymmetric();
  const CFuint lastST = _stComputers->size() - 1;

  DataHandle<CFreal> rhs = socket_rhs.getDataHandle();
  DataHandle<CFreal> diagMatrices = socket_diagMatrices.getDataHandle();
  DataHandle<CFint> upLocalIDsAll = socket_upLocalIDsAll.getDataHandle();
  
  SafePtr<TopologicalRegionSet> cells = MeshDataStack::getActive()->
    getTrs("InnerCells");
  const CFuint nbCells = cells->getLocalNbGeoEnts();
  
  CellTrsGeoBuilder::GeoData& geoData = _cellBuilder->getDataGE();
  
  // set on the analytical jacobian of the source term
  (*_stComputers)[ist]->setAnalyticalJacob(true);
  
  ///  @todo this way of computing the source term fails if
  ///        the source term depends on more than one state
  _flux = 0.0;
  
  for (CFuint iCell = 0; iCell < nbCells; ++iCell) {
    // build the GeometricEntity
    geoData.idx = iCell;
    GeometricEntity *const currCell = _cellBuilder->buildGE();
    State *const currState = currCell->getState(0);
    
    if (currState->isParUpdatable()) {
      const CFint stateID0 = upLocalIDsAll[currState->getLocalID()];
      cf_assert(stateID0 >= 0);
      
      _matIter0.resetPtr(&diagMatrices[stateID0*nbEqs2]);
      
      // compute the unperturbed source term and its analytical jacobian
      (*_stComputers)[ist]->computeSource(currCell, _flux);
      
      // update the RHS
      if (isAxisymm) {
	const RealVector& stateCoord = currState->getCoordinates();
	const CFreal ovR = 1./stateCoord[YY]; // 1/r
	cf_assert(stateCoord[YY] > 0.0);
	
	for (CFuint iEq = 0; iEq < nbEqs; ++iEq) {
	  rhs(iCell, iEq, nbEqs) += getResFactor()*_flux[iEq];
	  // all the rhs has to be divided by r
	  // if there are more than 1 source term, this does
	  // this operation has to be done only ONCE after the
	  // last source term contribution has been added,
	  if (ist == lastST) {
	    rhs(iCell, iEq, nbEqs) *= ovR;
	  }
	}
      }
      else {
	for (CFuint iEq = 0; iEq < nbEqs; ++iEq) {
	  rhs(iCell, iEq, nbEqs) += getResFactor()*_flux[iEq];
	}
      }
      
      // get the jacobian ,atrix for the source term
      RealMatrix& sourceJacob = (*_stComputers)[ist]->getJacobMatrix();
      
      if (isAxisymm) {
	const RealVector& stateCoord = currState->getCoordinates();
	cf_assert(stateCoord[YY] > 0.0);
	sourceJacob *= -getResFactor()/stateCoord[YY];
      }
      else {
	sourceJacob *= -getResFactor();
      }
      
      for (CFuint m = 0; m < nbEqs2; ++m) {
	_matIter0[m] = sourceJacob[m];
      }
    }
    
    // release the cell
    _cellBuilder->releaseGE();
  }
  
  // set off the analytical jacobian of the source term
  (*_stComputers)[ist]->setAnalyticalJacob(false);
}

//////////////////////////////////////////////////////////////////////////////

void FVMCC_ComputeRhsJacobDiag::computeSourceTermNumJacob(CFuint ist)
{
  SafePtr<LSSMatrix> jacobMatrix = _lss->getMatrix();

  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  const CFuint nbEqs2 = nbEqs*nbEqs;
  const bool isAxisymm = getMethodData().isAxisymmetric();
  const CFuint lastST = _stComputers->size() - 1;
  
  DataHandle<CFreal> rhs = socket_rhs.getDataHandle();
  DataHandle<CFreal> diagMatrices = socket_diagMatrices.getDataHandle();
  DataHandle<CFint> upLocalIDsAll = socket_upLocalIDsAll.getDataHandle();
  
  SafePtr<TopologicalRegionSet> cells = MeshDataStack::getActive()->
    getTrs("InnerCells");
  const CFuint nbCells = cells->getLocalNbGeoEnts();

  CellTrsGeoBuilder::GeoData& geoData = _cellBuilder->getDataGE();

  ///  @todo this way of computing the source term fails if
  ///        the source term depends on more than one state
  // add a diagonal block for the source term
    
  _flux = 0.;
  _pertFlux = 0.0;

  for (CFuint iCell = 0; iCell < nbCells; ++iCell) {
    // build the GeometricEntity
    geoData.idx = iCell;
    GeometricEntity *const currCell = _cellBuilder->buildGE();
    State *const currState = currCell->getState(0);
    
    if (currState->isParUpdatable()) {
      const CFint stateID0 = upLocalIDsAll[currState->getLocalID()];
      cf_assert(stateID0 >= 0);
      
      _matIter0.resetPtr(&diagMatrices[stateID0*nbEqs2]);
      
      // compute the unperturbed source term
      (*_stComputers)[ist]->computeSource(currCell, _flux);
      
      const vector<Node*>* const cellNodes = currCell->getNodes();
      const CFuint nbNodesInCell = cellNodes->size();
      
      // update the RHS
      if (isAxisymm) {
	const RealVector& stateCoord = currState->getCoordinates();
	const CFreal ovR = 1./stateCoord[YY]; // 1/r
	cf_assert(stateCoord[YY] > 0.0);
	
	for (CFuint iEq = 0; iEq < nbEqs; ++iEq) {
	  rhs(iCell, iEq, nbEqs) += getResFactor()*_flux[iEq];
	  // all the rhs has to be divided by r
	  // if there are more than 1 source term, this does
	  // this operation has to be done only ONCE after the
	  // last source term contribution has been added,
	  if (ist == lastST) {
	    rhs(iCell, iEq, nbEqs) *= ovR;
	  }
	}
      }
      else {
	for (CFuint iEq = 0; iEq < nbEqs; ++iEq) {
	  rhs(iCell, iEq, nbEqs) += getResFactor()*_flux[iEq];
	}
      }
      
      _fluxData->isPerturb = true;
      
      // compute the source term corresponding to the perturbed state
      for (CFuint iVar = 0; iVar < nbEqs; ++iVar) {
	CFLogDebugMax( "Perturbing iVar = " << iVar << "\n");
	// perturb the given component of the state vector
	_numericalJacob->perturb(iVar, (*currState)[iVar]);
	
	// if it has diffusive flux, nodal reconstruction is done
	// starting from the perturbed state
	if (_fluxData->hasDiffusiveTerm) {
	  _cellNodes.clear();
	  // store the cells
	  for (CFuint in = 0; in < nbNodesInCell; ++in) {
	    _cellNodes.push_back((*cellNodes)[in]);
	  }
	  //_nodalExtrapolator->extrapolateInNodes(_cellNodes);
	}
	
	(*_stComputers)[ist]->computeSource(currCell, _pertFlux);
	// compute the finite difference derivative of the source term
	_numericalJacob->computeDerivative(_flux, _pertFlux, _fluxDiff);
	
	if (isAxisymm) {
	  const RealVector& stateCoord = currState->getCoordinates();
	  cf_assert(stateCoord[YY] > 0.0);
	  _fluxDiff *= -getResFactor()/stateCoord[YY];
	}
	else {
	  _fluxDiff *= -getResFactor();
	}
	
	_matIter0.addColumn(_fluxDiff, iVar);
	
	// restore the unperturbed value
	_numericalJacob->restore((*currState)[iVar]);
      }
      
      _fluxData->isPerturb = false;
    }
    
    // release the cell
    _cellBuilder->releaseGE();
  }
}

//////////////////////////////////////////////////////////////////////////////

} // namespace FiniteVolume

} // namespace Numerics

} // namespace COOLFluiD
