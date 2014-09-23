#include "Framework/CFL.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/MeshData.hh"

#include "FiniteVolume/FiniteVolume.hh"
#include "FVMCC_ComputeRhsJacobTriGM.hh"
#include "FiniteVolume/FVMCC_BC.hh"

#include "Framework/LSSMatrix.hh"
#include "Framework/BlockAccumulator.hh"

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

MethodCommandProvider<FVMCC_ComputeRhsJacobTriGM,
                      CellCenterFVMData,
                      FiniteVolumeModule>
                      fvmcc_ComputeRhsJacobTriGM("NumJacobTriGM");

//////////////////////////////////////////////////////////////////////////////

FVMCC_ComputeRhsJacobTriGM::FVMCC_ComputeRhsJacobTriGM
(const std::string& name) :
  FVMCC_ComputeRhsJacob(name),
  socket_diagMatrices("diagMatrices"),
  socket_underDiagMatrices("underDiagMatrices"),
  socket_aboveDiagMatrices("aboveDiagMatrices"),
  socket_upLocalIDsAll("upLocalIDsAll"),
  _matIterDiag0(),
  _matIterDiag1(),
  _matIterUnderDiag0(),
  _matIterUnderDiag1(),
  _matIterAboveDiag0(),
  _matIterAboveDiag1()
{
}

//////////////////////////////////////////////////////////////////////////////

FVMCC_ComputeRhsJacobTriGM::~FVMCC_ComputeRhsJacobTriGM()
{
}

//////////////////////////////////////////////////////////////////////////////

void FVMCC_ComputeRhsJacobTriGM::defineConfigOptions(Config::OptionList& options)
{
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSource> >
FVMCC_ComputeRhsJacobTriGM::providesSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSource> > result = FVMCC_ComputeRhsJacob::providesSockets();

  result.push_back(&socket_diagMatrices);
  result.push_back(&socket_underDiagMatrices);
  result.push_back(&socket_aboveDiagMatrices);
  result.push_back(&socket_upLocalIDsAll);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

void FVMCC_ComputeRhsJacobTriGM::setup()
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

  DataHandle<CFreal> underDiagMatrices = socket_underDiagMatrices.getDataHandle();
  underDiagMatrices.resize(countUpdatables*nbEqs*nbEqs);

  DataHandle<CFreal> aboveDiagMatrices = socket_aboveDiagMatrices.getDataHandle();
  aboveDiagMatrices.resize(countUpdatables*nbEqs*nbEqs);

  _matIterDiag0.resize(nbEqs,nbEqs,false);
  _matIterDiag1.resize(nbEqs,nbEqs,false);

  _matIterUnderDiag0.resize(nbEqs, nbEqs, false);
  _matIterUnderDiag1.resize(nbEqs, nbEqs, false);

  _matIterAboveDiag0.resize(nbEqs, nbEqs, false);
  _matIterAboveDiag1.resize(nbEqs, nbEqs, false);

}

//////////////////////////////////////////////////////////////////////////////

void FVMCC_ComputeRhsJacobTriGM::computeBothJacobTerms()
{
  // set on the perturbation flag
  _fluxData->isPerturb = true;

  // set the index of the block corresponding to the current
  // state in the jacobian matrix
  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  const CFuint nbEqs2 = nbEqs*nbEqs;

  State& state0 = *_currFace->getState(0);
  State& state1 = *_currFace->getState(1);

  DataHandle<CFreal> diagMatrices = socket_diagMatrices.getDataHandle();
  DataHandle<CFreal> underDiagMatrices = socket_underDiagMatrices.getDataHandle();
  DataHandle<CFreal> aboveDiagMatrices = socket_aboveDiagMatrices.getDataHandle();

  DataHandle<CFint> upLocalIDsAll = socket_upLocalIDsAll.getDataHandle();

  const CFint stateID0 = upLocalIDsAll[state0.getLocalID()];
  const CFint stateID1 = upLocalIDsAll[state1.getLocalID()];
  cf_assert(stateID0 >= 0);
  cf_assert(stateID1 >= 0);

  // setting the pointer for diagMatrices
  _matIterDiag0.resetPtr(&diagMatrices[stateID0*nbEqs2]);
  _matIterDiag1.resetPtr(&diagMatrices[stateID1*nbEqs2]);

  // setting the pointer for underDiagMatrices
  _matIterUnderDiag0.resetPtr(&underDiagMatrices[stateID0*nbEqs2]);
  _matIterUnderDiag1.resetPtr(&underDiagMatrices[stateID1*nbEqs2]);

  // setting the pointer for aboveDiagMatrices
  _matIterAboveDiag0.resetPtr(&aboveDiagMatrices[stateID0*nbEqs2]);
  _matIterAboveDiag1.resetPtr(&aboveDiagMatrices[stateID1*nbEqs2]);

  _acc->setRowColIndex(0, stateID0);
  _acc->setRowColIndex(1, stateID1);

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
    // multiply for the residual factor
    _fluxDiff *= getResFactor();

    if (!getMethodData().isAxisymmetric()) {
      _acc->addValues(0, 0, iVar, &_fluxDiff[0]);
      _matIterDiag0.addColumn(_fluxDiff, iVar);

      // flux is opposite in sign for the other state
      _fluxDiff *= -1.;
      _acc->addValues(1, 0, iVar, &_fluxDiff[0]);
      // if stateID0 < stateID1 then assign 1 0 contribution to the under diagonal jacobian
      if(stateID0 < stateID1) {
        _matIterUnderDiag1.addColumn(_fluxDiff, iVar);
      }
      // if stateID0 > stateID1 then assign 1 0 contribution to the above diagonal jacobian
      else {
        _matIterAboveDiag1.addColumn(_fluxDiff, iVar);
      }
    }
    else {
      const CFreal invR0 = _rMid/(state0.getCoordinates())[YY];
      const CFreal invR1 = _rMid/(state1.getCoordinates())[YY];

      _axiFlux = _fluxDiff*invR0;
      _acc->addValues(0, 0, iVar, &_axiFlux[0]);
      _matIterDiag0.addColumn(_axiFlux, iVar);

      _axiFlux = _fluxDiff*(-invR1);
      _acc->addValues(1, 0, iVar, &_axiFlux[0]);
      // if stateID0 < stateID1 then assign 1 0 contribution to the under diagonal jacobian
      if(stateID0 < stateID1) {
        // flux is opposite in sign for the other state
        _matIterUnderDiag1.addColumn(_axiFlux, iVar);
      }
      // if stateID0 > stateID1 then assign 1 0 contribution to the above diagonal jacobian
      else {
        // flux is opposite in sign for the other state
        _matIterAboveDiag1.addColumn(_axiFlux, iVar);
      }
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
      _acc->addValues(0, 1, iVar, &_fluxDiff[0]);
      // if stateID0 < stateID1 then assign 0 1 contribution to the under-diagonal jacobian
      if(stateID0 < stateID1) {
        _matIterAboveDiag0.addColumn(_fluxDiff, iVar);
      }
      // if stateID0 > stateID1 then assign 0 1 contribution to the above-diagonal jacobian
      else {
        _matIterUnderDiag0.addColumn(_fluxDiff, iVar);
      }

      // flux is opposite in sign for the other state
      _fluxDiff *= -1.;

      _acc->addValues(1, 1, iVar, &_fluxDiff[0]);
      _matIterDiag1.addColumn(_fluxDiff, iVar);
    }
    else {
      const CFreal invR0 = _rMid/(state0.getCoordinates())[YY];
      const CFreal invR1 = _rMid/(state1.getCoordinates())[YY];

      _axiFlux = _fluxDiff*invR0;
      _acc->addValues(0, 1, iVar, &_axiFlux[0]);
      // if stateID0 < stateID1 then assign 0 1 contribution to the above-diagonal jacobian
      if(stateID0 < stateID1) {
        _matIterAboveDiag0.addColumn(_axiFlux, iVar);
      }
      // if stateID0 > stateID1 then assign 0 1 contribution to the under-diagonal jacobian
      else {
        _matIterUnderDiag0.addColumn(_axiFlux, iVar);
      }

      _axiFlux = _fluxDiff*(-invR1);
      _acc->addValues(1, 1, iVar, &_axiFlux[0]);
      _matIterDiag1.addColumn(_axiFlux, iVar);
    }

    // restore the unperturbed value
    _numericalJacob->restore(state1[iVar]);
    _polyRec->restoreValues(iVar, RIGHT);
  }

  // add the values in the jacobian matrix
  _lss->getMatrix()->addValues(*_acc);

  // reset to zero the entries in the block accumulator
  _acc->reset();
}

//////////////////////////////////////////////////////////////////////////////

void FVMCC_ComputeRhsJacobTriGM::computeJacobTerm(const CFuint idx)
{
  cf_assert(_currFace->getState(idx)->isParUpdatable());
  
  // set on the perturbation flag
  _fluxData->isPerturb = true;
  
  // set the index of the block corresponding to the current
  // state in the jacobian matrix
  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  const CFreal factor = pow(-1.,static_cast<CFreal>(idx))*getResFactor();

  State& state0 = *_currFace->getState(0);
  State& state1 = *_currFace->getState(1);

  _acc->setRowColIndex(0, state0.getLocalID());
  _acc->setRowColIndex(1, state1.getLocalID());

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
      //_nodalExtrapolator->extrapolateInNodes(*_currFace->getNodes());
      _diffusiveFlux->computeFlux(_dFlux);
      _pertFlux -= _dFlux;
    }

    // compute the finite difference derivative of the flux
    _numericalJacob->computeDerivative(getJacobianFlux(),_pertFlux,_fluxDiff);

    //multiply for the residual factor
    _fluxDiff *= factor;

    if (!getMethodData().isAxisymmetric()) {
      _acc->addValues(idx, 0, iVar, &_fluxDiff[0]);
    }
    else {
      const CFreal invR  = _rMid/(_currFace->getState(idx)->getCoordinates())[YY];
      _axiFlux = _fluxDiff*invR;
      _acc->addValues(idx, 0, iVar, &_axiFlux[0]);
    }

    // restore the unperturbed value
    _numericalJacob->restore(state0[iVar]);
    _polyRec->restoreValues(iVar, LEFT);
  }

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
      _acc->addValues(idx, 1, iVar, &_fluxDiff[0]);
    }
    else { 
      const CFreal invR   = _rMid/(_currFace->getState(idx)->getCoordinates())[YY];
      _axiFlux = _fluxDiff*invR;
      _acc->addValues(idx, 1, iVar, &_axiFlux[0]);
    }

    // restore the unperturbed value
    _numericalJacob->restore(state1[iVar]);
    _polyRec->restoreValues(iVar, RIGHT);
  }

  // add the values in the jacobian matrix
  _lss->getMatrix()->addValues(*_acc);

  // reset to zero the entries in the block accumulator
  _acc->reset();
}


//////////////////////////////////////////////////////////////////////////////

void FVMCC_ComputeRhsJacobTriGM::computeBoundaryJacobianTerm()
{
  // set on the perturbation flag
  _fluxData->isPerturb = true;
  
  // set the index of the block corresponding to the current
  // state in the jacobian matrix
  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  State& currState = *_currFace->getState(0);
  State& ghostState = *_currFace->getState(1);
  cf_assert(ghostState.isGhost());
  
  if (currState.isParUpdatable()) {
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
        _bAcc->addValues(0, 0, iVar, &_fluxDiff[0]);
      }
      else {
	const CFreal invR0 = _rMid/(currState.getCoordinates())[YY];
	_axiFlux = _fluxDiff*invR0;
        _bAcc->addValues(0, 0, iVar, &_axiFlux[0]);
      }

      // restore the unperturbed value
      _numericalJacob->restore(currState[iVar]);

      // restore the original ghost state
      ghostState = _origState;
    }

    // add the values in the jacobian matrix
    _lss->getMatrix()->addValues(*_bAcc);

    // reset to zero the entries in the block accumulator
    _bAcc->reset();
  }
}

//////////////////////////////////////////////////////////////////////////////

void FVMCC_ComputeRhsJacobTriGM::computeSourceTermAnalytJacob(CFuint ist)
{
  SafePtr<LSSMatrix> jacobMatrix = _lss->getMatrix();

  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  const bool isAxisymm = getMethodData().isAxisymmetric();
  const CFuint lastST = _stComputers->size() - 1;

  DataHandle<CFreal> rhs = socket_rhs.getDataHandle();

  SafePtr<TopologicalRegionSet> cells = MeshDataStack::getActive()->
    getTrs("InnerCells");
  const CFuint nbCells = cells->getLocalNbGeoEnts();

  CellTrsGeoBuilder::GeoData& geoData = _cellBuilder->getDataGE();

  // set on the analytical jacobian of the source term
  (*_stComputers)[ist]->setAnalyticalJacob(true);

  ///  @todo this way of computing the source term fails if
  ///        the source term depends on more than one state
  _flux = 0.0;
  
  //
  // static CFuint count = 0;
  // static std::vector<RealMatrix> st(socket_volumes.getDataHandle().size());
  //
    
  for (CFuint iCell = 0; iCell < nbCells; ++iCell) {
    // build the GeometricEntity
    geoData.idx = iCell;
    GeometricEntity *const currCell = _cellBuilder->buildGE();
    State *const currState = currCell->getState(0);
    
    if (currState->isParUpdatable()) {
      _bAcc->setRowColIndex(0, currState->getLocalID());
      
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
      
      // get the jacobian matrix for the source term
      RealMatrix& sourceJacob = (*_stComputers)[ist]->getJacobMatrix();
      
      if (isAxisymm) {
	const RealVector& stateCoord = currState->getCoordinates();
	cf_assert(stateCoord[YY] > 0.0);
	sourceJacob *= -getResFactor()/stateCoord[YY];
      }
      else {
	sourceJacob *= -getResFactor();
      }
            
      _bAcc->addValues(0, 0, sourceJacob);
      
      // print the block !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      // _bAcc->print();

      // add the values in the jacobian matrix
      jacobMatrix->addValues(*_bAcc);

      // reset to zero the entries in the block accumulator
      _bAcc->reset();
    }

    // release the cell
    _cellBuilder->releaseGE();
  }

  // set off the analytical jacobian of the source term
  (*_stComputers)[ist]->setAnalyticalJacob(false);
}

//////////////////////////////////////////////////////////////////////////////

void FVMCC_ComputeRhsJacobTriGM::computeSourceTermNumJacob(CFuint ist)
{
  SafePtr<LSSMatrix> jacobMatrix = _lss->getMatrix();

  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  const bool isAxisymm = getMethodData().isAxisymmetric();
  const CFuint lastST = _stComputers->size() - 1;

  DataHandle<CFreal> rhs = socket_rhs.getDataHandle();

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
      _bAcc->setRowColIndex(0, currState->getLocalID());
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

	_bAcc->addValues(0, 0, iVar, &_fluxDiff[0]);

	// restore the unperturbed value
	_numericalJacob->restore((*currState)[iVar]);
      }

      // print the block !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      // _bAcc->print();

      // add the values in the jacobian matrix
      jacobMatrix->addValues(*_bAcc);
      
      // reset to zero the entries in the block accumulator
      _bAcc->reset();
      
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
