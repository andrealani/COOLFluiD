#include "Framework/CFL.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/MeshData.hh"

#include "FiniteVolume/FiniteVolume.hh"
#include "FVMCC_ComputeRhsJacobTridiag.hh"
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

MethodCommandProvider<FVMCC_ComputeRhsJacobTridiag,
                      CellCenterFVMData,
                      FiniteVolumeModule>
                      fvmcc_computeRhsJacobTridiag("NumJacobTridiag");

//////////////////////////////////////////////////////////////////////////////

FVMCC_ComputeRhsJacobTridiag::FVMCC_ComputeRhsJacobTridiag
(const std::string& name) :
  FVMCC_ComputeRhsJacob(name),
  socket_diagMatrices("diagMatrices"),
  socket_underDiagMatrices("underDiagMatrices"),
  socket_aboveDiagMatrices("aboveDiagMatrices"),
  socket_upLocalIDsAll("upLocalIDsAll"),
  socket_dplrIsFirstInLine("dplrIsFirstInLine"),
  socket_dplrToLocalIDs("dplrToLocalIDs"),
  socket_localToDplrIDs("localToDplrIDs"),
  socket_dplrCellInLine("dplrCellInLine"),
  _matIterDiag0(),
  _matIterDiag1(),
  _matIterUnderDiag0(),
  _matIterUnderDiag1(),
  _matIterAboveDiag0(),
  _matIterAboveDiag1()
{
}

//////////////////////////////////////////////////////////////////////////////

FVMCC_ComputeRhsJacobTridiag::~FVMCC_ComputeRhsJacobTridiag()
{
}

//////////////////////////////////////////////////////////////////////////////

void FVMCC_ComputeRhsJacobTridiag::defineConfigOptions(Config::OptionList& options)
{
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSource> >
FVMCC_ComputeRhsJacobTridiag::providesSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSource> > result = FVMCC_ComputeRhsJacob::providesSockets();

  result.push_back(&socket_diagMatrices);
  result.push_back(&socket_underDiagMatrices);
  result.push_back(&socket_aboveDiagMatrices);
  result.push_back(&socket_upLocalIDsAll);
  result.push_back(&socket_dplrIsFirstInLine);
  result.push_back(&socket_dplrToLocalIDs);
  result.push_back(&socket_localToDplrIDs);
  result.push_back(&socket_dplrCellInLine);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

void FVMCC_ComputeRhsJacobTridiag::setup()
{
  FVMCC_ComputeRhsJacob::setup();

  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();

  DataHandle<State*, GLOBAL> states = socket_states.getDataHandle();
  DataHandle<CFint> upLocalIDsAll = socket_upLocalIDsAll.getDataHandle();

  const CFuint nbStates = states.size();
  upLocalIDsAll.resize(nbStates);

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

  DataHandle<bool> dplrIsFirstInLine = socket_dplrIsFirstInLine.getDataHandle();
  DataHandle<CFint> dplrToLocalIDs = socket_dplrToLocalIDs.getDataHandle();
  DataHandle<CFint> localToDplrIDs = socket_localToDplrIDs.getDataHandle();
  DataHandle<CFint> dplrCellInLine = socket_dplrCellInLine.getDataHandle();

  dplrIsFirstInLine.resize(nbStates);
  dplrToLocalIDs.resize(nbStates);
  localToDplrIDs.resize(nbStates);
  dplrCellInLine.resize(nbStates);

}

//////////////////////////////////////////////////////////////////////////////

void FVMCC_ComputeRhsJacobTridiag::computeBothJacobTerms()
{
  // set the index of the block corresponding to the current
  // state in the jacobian matrix
  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  const CFuint nbEqs2 = nbEqs*nbEqs;

  State& state0 = *_currFace->getState(0);
  State& state1 = *_currFace->getState(1);

  DataHandle<CFreal> rhs = socket_rhs.getDataHandle();

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
      _matIterDiag0.addColumn(_fluxDiff, iVar);

      // flux is opposite in sign for the other state
      _fluxDiff *= -1.;
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
      _matIterDiag0.addColumn(_axiFlux, iVar);

      _axiFlux = _fluxDiff*(-invR1);
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
      _matIterDiag1.addColumn(_fluxDiff, iVar);
    }
    else {
      const CFreal invR0 = _rMid/(state0.getCoordinates())[YY];
      const CFreal invR1 = _rMid/(state1.getCoordinates())[YY];

      _axiFlux = _fluxDiff*invR0;
      // if stateID0 < stateID1 then assign 0 1 contribution to the above-diagonal jacobian
      if(stateID0 < stateID1) {
        _matIterAboveDiag0.addColumn(_axiFlux, iVar);
      }
      // if stateID0 > stateID1 then assign 0 1 contribution to the under-diagonal jacobian
      else {
        _matIterUnderDiag0.addColumn(_axiFlux, iVar);
      }

      _axiFlux = _fluxDiff*(-invR1);
      _matIterDiag1.addColumn(_axiFlux, iVar);
    }

    // restore the unperturbed value
    _numericalJacob->restore(state1[iVar]);
    _polyRec->restoreValues(iVar, RIGHT);
  }
}

//////////////////////////////////////////////////////////////////////////////

void FVMCC_ComputeRhsJacobTridiag::computeJacobTerm(const CFuint idx)
{
  cf_assert(_currFace->getState(idx)->isParUpdatable());

  DataHandle<CFreal> diagMatrices = socket_diagMatrices.getDataHandle();
  DataHandle<CFreal> underDiagMatrices = socket_underDiagMatrices.getDataHandle();
  DataHandle<CFreal> aboveDiagMatrices = socket_aboveDiagMatrices.getDataHandle();

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

    _matIterDiag0.resetPtr(&diagMatrices[stateID0*nbEqs2]);

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
        _matIterDiag0.addColumn(_fluxDiff, iVar);
      }
      else {
        const CFreal invR  = _rMid/(_currFace->getState(idx)->getCoordinates())[YY];
        _axiFlux = _fluxDiff*invR;
        _matIterDiag0.addColumn(_axiFlux, iVar);
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

    _matIterDiag1.resetPtr(&diagMatrices[stateID1*nbEqs2]);

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
      _numericalJacob->computeDerivative(getJacobianFlux(), _pertFlux, _fluxDiff);

      // multiply for the residual factor
      _fluxDiff *= factor;

      if (!getMethodData().isAxisymmetric()) {
        _matIterDiag1.addColumn(_fluxDiff, iVar);
      }
      else {
        const CFreal invR   = _rMid/(_currFace->getState(idx)->getCoordinates())[YY];
        _axiFlux = _fluxDiff*invR;
        _matIterDiag1.addColumn(_axiFlux, iVar);
      }

      // restore the unperturbed value
      _numericalJacob->restore(state1[iVar]);
      _polyRec->restoreValues(iVar, RIGHT);
    }
  }
}


//////////////////////////////////////////////////////////////////////////////

void FVMCC_ComputeRhsJacobTridiag::computeBoundaryJacobianTerm()
{
  // set the index of the block corresponding to the current
  // state in the jacobian matrix
  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  const CFuint nbEqs2 = nbEqs*nbEqs;

  State& currState = *_currFace->getState(0);
  State& ghostState = *_currFace->getState(1);
  cf_assert(ghostState.isGhost());

  DataHandle<CFreal> diagMatrices = socket_diagMatrices.getDataHandle();
  //DataHandle<CFreal> underDiagMatrices = socket_underDiagMatrices.getDataHandle();
  //DataHandle<CFreal> aboveDiagMatrices = socket_aboveDiagMatrices.getDataHandle();

  DataHandle<CFint> upLocalIDsAll = socket_upLocalIDsAll.getDataHandle();

  if (currState.isParUpdatable()) {
    const CFint stateID0 = upLocalIDsAll[currState.getLocalID()];
    cf_assert(stateID0 >= 0);

    _matIterDiag0.resetPtr(&diagMatrices[stateID0*nbEqs2]); 

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
        _matIterDiag0.addColumn(_fluxDiff, iVar);
      }
      else {
        const CFreal invR0 = _rMid/(currState.getCoordinates())[YY];
        _axiFlux = _fluxDiff*invR0;
        _matIterDiag0.addColumn(_axiFlux, iVar);
      }

      // restore the unperturbed value
      _numericalJacob->restore(currState[iVar]);

      // restore the original ghost state
      ghostState = _origState;
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void FVMCC_ComputeRhsJacobTridiag::computeSourceTermAnalytJacob(CFuint ist)
{
  SafePtr<LSSMatrix> jacobMatrix = _lss->getMatrix();

  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  const CFuint nbEqs2 = nbEqs*nbEqs;

  const bool isAxisymm = getMethodData().isAxisymmetric();
  const CFuint lastST = _stComputers->size() - 1;

  DataHandle<CFreal> rhs = socket_rhs.getDataHandle();

  DataHandle<CFreal> diagMatrices = socket_diagMatrices.getDataHandle();
  //DataHandle<CFreal> underDiagMatrices = socket_underDiagMatrices.getDataHandle();
  //DataHandle<CFreal> aboveDiagMatrices = socket_aboveDiagMatrices.getDataHandle();

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

      _matIterDiag0.resetPtr(&diagMatrices[stateID0*nbEqs2]);

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
        _matIterDiag0[m] = sourceJacob[m];
      }
    }

    // release the cell
    _cellBuilder->releaseGE();
  }

  // set off the analytical jacobian of the source term
  (*_stComputers)[ist]->setAnalyticalJacob(false);
}

//////////////////////////////////////////////////////////////////////////////

void FVMCC_ComputeRhsJacobTridiag::computeSourceTermNumJacob(CFuint ist)
{
  SafePtr<LSSMatrix> jacobMatrix = _lss->getMatrix();

  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  const CFuint nbEqs2 = nbEqs*nbEqs;
  const bool isAxisymm = getMethodData().isAxisymmetric();
  const CFuint lastST = _stComputers->size() - 1;

  DataHandle<CFreal> rhs = socket_rhs.getDataHandle();

  DataHandle<CFreal> diagMatrices = socket_diagMatrices.getDataHandle();
  DataHandle<CFreal> underDiagMatrices = socket_underDiagMatrices.getDataHandle();
  DataHandle<CFreal> aboveDiagMatrices = socket_aboveDiagMatrices.getDataHandle();

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

      _matIterDiag0.resetPtr(&diagMatrices[stateID0*nbEqs2]);

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
        _matIterDiag0.addColumn(_fluxDiff, iVar);

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

void FVMCC_ComputeRhsJacobTridiag::initializeComputationRHS()
{
  // here build lines
  SafePtr<TopologicalRegionSet> cells = MeshDataStack::getActive()->getTrs("InnerCells");
  const CFuint nbCells = cells->getLocalNbGeoEnts();

  cout << "nbCells = " << nbCells << endl;

  CellTrsGeoBuilder::GeoData& geoData = _cellBuilder->getDataGE();

	/*
	vector<CFint> lines(nbCells, -1);     // IDs of line to which the cell belongs
	vector<CFint> neighbor1(nbCells, -1); // neighbor with highest gradient
	vector<CFint> neighbor2(nbCells, -1); // neigbor with second highest gradient
	CFint firstCellInLine = -1;
	
	// loop over all cells in the CFcase - looking for the highest and second highest gradient neighbors for each cell
	for (CFuint iCell = 0; iCell < nbCells; ++iCell) {
		geoData.idx = iCell;
		GeometricEntity *const currCell = _cellBuilder->buildGE();
		
		const vector<GeometricEntity*>& cellFaces = *currCell->getNeighborGeos();
		const CFuint nbFaces = cellFaces.size();
		State* state = currCell->getState(0);
		const CFuint stateID = state->getLocalID();
		
		CFreal maxGrad1 = 0.;
		CFreal maxGrad2 = 0.;
		
		for (CFuint iFace = 0; iFace < nbFaces; ++iFace) {
			// getting pointers to states
			State* stateA = cellFaces[iFace]->getState(0);
			State* stateB = cellFaces[iFace]->getState(1);
			
			const CFuint stateIDa = stateA->getLocalID();
			const CFuint stateIDb = stateB->getLocalID();
			
			CFint stateID0;
			CFint stateID1;
			
			State* state0 = CFNULL;
			State* state1 = CFNULL;
			
			if(stateIDa == stateID) {
				state0 = stateA;
				stateID0 = stateIDa;
				state1 = stateB;
				stateID1 = stateIDb;
			}
			else {
				state0 = stateB;
				stateID0 = stateIDb;
				state1 = stateA;
				stateID1 = stateIDa;
			}
			
			// gradient computation
			RealVector grad = *state0 - *state1;
			CFreal gradNorm = grad.norm2();
			
			// cell centers' distance computation
			RealVector& cellNode0 = state0->getCoordinates();
			RealVector& cellNode1 = state1->getCoordinates();
			RealVector dist = cellNode0 - cellNode1;
			CFreal distNorm = dist.norm2();
			
			gradNorm /= distNorm;
			
			if(gradNorm >= maxGrad1) {
				maxGrad2 = maxGrad1;
				maxGrad1 = gradNorm;
				
				neighbor2[iCell] = neighbor1[iCell];
				if(!state1->isGhost()) neighbor1[iCell] = stateID1;
				else neighbor1[iCell] = -1;
			}
			else if(gradNorm >= maxGrad2) {
				maxGrad2 = gradNorm;
				if(!state1->isGhost()) neighbor2[iCell] = stateID1;
				else neighbor2[iCell] = -1;
			}
		}
		_cellBuilder->releaseGE();
	}
	// two highest gradient neighbors have been found
	cout << "Two highest gradients found" << endl;
	
	CFuint nbLines = 0;
	for (CFuint iCell = 0; iCell < nbCells; ++iCell) {
		CFuint firstCellInLine;
		
		if(lines[iCell] < 0) {
			lines[iCell] = nbLines;
			firstCellInLine = iCell;
			
			CFuint currCellID = iCell;
			// forward line creation
			bool flag = true;
			do {
				const CFint neighborCell1 = neighbor1[currCellID];
				if (neighborCell1 < 0) {
					cout << "break1" << endl;
					break;
				}
				
				if((neighbor1[neighborCell1] == currCellID || neighbor2[neighborCell1] == currCellID) && lines[neighborCell1] < 0) {
					lines[neighborCell1] = nbLines;
					currCellID = neighborCell1;
				}
				else flag = false;
			} while(flag);
			
			// backward line creation
			currCellID = firstCellInLine;
			const CFint neighborCell2 = neighbor2[currCellID];
			
			if(neighborCell2 > 0 && ((neighbor1[neighborCell2] == currCellID || neighbor2[neighborCell2] == currCellID) && lines[neighborCell2] < 0)) {
				lines[neighborCell2] = nbLines;
				currCellID = neighborCell2;
				flag = true;
				do {
					const CFint neighborCell1 = neighbor1[currCellID];
					if (neighborCell1 < 0) {
						cout << "break3" << endl;
						break;
					}
					
					if((neighbor1[neighborCell1] == currCellID || neighbor2[neighborCell1] == currCellID) && lines[neighborCell1] < 0) {
						lines[neighborCell1] = nbLines;
						currCellID = neighborCell1;
					}
					else flag = false;
				} while(flag);
			}
			nbLines++;
		}
	}
	cout << "nbLines = " << nbLines << endl;
	*/
	
	DataHandle<bool> dplrIsFirstInLine = socket_dplrIsFirstInLine.getDataHandle();
	DataHandle<CFint> dplrToLocalIDs = socket_dplrToLocalIDs.getDataHandle();
	DataHandle<CFint> localToDplrIDs = socket_localToDplrIDs.getDataHandle();
	DataHandle<CFint> dplrCellInLine = socket_dplrCellInLine.getDataHandle();
	
	dplrIsFirstInLine = false;
	dplrToLocalIDs = -1;
	localToDplrIDs = -1;
	dplrCellInLine = -1;
	
	CFint prevCellID = -1;
	CFint currCellID = -1;
	CFint nextCellID = -1;
	
	CFreal maxGrad = 0.0;
	
	CFint nbLines = 0; // number of lines created in the grid
	CFint nbCellsInLine = 0; // number of cells in line
	CFint nbCellsInLines = 0; // number of cells in line
	
	// loop over all cells in the CFcase
  for (CFuint iCell = 0; iCell < nbCells; ++iCell) {
		currCellID = iCell;
		// check if this cell has been already visited
		if (dplrCellInLine[currCellID] < 0) {
			dplrIsFirstInLine[nbCellsInLines] = true;
			do {
				// build the GeometricEntity
				geoData.idx = currCellID;
				GeometricEntity *const currCell = _cellBuilder->buildGE();
				const vector<GeometricEntity*>& cellFaces = *currCell->getNeighborGeos();
				const CFuint nbFaces = cellFaces.size();
				
				State* state = currCell->getState(0);
				const CFuint stateID = state->getLocalID();
				
				nextCellID = -1;
				maxGrad = 0.0;
				for (CFuint iFace = 0; iFace < nbFaces; ++iFace) {
					// getting pointers to states
					State* stateA = cellFaces[iFace]->getState(0);
					State* stateB = cellFaces[iFace]->getState(1);
					
					const CFuint stateIDa = stateA->getLocalID();
					const CFuint stateIDb = stateB->getLocalID();
					
					State* state0 = CFNULL;
					State* state1 = CFNULL;
					
					CFint stateID0;
					CFint stateID1;
					
					if(stateIDa == stateID) {
						state0 = stateA;
						stateID0 = stateIDa;
						state1 = stateB;
						stateID1 = stateIDb;
					}
					else {
						state0 = stateB;
						stateID0 = stateIDb;
						state1 = stateA;
						stateID1 = stateIDa;
					}
					
					// check if the cell which is going to be connected has only one neighbor from the same line
					if (stateID1 != prevCellID && dplrCellInLine[stateID1] == nbLines) {
						nextCellID = -1;
						currCellID = -1;
						//cout << "break, stateID0 = " << stateID0 << " stateID1 = " << stateID1 << endl;
						break;
					}
					
					// check if the face is boundary face - if it is boundary face don't care
					if (!state1->isGhost() && stateID1 != prevCellID && dplrCellInLine[stateID1] < 0) { 
						RealVector& cellNode0 = state0->getCoordinates();
						RealVector& cellNode1 = state1->getCoordinates();
						
						// gradient computation
						RealVector grad = *state0 - *state1;
						CFreal gradNorm = grad.norm2();
						
						RealVector dist = cellNode0 - cellNode1;
						CFreal distNorm = dist.norm2();
						
						gradNorm /= distNorm;
						
						if(gradNorm >= maxGrad) {
							maxGrad = gradNorm;
							nextCellID = stateID1;
						}
					}
				}
				
				_cellBuilder->releaseGE();
				
				if(currCellID >= 0) {
					//cout << "Line" << nbLines << /* " StateLocalID = " << stateID << */ " currCellID = " << currCellID << /* " nextCellID = " << nextCellID << " prevCellID = " << prevCellID <<*/ endl;
					
					dplrCellInLine[currCellID] = nbLines;
					dplrToLocalIDs[nbCellsInLines] = currCellID;
					localToDplrIDs[currCellID] = nbCellsInLines;
					
					prevCellID = currCellID;
					currCellID = nextCellID;
					
					nbCellsInLine++;
					nbCellsInLines++;
				}
			} while(nextCellID > 0);
			
			nbCellsInLine = 0;
			nbLines++;
    }
  }
	cout << "nbLines = " << nbLines << "\n" << endl;

	/*
	CFuint aaa = 0;
	for(int i = 0; i < nbCells; i++){
		if(dplrIsFirstInLine[i]) {
			cout << "Next Line " << aaa << endl;
			aaa++;
		}
		//cout << dplrToLocalIDs[i] << endl;
		cout << "dplrToLocalIDs[" << i << "] = " << dplrToLocalIDs[i] << " localToDplrIDs[" << i << "] = " << localToDplrIDs[i] << endl;
	}
	cout << "nbLines2 = " << aaa << "\n" << endl;
	*/
	
	//cout << "nbCellsInLines = " << nbCellsInLines << " nbCells = " << nbCells << endl;
  // MHD2DPowellSourceTerm
}

//////////////////////////////////////////////////////////////////////////////

} // namespace FiniteVolume

} // namespace Numerics

} // namespace COOLFluiD
