#include "FiniteVolume/FiniteVolume.hh"
#include "FVMCC_ComputeRhsJacobDiagSimple.hh"
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

MethodCommandProvider<FVMCC_ComputeRhsJacobDiagSimple,
		      CellCenterFVMData,
		      FiniteVolumeModule>
fvmcc_computeRhsJacobDiagSimple("NumJacobDiagSimple");

//////////////////////////////////////////////////////////////////////////////

FVMCC_ComputeRhsJacobDiagSimple::FVMCC_ComputeRhsJacobDiagSimple
(const std::string& name) :
  FVMCC_ComputeRhsJacobDiag(name),
  _tmpJacob(),
  _avState(CFNULL)
{
}

//////////////////////////////////////////////////////////////////////////////

FVMCC_ComputeRhsJacobDiagSimple::~FVMCC_ComputeRhsJacobDiagSimple()
{
  deletePtr(_avState);
}

//////////////////////////////////////////////////////////////////////////////

void FVMCC_ComputeRhsJacobDiagSimple::setup()
{
  FVMCC_ComputeRhsJacobDiag::setup();

  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  _tmpJacob.resize(nbEqs, nbEqs);
  _avState = new State();
}

//////////////////////////////////////////////////////////////////////////////

void FVMCC_ComputeRhsJacobDiagSimple::computeBothJacobTerms()
{
  // set the index of the block corresponding to the current
  // state in the jacobian matrix
  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  const CFuint nbEqs2 = nbEqs*nbEqs;
  DataHandle<CFint> upLocalIDsAll = socket_upLocalIDsAll.getDataHandle();

  State& state0 = *_currFace->getState(0);
  State& state1 = *_currFace->getState(1);

  DataHandle<CFreal> diagMatrices = socket_diagMatrices.getDataHandle();
  const CFint stateID0 = upLocalIDsAll[state0.getLocalID()];
  const CFint stateID1 = upLocalIDsAll[state1.getLocalID()];
  cf_assert(stateID0 >= 0);
  cf_assert(stateID1 >= 0);

  _matIter0.resetPtr(&diagMatrices[stateID0*nbEqs2]);
  _matIter1.resetPtr(&diagMatrices[stateID1*nbEqs2]);

  RealVector& pData = PhysicalModelStack::getActive()->
    getImplementor()->getConvectiveTerm()->getPhysicalData();
  SafePtr<ConvectiveVarSet> solutionVars = getMethodData().getSolutionVar();
  SafePtr<VarSetMatrixTransformer> vs = getMethodData().getUpdateToSolutionInUpdateMatTrans();
  cf_assert(_faceIdx < _fluxData->faceAreas->size());
  const CFreal factor = 0.5*getResFactor()*(*_fluxData->faceAreas)[_faceIdx];

  *_avState = 0.5*(state0 + state1);

  // dFdU (U_0)
  _reconstrVar->computePhysicalData(state0, pData);
  solutionVars->computeProjectedJacobian(_fluxData->unitNormal, _tmpJacob);
  vs->setMatrix(state0);
  const CFreal absLambda0 = std::abs(_reconstrVar->getMaxEigenValue(*_avState, _fluxData->unitNormal));
  for (CFuint i = 0; i < nbEqs; ++i) {
    _tmpJacob(i,i) += absLambda0;
  }
  const RealVector& dUdP0 = (*vs->getMatrix());
  _matIter0 += _tmpJacob*dUdP0;
  _matIter0 *= factor;

  // dFdU (U_1)
  _reconstrVar->computePhysicalData(state1, pData);
  solutionVars->computeProjectedJacobian(_fluxData->unitNormal, _tmpJacob);
  vs->setMatrix(state1);

  const CFreal absLambda1 = std::abs(_reconstrVar->getMaxEigenValue(*_avState, _fluxData->unitNormal));
  for (CFuint i = 0; i < nbEqs; ++i) {
    _tmpJacob(i,i) += absLambda1;
  }
  const RealVector& dUdP1 = (*vs->getMatrix());
  _matIter1 += _tmpJacob*dUdP1;
  _matIter1 *= (-factor);

  //  if (_fluxData->hasDiffusiveTerm) {
  //       //
  //       }

  //     if (getMethodData().isAxisymmetric()) {
  //       _matIter0 *= (_rMid/(state0.getCoordinates())[YY]);
  //     }
}

//////////////////////////////////////////////////////////////////////////////

void FVMCC_ComputeRhsJacobDiagSimple::computeJacobTerm(const CFuint idx)
{
  cf_assert(_currFace->getState(idx)->isParUpdatable());

  // set the index of the block corresponding to the current
  // state in the jacobian matrix
  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  const CFuint nbEqs2 = nbEqs*nbEqs;

  DataHandle<CFreal> diagMatrices = socket_diagMatrices.getDataHandle();
  DataHandle<CFint> upLocalIDsAll = socket_upLocalIDsAll.getDataHandle();

  RealVector& pData = PhysicalModelStack::getActive()->getImplementor()->
    getConvectiveTerm()->getPhysicalData();
  SafePtr<ConvectiveVarSet> solutionVars = getMethodData().getSolutionVar();
  SafePtr<VarSetMatrixTransformer> vs = getMethodData().getUpdateToSolutionInUpdateMatTrans();
  const CFreal factor = 0.5*pow(-1.,static_cast<CFreal>(idx))*getResFactor()*
    (*_fluxData->faceAreas)[_faceIdx];

  if (idx == 0) {
    // dFdU (U_0)
    State& state0 = *_currFace->getState(0);
    const CFint stateID0 = upLocalIDsAll[state0.getLocalID()];
    cf_assert(stateID0 >= 0);

    _matIter0.resetPtr(&diagMatrices[stateID0*nbEqs2]);

    _reconstrVar->computePhysicalData(state0, pData);
    solutionVars->computeProjectedJacobian(_fluxData->unitNormal, _tmpJacob);
    vs->setMatrix(state0);

    const CFreal absLambda0 = std::abs(_reconstrVar->getMaxEigenValue(state0, _fluxData->unitNormal));
    for (CFuint i = 0; i < nbEqs; ++i) {
      _tmpJacob(i,i) += absLambda0;
    }

    _matIter0 = _tmpJacob*(*vs->getMatrix());
    _matIter0 *= factor;
  }

  if (idx == 1) {
    // dFdU (U_1)
    State& state1 = *_currFace->getState(1);
    const CFint stateID1 = upLocalIDsAll[state1.getLocalID()];
    cf_assert(stateID1 >= 0);

    _matIter1.resetPtr(&diagMatrices[stateID1*nbEqs2]);

    _reconstrVar->computePhysicalData(state1, pData);
    solutionVars->computeProjectedJacobian(_fluxData->unitNormal, _tmpJacob);
    vs->setMatrix(state1);
    const CFreal absLambda1 = std::abs(_reconstrVar->getMaxEigenValue(state1, _fluxData->unitNormal));
    for (CFuint i = 0; i < nbEqs; ++i) {
      _tmpJacob(i,i) += absLambda1;
    }

    _matIter1 = _tmpJacob*(*vs->getMatrix());
    _matIter1 *= factor;
  }
}


//////////////////////////////////////////////////////////////////////////////

void FVMCC_ComputeRhsJacobDiagSimple::computeBoundaryJacobianTerm()
{
  // set the index of the block corresponding to the current
  // state in the jacobian matrix
//   State& state0 = *_currFace->getState(0);

  DataHandle<CFreal> diagMatrices = socket_diagMatrices.getDataHandle();
  DataHandle<CFint> upLocalIDsAll = socket_upLocalIDsAll.getDataHandle();

  // set the index of the block corresponding to the current
  // state in the jacobian matrix
  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  const CFuint nbEqs2 = nbEqs*nbEqs;

  State& currState = *_currFace->getState(0);
  State& ghostState = *_currFace->getState(1);
  cf_assert(ghostState.isGhost());

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

} // namespace FiniteVolume

} // namespace Numerics

} // namespace COOLFluiD
