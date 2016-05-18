#include "FiniteVolume/FiniteVolume.hh"
#include "FVMCC_ALETimeRhs.hh"
#include "Framework/CFL.hh"
#include "Framework/SubSystemStatus.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/BlockAccumulator.hh"
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

MethodCommandProvider<FVMCC_ALETimeRhs, CellCenterFVMData, FiniteVolumeModule> fvmccALETimeRhs("ALETimeRhs");

//////////////////////////////////////////////////////////////////////////////

FVMCC_ALETimeRhs::FVMCC_ALETimeRhs
(const std::string& name) :
  FVMCC_PseudoSteadyTimeRhs(name),
  socket_pastVolumes("pastVolumes"),
  _pastDiagValue(0.0)
{
}

//////////////////////////////////////////////////////////////////////////////

FVMCC_ALETimeRhs::~FVMCC_ALETimeRhs()
{
}

//////////////////////////////////////////////////////////////////////////////

void FVMCC_ALETimeRhs::setup()
{
  FVMCC_PseudoSteadyTimeRhs::setup();
}

//////////////////////////////////////////////////////////////////////////////

void FVMCC_ALETimeRhs::execute()
{
  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();
  DataHandle<CFreal> volumes = socket_volumes.getDataHandle();
  DataHandle<CFreal> pastVolumes = socket_pastVolumes.getDataHandle();
  DataHandle<CFreal> updateCoeff = socket_updateCoeff.getDataHandle();

  const CFuint nbStates = states.size();
  const CFreal dt = SubSystemStatusStack::getActive()->getDT();

//   CFreal minDt = 0.0;
//   if (_useGlobalDT) {
//     // select the minimum delta T
//     minDt = volumes[0]/updateCoeff[0];
//     for (CFuint iState = 0; iState < nbStates; ++iState) {
//       minDt = min(volumes[iState]/updateCoeff[iState], minDt);
//       minDt = min(pastVolumes[iState]/updateCoeff[iState], minDt);
//     }
//
//     //set DT to minDt!!!
//     SubSystemStatusStack::getActive()->setDT(minDt);
//   }

  // add the diagonal entries in the jacobian (updateCoeff/CFL)
  for (CFuint iState = 0; iState < nbStates; ++iState) {

    if (!_useGlobalDT) {
      cf_assert (dt > 0.0);
      _pastDiagValue = pastVolumes[iState]/dt;
      _diagValue = volumes[iState]/dt;
    }
    else {
      ///this is not allowed as we are doing unsteady computation
      cf_assert(false);
    }

    if (states[iState]->isParUpdatable()) {

      if (!_useAnalyticalMatrix) {
        // compute the transformation matrix numerically
        computeNumericalTransMatrix(iState);
      }
      else {
        computeAnalyticalTransMatrix(iState);
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void FVMCC_ALETimeRhs::computeNumericalTransMatrix(const CFuint iState)
{
  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();
  DataHandle<State*> pastStates = socket_pastStates.getDataHandle();
  DataHandle<CFreal> rhs = socket_rhs.getDataHandle();

  State *const currState = states[iState];

  // this first transformed state HAS TO BE stored,
  // since the returned pointer will change pointee after
  // the second call to transform()
  _tempState = static_cast<RealVector&>
    (*_updateToSolutionVecTrans->transform(currState));


  const State *const tPastState =
    _updateToSolutionVecTrans->transform(pastStates[iState]);

  const CFreal resFactor = getMethodData().getResFactor();

  // first the contribution to the rhs is computed
  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  for (CFuint iEq = 0; iEq < nbEqs; ++iEq) {
    const CFreal dU = _tempState[iEq]*_diagValue - (*tPastState)[iEq]*_pastDiagValue;
    rhs(iState, iEq, nbEqs) -= (!_zeroDiagValue[iEq]) ? dU*resFactor : 0.0;
  }

  if ((!getMethodData().isSysMatrixFrozen()) && getMethodData().doComputeJacobian()) {
    // now the contribution to the jacobian matrix is calculated
    _acc->setRowColIndex(0, currState->getLocalID());

    for (CFuint iVar = 0; iVar < nbEqs; ++iVar) {
      // perturb the given component of the state vector
      _numericalJacob->perturb(iVar, (*currState)[iVar]);

      const RealVector& tempPertState = static_cast<RealVector&>
        (*_updateToSolutionVecTrans->transform(currState));
      
      // compute the finite difference derivative of the flux
      _numericalJacob->computeDerivative(_tempState,
                                         tempPertState,
                                         _fluxDiff);
      
      for (CFuint iEq = 0; iEq < nbEqs; ++iEq) {
	_fluxDiff[iEq] *= (!_zeroDiagValue[iEq]) ? _diagValue*resFactor : 0.0;
      }
      
      _acc->addValues(0, 0, iVar, &_fluxDiff[0]);
      // restore the unperturbed value
      _numericalJacob->restore((*currState)[iVar]);
    }

    // add the values in the jacobian matrix
    _lss->getMatrix()->addValues(*_acc);

    // reset to zero the entries in the block accumulator
    _acc->reset();
  }
}

//////////////////////////////////////////////////////////////////////////////

void FVMCC_ALETimeRhs::computeAnalyticalTransMatrix(const CFuint iState)
{
  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();
  DataHandle<State*> pastStates = socket_pastStates.getDataHandle();
  DataHandle<CFreal> rhs = socket_rhs.getDataHandle();

  State *const currState = states[iState];

  const State& tempState =
    *_updateToSolutionVecTrans->transform(currState);

  const State& tPastState =
    *_updateToSolutionVecTrans->transform(pastStates[iState]);

  const CFreal resFactor = getMethodData().getResFactor();

  // first the contribution to the rhs is computed
  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  for (CFuint iEq = 0; iEq < nbEqs; ++iEq) {
    const CFreal dU = tempState[iEq] * _diagValue - tPastState[iEq] * _pastDiagValue;
    rhs(iState, iEq, nbEqs) -= (!_zeroDiagValue[iEq]) ? dU*resFactor : 0.0;
  }

  if ((!getMethodData().isSysMatrixFrozen()) && getMethodData().doComputeJacobian()) {
    // now the contribution to the jacobian matrix is calculated
    _acc->setRowColIndex(0, currState->getLocalID());

    _updateToSolutionInUpdateMatTrans->setMatrix(*states[iState]);
    const RealMatrix& matrix = *_updateToSolutionInUpdateMatTrans->getMatrix();

    for (CFuint iVar = 0; iVar < nbEqs; ++iVar) {
      for (CFuint jVar = 0; jVar < nbEqs; ++jVar) {
        if (!_zeroDiagValue[jVar]) {
	  const CFreal value = matrix(iVar,jVar)*_diagValue*resFactor;
	  _acc->addValue(0, 0, iVar, jVar, value);
	}
      }
    }

    // add the values in the jacobian matrix
    _lss->getMatrix()->addValues(*_acc);

    // reset to zero the entries in the block accumulator
    _acc->reset();
  }
}

//////////////////////////////////////////////////////////////////////////////

vector<SafePtr<BaseDataSocketSink> >
FVMCC_ALETimeRhs::needsSockets()
{
  vector<SafePtr<BaseDataSocketSink> > result =
    FVMCC_PseudoSteadyTimeRhs::needsSockets();

  result.push_back(&socket_pastVolumes);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD
