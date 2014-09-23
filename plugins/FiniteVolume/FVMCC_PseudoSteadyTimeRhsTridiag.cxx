#include "FiniteVolume/FiniteVolume.hh"
#include "FVMCC_PseudoSteadyTimeRhsTridiag.hh"
#include "Framework/MethodCommandProvider.hh"
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

MethodCommandProvider<FVMCC_PseudoSteadyTimeRhsTridiag,
		      CellCenterFVMData,
		      FiniteVolumeModule>
fvmccPseudoSteadyTimeRhsTridiagProvider("PseudoSteadyTimeRhsTridiag");

//////////////////////////////////////////////////////////////////////////////

FVMCC_PseudoSteadyTimeRhsTridiag::FVMCC_PseudoSteadyTimeRhsTridiag
(const std::string& name) :
  FVMCC_PseudoSteadyTimeRhs(name),
  socket_diagMatrices("diagMatrices"),
  socket_upLocalIDsAll("upLocalIDsAll"),
  _matIter()
{
}

//////////////////////////////////////////////////////////////////////////////

FVMCC_PseudoSteadyTimeRhsTridiag::~FVMCC_PseudoSteadyTimeRhsTridiag()
{
}

//////////////////////////////////////////////////////////////////////////////

void FVMCC_PseudoSteadyTimeRhsTridiag::setup()
{
  FVMCC_PseudoSteadyTimeRhs::setup();

  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  _matIter.resize(nbEqs,nbEqs,false);
}

//////////////////////////////////////////////////////////////////////////////

void FVMCC_PseudoSteadyTimeRhsTridiag::computeNumericalTransMatrix(const CFuint iState)
{
  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();
  DataHandle<State*> pastStates = socket_pastStates.getDataHandle();
  DataHandle< CFreal> rhs = socket_rhs.getDataHandle();
  DataHandle< CFreal> diagMatrices = socket_diagMatrices.getDataHandle();
  DataHandle<CFint> upLocalIDsAll = socket_upLocalIDsAll.getDataHandle();

  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  const CFuint nbEqs2 = nbEqs*nbEqs;

  State *const currState = states[iState];
  const CFint stateID = upLocalIDsAll[currState->getLocalID()];
  cf_assert(stateID >= 0);

  _matIter.wrap(nbEqs, nbEqs, &diagMatrices[stateID*nbEqs2]);
  
  // this 1first transformed state HAS TO BE stored,
  // since the returned pointer will change pointee after
  // the second call to transform()
  _tempState = static_cast<RealVector&>
    (*_updateToSolutionVecTrans->transform(currState));

  const State *const tPastState =
    _updateToSolutionVecTrans->transform(pastStates[iState]);

  const CFreal resFactor = getMethodData().getResFactor();

  // first the contribution to the rhs is computed
  for (CFuint iEq = 0; iEq < nbEqs; ++iEq) {
    const CFreal dU = _tempState[iEq] - (*tPastState)[iEq];
    rhs(iState, iEq, nbEqs) -= (!_zeroDiagValue[iEq]) ? dU*_diagValue*resFactor : 0.0;
  }

  if ((!getMethodData().isSysMatrixFrozen()) && getMethodData().doComputeJacobian()) {
    for (CFuint iVar = 0; iVar < nbEqs; ++iVar) {
      // perturb the given component of the state vector
      _numericalJacob->perturb(iVar, (*currState)[iVar]);

      const RealVector& tempPertState = static_cast<RealVector&>
        (*_updateToSolutionVecTrans->transform(currState));

      // compute the finite difference derivative of the flux
      _numericalJacob->computeDerivative(_tempState, tempPertState, _fluxDiff);

      // _fluxDiff corresponds to a column vector of the dU/dP matrix
      for (CFuint iEq = 0; iEq < nbEqs; ++iEq) {
        _fluxDiff[iEq] *= (!_zeroDiagValue[iEq]) ? _diagValue*resFactor : 0.0;
      }

      _matIter.addColumn(_fluxDiff, iVar);

      // restore the unperturbed value
      _numericalJacob->restore((*currState)[iVar]);
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void FVMCC_PseudoSteadyTimeRhsTridiag::computeAnalyticalTransMatrix
(const CFuint iState)
{
  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();
  DataHandle<State*> pastStates = socket_pastStates.getDataHandle();
  DataHandle< CFreal> rhs = socket_rhs.getDataHandle();
  DataHandle< CFreal> diagMatrices = socket_diagMatrices.getDataHandle();
  DataHandle<CFint> upLocalIDsAll = socket_upLocalIDsAll.getDataHandle();

  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  const CFuint nbEqs2 = nbEqs*nbEqs;

  State *const currState = states[iState];
  const CFuint stateID = upLocalIDsAll[currState->getLocalID()];
  
  _matIter.wrap(nbEqs, nbEqs, &diagMatrices[stateID*nbEqs2]);
  
  const State& tempState =
    *_updateToSolutionVecTrans->transform(currState);

  const State& tPastState =
    *_updateToSolutionVecTrans->transform(pastStates[iState]);

  const CFreal resFactor = getMethodData().getResFactor();

  // first the contribution to the rhs is computed
  for (CFuint iEq = 0; iEq < nbEqs; ++iEq) {
    const CFreal dU = tempState[iEq] - tPastState[iEq];
    rhs(iState, iEq, nbEqs) -= (!_zeroDiagValue[iEq]) ? dU*_diagValue*resFactor : 0.0;
  }

  if ((!getMethodData().isSysMatrixFrozen()) && getMethodData().doComputeJacobian()) {
    // now the contribution to the jacobian matrix is calculated

    _updateToSolutionInUpdateMatTrans->setMatrix(*states[iState]);
    const RealMatrix& matrix = *_updateToSolutionInUpdateMatTrans->getMatrix();

    for (CFuint iVar = 0; iVar < nbEqs; ++iVar) {
      for (CFuint jVar = 0; jVar < nbEqs; ++jVar) {
	if (!_zeroDiagValue[jVar]) {
	  const CFreal value = matrix(iVar,jVar)*_diagValue*resFactor;
	  _matIter(iVar,jVar) += value;
	}
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

vector<SafePtr<BaseDataSocketSink> > FVMCC_PseudoSteadyTimeRhsTridiag::needsSockets()
{
  vector<SafePtr<BaseDataSocketSink> > result =
    FVMCC_PseudoSteadyTimeRhs::needsSockets();

  result.push_back(&socket_diagMatrices);
  result.push_back(&socket_upLocalIDsAll);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD
