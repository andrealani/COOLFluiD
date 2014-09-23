#include "FiniteVolume/FiniteVolume.hh"
#include "FVMCC_PseudoSteadyTimeRhsTriGM.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/MeshData.hh"

#include "Framework/CFL.hh"
#include "Framework/SubSystemStatus.hh"
#include "Framework/BlockAccumulator.hh"
#include "Framework/LSSMatrix.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<FVMCC_PseudoSteadyTimeRhsTriGM, 
		      CellCenterFVMData, 
		      FiniteVolumeModule> 
fvmccPseudoSteadyTimeRhsTriGM("PseudoSteadyTimeRhsTriGM");

//////////////////////////////////////////////////////////////////////////////

FVMCC_PseudoSteadyTimeRhsTriGM::FVMCC_PseudoSteadyTimeRhsTriGM
(const std::string& name) :
  FVMCC_PseudoSteadyTimeRhs(name),
  socket_diagMatrices("diagMatrices"),
  socket_upLocalIDsAll("upLocalIDsAll"),
  _matIter()
{
}

//////////////////////////////////////////////////////////////////////////////

FVMCC_PseudoSteadyTimeRhsTriGM::~FVMCC_PseudoSteadyTimeRhsTriGM()
{
}

//////////////////////////////////////////////////////////////////////////////

void FVMCC_PseudoSteadyTimeRhsTriGM::setup()
{
  FVMCC_PseudoSteadyTimeRhs::setup();

  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  _matIter.resize(nbEqs,nbEqs,false);
}

//////////////////////////////////////////////////////////////////////////////

void FVMCC_PseudoSteadyTimeRhsTriGM::computeNumericalTransMatrix(const CFuint iState)
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
  
  // this first transformed state HAS TO BE stored,
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
    // now the contribution to the jacobian matrix is calculated
    _acc->setRowColIndex(0, currState->getLocalID());

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

      _acc->addValues(0, 0, iVar, &_fluxDiff[0]);
      _matIter.addColumn(_fluxDiff, iVar);

      // restore the unperturbed value
      _numericalJacob->restore((*currState)[iVar]);
    }

    // add the values in the jacobian matrix
    getMethodData().getLSSMatrix(0)->addValues(*_acc);

    // reset to zero the entries in the block accumulator
    _acc->reset();
  }
}

//////////////////////////////////////////////////////////////////////////////

void FVMCC_PseudoSteadyTimeRhsTriGM::computeAnalyticalTransMatrix
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
    getMethodData().getLSSMatrix(0)->addValues(*_acc);

    // reset to zero the entries in the block accumulator
    _acc->reset();
  }
}

//////////////////////////////////////////////////////////////////////////////

vector<SafePtr<BaseDataSocketSink> >
FVMCC_PseudoSteadyTimeRhsTriGM::needsSockets()
{
  vector<SafePtr<BaseDataSocketSink> > result =
    FVMCC_StdComputeTimeRhs::needsSockets();

  result.push_back(&socket_rhs);
  result.push_back(&socket_pastStates);
  result.push_back(&socket_volumes);

  result.push_back(&socket_diagMatrices);
  result.push_back(&socket_upLocalIDsAll);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD
