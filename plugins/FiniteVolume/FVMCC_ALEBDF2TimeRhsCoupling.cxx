#include "FiniteVolume/FiniteVolume.hh"
#include "FVMCC_ALEBDF2TimeRhsCoupling.hh"
#include "Framework/CFL.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/BlockAccumulator.hh"
#include "Framework/LSSMatrix.hh"
#include "Framework/MeshData.hh"
#include "Framework/SubSystemStatus.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<FVMCC_ALEBDF2TimeRhsCoupling, CellCenterFVMData, FiniteVolumeModule> fvmccALEBDF2TimeRhsCoupling("ALEBDF2TimeRhsCoupling");

//////////////////////////////////////////////////////////////////////////////

FVMCC_ALEBDF2TimeRhsCoupling::FVMCC_ALEBDF2TimeRhsCoupling
(const std::string& name) :
  FVMCC_ALETimeRhsCoupling(name),
  socket_pastTimeRhs("pastTimeRhs")
{
}

//////////////////////////////////////////////////////////////////////////////

FVMCC_ALEBDF2TimeRhsCoupling::~FVMCC_ALEBDF2TimeRhsCoupling()
{
}

//////////////////////////////////////////////////////////////////////////////

void FVMCC_ALEBDF2TimeRhsCoupling::setup()
{
  FVMCC_ALETimeRhsCoupling::setup();
}

//////////////////////////////////////////////////////////////////////////////

void FVMCC_ALEBDF2TimeRhsCoupling::computeNumericalTransMatrix(const CFuint iState, const CFuint iLSS)
{
  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();
  DataHandle<State*> pastStates = socket_pastStates.getDataHandle();
  DataHandle<CFreal> rhs = socket_rhs.getDataHandle();
  DataHandle<CFreal> pastTimeRhs = socket_pastTimeRhs.getDataHandle();
  State *const currState = states[iState];

  // this first transformed state HAS TO BE stored,
  // since the returned pointer will change pointee after
  // the second call to transform()
  _tempState.slice(_start, _nbEqs) =
    (static_cast<RealVector&>(*_updateToSolutionVecTrans->
			      transform(currState))).slice(_start, _nbEqs);
  
  // here there is no overwork because the transformer performs a
  // partial transformation and returns a pointer to a partially changed
  // State
  const State *const tPastState =
    _updateToSolutionVecTrans->transform(pastStates[iState]);
  
  const CFreal resFactor = getMethodData().getResFactor();
  const CFuint totalNbEqs = PhysicalModelStack::getActive()->getNbEq();
  const vector<CFuint>& currEqs = *_equations[iLSS];

  // first the contribution to the rhs is computed
  for (CFuint iEq = 0; iEq < _nbEqs; ++iEq) {
    const CFuint iVar = currEqs[iEq];
    const CFreal dU = _tempState[iVar]*_diagValue - (*tPastState)[iVar]*_pastDiagValue;
    rhs(iState, iVar, totalNbEqs) -= dU*(1.+ resFactor) - resFactor*pastTimeRhs(iState, iVar, totalNbEqs);
    if((SubSystemStatusStack::getActive()->isLastStep())
    && (SubSystemStatusStack::getActive()->isSubIterationLastStep()))
    {
     pastTimeRhs(iState, iVar, totalNbEqs) = dU;
    }
  }

  if(!getMethodData().isSysMatrixFrozen()) {
  // now the contribution to the jacobian matrix is calculated
  _acc[iLSS]->setRowColIndex(0, currState->getLocalID());

  for (CFuint iEq = 0; iEq < _nbEqs; ++iEq) {
    const CFuint iVar = currEqs[iEq];
    // perturb the given component of the state vector
    _numericalJacob->perturb(iVar, (*currState)[iVar]);

    // here there is no overwork because the transformer performs a
    // partial transformation and returns a pointer to a partially changed
    // RealVector
    RealVector& tempPertState = static_cast<RealVector&>
      (*_updateToSolutionVecTrans->transform(currState));
    
    // compute the finite difference derivative of the flux
    _numericalJacob->computeDerivative(_tempState.slice(_start, _nbEqs),
				       tempPertState.slice(_start, _nbEqs),
				       _fluxDiff.slice(_start, _nbEqs));

    _fluxDiff.slice(_start, _nbEqs) *= _diagValue*(1.+ resFactor);

    _acc[iLSS]->addValues(0, 0, iEq, &_fluxDiff[_start]);

    // restore the unperturbed value
    _numericalJacob->restore((*currState)[iVar]);
  }
}
}

//////////////////////////////////////////////////////////////////////////////

void FVMCC_ALEBDF2TimeRhsCoupling::computeAnalyticalTransMatrix
(const CFuint iState, const CFuint iLSS)
{
  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();
  DataHandle<State*> pastStates = socket_pastStates.getDataHandle();
  DataHandle<CFreal> rhs = socket_rhs.getDataHandle();
  DataHandle<CFreal> pastTimeRhs = socket_pastTimeRhs.getDataHandle();
  State *const currState = states[iState];

  const State& tempState =
    *_updateToSolutionVecTrans->transform(currState);

  const State& tPastState =
    *_updateToSolutionVecTrans->transform(pastStates[iState]);

  const CFreal resFactor = getMethodData().getResFactor();
  const CFuint totalNbEqs = PhysicalModelStack::getActive()->getNbEq();
  const vector<CFuint>& currEqs = *_equations[iLSS];

  // first the contribution to the rhs is computed
  for (CFuint iEq = 0; iEq < _nbEqs; ++iEq) {
    const CFuint iVar = currEqs[iEq];
    const CFreal dU = tempState[iVar]*_diagValue - (tPastState)[iVar]*_pastDiagValue;
    rhs(iState, iVar, totalNbEqs) -= dU*(1.+ resFactor) - resFactor*pastTimeRhs(iState, iVar, totalNbEqs);
    if((SubSystemStatusStack::getActive()->isLastStep())
    && (SubSystemStatusStack::getActive()->isSubIterationLastStep()))
    {
     pastTimeRhs(iState, iVar, totalNbEqs) = dU;
    }
  }

  if(!getMethodData().isSysMatrixFrozen()) {
  // now the contribution to the jacobian matrix is calculated
  _acc[iLSS]->setRowColIndex(0, currState->getLocalID());

  _updateToSolutionInUpdateMatTrans->setMatrix(*states[iState]);
  const RealMatrix& matrix = *_updateToSolutionInUpdateMatTrans->getMatrix();

  for (CFuint iEq = 0; iEq < _nbEqs; ++iEq) {
    const CFuint iVar = currEqs[iEq];
    for (CFuint jEq = 0; jEq < _nbEqs; ++jEq) {
      const CFuint jVar = currEqs[jEq];
        const CFreal value = matrix(iVar,jVar)*_diagValue*resFactor;
      _acc[iLSS]->addValue(0, 0, iEq, jEq, value);
    }
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

vector<SafePtr<BaseDataSocketSink> >
FVMCC_ALEBDF2TimeRhsCoupling::needsSockets()
{
  vector<SafePtr<BaseDataSocketSink> > result =
    FVMCC_ALETimeRhsCoupling::needsSockets();

  result.push_back(&socket_pastTimeRhs);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD
