#include "FiniteVolume/FiniteVolume.hh"
#include "FVMCC_BDF2TimeRhs.hh"
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

MethodCommandProvider<FVMCC_BDF2TimeRhs,
                      CellCenterFVMData,
                      FiniteVolumeModule>
                      fvmccBDF2TimeRhs("BDF2TimeRhs");

//////////////////////////////////////////////////////////////////////////////

FVMCC_BDF2TimeRhs::FVMCC_BDF2TimeRhs(const std::string& name) :
  FVMCC_PseudoSteadyTimeRhs(name),
  socket_pastTimeRhs("pastTimeRhs")
{
}

//////////////////////////////////////////////////////////////////////////////

FVMCC_BDF2TimeRhs::~FVMCC_BDF2TimeRhs()
{
}

//////////////////////////////////////////////////////////////////////////////

void FVMCC_BDF2TimeRhs::setup()
{
  FVMCC_PseudoSteadyTimeRhs::setup();
}

//////////////////////////////////////////////////////////////////////////////

void FVMCC_BDF2TimeRhs::computeNumericalTransMatrix(const CFuint iState)
{
  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();
  DataHandle<State*> pastStates = socket_pastStates.getDataHandle();
  DataHandle<CFreal> rhs = socket_rhs.getDataHandle();
  DataHandle<CFreal> pastTimeRhs = socket_pastTimeRhs.getDataHandle();
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
    const CFreal dU = _tempState[iEq]*_diagValue - (*tPastState)[iEq]*_diagValue;
    rhs(iState, iEq, nbEqs) -= (!_zeroDiagValue[iEq]) ?
      dU*(1.+ resFactor) - resFactor*pastTimeRhs(iState, iEq, nbEqs) : 0.0;

    if((SubSystemStatusStack::getActive()->isLastStep())
       && (SubSystemStatusStack::getActive()->isSubIterationLastStep()))
      {
	pastTimeRhs(iState, iEq, nbEqs) = dU;
      }
  }
  
  if (getMethodData().doComputeJacobian()) {
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
      
      // _fluxDiff corresponds to a column vector of the dU/dP matrix
      for (CFuint iEq = 0; iEq < nbEqs; ++iEq) {
	_fluxDiff[iEq] *= (!_zeroDiagValue[iEq]) ? _diagValue*(1.+ resFactor) : 0.0;
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

void FVMCC_BDF2TimeRhs::computeAnalyticalTransMatrix(const CFuint iState)
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

  // first the contribution to the rhs is computed
  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  for (CFuint iEq = 0; iEq < nbEqs; ++iEq) {
    const CFreal dU = tempState[iEq] * _diagValue - tPastState[iEq] * _diagValue;
    rhs(iState, iEq, nbEqs) -= dU*(1.+ resFactor) - resFactor*pastTimeRhs(iState, iEq, nbEqs);
    if(SubSystemStatusStack::getActive()->isLastStep()) pastTimeRhs(iState, iEq, nbEqs) = dU;
  }
  
  if (getMethodData().doComputeJacobian()) {
    // now the contribution to the jacobian matrix is calculated
    _acc->setRowColIndex(0, currState->getLocalID());
    
    _updateToSolutionInUpdateMatTrans->setMatrix(*states[iState]);
    const RealMatrix& matrix = *_updateToSolutionInUpdateMatTrans->getMatrix();
    
    for (CFuint iVar = 0; iVar < nbEqs; ++iVar) {
      for (CFuint jVar = 0; jVar < nbEqs; ++jVar) {
	if (!_zeroDiagValue[jVar]) {
	  const CFreal value = matrix(iVar,jVar)*_diagValue*(1. + resFactor);
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
FVMCC_BDF2TimeRhs::needsSockets()
{
  vector<SafePtr<BaseDataSocketSink> > result =
    FVMCC_PseudoSteadyTimeRhs::needsSockets();

  result.push_back(&socket_pastTimeRhs);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD
