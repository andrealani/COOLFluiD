#include "FiniteVolume/FiniteVolume.hh"
#include "FVMCC_ALETimeRhsCoupling.hh"
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

MethodCommandProvider<FVMCC_ALETimeRhsCoupling, CellCenterFVMData, FiniteVolumeModule> fvmccALETimeRhsCoupling("ALETimeRhsCoupling");

//////////////////////////////////////////////////////////////////////////////

FVMCC_ALETimeRhsCoupling::FVMCC_ALETimeRhsCoupling
(const std::string& name) :
  FVMCC_PseudoSteadyTimeRhsCoupling(name),
  socket_pastVolumes("pastVolumes"),
  _pastDiagValue(0.0)
{
}

//////////////////////////////////////////////////////////////////////////////

FVMCC_ALETimeRhsCoupling::~FVMCC_ALETimeRhsCoupling()
{
}

//////////////////////////////////////////////////////////////////////////////

void FVMCC_ALETimeRhsCoupling::setup()
{
  FVMCC_PseudoSteadyTimeRhsCoupling::setup();
}

//////////////////////////////////////////////////////////////////////////////

void FVMCC_ALETimeRhsCoupling::execute()
{
  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();
  DataHandle<CFreal> volumes = socket_volumes.getDataHandle();
  DataHandle<CFreal> pastVolumes = socket_pastVolumes.getDataHandle();
  DataHandle<CFreal> updateCoeff = socket_updateCoeff.getDataHandle();

  const CFuint nbStates = states.size();
  const CFreal dt = SubSystemStatusStack::getActive()->getDT();
  const CFuint nbLSS = _lss.size();

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

  // reset the equation subsystem descriptor
  PhysicalModelStack::getActive()->resetEquationSubSysDescriptor();

  // add the diagonal entries in the jacobian (updateCoeff/CFL)
  for (CFuint iState = 0; iState < nbStates; ++iState) {

   if (states[iState]->isParUpdatable()) {

      // loop over the LSSs
      for (CFuint iLSS = 0; iLSS < nbLSS; ++iLSS) {

	const vector<CFuint>& currEqs = *_equations[iLSS];
	_nbEqs = currEqs.size();
	_start = currEqs[0];

	// set the equation subsystem descriptor
	PhysicalModelStack::getActive()->setEquationSubSysDescriptor
	  (_start, _nbEqs, iLSS);

 	if (!_useGlobalDT) {
          cf_assert (dt > 0.0);
          _pastDiagValue = pastVolumes[iState]/dt;
          _diagValue = volumes[iState]/dt;
	}
	else {
	  ///this is not allowed as we are doing unsteady computation
          cf_assert(false);
	}

	if (_annullDiagValue[iLSS]) {
	  _diagValue = 0.0;
	}

	if (!_useAnalyticalMatrix) {
	  // compute the transformation matrix numerically
	  computeNumericalTransMatrix(iState, iLSS);
	}
	else {
	  computeAnalyticalTransMatrix(iState, iLSS);
	}
      }
   }
  }
  
  PhysicalModelStack::getActive()->resetEquationSubSysDescriptor();
}

//////////////////////////////////////////////////////////////////////////////

void FVMCC_ALETimeRhsCoupling::computeNumericalTransMatrix
(const CFuint iState, const CFuint iLSS)
{
  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();
  DataHandle<State*> pastStates = socket_pastStates.getDataHandle();
  DataHandle<CFreal> rhs = socket_rhs.getDataHandle();
  
  State *const currState = states[iState];
  
  // this first transformed state HAS TO BE stored,
  // since the returned pointer will change pointee after
  // the second call to transform()
  _tempState.slice(_start, _nbEqs) = 
    (static_cast<RealVector&>(*_updateToSolutionVecTrans->transform(currState))).slice(_start, _nbEqs);
  
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
    rhs(iState, iVar, totalNbEqs) -= dU*resFactor;
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
      
      _fluxDiff.slice(_start, _nbEqs) *= _diagValue*resFactor;
      
      _acc[iLSS]->addValues(0, 0, iEq, &_fluxDiff[_start]);
      
      // restore the unperturbed value
      _numericalJacob->restore((*currState)[iVar]);
    }
    
    // add the values in the jacobian matrix
    _jacobMatrix[iLSS]->addValues(*_acc[iLSS]);
    
    // reset to zero the entries in the block accumulator
    _acc[iLSS]->reset();
  }
}

//////////////////////////////////////////////////////////////////////////////

void FVMCC_ALETimeRhsCoupling::computeAnalyticalTransMatrix
(const CFuint iState, const CFuint iLSS)
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
  const CFuint totalNbEqs = PhysicalModelStack::getActive()->getNbEq();
  const vector<CFuint>& currEqs = *_equations[iLSS];

  // first the contribution to the rhs is computed
  for (CFuint iEq = 0; iEq < _nbEqs; ++iEq) {
    const CFuint iVar = currEqs[iEq];
    const CFreal dU = tempState[iVar]*_diagValue - (tPastState)[iVar]*_pastDiagValue;
    rhs(iState, iVar, totalNbEqs) -= dU*resFactor;
  }

  if (!getMethodData().isSysMatrixFrozen()) {
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
    
    // add the values in the jacobian matrix
    _jacobMatrix[iLSS]->addValues(*_acc[iLSS]);
    
    // reset to zero the entries in the block accumulator
    _acc[iLSS]->reset();
  }
}

//////////////////////////////////////////////////////////////////////////////

vector<SafePtr<BaseDataSocketSink> >
FVMCC_ALETimeRhsCoupling::needsSockets()
{
  vector<SafePtr<BaseDataSocketSink> > result =
    FVMCC_PseudoSteadyTimeRhsCoupling::needsSockets();

  result.push_back(&socket_pastVolumes);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD
