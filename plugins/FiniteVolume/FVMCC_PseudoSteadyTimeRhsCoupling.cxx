#include "FiniteVolume/FiniteVolume.hh"
#include "FVMCC_PseudoSteadyTimeRhsCoupling.hh"
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

MethodCommandProvider<FVMCC_PseudoSteadyTimeRhsCoupling, CellCenterFVMData, FiniteVolumeModule>
fvmccPseudoSteadyTimeRhsCoupling("PseudoSteadyTimeRhsCoupling");

//////////////////////////////////////////////////////////////////////////////

FVMCC_PseudoSteadyTimeRhsCoupling::FVMCC_PseudoSteadyTimeRhsCoupling
(const std::string& name) :
  FVMCC_StdComputeTimeRhsCoupling(name),
  _numericalJacob(CFNULL),
  socket_rhs("rhs"),
  socket_pastStates("pastStates"),
  socket_volumes("volumes"),
  _updateToSolutionVecTrans(CFNULL),
  _updateToSolutionInUpdateMatTrans(CFNULL),
  _fluxDiff(),
  _tempState(),
  _tempPertState(),
  _acc(),
  _diagValue(0.0),
  _nbEqs(0),
  _start(0)
{
  
}

//////////////////////////////////////////////////////////////////////////////

FVMCC_PseudoSteadyTimeRhsCoupling::~FVMCC_PseudoSteadyTimeRhsCoupling()
{
  for (CFuint i = 0; i < _acc.size(); ++i) {
    deletePtr(_acc[i]);
  }
}

//////////////////////////////////////////////////////////////////////////////

void FVMCC_PseudoSteadyTimeRhsCoupling::setup()
{
  FVMCC_StdComputeTimeRhsCoupling::setup();

  _fluxDiff.resize(PhysicalModelStack::getActive()->getNbEq());
  _tempState.resize(PhysicalModelStack::getActive()->getNbEq());
  _tempPertState.resize(PhysicalModelStack::getActive()->getNbEq());

  // numerical jacobian
  _numericalJacob = &getMethodData().getNumericalJacobian();

  _updateToSolutionVecTrans =
    getMethodData().getUpdateToSolutionVecTrans();

  _updateToSolutionInUpdateMatTrans =
    getMethodData().getUpdateToSolutionInUpdateMatTrans();

  const CFuint nbLSS = getMethodData().getLinearSystemSolver().size();
  _acc.resize(nbLSS);
  for (CFuint i = 0; i < nbLSS; ++i) {
    const CFuint nbSysEq = _lss[i]->getNbSysEqs();
    _acc[i] = _lss[i]->createBlockAccumulator(1, 1, nbSysEq);
  }

  if (_annullDiagValue.size() == 0) {
    _annullDiagValue.resize(nbLSS);
    _annullDiagValue.assign(nbLSS, false);
  }
  else if (_annullDiagValue.size() != nbLSS) {
    CFout << "WATCH OUT: FVMCC_StdComputeTimeRhsCoupling::setup() => _annullDiagValue.size() != nbLSS "
	  << "_annullDiagValue assigned to FALSE \n";
    _annullDiagValue.resize(nbLSS);
    _annullDiagValue.assign(nbLSS, false);
  }
}

//////////////////////////////////////////////////////////////////////////////

void FVMCC_PseudoSteadyTimeRhsCoupling::execute()
{
  CFLog(VERBOSE, "FVMCC_PseudoSteadyTimeRhsCoupling::execute() => START\n");
  
  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();
  DataHandle< CFreal> updateCoeff = socket_updateCoeff.getDataHandle();

  const CFuint nbStates = states.size();
  const CFreal cfl = getMethodData().getCFL()->getCFLValue();
  DataHandle<CFreal> volumes = socket_volumes.getDataHandle();
  const CFreal dt = SubSystemStatusStack::getActive()->getDT();
  const CFuint nbLSS = _lss.size();

  CFreal minDt = 0.0;
  if (_useGlobalDT) {
    // select the minimum delta T
    minDt = volumes[0]/updateCoeff[0];
    for (CFuint iState = 0; iState < nbStates; ++iState) {
      // loop over the LSSs
      for (CFuint iLSS = 0; iLSS < nbLSS; ++iLSS) {
	minDt = min(volumes[iState]/updateCoeff[iState*nbLSS + iLSS], minDt);
      }
    }
  }

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
	  _diagValue = (dt > 0.0) ? volumes[iState]/dt :
	    updateCoeff[iState*nbLSS + iLSS]/cfl;
	}
	else {
	  _diagValue = volumes[iState]/(minDt*cfl);
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
	
	if (getMethodData().doComputeJacobian()) {
	  // add the values in the jacobian matrix
	  _jacobMatrix[iLSS]->addValues(*_acc[iLSS]);
	  
	  // _acc[iLSS]->print(); 
	  
	  // reset to zero the entries in the block accumulator
	  _acc[iLSS]->reset();
	}
      }
    }
  }
  
  PhysicalModelStack::getActive()->resetEquationSubSysDescriptor();
  
  CFLog(VERBOSE, "FVMCC_PseudoSteadyTimeRhsCoupling::execute() => END\n");
}
      
//////////////////////////////////////////////////////////////////////////////

void FVMCC_PseudoSteadyTimeRhsCoupling::computeNumericalTransMatrix
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
    const CFreal dU = _tempState[iVar] - (*tPastState)[iVar];
    rhs(iState, iVar, totalNbEqs) -= dU*_diagValue*resFactor;
  }
  
  if (getMethodData().doComputeJacobian()) {
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
  }
}
      
//////////////////////////////////////////////////////////////////////////////

void FVMCC_PseudoSteadyTimeRhsCoupling::computeAnalyticalTransMatrix
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
    const CFreal dU = tempState[iVar] - tPastState[iVar];
    rhs(iState, iVar, totalNbEqs) -= dU*_diagValue*resFactor;
  }
  
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

//////////////////////////////////////////////////////////////////////////////

vector<SafePtr<BaseDataSocketSink> >
FVMCC_PseudoSteadyTimeRhsCoupling::needsSockets()
{
  vector<SafePtr<BaseDataSocketSink> > result =
    FVMCC_StdComputeTimeRhsCoupling::needsSockets();

  result.push_back(&socket_rhs);
  result.push_back(&socket_pastStates);
  result.push_back(&socket_volumes);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD
