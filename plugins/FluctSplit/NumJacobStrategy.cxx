#include "FluctSplit/NumJacobStrategy.hh"
#include "FluctSplit/ComputeDiffusiveTerm.hh"

#include "Framework/MeshData.hh"
#include "Framework/MethodStrategyProvider.hh"
#include "Framework/BlockAccumulator.hh"
#include "FluctSplit/FluctSplit.hh"
#include "FluctSplit/FluctuationSplitData.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

MethodStrategyProvider<NumJacobStrategy,
                       FluctuationSplitData,
                       ComputeJacobStrategy,
                       FluctSplitModule>
                       numJacobStrategyProvider("Numerical");

//////////////////////////////////////////////////////////////////////////////

NumJacobStrategy::NumJacobStrategy(const std::string& name) :
  ComputeJacobStrategy(name),
  _tempRes(),
  _otherResidual(0),
  _tResidual(0)
{
}

//////////////////////////////////////////////////////////////////////////////

NumJacobStrategy::~NumJacobStrategy()
{
}

//////////////////////////////////////////////////////////////////////////////

void NumJacobStrategy::setup()
{
  // first call parent method
  ComputeJacobStrategy::setup();

  _tempRes.resize(PhysicalModelStack::getActive()->getNbEq());
  _otherResidual.resize(MeshDataStack::getActive()->Statistics().getMaxNbStatesInCell());
  for (CFuint i = 0; i < _otherResidual.size(); ++i) {
    _otherResidual[i].resize(PhysicalModelStack::getActive()->getNbEq());
    _otherResidual[i] = 0.;
  }

  _tResidual.resize(MeshDataStack::getActive()->Statistics().getMaxNbStatesInCell());
  for (CFuint i = 0; i < _tResidual.size(); ++i) {
    _tResidual[i].resize(PhysicalModelStack::getActive()->getNbEq());
    _tResidual[i] = 0.;
  }
}

//////////////////////////////////////////////////////////////////////////////

void NumJacobStrategy::computeJacobianTerm
(GeometricEntity *const cell,
 const vector<RealVector>& residual,
 BlockAccumulator *const acc,
 const std::vector<CFuint>&  equationIDs)
{
  getMethodData().switchToImplicit();

  vector<State*> *const states = cell->getStates();
  
  NumericalJacobian& numericalJacob = getMethodData().getNumericalJacobian();
  
  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  const CFuint nbStatesInCell = states->size();

  SafePtr<ConvectiveVarSet> distribVar =
    getMethodData().getDistribVar();

  SafePtr<ComputeDiffusiveTerm> diffTermComputer =
    getMethodData().getDiffusiveTermComputer();

  // here we store the residual (referred to solution variables)
  _tResidual = residual;

  cf_assert( _fsStrategy.isNotNull() );
  
  DistributionData& dd = getMethodData().getDistributionData();
  dd.cell   = const_cast<GeometricEntity*>(cell);
  dd.cellID = cell->getID();
  dd.states = states;
  
  // compute the perturbed states for the evaluation of the
  // jacobian matrix looping over the state vectors in this cell
  for (CFuint iState = 0; iState < nbStatesInCell; ++iState) {
    State *const currState = (*states)[iState];

    CFLogDebugMed( "Perturbing iState = "
                   << iState
                   << " with stateID = "
                   << currState->getLocalID() << "\n");
    
    // set the index of the block corresponding to the current
    // state in the jacobian matrix
    acc->setRowColIndex(iState, currState->getLocalID());

    // loop over the variables in the state vector to perturb one
    // component at a time
    for (CFuint iVar = 0; iVar < nbEqs; ++iVar)
    {
      CFLogDebugMax( "Perturbing iVar = " << iVar << "\n");
      
      dd.iVar = iVar;
      
      // reset the residual to 0.
      cleanOtherResidual();
      
      // perturb the given component of the state vector
      numericalJacob.perturb(iVar, (*currState)[iVar]);
      
      _fsStrategy->computeFluctuation(_otherResidual);
      
      vector<RealVector> *const tBackResidual =
        _distToSolutionMatTrans->transformFromRef(&_otherResidual);
      
      if (_hasDiffusiveTerm) {
	// reset the residual to 0.
	// cleanDiffResidual();
	
	diffTermComputer->computeDiffusiveTerm(cell, _diffResidual, false);
	
	for (CFuint i = 0; i < nbStatesInCell; ++i) {
	  // CFLog(INFO, "tBackResidual[" << i << "] = " << (*tBackResidual)[i] << "\n");
	  // CFLog(INFO, "_diffResidua [" << i << "] = " << _diffResidual[i] << "\n");
	  (*tBackResidual)[i] -= _diffResidual[i];
	}
      }
      
      // compute and distribute the jacobian contributions
      for (CFuint jState = 0; jState < nbStatesInCell; ++jState)
      {
        CFLogDebugMax( "Perturbing jState = " << jState << "\n");

        if ((*states)[jState]->isParUpdatable()) {
          // jacobian contribution (dR_jState/dU_iState)_k
          // compute (R_jState[U_iState + dU_k] - R_jState[U_iState])/eps
          numericalJacob.computeDerivative(_tResidual[jState],
					   (*tBackResidual)[jState],
					   _tempRes);
	  
	  //    cout << "_tempRes = " << _tempRes << endl;
	  
          acc->addValues(jState, iState, iVar, &_tempRes[0]);
        }
      }
      
      // restore the unperturbed value
      numericalJacob.restore((*currState)[iVar]);
    }
  }

  getMethodData().switchToExplicit();
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD
