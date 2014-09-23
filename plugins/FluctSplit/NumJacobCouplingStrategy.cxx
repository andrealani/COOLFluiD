#include "FluctSplit/NumJacobCouplingStrategy.hh"
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

MethodStrategyProvider<NumJacobCouplingStrategy,
                       FluctuationSplitData,
                       ComputeJacobStrategy,
                       FluctSplitModule>
                       numJacobCouoplingStrategyProvider("NumericalCoupling");

//////////////////////////////////////////////////////////////////////////////

NumJacobCouplingStrategy::NumJacobCouplingStrategy(const std::string& name) :
  ComputeJacobStrategy(name),
  _tempRes(),
  _otherResidual(0),
  _tResidual(0)
{
}

//////////////////////////////////////////////////////////////////////////////

NumJacobCouplingStrategy::~NumJacobCouplingStrategy()
{
}

//////////////////////////////////////////////////////////////////////////////

void NumJacobCouplingStrategy::setup()
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

void NumJacobCouplingStrategy::computeJacobianTerm
 (GeometricEntity *const cell,
  const vector<RealVector>& residual,
  BlockAccumulator *const acc,
  const std::vector<CFuint>&  equationIDs)
{
  getMethodData().switchToImplicit();
  
  vector<State*> *const states = cell->getStates();
  
  NumericalJacobian& numericalJacob = getMethodData().getNumericalJacobian();
  
  const CFuint nbEqs = equationIDs.size();
  const CFuint start = equationIDs[0];
  const CFuint nbStatesInCell = states->size();
  
  SafePtr<ConvectiveVarSet> distribVar = getMethodData().getDistribVar();
  SafePtr<ComputeDiffusiveTerm> diffTermComputer = getMethodData().getDiffusiveTermComputer();
  
  // here we store the residual (referred to solution variables)
  const CFuint end = start + nbEqs;
  for (CFuint i = 0; i < _tResidual.size(); ++i) {
    for (CFuint j = start; j < end; ++j) {
      _tResidual[i][j] = residual[i][j];
    }
  }
  
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
      
      // think about a way of avoiding to change all schemes to split between the two systems  
      _fsStrategy->computeFluctuation(_otherResidual);
      
      // if you distribute in conservative variables this is Identity
      vector<RealVector> *const tBackResidual =
        _distToSolutionMatTrans->transformFromRef(&_otherResidual); 

      if (_hasDiffusiveTerm)
      {
	// modify the diffusive term to split between the subsystems
	diffTermComputer->computeDiffusiveTerm(cell, _diffResidual, false);

	for (CFuint i = 0; i < nbStatesInCell; ++i)
	{
	  (*tBackResidual)[i].slice(start, nbEqs) -= _diffResidual[i].slice(start, nbEqs);
	}
      }

      // compute and distribute the jacobian contributions
      for (CFuint jState = 0; jState < nbStatesInCell; ++jState)
      {

        CFLogDebugMax( "Perturbing jState = " << jState << "\n");

        if ((*states)[jState]->isParUpdatable()) {
          // jacobian contribution (dR_jState/dU_iState)_k
          // compute (R_jState[U_iState + dU_k] - R_jState[U_iState])/eps
          numericalJacob.computeDerivative(_tResidual[jState].slice(start, nbEqs),
					   (*tBackResidual)[jState].slice(start, nbEqs),
					   _tempRes.slice(start, nbEqs));
	  
	  
	  //    cout << "_tempRes = " << _tempRes << endl;
	  
          acc->addValues(jState, iState, iVar, &_tempRes[start]);
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
