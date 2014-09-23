#include "FluctSplit/TwoLayerNumJacobStrategy.hh"

#include "Framework/MeshData.hh"
#include "Framework/MethodStrategyProvider.hh"
#include "Framework/BlockAccumulator.hh"
#include "FluctSplit/FluctSplitSpaceTime.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

MethodStrategyProvider<TwoLayerNumJacobStrategy,
                       FluctuationSplitData,
                       ComputeJacobStrategy,
                       FluctSplitSpaceTimeModule>
                       twoLayerNumJacobStrategyProvider("Numerical2Layer");

//////////////////////////////////////////////////////////////////////////////

TwoLayerNumJacobStrategy::TwoLayerNumJacobStrategy(const std::string& name) :
  ComputeJacobStrategy(name),
  socket_interStates("interStates"),
  _tempRes(),
  _otherResidual(0),
  _tResidual(0)
{
}

//////////////////////////////////////////////////////////////////////////////

TwoLayerNumJacobStrategy::~TwoLayerNumJacobStrategy()
{
}

//////////////////////////////////////////////////////////////////////////////

void TwoLayerNumJacobStrategy::setup()
{
  CFAUTOTRACE;

  cf_assert(socket_interStates.isConnected());
  // first call parent method
  ComputeJacobStrategy::setup();

  _tempRes.resize(PhysicalModelStack::getActive()->getNbEq()*2);

  _otherResidual.resize(MeshDataStack::getActive()->Statistics().getMaxNbStatesInCell());
  for (CFuint i = 0; i < _otherResidual.size(); ++i) {
    _otherResidual[i].resize(PhysicalModelStack::getActive()->getNbEq()*2);
    _otherResidual[i] = 0.;
  }

  _tResidual.resize(MeshDataStack::getActive()->Statistics().getMaxNbStatesInCell());
  for (CFuint i = 0; i < _tResidual.size(); ++i) {
    _tResidual[i].resize(PhysicalModelStack::getActive()->getNbEq()*2);
    _tResidual[i] = 0.;
  }
}

//////////////////////////////////////////////////////////////////////////////

void TwoLayerNumJacobStrategy::computeJacobianTerm(GeometricEntity *const cell,
						   const vector<RealVector>& residual,
						   BlockAccumulator *const acc,
						   const std::vector<CFuint>&  equationIDs)
{
  NumericalJacobian& numericalJacob = getMethodData().getNumericalJacobian();

  vector<State*> *const states = cell->getStates();

  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  const CFuint nbStatesInCell = states->size();

  DataHandle<State*> interStates = socket_interStates.getDataHandle();

  SafePtr<ConvectiveVarSet> distribVar =
    getMethodData().getDistribVar();

  // here we do the same transformation a second time
  // (try to avoid this !!!!!) after ComputeRhsJacobImpl
  // transform back the residual to solution variables
  _tResidual = *_distToSolutionMatTrans->transformMultiFromRef
    (const_cast<vector<RealVector>*>(&residual),2);

  // compute the perturbed states for the evaluation of the
  // jacobian matrix looping over the state vectors in this cell
  for (CFuint iState = 0; iState < nbStatesInCell; ++iState) {
    State *const currState = (*states)[iState];
    const CFuint stateID = currState->getLocalID();
    State *const intState = interStates[stateID];

    CFLogDebugMin( "Perturbing iState = " << iState <<
          " with stateID = " << stateID << "\n");

    // set the index of the block corresponding to the current
    // state in the jacobian matrix
    ///@todo ??????????? CHANGE THIS !!!!!!!!!!!!!!!!!
    acc->setRowColIndex(iState, currState->getLocalID());
    //CFout << "GlobalID: " << currState->getGlobalID() << "\n";
    //CFout << "LocalID: " << currState->getLocalID() << "\n";
    //CFout << "iState: " << iState << "\n";
    
    // loop over the variables in the state vector to perturb one
    // component at a time
    for (CFuint iVar = 0; iVar < nbEqs*2; ++iVar) {
      CFLogDebugMax( "Perturbing iVar = " << iVar << "\n");

      // reset the residual to 0.
      cleanOtherResidual();

      // perturb the given component of the state vector
      if (iVar < nbEqs) numericalJacob.perturb(iVar, (*intState)[iVar]);
      else numericalJacob.perturb((iVar-nbEqs), (*currState)[iVar-nbEqs]);

      _fsStrategy->computeFluctuation(_otherResidual);

      const vector<RealVector> *const tBackResidual =
        _distToSolutionMatTrans->transformMultiFromRef(&_otherResidual,2);

      // compute and distribute the jacobian contributions
      for (CFuint jState = 0; jState < nbStatesInCell; ++jState) {
        CFLogDebugMax( "Perturbing jState = " << jState << "\n");

        if ((*states)[jState]->isParUpdatable()) {
          // jacobian contribution (dR_jState/dU_iState)_k
          // compute (R_jState[U_iState + dU_k] - R_jState[U_iState])/eps
          numericalJacob.computeDerivative(_tResidual[jState],
                                            (*tBackResidual)[jState],
                                            _tempRes);
          acc->addValues(jState, iState, iVar, &_tempRes[0]);

//CFout << "(" << jState << "," << iState <<"," << iVar << ")" << "\n";
//CFout << _tResidual[jState] << "\n";
//CFout << (*tBackResidual)[jState] << "\n";
        }
      }

      // restore the unperturbed value
      if (iVar < nbEqs) numericalJacob.restore((*intState)[iVar]);
      else numericalJacob.restore((*currState)[iVar-nbEqs]);
    }
  }

}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD
