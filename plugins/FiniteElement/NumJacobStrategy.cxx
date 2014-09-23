#include "FiniteElement/FiniteElement.hh"
#include "NumJacobStrategy.hh"
#include "Framework/MeshData.hh"
#include "Framework/MethodStrategyProvider.hh"
#include "Framework/BlockAccumulator.hh"
#include "Framework/NumericalJacobian.hh"
#include "ComputeResidualStrategy.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteElement {

//////////////////////////////////////////////////////////////////////////////

MethodStrategyProvider<NumJacobStrategy, FiniteElementMethodData,ComputeJacobStrategy, FiniteElementModule> numJacobStrategyProvider("Numerical");

//////////////////////////////////////////////////////////////////////////////

NumJacobStrategy::NumJacobStrategy(const std::string& name) :
  ComputeJacobStrategy(name),
  _tempRes(),
  _otherResidual(0)
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
}

//////////////////////////////////////////////////////////////////////////////

void NumJacobStrategy::computeJacobianTerm()
{
  FiniteElementMethodData& femdata  = getMethodData();
  LocalElementData& local_elem_data = femdata.getLocalElementData();
  
  BlockAccumulator& acc = *local_elem_data.blockacc;
//   RealMatrix& elemMat   = *local_elem_data.stiff_mat;
//   RealVector& elemVec   = *local_elem_data.load_vec;
  const CFuint nbEqs    =  local_elem_data.nbEqs;
  const CFuint nbStates =  local_elem_data.nbStates;
  GeometricEntity& cell = *local_elem_data.cell;
  std::vector<RealVector>& residual = *local_elem_data.residual;

  NumericalJacobian& numericalJacob = femdata.getNumericalJacobian();
  SafePtr<ComputeResidualStrategy> feStrategy = femdata.getResidualStrategy();

  vector<State*>& states = *cell.getStates();

  // compute the perturbed states for the evaluation of the
  // jacobian matrix looping over the state vectors in this cell
  for (CFuint iState = 0; iState < nbStates; ++iState) {

    State& currState = *states[iState];

    CFLogDebugMax( "Perturbing iState = " << iState << " with stateID = " << currState.getLocalID() << "\n");

    // set the index of the block corresponding to the current
    // state in the jacobian matrix
    acc.setRowColIndex(iState, currState.getLocalID());
    // loop over the variables in the state vector to perturb one
    // component at a time
    for (CFuint iVar = 0; iVar < nbEqs; ++iVar) {
      CFLogDebugMax( "Perturbing iVar = " << iVar << "\n");

      // reset the residual to 0.
      cleanOtherResidual();

      // perturb the given component of the state vector
      numericalJacob.perturb(iVar, currState[iVar]);

      // compute element matrix and vector and put in _otherResidual the perturbed residual
      feStrategy->computeElementResidual(_otherResidual);

      // compute and distribute the jacobian contributions
      for (CFuint jState = 0; jState < nbStates; ++jState) {
        CFLogDebugMax( "Perturbing jState = " << jState << "\n");

        if (states[jState]->isParUpdatable()) {

          // jacobian contribution (dR_jState/dU_iState)_k
          // compute (R_jState[U_iState + dU_k] - R_jState[U_iState])/eps
          numericalJacob.computeDerivative(residual[jState],
                                           _otherResidual[jState],
                                           _tempRes);
          acc.addValues(jState, iState, iVar, &_tempRes[0]);
        }
      }
      // restore the unperturbed value
      numericalJacob.restore(currState[iVar]);
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteElement

  } // namespace Numerics

} // namespace COOLFluiD
