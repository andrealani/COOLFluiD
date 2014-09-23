#include "FluctSplit/FluctSplit.hh"
#include "PseudoSteadyTimeRhs.hh"
#include "Framework/CFL.hh"
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



    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<PseudoSteadyTimeRhs, FluctuationSplitData, FluctSplitModule> fvmccPseudoSteadyTimeRhs("PseudoSteadyTimeRhs");

//////////////////////////////////////////////////////////////////////////////

PseudoSteadyTimeRhs::PseudoSteadyTimeRhs
(const std::string& name) :
  StdComputeTimeRhs(name),
  _numericalJacob(CFNULL),
  socket_pastStates("pastStates"),
  socket_rhs("rhs"),
  _fluxDiff(),
  _tempState(),
  _tempPertState(),
  _acc(CFNULL)
{
}

//////////////////////////////////////////////////////////////////////////////

PseudoSteadyTimeRhs::~PseudoSteadyTimeRhs()
{
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
PseudoSteadyTimeRhs::needsSockets()
{

  std::vector<Common::SafePtr<BaseDataSocketSink> > result = StdComputeTimeRhs::needsSockets();

  result.push_back(&socket_rhs);
  result.push_back(&socket_pastStates);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

void PseudoSteadyTimeRhs::setup()
{
  StdComputeTimeRhs::setup();

  _fluxDiff.resize(PhysicalModelStack::getActive()->getNbEq());
  _tempState.resize(PhysicalModelStack::getActive()->getNbEq());
  _tempPertState.resize(PhysicalModelStack::getActive()->getNbEq());

  // numerical jacobian
  _numericalJacob = &getMethodData().getNumericalJacobian();

  _acc.reset(getMethodData().getLinearSystemSolver()[0]->
	     createBlockAccumulator(1, 1, PhysicalModelStack::getActive()->getNbEq()));
}

//////////////////////////////////////////////////////////////////////////////

void PseudoSteadyTimeRhs::execute()
{
  DataHandle<State*,Framework::GLOBAL> pastStates = socket_states.getDataHandle();
  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();
  DataHandle< CFreal> rhs = socket_rhs.getDataHandle();
  DataHandle< CFreal> updateCoeff = socket_updateCoeff.getDataHandle();
  DataHandle< vector<bool> > discardTimeJacob =
    socket_discardTimeJacob.getDataHandle();

  SafePtr<LSSMatrix> jacobMatrix = _lss->getMatrix();

  SafePtr<VarSetTransformer> updateToSolutionVecTrans =
    getMethodData().getUpdateToSolutionVecTrans();

  const CFuint nbStates = states.size();
  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  const CFreal cfl = getMethodData().getCFL()->getCFLValue();

  RealVector coeffAddRhs(0.0, nbEqs);

  // add the diagonal entries in the jacobian (updateCoeff/CFL)
  for (CFuint iState = 0; iState < nbStates; ++iState) {
    const CFreal diagValue = updateCoeff[iState]/cfl;

    State *const currState = states[iState];

    if (currState->isParUpdatable()) {

      // this first transformed state HAS TO BE stored,
      // since the returned pointer will change pointee after
      // the second call to transform()
      _tempState = static_cast<RealVector&>
	(*updateToSolutionVecTrans->transform(currState));

      const State *const tPastState =
	updateToSolutionVecTrans->transform(pastStates[iState]);

      // first the contribution to the rhs is computed
      for (CFuint iEq = 0; iEq < nbEqs; ++iEq) {
	const CFreal dU = _tempState[iEq] - (*tPastState)[iEq];

	coeffAddRhs[iEq] = 1.0;
	if (discardTimeJacob[iState].size() > 0) {
	  cf_assert(discardTimeJacob[iState].size() == nbEqs);
	  if(discardTimeJacob[iState][iEq]) {
	    coeffAddRhs[iEq] = 0.0;
	  }
	}

	rhs(iState, iEq, nbEqs) -= coeffAddRhs[iEq]*dU*diagValue;
      }

      // now the contribution to the jacobian matrix is calculated
      _acc->setRowColIndex(0, currState->getLocalID());
      for (CFuint iVar = 0; iVar < nbEqs; ++iVar) {
	// perturb the given component of the state vector
	_numericalJacob->perturb(iVar, (*currState)[iVar]);

	// this first transformed state HAS TO BE stored,
	// since the returned pointer will chenge pointee after
	// the second call to transform()
	_tempPertState = static_cast<RealVector&>
	  (*updateToSolutionVecTrans->transform(currState));

	// compute the finite difference derivative of the flux
	_numericalJacob->computeDerivative(_tempState,
					   _tempPertState,
					   _fluxDiff);

	// remove contribution from components whose time jacobian
	// has to be discarded
	_fluxDiff *= diagValue*coeffAddRhs;

	_acc->addValues(0, 0, iVar, &_fluxDiff[0]);
	// restore the unperturbed value
	_numericalJacob->restore((*currState)[iVar]);
      }

      // add the values in the jacobian matrix
      jacobMatrix->addValues(*_acc);

      // reset to zero the entries in the block accumulator
      _acc->reset();
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD
