#include "FluctSplit/FluctSplit.hh"
#include "PseudoSteadyTimeAnalytMatRhs.hh"
#include "Framework/CFL.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/BlockAccumulator.hh"
#include "Framework/LSSMatrix.hh"
#include "MathTools/MatrixInverter.hh"
#include "Framework/MeshData.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::MathTools;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<PseudoSteadyTimeAnalytMatRhs, FluctuationSplitData, FluctSplitModule> fvmccPseudoSteadyTimeAnalytMatRhs("PseudoSteadyTimeAnalytMatRhs");

//////////////////////////////////////////////////////////////////////////////

PseudoSteadyTimeAnalytMatRhs::PseudoSteadyTimeAnalytMatRhs
(const std::string& name) :
  StdComputeTimeRhs(name),
  _inverter(CFNULL),
  socket_rhs("rhs"),
  socket_pastStates("pastStates"),
  _acc(CFNULL)
{
  _inverter = MatrixInverter::create(PhysicalModelStack::getActive()->getNbEq(), false);
}

//////////////////////////////////////////////////////////////////////////////

PseudoSteadyTimeAnalytMatRhs::~PseudoSteadyTimeAnalytMatRhs()
{
  deletePtr(_inverter);
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
PseudoSteadyTimeAnalytMatRhs::needsSockets()
{

  std::vector<Common::SafePtr<BaseDataSocketSink> > result = StdComputeTimeRhs::needsSockets();

  result.push_back(&socket_rhs);
  result.push_back(&socket_pastStates);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

void PseudoSteadyTimeAnalytMatRhs::setup()
{
  StdComputeTimeRhs::setup();

  _acc.reset(getMethodData().getLinearSystemSolver()[0]->
	     createBlockAccumulator(1, 1, PhysicalModelStack::getActive()->getNbEq()));
}

//////////////////////////////////////////////////////////////////////////////

void PseudoSteadyTimeAnalytMatRhs::execute()
{
  DataHandle< CFreal> rhs = socket_rhs.getDataHandle();

  DataHandle< CFreal> updateCoeff = socket_updateCoeff.getDataHandle();

  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();

  DataHandle<State*> pastStates = socket_pastStates.getDataHandle();

  SafePtr<LSSMatrix> jacobMatrix = _lss->getMatrix();

  SafePtr<VarSetTransformer> updateToSolutionVecTrans =
    getMethodData().getUpdateToSolutionVecTrans();

  SafePtr<VarSetMatrixTransformer> solutionToUpdateInUpdateMatTrans =
    getMethodData().getSolToUpdateInUpdateMatTrans();

  const CFuint nbStates = states.size();
  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  const CFreal cfl = getMethodData().getCFL()->getCFLValue();

  RealVector tempState(nbEqs);
  RealMatrix inverseMat(nbEqs, nbEqs);

  // add the diagonal entries in the jacobian (updateCoeff/CFL)
  for (CFuint iState = 0; iState < nbStates; ++iState) {
    const CFreal diagValue = updateCoeff[iState]/cfl;

    State *const currState = states[iState];

    if (currState->isParUpdatable()) {

      // this first transformed state HAS TO BE stored,
      // since the returned pointer will change pointee after
      // the second call to transform()
      tempState = static_cast<RealVector&>
	(*updateToSolutionVecTrans->transform(currState));

      const State *const tPastState =
	updateToSolutionVecTrans->transform(pastStates[iState]);

      // first the contribution to the rhs is computed
      for (CFuint iEq = 0; iEq < nbEqs; ++iEq) {
	const CFreal dU = tempState[iEq] - (*tPastState)[iEq];
	rhs(iState, iEq, nbEqs) -= dU*diagValue;
      }

      // now the contribution to the jacobian matrix is calculated
      _acc->setRowColIndex(0, currState->getLocalID());

      // set and get the transformation matrix in the update variables
      solutionToUpdateInUpdateMatTrans->setMatrix(*currState);
      const RealMatrix& matrix = *solutionToUpdateInUpdateMatTrans->getMatrix();

      // invert the matrix
      _inverter->invert(matrix, inverseMat);

      inverseMat *= diagValue;
      for (CFuint iVar = 0; iVar < nbEqs; ++iVar) {
	for (CFuint jVar = 0; jVar < nbEqs; ++jVar) {
	  _acc->addValue(0, 0, iVar, jVar,
			 inverseMat(iVar,jVar));
	}
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
