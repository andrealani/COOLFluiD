#include "FluctSplit/PicardJacobStrategy.hh"
#include "FluctSplit/ComputeDiffusiveTerm.hh"

#include "Framework/MeshData.hh"
#include "Framework/MethodStrategyProvider.hh"
#include "Framework/BlockAccumulator.hh"
#include "FluctSplit/FluctSplit.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

MethodStrategyProvider<PicardJacobStrategy,
                       FluctuationSplitData,
                       ComputeJacobStrategy,
                       FluctSplitModule>
                       picardJacobStrategyProvider("Picard");

//////////////////////////////////////////////////////////////////////////////

PicardJacobStrategy::PicardJacobStrategy(const std::string& name) :
  ComputeJacobStrategy(name),
  _splitter(CFNULL),
  _sysSplitter(CFNULL),
  _scalarSplitter(CFNULL),
  _jacob(),
  _temp0(),
  _temp1(),
  _temp2(),
  _jacobBlocks()
{
}

//////////////////////////////////////////////////////////////////////////////

PicardJacobStrategy::~PicardJacobStrategy()
{
  for (CFuint i = 0; i < _jacobBlocks.size(); ++i) {
    deletePtr(_jacobBlocks[i]);
  }
}

//////////////////////////////////////////////////////////////////////////////

void PicardJacobStrategy::setup()
{
  // first call parent method
  ComputeJacobStrategy::setup();

  _jacob.resize(PhysicalModelStack::getActive()->getNbEq(),PhysicalModelStack::getActive()->getNbEq());
  _temp0.resize(PhysicalModelStack::getActive()->getNbEq(),PhysicalModelStack::getActive()->getNbEq());
  _temp1.resize(PhysicalModelStack::getActive()->getNbEq(),PhysicalModelStack::getActive()->getNbEq());
  _temp2.resize(PhysicalModelStack::getActive()->getNbEq(),PhysicalModelStack::getActive()->getNbEq());

  // get the splitter
  if (getMethodData().isMultipleSplitter()) {
    // it is assumed that splitters have been already set up before
    // @see concrete FluctuationSplitStrategy's
    _scalarSplitter = getMethodData().getScalarSplitter();
    _sysSplitter = getMethodData().getSysSplitter();
  }
  else {
    // it is assumed that splitters have been already set up before
    // @see concrete FluctuationSplitStrategy's
    _splitter = getMethodData().getSplitter();
  }

  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  const CFuint maxNbStatesInCell = MeshDataStack::getActive()->Statistics().getMaxNbStatesInCell();
  _jacobBlocks.resize(maxNbStatesInCell*maxNbStatesInCell);
  _diffjacob.resize(maxNbStatesInCell*maxNbStatesInCell);

  for (CFuint i = 0; i < _jacobBlocks.size(); ++i) {
    _jacobBlocks[i] = new RealMatrix(nbEqs, nbEqs);
    _diffjacob[i] = new RealMatrix(nbEqs, nbEqs);

  }
}

//////////////////////////////////////////////////////////////////////////////

void PicardJacobStrategy::computeJacobianTerm(GeometricEntity *const cell,
					      const vector<RealVector>& residual,
					      BlockAccumulator *const acc,
					      const std::vector<CFuint>&  equationIDs)
{
  SafePtr<ComputeDiffusiveTerm> diffTermComputer =
    getMethodData().getDiffusiveTermComputer();

  if (getMethodData().isMultipleSplitter()) {
    // compute part of the jacobian contribution
    _sysSplitter->computePicardJacobPart(_jacobBlocks);
    _scalarSplitter->computePicardJacobPart(_jacobBlocks);
  }
  else {
    // compute the jacobian contribution
    _splitter->computePicardJacob(_jacobBlocks);
  }

  //  if (_hasDiffusiveTerm)
  //    {
  //	diffTermComputer->computePicardDiffJacob(cell,_diffjacob);
  //   }
  // _jacobBlocks = _jacobBlocks + _diffjacob;

  vector<State*> *const states = cell->getStates();

  // insert the blocks of values in the acccumulator
  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  const CFuint nbStatesInCell = states->size();

  // matrix of the transformation distribution-solution variables
  const RealMatrix *const dUdW = _distToSolutionMatTrans->getMatrix();

  // matrix of the transformation linearization-solution variables
  const RealMatrix *const dWdZ = _linearToDistMatTrans->getMatrix();

  const bool isResTransformationNeeded = getMethodData().isResidualTransformationNeeded();

  for (CFuint iState = 0; iState < nbStatesInCell; ++iState)
  {
    const CFuint nStart = iState*nbStatesInCell;
    acc->setRowColIndex(iState, (*states)[iState]->getLocalID());

    const State& currState = *(*states)[iState];

    // set the transformation matrix from linear to solution with
    // solution variables as input
    _solutionToLinearInUpdateMatTrans->setMatrix(currState);
    const RealMatrix *const dZdU =
      _solutionToLinearInUpdateMatTrans->getMatrix();

    for (CFuint jState = 0; jState < nbStatesInCell; ++jState) {
      if ((*states)[jState]->isParUpdatable()) {
	const RealMatrix *const block = _jacobBlocks[nStart + jState];

	if (!isResTransformationNeeded) {
	  _temp1 = (*dWdZ)*(*dZdU);
	}
	else {
	  _updateToSolutionInUpdateMatTrans->setMatrix(currState);
	  const RealMatrix& dUdP =
	    *_updateToSolutionInUpdateMatTrans->getMatrix();


	  _temp0 = (*dZdU)*dUdP;
	  _temp1 = (*dWdZ)*_temp0;
	}

	_temp2 = (*block)*_temp1;
	_jacob = (*dUdW)*_temp2;

	for (CFuint iVar = 0; iVar < nbEqs; ++iVar) {
	  for (CFuint jVar = 0; jVar < nbEqs; ++jVar) {
	    acc->addValue(iState, jState, iVar, jVar,
			  _jacob(iVar,jVar));
	  }
	}
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD
