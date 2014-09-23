#include "FluctSplit/ComputeJacobStrategy.hh"
#include "FluctSplit/ComputeDiffusiveTerm.hh"
#include "Framework/MeshData.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

ComputeJacobStrategy::ComputeJacobStrategy(const std::string& name) :
  MethodStrategy<FluctuationSplitData>(name),
  _fsStrategy(CFNULL),
  _solutionToDistMatTrans(CFNULL),
  _distToSolutionMatTrans(CFNULL),
  _linearToDistMatTrans(CFNULL),
  _solutionToLinearInUpdateMatTrans(CFNULL),
  _solutionToLinearMatTrans(CFNULL),
  _updateToLinearVecTrans(CFNULL),
  _updateToSolutionInUpdateMatTrans(CFNULL),
  _diffResidual(0),
  _hasDiffusiveTerm(false)
{
}

//////////////////////////////////////////////////////////////////////////////

ComputeJacobStrategy::~ComputeJacobStrategy()
{
}

//////////////////////////////////////////////////////////////////////////////

void ComputeJacobStrategy::setup()
{
  _fsStrategy = getMethodData().getFluctSplitStrategy();
//  _adStrategy = getMethodData().getArtificialDiffusionStrategy();

  _solutionToDistMatTrans = getMethodData().getSolutionToDistribMatTrans();
  _distToSolutionMatTrans = getMethodData().getDistribToSolutionMatTrans();
  _linearToDistMatTrans = getMethodData().getLinearToDistribMatTrans();
  _solutionToLinearInUpdateMatTrans =
    getMethodData().getSolutionToLinearInUpdateMatTrans();
  _solutionToLinearMatTrans = getMethodData().getSolutionToLinearMatTrans();
  _updateToLinearVecTrans = getMethodData().getUpdateToLinearVecTrans();

  _updateToSolutionInUpdateMatTrans = getMethodData().getUpdateToSolutionInUpdateMatTrans();

  const CFuint maxNbStatesInCell = MeshDataStack::getActive()->Statistics().getMaxNbStatesInCell();

  _diffResidual.resize(maxNbStatesInCell);
  // allocating data for the temporary local residual
  for (CFuint i = 0; i < maxNbStatesInCell; ++i) {
    _diffResidual[i].resize(PhysicalModelStack::getActive()->getNbEq());
    _diffResidual[i] = 0.0;
  }

  // flag telling if a diffusive term has to be computed
  _hasDiffusiveTerm = !getMethodData().getArtificialDiffusionStrategy()->isNull();
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD
