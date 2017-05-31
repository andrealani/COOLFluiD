#include "FluctSplit/NullJacobStrategy.hh"

#include "Framework/MeshData.hh"
#include "Framework/MethodStrategyProvider.hh"
#include "FluctSplit/FluctSplit.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

MethodStrategyProvider<NullJacobStrategy,
                       FluctuationSplitData,
                       ComputeJacobStrategy,
                       FluctSplitModule>
                       nullJacobStrategyProvider("Null");

//////////////////////////////////////////////////////////////////////////////

NullJacobStrategy::NullJacobStrategy(const std::string& name) :
  ComputeJacobStrategy(name)
{
}

//////////////////////////////////////////////////////////////////////////////

NullJacobStrategy::~NullJacobStrategy()
{
}

//////////////////////////////////////////////////////////////////////////////

void NullJacobStrategy::setup()
{
  ComputeJacobStrategy::setup();
}

//////////////////////////////////////////////////////////////////////////////

void NullJacobStrategy::computeJacobianTerm (GeometricEntity *const cell,
					     const vector<RealVector>& residual,
					     BlockAccumulator *const acc,
					     const std::vector<CFuint>&  equationIDs)
{
  CFLog(VERBOSE, "Calling NullJacobStrategy::computeJacobianTerm()" << "\n");
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD
