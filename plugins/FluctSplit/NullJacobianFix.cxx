#include "Framework/MethodStrategyProvider.hh"
#include "FluctSplit/NullJacobianFix.hh"
#include "FluctSplit/FluctuationSplitData.hh"
#include "FluctSplit/FluctSplit.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

MethodStrategyProvider<NullJacobianFix,
                       FluctuationSplitData,
                       ComputeJacobianFix,
                       FluctSplitModule>
nullJacobianFixProvider("Null");

//////////////////////////////////////////////////////////////////////////////

NullJacobianFix::NullJacobianFix(const std::string& name) :
  ComputeJacobianFix(name)
{
}

//////////////////////////////////////////////////////////////////////////////

NullJacobianFix::~NullJacobianFix()
{
}

//////////////////////////////////////////////////////////////////////////////

void NullJacobianFix::setup()
{
  ComputeJacobianFix::setup();
}

//////////////////////////////////////////////////////////////////////////////

void NullJacobianFix::computeFix(const InwardNormalsData& normalsData,
				 RealVector& delta)
{
  CFLog(DEBUG_MIN,"Calling NullJacobianFix::computeFix() : this is a Null JacobianFix.\n");
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD
