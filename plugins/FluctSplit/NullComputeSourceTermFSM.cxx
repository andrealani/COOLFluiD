#include "NullComputeSourceTermFSM.hh"
#include "Common/CFLog.hh"
#include "Framework/State.hh"
#include "FluctSplit/FluctSplit.hh"
#include "Framework/MethodStrategyProvider.hh"
#include "FluctSplit/FluctuationSplitData.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

MethodStrategyProvider<NullComputeSourceTermFSM,
		       FluctuationSplitData,
		       ComputeSourceTermFSM,
		       FluctSplitModule>
nullSTProvider("Null");

//////////////////////////////////////////////////////////////////////////////

NullComputeSourceTermFSM::NullComputeSourceTermFSM(const std::string& name) :
  ComputeSourceTermFSM(name)
{
}

//////////////////////////////////////////////////////////////////////////////

NullComputeSourceTermFSM::~NullComputeSourceTermFSM()
{
}

//////////////////////////////////////////////////////////////////////////////

void NullComputeSourceTermFSM::setup()
{
  CFLog(VERBOSE, "Calling NullComputeSourceTermFSM::setup()\n");
}

//////////////////////////////////////////////////////////////////////////////

void NullComputeSourceTermFSM::computeSourceFSM(Framework::GeometricEntity *const cell,
			                  RealVector& source,
			                  const FluctSplit::InwardNormalsData& normalsData)
{
  CFLog(VERBOSE, "Calling NullComputeSourceTermFSM::computeSource()\n");
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
