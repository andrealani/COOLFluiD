#include "NullComputeDiffusiveTerm.hh"
#include "Environment/ObjectProvider.hh"
#include "FluctSplit/FluctSplit.hh"
#include "Framework/MethodStrategyProvider.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

MethodStrategyProvider<NullComputeDiffusiveTerm,
		       FluctuationSplitData,
		       ComputeDiffusiveTerm,
		       FluctSplitModule>
nullDiffusiveTermProvider("Null");

//////////////////////////////////////////////////////////////////////////////

NullComputeDiffusiveTerm::NullComputeDiffusiveTerm(const std::string& name) :
  ComputeDiffusiveTerm(name)
{
}

//////////////////////////////////////////////////////////////////////////////

NullComputeDiffusiveTerm::~NullComputeDiffusiveTerm()
{
}

//////////////////////////////////////////////////////////////////////////////

void NullComputeDiffusiveTerm::setDiffusiveVarSet(SafePtr<DiffusiveVarSet> diffVar)
{
  CFLog(VERBOSE, "Calling NullComputeDiffusiveTerm::setDiffusiveVarSet()" << "\n");
}

//////////////////////////////////////////////////////////////////////////////

void NullComputeDiffusiveTerm::setUpdateVarSet(SafePtr<ConvectiveVarSet> updateVar)
{
  CFLog(VERBOSE, "Calling NullComputeDiffusiveTerm::setUpdateVarSet()" << "\n");
}

//////////////////////////////////////////////////////////////////////////////

void NullComputeDiffusiveTerm::computeDiffusiveTerm
(GeometricEntity *const geo, vector<RealVector>& result,
 const bool updateCoeff)
{
  CFLog(VERBOSE, "Calling NullComputeDiffusiveTerm::computeTerm()" <<
	"\n");
}

//////////////////////////////////////////////////////////////////////////////

void NullComputeDiffusiveTerm::setup()
{
  CFLog(VERBOSE, "Calling NullComputeDiffusiveTerm::setup()" << "\n");
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
