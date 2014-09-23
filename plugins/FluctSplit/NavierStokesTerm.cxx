#include "Environment/ObjectProvider.hh"
#include "Framework/MethodStrategyProvider.hh"
#include "FluctSplit/FluctSplitNavierStokes.hh"
#include "FluctSplit/DiffTermCarbuncleFixEuler.hh"
#include "FluctSplit/NavierStokesTerm.hh"
#include "NavierStokes/EulerVarSet.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Physics::NavierStokes;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

MethodStrategyProvider<NavierStokesTerm<ComputeDiffusiveTerm>,
		       FluctuationSplitData,
		       ComputeDiffusiveTerm,
		       FluctSplitNavierStokesModule>
navierStokesDiffusiveTermProvider("NavierStokes");

//////////////////////////////////////////////////////////////////////////////

MethodStrategyProvider<NavierStokesTerm<DiffTermCarbuncleFixEuler<EulerVarSet> >,
		       FluctuationSplitData,
		       ComputeDiffusiveTerm,
		       FluctSplitNavierStokesModule>
carbuncleFixNavierStokesProvider("CarbuncleFixNavierStokes");

//////////////////////////////////////////////////////////////////////////////

} // namespace FluctSplit

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
