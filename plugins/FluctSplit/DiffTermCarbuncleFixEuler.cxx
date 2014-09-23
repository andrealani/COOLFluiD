#include "FluctSplit/DiffTermCarbuncleFixEuler.hh"
//#include "FluctSplit/DiffTermCarbuncleFixEuler2.hh"
#include "Framework/MethodStrategyProvider.hh"
#include "NavierStokes/EulerVarSet.hh"
#include "FluctSplit/FluctSplitNavierStokes.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Physics::NavierStokes;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

MethodStrategyProvider<DiffTermCarbuncleFixEuler<EulerVarSet>,
		       FluctuationSplitData,
		       ComputeDiffusiveTerm,
		       FluctSplitNavierStokesModule>
carbuncleFixEulerProvider("CarbuncleFixEuler");

// MethodStrategyProvider<DiffTermCarbuncleFixEuler<EulerVarSet>,
//            FluctuationSplitData,
//            ComputeDiffusiveTerm,
//            FluctSplitNavierStokesModule>
// carbuncleFixEuler2Provider("CarbuncleFixEuler2");

//MethodStrategyProvider<DiffTermCarbuncleFixEuler2<EulerVarSet>,
//           FluctuationSplitData,
//            ComputeDiffusiveTerm,
//            FluctSplitNavierStokesModule>
// carbuncleFixEuler2Provider("CarbuncleFixEuler2");

//////////////////////////////////////////////////////////////////////////////

  } // namespace FluctSplit

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
