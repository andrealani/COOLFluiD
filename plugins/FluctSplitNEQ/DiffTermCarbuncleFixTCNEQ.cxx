#include "FluctSplit/DiffTermCarbuncleFixEuler.hh"
#include "Framework/MethodStrategyProvider.hh"
#include "FluctSplitNEQ/FluctSplitNEQ.hh"
#include "NavierStokes/MultiScalarVarSet.hh"
#include "NavierStokes/Euler2DVarSet.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Physics::NavierStokes;
using namespace COOLFluiD::FluctSplit;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace FluctSplitNEQ {

//////////////////////////////////////////////////////////////////////////////

MethodStrategyProvider<DiffTermCarbuncleFixEuler<MultiScalarVarSet<Euler2DVarSet> >,
		       FluctuationSplitData,
		       ComputeDiffusiveTerm,
		       FluctSplitNEQModule>
carbuncleFixTCNEQProvider("CarbuncleFixTCNEQ");

//////////////////////////////////////////////////////////////////////////////

  } // namespace FluctSplitNEQ

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
