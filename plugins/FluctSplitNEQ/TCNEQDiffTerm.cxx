#include "TCNEQDiffTerm.hh"
#include "Framework/MethodStrategyProvider.hh"
#include "NavierStokes/Euler2DVarSet.hh"
#include "NavierStokes/Euler3DVarSet.hh"
#include "NavierStokes/MultiScalarVarSet.hh"
#include "FluctSplitNEQ/FluctSplitNEQ.hh"
#include "FluctSplit/DiffTermCarbuncleFixEuler.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Physics::NavierStokes;
using namespace COOLFluiD::FluctSplit;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace FluctSplitNEQ {

//////////////////////////////////////////////////////////////////////////////

MethodStrategyProvider<TCNEQDiffTerm<ComputeDiffusiveTerm, MultiScalarVarSet<Euler2DVarSet> >,
		       FluctuationSplitData,
		       ComputeDiffusiveTerm,
		       FluctSplitNEQModule>
tcneqDiffTerm2DProvider("TCNEQ2D");

MethodStrategyProvider<TCNEQDiffTerm<ComputeDiffusiveTerm, MultiScalarVarSet<Euler3DVarSet> >,
		       FluctuationSplitData,
		       ComputeDiffusiveTerm,
		       FluctSplitNEQModule>
tcneqDiffTerm3DProvider("TCNEQ3D");

MethodStrategyProvider<TCNEQDiffTerm<DiffTermCarbuncleFixEuler<MultiScalarVarSet<Euler2DVarSet> >, MultiScalarVarSet<Euler2DVarSet> >,
		       FluctuationSplitData,
		       ComputeDiffusiveTerm,
		       FluctSplitNEQModule>
carbuncleFixTcneqDiffTerm2DProvider("CarbuncleFixTCNEQ2D");

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplitNEQ

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
