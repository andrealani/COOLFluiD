#include "NavierStokes2DSourceTerm.hh"
#include "FluctSplit/FluctSplitNavierStokes.hh"
#include "NavierStokes/Euler2DVarSet.hh"
#include "NavierStokes/MultiScalarVarSet.hh"
#include "Framework/MethodStrategyProvider.hh"
#include "FluctSplit/FluctSplitNavierStokes.hh"
#include "FluctSplit/FluctuationSplitData.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Physics::NavierStokes;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

MethodStrategyProvider<NavierStokes2DSourceTerm<Euler2DVarSet>,
		       FluctuationSplitData,
		       ComputeSourceTermFSM,
		       FluctSplitNavierStokesModule>
ns2DSTProvider("NavierStokes2DAxiST");

MethodStrategyProvider<NavierStokes2DSourceTerm<MultiScalarVarSet<Euler2DVarSet> >,
		       FluctuationSplitData,
		       ComputeSourceTermFSM,
		       FluctSplitNavierStokesModule>
ns2DMultiScalarSTProvider("NavierStokes2DMSAxiST");

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
