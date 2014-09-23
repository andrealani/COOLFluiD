#include "TCNEQSourceTerm.hh"
#include "TCNEQAxiSourceTerm.hh"
#include "Framework/MethodStrategyProvider.hh"
#include "FluctSplitNEQ/FluctSplitNEQ.hh"
#include "NavierStokes/Euler2DVarSet.hh"
#include "NavierStokes/MultiScalarVarSet.hh"
#include "FluctSplit/FluctuationSplitData.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Physics::NavierStokes;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



    namespace FluctSplitNEQ {

//////////////////////////////////////////////////////////////////////////////

MethodStrategyProvider<TCNEQSourceTerm<MultiScalarVarSet<Euler2DVarSet> >,
		       FluctSplit::FluctuationSplitData,
		       FluctSplit::ComputeSourceTermFSM,
		       FluctSplitNEQModule>
tcneqEuler2DSTProvider("Euler2DNEQST");

MethodStrategyProvider<TCNEQAxiSourceTerm<MultiScalarVarSet<Euler2DVarSet> >,
		       FluctSplit::FluctuationSplitData,
		       FluctSplit::ComputeSourceTermFSM,
		       FluctSplitNEQModule>
tcneqAxiEuler2DSTProvider("Euler2DNEQAxiST");

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplitNEQ



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
