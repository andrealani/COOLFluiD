#include "NavierStokes2DKOmegaSourceTerm.hh"
#include "FluctSplit/FluctSplitNavierStokes.hh"
#include "NavierStokes/Euler2DVarSet.hh"
#include "NavierStokes/MultiScalarVarSet.hh"
#include "Framework/MethodStrategyProvider.hh"
#include "FluctSplit/FluctuationSplitData.hh"
#include "FluctSplitKOmega.hh"
#include "KOmega/NavierStokesKOmegaVarSetTypes.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Physics::NavierStokes;
using namespace COOLFluiD::Physics::KOmega;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

//   namespace Numerics {

    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

MethodStrategyProvider<NavierStokes2DKOmegaSourceTerm<NavierStokes2DKOmega>,
		       FluctuationSplitData,
		       ComputeSourceTermFSM,
		       FluctSplitNavierStokesModule>
ns2DKOProvider("NavierStokes2DKOmegaSourceTerm");

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit

//   } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
