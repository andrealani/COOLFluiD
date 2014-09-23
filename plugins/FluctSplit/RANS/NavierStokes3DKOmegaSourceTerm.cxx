#include "NavierStokes3DKOmegaSourceTerm.hh"
#include "FluctSplit/FluctSplitNavierStokes.hh"
#include "NavierStokes/Euler3DVarSet.hh"
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

MethodStrategyProvider<NavierStokes3DKOmegaSourceTerm<NavierStokes3DKOmega>,
		       FluctuationSplitData,
		       ComputeSourceTermFSM,
		       FluctSplitNavierStokesModule>
ns3DKOProvider("NavierStokes3DKOmegaSourceTerm");

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit

//   } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
