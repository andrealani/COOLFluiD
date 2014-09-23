#include "FiniteVolumeTurb/NavierStokesKOmega2DSourceTerm.hh"
#include "FiniteVolumeCombustion/FiniteVolumeCombustion.hh"
#include "Framework/MethodStrategyProvider.hh"
#include "FiniteVolume/CellCenterFVMData.hh"
#include "NEQKOmega/NavierStokesNEQKOmegaRhoivt.hh"
#include "NEQ/NavierStokesCNEQVarSet.hh"
#include "KOmega/NavierStokesKOmegaVarSet.hh"
#include "NavierStokes/NavierStokes2DVarSet.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Physics::NavierStokes;
using namespace COOLFluiD::Physics::NEQKOmega;
using namespace COOLFluiD::Physics::NEQ;
using namespace COOLFluiD::Physics::KOmega;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

MethodStrategyProvider<NavierStokesKOmega2DSourceTerm
		       <NavierStokesNEQKOmegaRhoivt
			<NavierStokesKOmegaVarSet
			 <NavierStokesCNEQVarSet
			  <NavierStokes2DVarSet>, 2> > >,
		       CellCenterFVMData,
		       ComputeSourceTerm<CellCenterFVMData>,
		       FiniteVolumeKOmegaModule>
navierStokesNEQKOmega2DSourceTermFVMCCProvider("NavierStokesNEQKOmega2DSourceTerm");
      
//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
