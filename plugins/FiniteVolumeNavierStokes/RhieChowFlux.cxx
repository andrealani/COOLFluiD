#include "Framework/MethodStrategyProvider.hh"
#include "FiniteVolumeNavierStokes/FiniteVolumeNavierStokes.hh"
#include "FiniteVolumeNavierStokes/RhieChowFlux.hh"
#include "FiniteVolumeNavierStokes/RhieChowFluxBlended.hh"
#include "NavierStokes/Euler2DVarSet.hh"
#include "NavierStokes/Euler3DVarSet.hh"
#include "NavierStokes/MultiScalarVarSet.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Physics::NavierStokes;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

MethodStrategyProvider<RhieChowFlux<Euler2DVarSet>,
		       CellCenterFVMData,
		       FluxSplitter<CellCenterFVMData>,
		       FiniteVolumeNavierStokesModule>
rhieChow2DFluxProvider("RhieChow2D");

MethodStrategyProvider<RhieChowFlux<Euler3DVarSet>,
		       CellCenterFVMData,
		       FluxSplitter<CellCenterFVMData>,
		       FiniteVolumeNavierStokesModule>
rhieChowFlux3DProvider("RhieChow3D");

// 2D and 1/2
MethodStrategyProvider<RhieChowFlux<Euler3DVarSet>,
		       CellCenterFVMData,
		       FluxSplitter<CellCenterFVMData>,
		       FiniteVolumeNavierStokesModule>
rhieChowFlux2DHalfProvider("RhieChow2DHalf");
      
MethodStrategyProvider<RhieChowFlux<MultiScalarVarSet<Euler2DVarSet> >,
		       CellCenterFVMData,
		       FluxSplitter<CellCenterFVMData>,
		       FiniteVolumeNavierStokesModule>
rhieChowMS2DFluxProvider("RhieChowMS2D");

MethodStrategyProvider<RhieChowFlux<MultiScalarVarSet<Euler3DVarSet> >,
		       CellCenterFVMData,
		       FluxSplitter<CellCenterFVMData>,
		       FiniteVolumeNavierStokesModule>
rhieChowFluxMS3DProvider("RhieChowMS3D");

MethodStrategyProvider<RhieChowFluxBlended<Euler2DVarSet>,
		       CellCenterFVMData,
		       FluxSplitter<CellCenterFVMData>,
		       FiniteVolumeNavierStokesModule>
rhieChowBlended2DFluxProvider("RhieChowBlended2D");

MethodStrategyProvider<RhieChowFluxBlended<Euler3DVarSet>,
		       CellCenterFVMData,
		       FluxSplitter<CellCenterFVMData>,
		       FiniteVolumeNavierStokesModule>
rhieChowBlendedFlux3DProvider("RhieChowBlended3D");
      
MethodStrategyProvider<RhieChowFluxBlended<MultiScalarVarSet<Euler2DVarSet> >,
		       CellCenterFVMData,
		       FluxSplitter<CellCenterFVMData>,
		       FiniteVolumeNavierStokesModule>
rhieChowBlendedMS2DFluxProvider("RhieChowBlendedMS2D");

MethodStrategyProvider<RhieChowFluxBlended<MultiScalarVarSet<Euler3DVarSet> >,
		       CellCenterFVMData,
		       FluxSplitter<CellCenterFVMData>,
		       FiniteVolumeNavierStokesModule>
rhieChowBlendedFluxMS3DProvider("RhieChowBlendedMS3D");

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
