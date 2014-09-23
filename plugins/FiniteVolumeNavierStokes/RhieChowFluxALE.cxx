#include "Framework/MethodStrategyProvider.hh"
#include "FiniteVolumeNavierStokes/FiniteVolumeNavierStokes.hh"
#include "FiniteVolumeNavierStokes/RhieChowFlux.hh"
#include "FiniteVolumeNavierStokes/RhieChowFluxBlended.hh"
#include "NavierStokes/MultiScalarVarSet.hh"
#include "NavierStokes/Euler2DVarSet.hh"
#include "NavierStokes/Euler3DVarSet.hh"
#include "FiniteVolumeNavierStokes/RhieChowFluxALE.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Physics::NavierStokes;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

MethodStrategyProvider<RhieChowFluxALE<RhieChowFlux<Euler2DVarSet> >,
		       CellCenterFVMData,
		       FluxSplitter<CellCenterFVMData>,
		       FiniteVolumeNavierStokesModule>
rhieChow2DFluxALEProvider("RhieChow2DALE");

MethodStrategyProvider<RhieChowFluxALE<RhieChowFlux<Euler3DVarSet> >,
		       CellCenterFVMData,
		       FluxSplitter<CellCenterFVMData>,
		       FiniteVolumeNavierStokesModule>
rhieChowFluxALE3DProvider("RhieChow3DALE");
      
MethodStrategyProvider<RhieChowFluxALE<RhieChowFlux<MultiScalarVarSet<Euler2DVarSet> > >,
		       CellCenterFVMData,
		       FluxSplitter<CellCenterFVMData>,
		       FiniteVolumeNavierStokesModule>
rhieChowMS2DFluxALEProvider("RhieChowMS2DALE");

MethodStrategyProvider<RhieChowFluxALE<RhieChowFlux<MultiScalarVarSet<Euler3DVarSet> > >,
		       CellCenterFVMData,
		       FluxSplitter<CellCenterFVMData>,
		       FiniteVolumeNavierStokesModule>
rhieChowFluxALEMS3DProvider("RhieChowMS3DALE");

MethodStrategyProvider<RhieChowFluxALE<RhieChowFluxBlended<Euler2DVarSet> >,
		       CellCenterFVMData,
		       FluxSplitter<CellCenterFVMData>,
		       FiniteVolumeNavierStokesModule>
rhieChowBlendedALE2DFluxProvider("RhieChowBlended2DALE");

MethodStrategyProvider<RhieChowFluxALE<RhieChowFluxBlended<Euler3DVarSet> >,
		       CellCenterFVMData,
		       FluxSplitter<CellCenterFVMData>,
		       FiniteVolumeNavierStokesModule>
rhieChowBlendedALEFlux3DProvider("RhieChowBlended3DALE");
      
MethodStrategyProvider<RhieChowFluxALE<RhieChowFluxBlended<MultiScalarVarSet<Euler2DVarSet> > >,
		       CellCenterFVMData,
		       FluxSplitter<CellCenterFVMData>,
		       FiniteVolumeNavierStokesModule>
rhieChowBlendedALEMS2DFluxProvider("RhieChowBlendedMS2DALE");

MethodStrategyProvider<RhieChowFluxALE<RhieChowFluxBlended<MultiScalarVarSet<Euler3DVarSet> > >,
		       CellCenterFVMData,
		       FluxSplitter<CellCenterFVMData>,
		       FiniteVolumeNavierStokesModule>
rhieChowBlendedALEFluxMS3DProvider("RhieChowBlendedMS3DALE");

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
