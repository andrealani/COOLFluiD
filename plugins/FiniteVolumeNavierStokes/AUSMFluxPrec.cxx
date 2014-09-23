#include "FiniteVolumeNavierStokes/AUSMPlusFluxPrec.hh"
#include "NavierStokes/MultiScalarVarSet.hh"
#include "NavierStokes/Euler1DVarSet.hh"
#include "NavierStokes/Euler2DVarSet.hh"
#include "NavierStokes/Euler3DVarSet.hh"
#include "Framework/MethodStrategyProvider.hh"
#include "FiniteVolumeNavierStokes/FiniteVolumeNavierStokes.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Physics::NavierStokes;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

/// AUSM+P Comp
MethodStrategyProvider<AUSMPlusFluxPrec<Euler1DVarSet>,
		       CellCenterFVMData,
		       FluxSplitter<CellCenterFVMData>,
		       FiniteVolumeNavierStokesModule>
ausmPlusPrec1dProvider("AUSMPlusP1D");

MethodStrategyProvider<AUSMPlusFluxPrec<Euler2DVarSet>,
		       CellCenterFVMData,
		       FluxSplitter<CellCenterFVMData>,
		       FiniteVolumeNavierStokesModule>
ausmPlusPrec2dProvider("AUSMPlusP2D");

MethodStrategyProvider<AUSMPlusFluxPrec<Euler3DVarSet>,
		       CellCenterFVMData,
		       FluxSplitter<CellCenterFVMData>,
		       FiniteVolumeNavierStokesModule>
ausmPlusPrec3dProvider("AUSMPlusP3D");
      
MethodStrategyProvider<AUSMPlusFluxPrec<MultiScalarVarSet<Euler1DVarSet> >,
		       CellCenterFVMData,
		       FluxSplitter<CellCenterFVMData>,
		       FiniteVolumeNavierStokesModule>
ausmPlusPrecMS1dProvider("AUSMPlusPMS1D");

MethodStrategyProvider<AUSMPlusFluxPrec<MultiScalarVarSet<Euler2DVarSet> >,
		       CellCenterFVMData,
		       FluxSplitter<CellCenterFVMData>,
		       FiniteVolumeNavierStokesModule>
ausmPlusPrecMS2dProvider("AUSMPlusPMS2D");
      
MethodStrategyProvider<AUSMPlusFluxPrec<MultiScalarVarSet<Euler3DVarSet> >,
		       CellCenterFVMData,
		       FluxSplitter<CellCenterFVMData>,
		       FiniteVolumeNavierStokesModule>
ausmPlusPrecMS3dProvider("AUSMPlusPMS3D");
      
//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
