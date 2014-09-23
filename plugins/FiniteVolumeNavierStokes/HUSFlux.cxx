#include "Framework/MethodStrategyProvider.hh"
#include "FiniteVolumeNavierStokes/FiniteVolumeNavierStokes.hh"
#include "NavierStokes/MultiScalarVarSet.hh"
#include "NavierStokes/Euler2DVarSet.hh"
#include "NavierStokes/Euler3DVarSet.hh"
#include "FiniteVolumeNavierStokes/HUSFlux2D.hh"
#include "FiniteVolumeNavierStokes/HUSFlux3D.hh"

//////////////////////////////////////////////////////////////////////

using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Physics::NavierStokes;

//////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////
      
MethodStrategyProvider<HUSFlux2D<Euler2DVarSet>,
		       CellCenterFVMData,
		       FluxSplitter<CellCenterFVMData>,
		       FiniteVolumeNavierStokesModule>
HUSFlux2DProvider("HUS2D");
  
MethodStrategyProvider<HUSFlux2D
		       <MultiScalarVarSet<Euler2DVarSet> >,
			CellCenterFVMData,
			FluxSplitter<CellCenterFVMData>,
			FiniteVolumeNavierStokesModule>
HUSFluxMS2DProvider("HUSMS2D");
      
MethodStrategyProvider<HUSFlux3D<Euler3DVarSet>,
 		       CellCenterFVMData,
		       FluxSplitter<CellCenterFVMData>,
 		       FiniteVolumeNavierStokesModule>
HUSFlux3DProvider("HUS3D");
  
MethodStrategyProvider<HUSFlux3D 
 		       <MultiScalarVarSet<Euler3DVarSet> >,
 			CellCenterFVMData,
			FluxSplitter<CellCenterFVMData>,
 			FiniteVolumeNavierStokesModule>
HUSFluxMS3DProvider("HUSMS3D");
  
//////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////
