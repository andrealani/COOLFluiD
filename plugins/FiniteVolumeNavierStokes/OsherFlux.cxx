#include "FiniteVolumeNavierStokes/OsherFlux2D.hh"
#include "Framework/MethodStrategyProvider.hh"
#include "FiniteVolumeNavierStokes/FiniteVolumeNavierStokes.hh"
#include "NavierStokes/MultiScalarVarSet.hh"
#include "NavierStokes/Euler2DVarSet.hh"
#include "NavierStokes/Euler3DVarSet.hh"

//////////////////////////////////////////////////////////////////////

using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Physics::NavierStokes;

//////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////

MethodStrategyProvider<OsherFlux2D<Euler2DVarSet>,
			CellCenterFVMData,
			FluxSplitter<CellCenterFVMData>,
			FiniteVolumeNavierStokesModule>
OsherFlux2DProvider("OsherFlux2D");
  
MethodStrategyProvider<OsherFlux2D
		       <MultiScalarVarSet<Euler2DVarSet> >,
			CellCenterFVMData,
			FluxSplitter<CellCenterFVMData>,
			FiniteVolumeNavierStokesModule>
OsherFluxMS2DProvider("OsherMS2D");
  
//////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////
