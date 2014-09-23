#include "Framework/MethodStrategyProvider.hh"
#include "NavierStokes/NavierStokesVarSet.hh"
#include "FiniteVolumeNavierStokes/NSFluxCoupling.hh"
#include "FiniteVolumeNavierStokes/FiniteVolumeNavierStokes.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Physics::NavierStokes;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

MethodStrategyProvider<NSFluxCoupling<NavierStokesVarSet>,
                       CellCenterFVMData,
		       ComputeDiffusiveFlux,
                       FiniteVolumeNavierStokesModule>
nsFluxCouplingProvider("NavierStokesCoupling");

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
