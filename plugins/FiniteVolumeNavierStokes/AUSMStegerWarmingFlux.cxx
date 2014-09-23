#include "FiniteVolumeNavierStokes/AUSMPlusUpFlux.hh"
#include "NavierStokes/Euler2DVarSet.hh"
#include "NavierStokes/Euler3DVarSet.hh"
#include "FiniteVolume/StegerWarmingFluxT.hh"
#include "Framework/MethodStrategyProvider.hh"
#include "FiniteVolumeNavierStokes/FiniteVolumeNavierStokes.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Physics::NavierStokes;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

MethodStrategyProvider<StegerWarmingFluxT< AUSMPlusUpFlux<Euler2DVarSet> >,
                       CellCenterFVMData,
                       FluxSplitter<CellCenterFVMData>,
                       FiniteVolumeNavierStokesModule>
ausmPlusUpSW2DFluxProvider("AUSMPlusUpStegerWarming2D");

MethodStrategyProvider<StegerWarmingFluxT< AUSMPlusUpFlux<Euler3DVarSet> >,
                       CellCenterFVMData,
                       FluxSplitter<CellCenterFVMData>,
                       FiniteVolumeNavierStokesModule>
ausmPlusUpSW3DFluxProvider("AUSMPlusUpStegerWarming3D");

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
