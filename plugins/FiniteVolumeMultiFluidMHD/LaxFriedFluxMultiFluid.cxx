#include "FiniteVolumeMultiFluidMHD/LaxFriedFluxMultiFluid.hh"
#include "MultiFluidMHD/MultiFluidMHDVarSet.hh"
#include "Maxwell/Maxwell2DProjectionVarSet.hh"
#include "Maxwell/Maxwell3DProjectionVarSet.hh"
#include "Framework/MethodStrategyProvider.hh"
#include "FiniteVolumeMultiFluidMHD/FiniteVolumeMultiFluidMHD.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Physics::Maxwell;
using namespace COOLFluiD::Physics::MultiFluidMHD;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

MethodStrategyProvider<LaxFriedFluxMultiFluid<MultiFluidMHDVarSet<Maxwell2DProjectionVarSet> >,
		       CellCenterFVMData,
                       FluxSplitter<CellCenterFVMData>,
		       FiniteVolumeMultiFluidMHDModule>
LaxFriedFlux2DMultiFluidProvider("LaxFriedFluxMultiFluid2D");

MethodStrategyProvider<LaxFriedFluxMultiFluid<MultiFluidMHDVarSet<Maxwell3DProjectionVarSet> >,
		       CellCenterFVMData,
                       FluxSplitter<CellCenterFVMData>,
		       FiniteVolumeMultiFluidMHDModule>
LaxFriedFlux3DMultiFluidProvider("LaxFriedFluxMultiFluid3D");

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
