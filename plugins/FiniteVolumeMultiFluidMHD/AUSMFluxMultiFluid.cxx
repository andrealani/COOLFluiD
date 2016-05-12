#include "FiniteVolumeMultiFluidMHD/AUSMPlusUpFluxMultiFluid.hh"
#include "MultiFluidMHD/MultiFluidMHDVarSet.hh"
#include "Maxwell/Maxwell2DProjectionVarSet.hh"
#include "Maxwell/Maxwell3DProjectionVarSet.hh"
#include "Framework/MethodStrategyProvider.hh"
#include "FiniteVolumeMultiFluidMHD/FiniteVolumeMultiFluidMHD.hh"
#include "FiniteVolumeMultiFluidMHD/AUSMFluxMultiFluidALE.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Physics::Maxwell;
using namespace COOLFluiD::Physics::MultiFluidMHD;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

/// AUSM+up

MethodStrategyProvider<AUSMPlusUpFluxMultiFluid<MultiFluidMHDVarSet<Maxwell2DProjectionVarSet> >,
		       CellCenterFVMData,
                       FluxSplitter<CellCenterFVMData>,
		       FiniteVolumeMultiFluidMHDModule>
ausmPlusUp2DMultiFluidProvider("AUSMPlusUpMultiFluid2D");

MethodStrategyProvider<AUSMPlusUpFluxMultiFluid<MultiFluidMHDVarSet<Maxwell3DProjectionVarSet> >,
		       CellCenterFVMData,
                       FluxSplitter<CellCenterFVMData>,
		       FiniteVolumeMultiFluidMHDModule>
ausmPlusUp3DMultiFluidProvider("AUSMPlusUpMultiFluid3D");

/// AUSM+up ALE

MethodStrategyProvider< AUSMFluxMultiFluidALE<AUSMPlusUpFluxMultiFluid<MultiFluidMHDVarSet<Maxwell2DProjectionVarSet> > >,
		       CellCenterFVMData,
                       FluxSplitter<CellCenterFVMData>,
		       FiniteVolumeMultiFluidMHDModule>
ausmPlusUp2DALEMultiFluidProvider("AUSMPlusUpMultiFluid2DALE");

MethodStrategyProvider< AUSMFluxMultiFluidALE<AUSMPlusUpFluxMultiFluid<MultiFluidMHDVarSet<Maxwell3DProjectionVarSet> > >,
		       CellCenterFVMData,
                       FluxSplitter<CellCenterFVMData>,
		       FiniteVolumeMultiFluidMHDModule>
ausmPlusUp3DALEMultiFluidProvider("AUSMPlusUpMultiFluid3DALE");

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
