#include "FiniteVolumeMultiFluidMHD/DriftWaves2DHalfThreeFluid.hh"
#include "FiniteVolumeMultiFluidMHD/FiniteVolumeMultiFluidMHD.hh"
#include "Framework/MethodStrategyProvider.hh"
#include "Maxwell/Maxwell2DProjectionVarSet.hh"
#include "MultiFluidMHD/MultiFluidMHDVarSet.hh"
#include "Framework/MethodStrategyProvider.hh"
#include "FiniteVolume/CellCenterFVMData.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Physics::Maxwell;
using namespace COOLFluiD::Physics::MultiFluidMHD;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

MethodStrategyProvider<DriftWaves2DHalfThreeFluid<MultiFluidMHDVarSet<Maxwell2DProjectionVarSet> >,
		       CellCenterFVMData,
		       Framework::ComputeSourceTerm<CellCenterFVMData>,
		       FiniteVolumeMultiFluidMHDModule>
DriftWaves2DHalfThreeFluidProvider("DriftWaves2DHalfThreeFluid"); //This word is used for the CFcase (check the source term)
//example: Simulator.SubSystem.CellCenterFVM.Data.SourceTerm = DriftWaves2DHalfThreeFluid

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
