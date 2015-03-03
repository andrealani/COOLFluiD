#include "FiniteVolumeMultiFluidMHD/GEMMHDST3DTwoFluid.hh"
#include "FiniteVolumeMultiFluidMHD/FiniteVolumeMultiFluidMHD.hh"
#include "Framework/MethodStrategyProvider.hh"
#include "Maxwell/Maxwell3DProjectionVarSet.hh"
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

MethodStrategyProvider<GEMMHDST3DTwoFluid<MultiFluidMHDVarSet<Maxwell3DProjectionVarSet> >,
		       CellCenterFVMData,
		       Framework::ComputeSourceTerm<CellCenterFVMData>,
		       FiniteVolumeMultiFluidMHDModule>
GEMMHDST3DTwoFluidProvider("GEMMHDST3DTwoFluid"); //This word is used for the CFcase (check the source term)
//example: Simulator.SubSystem.CellCenterFVM.Data.SourceTerm = GEMMHDST3DTwoFluid

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
