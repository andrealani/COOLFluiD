#include "FiniteVolumeMultiFluidMHD/ThreeFluidMHDST2D.hh"
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

MethodStrategyProvider<ThreeFluidMHDST2D<MultiFluidMHDVarSet<Maxwell2DProjectionVarSet> >,
		       CellCenterFVMData,
		       Framework::ComputeSourceTerm<CellCenterFVMData>,
		       FiniteVolumeMultiFluidMHDModule>
ThreeFluidMHDST2DProvider("ThreeFluidMHDST2D"); //This word is used for the CFcase (check the source term)
//example: Simulator.SubSystem.CellCenterFVM.Data.SourceTerm = ThreeFluidMHDST2D

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
