#include "FiniteVolumeMultiFluidMHD/MultiFluidMHDSTNoRadiation.hh"
#include "FiniteVolumeMultiFluidMHD/FiniteVolumeMultiFluidMHD.hh"
#include "Framework/MethodStrategyProvider.hh"
#include "Maxwell/Maxwell2DProjectionVarSet.hh"
//#include "Maxwell/Maxwell3DProjectionVarSet.hh"
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

MethodStrategyProvider<MultiFluidMHDSTNoRadiation<MultiFluidMHDVarSet<Maxwell2DProjectionVarSet> >,
		       CellCenterFVMData,
		       Framework::ComputeSourceTerm<CellCenterFVMData>,
		       FiniteVolumeMultiFluidMHDModule>
multiFluidMHDSTNoRadiation2DProvider("MultiFluidMHDSTNoRadiation2D");

// MethodStrategyProvider<MultiFluidMHDSTNoRadiation<MultiFluidMHDVarSet<Maxwell3DProjectionVarSet> >,
// 		       CellCenterFVMData,
// 		       Framework::ComputeSourceTerm<CellCenterFVMData>,
// 		       FiniteVolumeMultiFluidMHDModule>
// multiFluidMHDSTNoRadiation3DProvider("MultiFluidMHDST3D");

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
