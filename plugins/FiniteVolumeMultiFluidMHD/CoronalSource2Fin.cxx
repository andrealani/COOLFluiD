#include "FiniteVolumeMultiFluidMHD/CoronalSource2Fin.hh"
#include "FiniteVolumeMultiFluidMHD/FiniteVolumeMultiFluidMHD.hh"
#include "Framework/MethodStrategyProvider.hh"
//#include "Maxwell/Maxwell2DProjectionVarSet.hh"
#include "Maxwell/Maxwell3DProjectionVarSet.hh"
#include "MultiFluidMHD/MultiFluidMHDVarSet.hh"
#include "Framework/MethodStrategyProvider.hh"
#include "FiniteVolume/CellCenterFVMData.hh"
#include "MultiFluidMHD/DiffMFMHDVarSet.hh"

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

// MethodStrategyProvider<CoronalSource2Fin<MultiFluidMHDVarSet<Maxwell2DProjectionVarSet> >,
//		       CellCenterFVMData,
//		       Framework::ComputeSourceTerm<CellCenterFVMData>,
//		       FiniteVolumeMultiFluidMHDModule>
// twoFluidGravMHDST2DProvider("CoronalSource2Fin");

MethodStrategyProvider<CoronalSource2Fin<MultiFluidMHDVarSet<Maxwell3DProjectionVarSet> >,
		       CellCenterFVMData,
 		       Framework::ComputeSourceTerm<CellCenterFVMData>,
 		       FiniteVolumeMultiFluidMHDModule>
coronalSource2Fin("CoronalSource2Fin");

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
