#include "FiniteVolumeMaxwell/StegerWarmingMaxwellProjection2D.hh"
#include "FiniteVolumeMaxwell/StegerWarmingMaxwellProjection3D.hh"
#include "FiniteVolumeMaxwell/FiniteVolumeMaxwell.hh"
#include "Framework/MethodStrategyProvider.hh"
#include "Maxwell/Maxwell2DProjectionVarSet.hh"
#include "Maxwell/Maxwell3DProjectionVarSet.hh"
//#include "MultiFluidMHD/MultiFluidMHDVarSet.hh"
#include "Framework/MethodStrategyProvider.hh"
#include "FiniteVolume/CellCenterFVMData.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Physics::Maxwell;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

MethodStrategyProvider<StegerWarmingMaxwellProjection2D<Maxwell2DProjectionVarSet>,
		       CellCenterFVMData,
		       FluxSplitter<CellCenterFVMData>,
		       FiniteVolumeMaxwellModule>
stegerWarmingMaxwellProjection2DProvider("StegerWarmingMaxwellProjection2D");

MethodStrategyProvider<StegerWarmingMaxwellProjection3D<Maxwell3DProjectionVarSet>,
		       CellCenterFVMData,
		       FluxSplitter<CellCenterFVMData>,
		       FiniteVolumeMaxwellModule>
stegerWarmingMaxwellProjection3DProvider("StegerWarmingMaxwellProjection3D");

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
