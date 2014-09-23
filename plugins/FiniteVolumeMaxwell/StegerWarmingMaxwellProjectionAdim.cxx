#include "FiniteVolumeMaxwell/StegerWarmingMaxwellProjectionAdim2D.hh"
#include "FiniteVolumeMaxwell/StegerWarmingMaxwellProjectionAdim3D.hh"
#include "FiniteVolumeMaxwell/FiniteVolumeMaxwell.hh"
#include "Framework/MethodStrategyProvider.hh"
#include "Maxwell/Maxwell2DProjectionAdimVarSet.hh"
#include "Maxwell/Maxwell3DProjectionAdimVarSet.hh"
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

MethodStrategyProvider<StegerWarmingMaxwellProjectionAdim2D<Maxwell2DProjectionAdimVarSet>,
		       CellCenterFVMData,
		       FluxSplitter<CellCenterFVMData>,
		       FiniteVolumeMaxwellModule>
stegerWarmingMaxwellProjectionAdim2DProvider("StegerWarmingMaxwellProjectionAdim2D");

MethodStrategyProvider<StegerWarmingMaxwellProjectionAdim3D<Maxwell3DProjectionAdimVarSet>,
		       CellCenterFVMData,
		       FluxSplitter<CellCenterFVMData>,
		       FiniteVolumeMaxwellModule>
stegerWarmingMaxwellProjectionAdim3DProvider("StegerWarmingMaxwellProjectionAdim3D");

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
