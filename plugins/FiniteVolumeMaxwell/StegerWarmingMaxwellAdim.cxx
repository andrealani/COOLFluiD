#include "FiniteVolumeMaxwell/StegerWarmingMaxwellAdim2D.hh"
#include "FiniteVolumeMaxwell/StegerWarmingMaxwellAdim3D.hh"
#include "FiniteVolumeMaxwell/FiniteVolumeMaxwell.hh"
#include "Framework/MethodStrategyProvider.hh"
#include "Maxwell/Maxwell2DVarSetAdim.hh"
#include "Maxwell/Maxwell3DVarSetAdim.hh"
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

MethodStrategyProvider<StegerWarmingMaxwellAdim2D<Maxwell2DVarSetAdim>,
		       CellCenterFVMData,
		       FluxSplitter<CellCenterFVMData>,
		       FiniteVolumeMaxwellModule>
stegerWarmingMaxwellAdim2DProvider("StegerWarmingMaxwellAdim2D");

MethodStrategyProvider<StegerWarmingMaxwellAdim3D<Maxwell3DVarSetAdim>,
		       CellCenterFVMData,
		       FluxSplitter<CellCenterFVMData>,
		       FiniteVolumeMaxwellModule>
stegerWarmingMaxwellAdim3DProvider("StegerWarmingMaxwellAdim3D");

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
