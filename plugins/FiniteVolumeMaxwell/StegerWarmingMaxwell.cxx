#include "FiniteVolumeMaxwell/StegerWarmingMaxwell2D.hh"
#include "FiniteVolumeMaxwell/StegerWarmingMaxwell3D.hh"
#include "FiniteVolumeMaxwell/FiniteVolumeMaxwell.hh"
#include "Framework/MethodStrategyProvider.hh"
#include "Maxwell/Maxwell2DVarSet.hh"
#include "Maxwell/Maxwell3DVarSet.hh"
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

MethodStrategyProvider<StegerWarmingMaxwell2D<Maxwell2DVarSet>,
		       CellCenterFVMData,
		       FluxSplitter<CellCenterFVMData>,
		       FiniteVolumeMaxwellModule>
StegerWarmingMaxwell2DProvider("StegerWarmingMaxwell2D");

MethodStrategyProvider<StegerWarmingMaxwell3D<Maxwell3DVarSet>,
		       CellCenterFVMData,
		       FluxSplitter<CellCenterFVMData>,
		       FiniteVolumeMaxwellModule>
StegerWarmingMaxwell3DProvider("StegerWarmingMaxwell3D");

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
