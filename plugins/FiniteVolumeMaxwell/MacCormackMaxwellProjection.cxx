#include "FiniteVolumeMaxwell/MacCormackMaxwell2DProjection.hh"
#include "FiniteVolumeMaxwell/MacCormackMaxwell3DProjection.hh"
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

MethodStrategyProvider<MacCormackMaxwell2DProjection<Maxwell2DProjectionVarSet>,
		       CellCenterFVMData,
		       FluxSplitter<CellCenterFVMData>,
		       FiniteVolumeMaxwellModule>
macCormackMaxwell2DProjectionProvider("MacCormackMaxwell2DProjection");

MethodStrategyProvider<MacCormackMaxwell3DProjection<Maxwell3DProjectionVarSet>,
		       CellCenterFVMData,
		       FluxSplitter<CellCenterFVMData>,
		       FiniteVolumeMaxwellModule>
macCormackMaxwell3DProjectionProvider("MacCormackMaxwell3DProjection");

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
