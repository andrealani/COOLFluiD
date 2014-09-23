#include "FiniteVolumeMaxwell/MunzFluxMaxwell2D.hh"
#include "FiniteVolumeMaxwell/MunzFluxMaxwell3D.hh"
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

MethodStrategyProvider<MunzFluxMaxwell2D<Maxwell2DProjectionVarSet>,
		       CellCenterFVMData,
		       FluxSplitter<CellCenterFVMData>,
		       FiniteVolumeMaxwellModule>
MunzFluxMaxwell2DProvider("MunzFluxMaxwell2D");

MethodStrategyProvider<MunzFluxMaxwell3D<Maxwell3DProjectionVarSet>,
		       CellCenterFVMData,
		       FluxSplitter<CellCenterFVMData>,
		       FiniteVolumeMaxwellModule>
MunzFluxMaxwell3DProvider("MunzFluxMaxwell3D");

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
