#include "FiniteVolumeNavierStokes/FiniteVolumeNavierStokes.hh"
#include "FiniteVolumeNavierStokes/CoupledNoSlipWallIsothermalNSvt_Ghost.hh"
#include "Framework/MethodCommandProvider.hh"
#include "NavierStokes/EulerTerm.hh"
#include "Framework/MeshData.hh"
#include "Framework/SubSystemStatus.hh"
#include "Framework/NamespaceSwitcher.hh"
#include "Framework/MapGeoToTrsAndIdx.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Physics::NavierStokes;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<CoupledNoSlipWallIsothermalNSvt_Ghost<EulerTerm>,
		      CellCenterFVMData,
		      FiniteVolumeNavierStokesModule>
coupledNoSlipWallIsothermalNSvtGhostFVMCCProvider("CoupledNoSlipWallIsothermalNSvt_GhostFVMCC");

MethodCommandProvider<CoupledNoSlipWallIsothermalNSvt_Ghost<EulerTerm>,
		      CellCenterFVMData,
		      FiniteVolumeNavierStokesModule>
coupledNoSlipWallIsothermalNSPvtGhostFVMCCProvider("CoupledNoSlipWallIsothermalNSPvt_GhostFVMCC");

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD
