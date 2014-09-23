#include "FiniteVolumeNEQ/CoupledNoSlipWallIsothermalNSrvtLTE_Nodes.hh"
#include "FiniteVolumeNEQ/FiniteVolumeNEQ.hh"
#include "NavierStokes/EulerTerm.hh"
#include "Framework/MultiScalarTerm.hh"

#include "Framework/MethodCommandProvider.hh"
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

MethodCommandProvider<CoupledNoSlipWallIsothermalNSrvtLTE_Nodes<MultiScalarTerm<EulerTerm> >,
		      CellCenterFVMData,
		      FiniteVolumeNEQModule>
coupledNoSlipWallIsothermalNSrvtLTE_NodesFVMCCProvider("CoupledNoSlipWallIsothermalNSrvtLTE_NodesFVMCC");

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD
