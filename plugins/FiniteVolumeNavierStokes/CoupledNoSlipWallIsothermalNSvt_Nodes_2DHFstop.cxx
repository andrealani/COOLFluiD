#include "FiniteVolumeNavierStokes/FiniteVolumeNavierStokes.hh"
#include "FiniteVolumeNavierStokes/CoupledNoSlipWallIsothermalNSvt_Nodes_2DHF.hh"
#include "FiniteVolumeNavierStokes/CoupledNoSlipWallIsothermalNSrvt_Nodes.hh"
#include "Framework/MethodCommandProvider.hh"
#include "NavierStokes/EulerTerm.hh"
#include "Framework/MultiScalarTerm.hh"

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

MethodCommandProvider<CoupledNoSlipWallIsothermalNSvt_Nodes_2DHF<EulerTerm>,
		      CellCenterFVMData,
		      FiniteVolumeNavierStokesModule>
coupledNoSlipWallIsothermalNSPvtGhostFVMCCProvider("CoupledNoSlipWallIsothermalNSvt_NodesFVMCC");

MethodCommandProvider<CoupledNoSlipWallIsothermalNSvt_Nodes_2DHF<EulerTerm>,
		      CellCenterFVMData,
		      FiniteVolumeNavierStokesModule>
coupledNoSlipWallIsothermalNSPvtGhostFVMCCProvider("CoupledNoSlipWallIsothermalNSPvt_NodesFVMCC");

MethodCommandProvider<CoupledNoSlipWallIsothermalNSvt_Nodes_2DHF
		      <MultiScalarTerm<EulerTerm> >,
		      CellCenterFVMData,
		      FiniteVolumeNavierStokesModule>
coupledNoSlipWallEulerMultiIsotNSvtFVMCCProvider
("CoupledNoSlipWallIsothermalNSvtMultiFVMCC_Nodes");

MethodCommandProvider<CoupledNoSlipWallIsothermalNSrvt_Nodes
		      <MultiScalarTerm<EulerTerm> >,
		      CellCenterFVMData,
		      FiniteVolumeNavierStokesModule>
CoupledNoSlipWallIsothermalNSvt_NodesFVMCCProvider
("CoupledNoSlipWallIsothermalNSrvt_NodesFVMCC");


//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD
