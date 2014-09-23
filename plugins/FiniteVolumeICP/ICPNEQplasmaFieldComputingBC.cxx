#include "FiniteVolumeICP/ICPplasmaFieldComputingBC.hh"
#include "FiniteVolumeICP/FiniteVolumeICPNEQ.hh"
#include "FiniteVolumeNavierStokes/NoSlipWallIsothermalNSrvt.hh"
#include "FiniteVolumeNEQ/SubOutletEulerP.hh"
#include "FiniteVolumeNEQ/SubInletEulerMassFlowT.hh"
#include "FiniteVolume/MirrorVelocity.hh"
#include "FiniteVolume/SuperOutlet.hh"
#include "Framework/MethodCommandProvider.hh"
#include "NavierStokes/Euler2DVarSet.hh"
#include "Framework/MultiScalarTerm.hh"
#include "ICP/ICPReactionTerm.hh"
#include "NEQ/NEQReactionTerm.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Numerics::FiniteVolume;
using namespace COOLFluiD::Physics::NavierStokes;
using namespace COOLFluiD::Physics::ICP;
using namespace COOLFluiD::Physics::NEQ;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolumeICP {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<ICPplasmaFieldComputingBC<MirrorVelocity, ICPReactionTerm<NEQReactionTerm> >, 
		      CellCenterFVMData, 
		      FiniteVolumeICPNEQModule> 
EpComputingMirrorICPNEQ2DPuvtFVMCCProvider("EpComputingMirrorICPNEQ2DPuvtFVMCC");

MethodCommandProvider<ICPplasmaFieldComputingBC<NoSlipWallIsothermalNSrvt<MultiScalarTerm<EulerTerm> >, ICPReactionTerm<NEQReactionTerm> >, 
		      CellCenterFVMData, 
		      FiniteVolumeICPNEQModule> 
EpComputingNoSlipWallIsothermalICPNEQPvtFVMCCProvider("EpComputingNoSlipWallIsothermalICPNEQPvtFVMCC");

MethodCommandProvider<ICPplasmaFieldComputingBC<SubOutletEulerP<Euler2DVarSet>, ICPReactionTerm<NEQReactionTerm> >,
		      CellCenterFVMData, 
		      FiniteVolumeICPNEQModule> 
EpComputingSubOutletICPNEQ2DPuvtFVMCCProvider("EpComputingSubOutletICPNEQ2DPuvtFVMCC");

MethodCommandProvider<ICPplasmaFieldComputingBC<SuperOutlet, ICPReactionTerm<NEQReactionTerm> >,
		      CellCenterFVMData, 
		      FiniteVolumeICPNEQModule> 
EpComputingSuperOutletICPNEQ2DPuvtFVMCCProvider("EpComputingSuperOutletICPNEQ2DPuvtFVMCC");

MethodCommandProvider<ICPplasmaFieldComputingBC<SubInletEulerMassFlowT, ICPReactionTerm<NEQReactionTerm> >,
		      CellCenterFVMData, 
		      FiniteVolumeICPNEQModule> 
EpComputingSubInletICPNEQ2DPuvtUVTFVMCCProvider("EpComputingSubInletICPNEQ2DPuvtUVTFVMCC");

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD
