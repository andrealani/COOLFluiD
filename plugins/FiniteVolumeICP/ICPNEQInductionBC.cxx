#include "FiniteVolumeICP/ICPInductionBC.hh"
#include "FiniteVolumeICP/FiniteVolumeICPNEQ.hh"
#include "FiniteVolumeNavierStokes/NoSlipWallIsothermalNSrvt.hh"
#include "FiniteVolumeNEQ/SubOutletEulerP.hh"
#include "FiniteVolumeNEQ/SubInletEulerMassFlowT.hh"
#include "FiniteVolume/MirrorVelocity.hh"
#include "FiniteVolume/SuperOutlet.hh"
#include "Framework/MethodCommandProvider.hh"
#include "NavierStokes/Euler2DVarSet.hh"
#include "Framework/MultiScalarTerm.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Numerics::FiniteVolume;
using namespace COOLFluiD::Physics::NavierStokes;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolumeICP {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<ICPInductionBC<MirrorVelocity>, 
		      CellCenterFVMData, 
		      FiniteVolumeICPNEQModule> 
mirrorICPNEQ2DPuvtFVMCCProvider("MirrorICPNEQ2DPuvtFVMCC");

MethodCommandProvider<ICPInductionBC<NoSlipWallIsothermalNSrvt<MultiScalarTerm<EulerTerm> > >, 
		      CellCenterFVMData, 
		      FiniteVolumeICPNEQModule> 
noSlipWallIsothermalICPNEQPvtFVMCCProvider("NoSlipWallIsothermalICPNEQPvtFVMCC");
        
MethodCommandProvider<ICPInductionBC<SubOutletEulerP<Euler2DVarSet> >,
		      CellCenterFVMData, 
		      FiniteVolumeICPNEQModule> 
subOutletICPNEQ2DPuvtFVMCCProvider("SubOutletICPNEQ2DPuvtFVMCC");

MethodCommandProvider<ICPInductionBC<SuperOutlet>,
		      CellCenterFVMData, 
		      FiniteVolumeICPNEQModule> 
superOutletICPNEQ2DPuvtFVMCCProvider("SuperOutletICPNEQ2DPuvtFVMCC");

MethodCommandProvider<ICPInductionBC<SubInletEulerMassFlowT>,
		      CellCenterFVMData, 
		      FiniteVolumeICPNEQModule> 
subInletICPNEQ2DPuvtUVTFVMCCProvider("SubInletICPNEQ2DPuvtUVTFVMCC");

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD
