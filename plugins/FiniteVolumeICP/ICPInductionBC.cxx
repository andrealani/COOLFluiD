#include "FiniteVolumeICP/ICPInductionBC.hh"
#include "FiniteVolumeICP/FiniteVolumeICP.hh"
#include "FiniteVolumeNavierStokes/NoSlipWallIsothermalNSvt.hh"
#include "FiniteVolumeNavierStokes/NoSlipWallIsothermalNSPvt.hh"
#include "FiniteVolume/MirrorVelocity.hh"
#include "FiniteVolumeNavierStokes/SubInletEuler2DPuvtUVTLTE.hh"
#include "FiniteVolumeNavierStokes/SubOutletEulerPvt.hh"
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
		      FiniteVolumeICPModule> 
mirrorICP2DPuvtFVMCCProvider("MirrorICP2DPuvtFVMCC");

MethodCommandProvider<ICPInductionBC<NoSlipWallIsothermalNSvt<MultiScalarTerm<EulerTerm> > >, 
		      CellCenterFVMData, 
		      FiniteVolumeICPModule> 
noSlipWallIsothermalICPPvtFVMCCProvider("NoSlipWallIsothermalICPPvtFVMCC");
        
MethodCommandProvider<ICPInductionBC<SubOutletEulerPvt<Euler2DVarSet> >,
		      CellCenterFVMData, 
		      FiniteVolumeICPModule> 
subOutletICP2DPuvtFVMCCProvider("SubOutletICP2DPuvtFVMCC");

MethodCommandProvider<ICPInductionBC<SuperOutlet>,
		      CellCenterFVMData, 
		      FiniteVolumeICPModule> 
superOutletICP2DPuvtFVMCCProvider("SuperOutletICP2DPuvtFVMCC");

MethodCommandProvider<ICPInductionBC<SubInletEuler2DPuvtUVTLTE>,
		      CellCenterFVMData, 
		      FiniteVolumeICPModule> 
subInletICP2DPuvtUVTFVMCCProvider("SubInletICP2DPuvtUVTFVMCC");

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD
