#include "FiniteVolumeICP/ICPplasmaFieldComputingBC.hh"
#include "FiniteVolumeICP/FiniteVolumeICP.hh"
#include "FiniteVolumeNavierStokes/NoSlipWallIsothermalNSvt.hh"
#include "FiniteVolumeNavierStokes/NoSlipWallIsothermalNSPvt.hh"
#include "FiniteVolume/MirrorVelocity.hh"
#include "FiniteVolume/SuperOutlet.hh"
#include "FiniteVolumeNavierStokes/SubInletEulerPvtVTLTE.hh"
#include "FiniteVolumeNavierStokes/SubOutletEulerPvt.hh"
#include "Framework/MethodCommandProvider.hh"
#include "NavierStokes/Euler2DVarSet.hh"
#include "Framework/MultiScalarTerm.hh"
#include "ICP/ICPReactionTerm.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Numerics::FiniteVolume;
using namespace COOLFluiD::Physics::NavierStokes;
using namespace COOLFluiD::Physics::ICP;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolumeICP {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<ICPplasmaFieldComputingBC<MirrorVelocity, ICPReactionTerm<BaseTerm> >, 
		      CellCenterFVMData, 
		      FiniteVolumeICPModule> 
EpComputingMirrorICP2DPuvtFVMCCProvider("EpComputingMirrorICP2DPuvtFVMCC");

MethodCommandProvider<ICPplasmaFieldComputingBC<NoSlipWallIsothermalNSvt<MultiScalarTerm<EulerTerm> >, ICPReactionTerm<BaseTerm> >, 
		      CellCenterFVMData, 
		      FiniteVolumeICPModule> 
EpComputingNoSlipWallIsothermalICPPvtFVMCCProvider("EpComputingNoSlipWallIsothermalICPPvtFVMCC");
        
MethodCommandProvider<ICPplasmaFieldComputingBC<SubOutletEulerPvt<Euler2DVarSet>, ICPReactionTerm<BaseTerm> >,
		      CellCenterFVMData, 
		      FiniteVolumeICPModule> 
EpComputingSubOutletICP2DPuvtFVMCCProvider("EpComputingSubOutletICP2DPuvtFVMCC");

MethodCommandProvider<ICPplasmaFieldComputingBC<SuperOutlet, ICPReactionTerm<BaseTerm> >,
		      CellCenterFVMData, 
		      FiniteVolumeICPModule> 
EpComputingSuperOutletICP2DPuvtFVMCCProvider("EpComputingSuperOutletICP2DPuvtFVMCC");

MethodCommandProvider<ICPplasmaFieldComputingBC<SubInletEulerPvtVTLTE, ICPReactionTerm<BaseTerm> >,
		      CellCenterFVMData, 
		      FiniteVolumeICPModule> 
EpComputingSubInletICP2DPuvtUVTFVMCCProvider("EpComputingSubInletICP2DPuvtUVTFVMCC");

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD
