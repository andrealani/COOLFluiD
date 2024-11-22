#include "FiniteVolumePoissonNEQ/PoissonNEQBC.hh"
#include "FiniteVolumePoissonNEQ/FiniteVolumePoissonNEQ.hh"
#include "FiniteVolumeNavierStokes/NoSlipWallIsothermalNSrvt.hh"
#include "FiniteVolume/MirrorVelocity.hh"
#include "FiniteVolume/SuperInlet.hh"
#include "FiniteVolume/SuperOutlet.hh"
#include "Framework/MultiScalarTerm.hh"
#include "NavierStokes/EulerTerm.hh"
#include "NavierStokes/EulerVarSet.hh"
#include "Framework/MethodCommandProvider.hh"
//Vatsalya: new headers
#include "FiniteVolumeNEQ/SubOutletEulerP.hh"
#include "NavierStokes/Euler3DVarSet.hh" //since we do 3D only
#include "FiniteVolumeNEQ/SubInletInterpYiVTTv.hh"
//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Numerics::FiniteVolume;
using namespace COOLFluiD::Physics::NavierStokes;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolumePoissonNEQ {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<PoissonNEQBC<MirrorVelocity>, 
		      CellCenterFVMData, 
		      FiniteVolumePoissonNEQModule> 
mirrorVelocityPoissonNEQFVMCCProvider("MirrorVelocityPoissonNEQFVMCC");

MethodCommandProvider<PoissonNEQBC<SuperInlet>, 
		      CellCenterFVMData, 
		      FiniteVolumePoissonNEQModule> 
superInletPoissonNEQFVMCCProvider("SuperInletPoissonNEQFVMCC");
      
MethodCommandProvider<PoissonNEQBC<SuperOutlet>, 
		      CellCenterFVMData, 
		      FiniteVolumePoissonNEQModule> 
superOutletPoissonNEQFVMCCProvider("SuperOutletPoissonNEQFVMCC");
      
MethodCommandProvider<PoissonNEQBC<NoSlipWallIsothermalNSrvt<MultiScalarTerm<EulerTerm> > >, 
		      CellCenterFVMData, 
		      FiniteVolumePoissonNEQModule> 
noSlipWallIsothermalNSrvtMultiPoissonNEQFVMCCProvider("NoSlipWallIsothermalNSrvtMultiPoissonNEQFVMCC");

//vatsalya:new bc for ICP
MethodCommandProvider<PoissonNEQBC<SubOutletEulerP<Euler3DVarSet> >,
		      CellCenterFVMData, 
		      FiniteVolumePoissonNEQModule> 
subOutletEuler3DPPoissonNEQFVMCCProvider("SubOutletEuler3DPPoissonNEQFVMCC");

MethodCommandProvider<PoissonNEQBC<SubInletInterpYiVTTv>, 
		      CellCenterFVMData, 
		      FiniteVolumePoissonNEQModule> 
subInletInterpYiVTTvPoissonNEQFVMCCProvider("SubInletInterpYiVTTvPoissonNEQFVMCC");


//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolumePoissonNEQ

  } // namespace Numerics

} // namespace COOLFluiD
