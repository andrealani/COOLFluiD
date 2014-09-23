#include "FiniteVolumeNavierStokes/NoSlipWallIsothermalNSrvt.hh"
#include "FiniteVolumeNavierStokes/FiniteVolumeNavierStokes.hh"
#include "Framework/MethodCommandProvider.hh"
#include "NavierStokes/EulerTerm.hh"
#include "Framework/MultiScalarTerm.hh"
  
//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Physics::NavierStokes;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<NoSlipWallIsothermalNSvt<EulerTerm>, 
		      CellCenterFVMData, 
		      FiniteVolumeNavierStokesModule> 
noSlipWallEulerIsotNSvtFVMCCProvider
("NoSlipWallIsothermalNSvtFVMCC");

MethodCommandProvider<NoSlipWallIsothermalNSvt
		      <MultiScalarTerm<EulerTerm> >, 
		      CellCenterFVMData, 
		      FiniteVolumeNavierStokesModule> 
noSlipWallEulerMultiIsotNSvtFVMCCProvider
("NoSlipWallIsothermalNSvtMultiFVMCC");

MethodCommandProvider<NoSlipWallIsothermalNSrvt
		      <MultiScalarTerm<EulerTerm> >, 
		      CellCenterFVMData, 
		      FiniteVolumeNavierStokesModule> 
noSlipWallEulerMultiIsotNSrvtFVMCCProvider
("NoSlipWallIsothermalNSrvtMultiFVMCC");
                                        
//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD
