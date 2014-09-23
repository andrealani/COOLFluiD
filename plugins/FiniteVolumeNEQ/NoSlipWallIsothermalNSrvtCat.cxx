#include "FiniteVolumeNEQ/NoSlipWallIsothermalNSrvtCat.hh"
#include "FiniteVolumeNEQ/NoSlipWallIsothermalNSrvtCatT.hh"
#include "FiniteVolumeNEQ/NoSlipWallIsothermalNSrvtCatR.hh"
#include "FiniteVolumeNEQ/FiniteVolumeNEQ.hh"
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

MethodCommandProvider<NoSlipWallIsothermalNSrvtCat
		      <MultiScalarTerm<EulerTerm> >, 
		      CellCenterFVMData, 
		      FiniteVolumeNEQModule> 
noSlipWallIsothermalNSrvtCatMultiFVMCCProvider
("NoSlipWallIsothermalNSrvtCatMultiFVMCC");
  
MethodCommandProvider<NoSlipWallIsothermalNSrvtCatT
		      <MultiScalarTerm<EulerTerm> >, 
		      CellCenterFVMData, 
		      FiniteVolumeNEQModule> 
noSlipWallIsothermalNSrvtCatTMultiFVMCCProvider
("NoSlipWallIsothermalNSrvtCatFVMCC");
      
MethodCommandProvider<NoSlipWallIsothermalNSrvtCatR
		      <MultiScalarTerm<EulerTerm> >, 
		      CellCenterFVMData, 
		      FiniteVolumeNEQModule> 
noSlipWallIsothermalNSrvtCatRMultiFVMCCProvider
("NoSlipWallIsothermalNSrvtCatRFVMCC");

/*MethodCommandProvider<NoSlipWallIsothermalNSrvtCat_oldT
		      <MultiScalarTerm<EulerTerm> >, 
		      CellCenterFVMData, 
		      FiniteVolumeNEQModule> 
noSlipWallIsothermalNSrvtCatOLDTMultiFVMCCProvider
("NoSlipWallIsothermalNSrvtCat_oldTFVMCC");*/

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD
