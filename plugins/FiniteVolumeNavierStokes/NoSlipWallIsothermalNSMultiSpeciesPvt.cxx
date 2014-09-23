#include "FiniteVolumeNavierStokes/NoSlipWallIsothermalNSMultiSpeciesPvt.hh"
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

MethodCommandProvider<NoSlipWallIsothermalNSMultiSpeciesPvt
		      <MultiScalarTerm<EulerTerm> >, 
		      CellCenterFVMData, 
		      FiniteVolumeNavierStokesModule> 
noSlipWallIsotNSLTEDemixPvtFVMCCProvider
("NoSlipWallIsothermalNSLTEDemixPvtFVMCC");
      
MethodCommandProvider<NoSlipWallIsothermalNSMultiSpeciesPvt
		      <MultiScalarTerm<EulerTerm> >,
		      CellCenterFVMData, 
		      FiniteVolumeNavierStokesModule> 
noSlipWallIsotNSChemNEQPvtFVMCCProvider
("NoSlipWallIsothermalNSChemNEQPvtFVMCC");
                                        
//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD
