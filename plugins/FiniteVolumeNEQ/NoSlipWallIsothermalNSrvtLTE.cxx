#include "FiniteVolumeNEQ/NoSlipWallIsothermalNSrvtLTE.hh"
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

MethodCommandProvider<NoSlipWallIsothermalNSrvtLTE
		      <MultiScalarTerm<EulerTerm> >, 
		      CellCenterFVMData, 
		      FiniteVolumeNEQModule> 
noSlipWallIsothermalNSrvtLTEMultiFVMCCProvider
("NoSlipWallIsothermalNSrvtLTEMultiFVMCC");
  
//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD
