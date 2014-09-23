#include "FiniteVolumeGReKO/FiniteVolumeGReKO.hh"
#include "FiniteVolumeGReKO/NoSlipWallIsothermalNSGReKOPvt.hh"
#include "Framework/MethodCommandProvider.hh"
#include "NavierStokes/NSTurbTerm.hh"
#include "NavierStokes/Euler2DVarSet.hh"
#include "NavierStokes/Euler3DVarSet.hh"
#include "NavierStokes/MultiScalarVarSet.hh"
#include "NavierStokes/NavierStokesTurbVarSetTypes.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Physics::NavierStokes;
using namespace COOLFluiD::MathTools;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<NoSlipWallIsothermalNSGReKOPvt
		      <MultiScalarVarSet<Euler2DVarSet>,
		       NavierStokesTurb2DVarSet>,
		      CellCenterFVMData,
		      FiniteVolumeGReKOModule>
noSlipWallIsothermalNSGReKOPvt2DFVMCCProvider("NoSlipWallIsothermalNSGReKOPvt2DFVMCC");

//MethodCommandProvider<NoSlipWallIsothermalNSGReKOPvt
////		      <MultiScalarVarSet<Euler3DVarSet>,
//		       NavierStokesTurbVarSet>,
//		      CellCenterFVMData,
//		      FiniteVolumeGReKOModule>
//noSlipWallIsothermalNSGReKOPvt3DFVMCCProvider("NoSlipWallIsothermalNSGReKOPvt3DFVMCC");

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD
