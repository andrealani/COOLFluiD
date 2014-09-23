#include "FiniteVolumeTurb/FiniteVolumeKOmega.hh"
#include "FiniteVolumeTurb/NoSlipWallIsothermalNSKOmPvt.hh"
#include "Framework/MethodCommandProvider.hh"
#include "NavierStokes/NSTurbTerm.hh"
#include "NavierStokes/EulerVarSet.hh"
#include "NavierStokes/NavierStokesTurbVarSet.hh"
#include "NavierStokes/NavierStokes2DVarSet.hh"
#include "NavierStokes/NavierStokes3DVarSet.hh"

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

MethodCommandProvider<NoSlipWallIsothermalNSKOmPvt
		      <EulerVarSet, NavierStokesTurbVarSet<NavierStokes2DVarSet, 0> >,
		      CellCenterFVMData,
		      FiniteVolumeKOmegaModule>
noSlipWallIsothermalNSKOmPvt2DFVMCCProvider("NoSlipWallIsothermalNSKOmPvt2DFVMCC");

MethodCommandProvider<NoSlipWallIsothermalNSKOmPvt
		      <EulerVarSet, NavierStokesTurbVarSet<NavierStokes3DVarSet, 0> >,
		      CellCenterFVMData,
		      FiniteVolumeKOmegaModule>
noSlipWallIsothermalNSKOmPvt3DFVMCCProvider("NoSlipWallIsothermalNSKOmPvt3DFVMCC");

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD
