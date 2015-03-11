#include "FiniteVolumeNavierStokes/FiniteVolumeNavierStokes.hh"
#include "FiniteVolumeNavierStokes/NoSlipWallAdiabaticNSTurb.hh"
#include "Framework/MethodCommandProvider.hh"
#include "NavierStokes/MultiScalarVarSet.hh"
#include "NavierStokes/NavierStokes2DVarSet.hh"
#include "NavierStokes/NavierStokes3DVarSet.hh"
#include "NavierStokes/NavierStokesTurbVarSet.hh"
#include "NavierStokes/Euler2DPuvt.hh"
#include "NavierStokes/Euler3DPvt.hh"
#include "NavierStokes/Euler3DVarSet.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Physics::NavierStokes;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<NoSlipWallAdiabaticNSTurb
		      <MultiScalarVarSet<Euler2DPuvt>, 
		       NavierStokesTurbVarSet<NavierStokes2DVarSet, 0> >, 
		      CellCenterFVMData,
		      FiniteVolumeNavierStokesModule>
NoSlipWallAdiabaticNSTurb2DFVMCCProvider("NoSlipWallAdiabaticNSTurb2DFVMCC");

MethodCommandProvider<NoSlipWallAdiabaticNSTurb 
		      <MultiScalarVarSet<Euler3DPvt<Euler3DVarSet> >, 
		       NavierStokesTurbVarSet<NavierStokes3DVarSet, 0> >,
		      CellCenterFVMData,
		      FiniteVolumeNavierStokesModule>
NoSlipWallAdiabaticNSTurb3DFVMCCProvider("NoSlipWallAdiabaticNSTurb3DFVMCC");

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD
