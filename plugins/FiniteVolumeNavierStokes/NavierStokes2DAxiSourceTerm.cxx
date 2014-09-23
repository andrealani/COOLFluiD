#include "NavierStokes2DAxiSourceTerm.hh"
#include "NavierStokes/Euler2DVarSet.hh"
#include "NavierStokes/NavierStokesVarSet.hh"
#include "NavierStokes/NavierStokes2DVarSet.hh"
#include "NavierStokes/MultiScalarVarSet.hh"
#include "FiniteVolumeNavierStokes/FiniteVolumeNavierStokes.hh"
#include "Framework/MethodStrategyProvider.hh"
#include "FiniteVolume/CellCenterFVMData.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Physics::NavierStokes;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

MethodStrategyProvider<NavierStokes2DAxiSourceTerm
		       <Euler2DVarSet,NavierStokes2DVarSet>,
		       CellCenterFVMData, 
		       ComputeSourceTerm<CellCenterFVMData>,
		       FiniteVolumeNavierStokesModule>
navierStokes2DAxiSTFVMCCProvider("NavierStokes2DAxiST");

MethodStrategyProvider<NavierStokes2DAxiSourceTerm
		       <MultiScalarVarSet<Euler2DVarSet>,NavierStokesVarSet>,
		       CellCenterFVMData, 
		       ComputeSourceTerm<CellCenterFVMData>,
		       FiniteVolumeNavierStokesModule>
navierStokes2DMSAxiSTFVMCCProvider("NavierStokes2DMSAxiST");

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
