#include "FiniteVolumeNavierStokes/NavierStokes2DAxiSourceTerm.hh"
#include "NavierStokes2DNEQSourceTerm.hh"
#include "NavierStokes/Euler2DVarSet.hh"
#include "NavierStokes/MultiScalarVarSet.hh"
#include "NEQ/NavierStokesCNEQVarSet.hh"
#include "NEQ/NavierStokesTCNEQVarSet.hh"
#include "NavierStokes/NavierStokes2DVarSet.hh"
#include "Framework/MethodStrategyProvider.hh"
#include "FiniteVolumeNEQ/FiniteVolumeNEQ.hh"
#include "FiniteVolume/CellCenterFVMData.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Physics::NavierStokes;
using namespace COOLFluiD::Physics::NEQ;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

MethodStrategyProvider<NavierStokes2DAxiSourceTerm
		       <MultiScalarVarSet<Euler2DVarSet>,
			NavierStokesCNEQVarSet<NavierStokes2DVarSet> >,
		       CellCenterFVMData,
		       ComputeSourceTerm<CellCenterFVMData>,
		       FiniteVolumeNEQModule>
navierStokes2DNEQAxiSTFVMCCProvider("NavierStokes2DNEQAxiST");

MethodStrategyProvider<NavierStokes2DAxiSourceTerm
		       <MultiScalarVarSet<Euler2DVarSet>,
			NavierStokesTCNEQVarSet<NavierStokes2DVarSet> >,
		       CellCenterFVMData,
		       ComputeSourceTerm<CellCenterFVMData>,
		       FiniteVolumeNEQModule>
navierStokes2DTCNEQAxiSTFVMCCProvider("NavierStokes2DTCNEQAxiST");

MethodStrategyProvider<NavierStokes2DNEQSourceTerm
		       <MultiScalarVarSet<Euler2DVarSet>,
			NavierStokesCNEQVarSet<NavierStokes2DVarSet> >,
		       CellCenterFVMData,
		       ComputeSourceTerm<CellCenterFVMData>,
		       FiniteVolumeNEQModule>
navierStokes2DNEQSTFVMCCProvider("NavierStokes2DNEQST");

MethodStrategyProvider<NavierStokes2DNEQSourceTerm
		       <MultiScalarVarSet<Euler2DVarSet>,
			NavierStokesTCNEQVarSet<NavierStokes2DVarSet> >,
		       CellCenterFVMData,
		       ComputeSourceTerm<CellCenterFVMData>,
		       FiniteVolumeNEQModule>
navierStokes2DTCNEQSTFVMCCProvider("NavierStokes2DTCNEQST");

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
