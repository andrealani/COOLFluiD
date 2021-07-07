#include "FiniteVolumeNEQ/PoissonNEQSourceTerm.hh"
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

MethodStrategyProvider<PoissonNEQSourceTerm
		       <MultiScalarVarSet<Euler2DVarSet>,
			NavierStokesCNEQVarSet<NavierStokes2DVarSet> >,
		       CellCenterFVMData,
		       ComputeSourceTerm<CellCenterFVMData>,
		       FiniteVolumeNEQModule>
poisson2DCNEQSTFVMCCProvider("PoissonCNEQSourceTerm");
      
MethodStrategyProvider<PoissonNEQSourceTerm
		       <MultiScalarVarSet<Euler2DVarSet>,
			NavierStokesTCNEQVarSet<NavierStokes2DVarSet> >,
		       CellCenterFVMData,
		       ComputeSourceTerm<CellCenterFVMData>,
		       FiniteVolumeNEQModule>
poisson2DTCNEQSTFVMCCProvider("PoissonTCNEQSourceTerm");
      
//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
