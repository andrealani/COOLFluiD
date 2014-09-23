#include "NavierStokes/Euler2DVarSet.hh"
#include "NavierStokes/MultiScalarVarSet.hh"
#include "FiniteVolumeICP/FiniteVolumeICPNEQ.hh"
#include "FiniteVolumeICP/ICPInductionEquationSourceTerm.hh"
#include "Framework/MethodStrategyProvider.hh"
#include "FiniteVolume/CellCenterFVMData.hh"
#include "ICP/ICPReactionTerm.hh"
#include "NEQ/NEQReactionTerm.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace COOLFluiD::Environment;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Physics::NavierStokes;
using namespace COOLFluiD::Physics::ICP;
using namespace COOLFluiD::Physics::NEQ;
using namespace COOLFluiD::Numerics::FiniteVolume;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {
 
    namespace FiniteVolumeICP {

//////////////////////////////////////////////////////////////////////////////

MethodStrategyProvider<ICPInductionEquationSourceTerm
		       <MultiScalarVarSet<Euler2DVarSet>,
			ICPReactionTerm<NEQReactionTerm> >, 
		       CellCenterFVMData, 
		       ComputeSourceTerm<CellCenterFVMData>, 
		       FiniteVolumeICPNEQModule> 
icpneqInductionEquationSTFVMCCProvider("ICPNEQInductionEquationST");

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolumeICP

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
