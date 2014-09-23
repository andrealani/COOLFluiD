#include "NavierStokes/Euler2DVarSet.hh"
#include "FiniteVolumeICP/FiniteVolumeICP.hh"
#include "FiniteVolumeICP/ICPInductionEquationSourceTerm.hh"
#include "Framework/MethodStrategyProvider.hh"
#include "FiniteVolume/CellCenterFVMData.hh"
#include "ICP/ICPReactionTerm.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace COOLFluiD::Environment;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Physics::NavierStokes;
using namespace COOLFluiD::Physics::ICP;
using namespace COOLFluiD::Numerics::FiniteVolume;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolumeICP {

//////////////////////////////////////////////////////////////////////////////

MethodStrategyProvider<ICPInductionEquationSourceTerm
		       <Euler2DVarSet,ICPReactionTerm<BaseTerm> >, 
		       CellCenterFVMData, 
		       ComputeSourceTerm<CellCenterFVMData>, 
		       FiniteVolumeICPModule> 
icpInductionEquationSTFVMCCProvider("ICPInductionEquationST");
    
//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolumeICP

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
