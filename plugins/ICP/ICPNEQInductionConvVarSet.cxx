#include "ICP/ICPNEQ.hh"
#include "ICP/ICPInductionConvVarSet.hh"
#include "NEQ/Euler2DNEQCons.hh"
#include "NEQ/Euler2DNEQPivtTv.hh"
#include "NEQ/Euler2DNEQRhoivtTv.hh"
#include "Environment/ObjectProvider.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Physics::NavierStokes;
using namespace COOLFluiD::Physics::NEQ;
using namespace COOLFluiD::Physics::ICP;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace ICP {

//////////////////////////////////////////////////////////////////////////////

/// Cons VarSet does nothing but is needed for giving variables so it 
/// is parametrized with the non-LTE varset.
Environment::ObjectProvider<ICPInductionConvVarSet<Euler2DNEQCons>, 
			    ConvectiveVarSet, ICPNEQModule, 1> 
icpNEQ2DConsConvProvider("ICPNEQ2DCons");
    
// this choice is ill conditioned if dp=0
Environment::ObjectProvider<ICPInductionConvVarSet<Euler2DNEQPivt>, 
			    ConvectiveVarSet, ICPNEQModule, 1> 
icpNEQ2DPivtConvProvider("ICPNEQ2DPivt");

// this choice is ill conditioned if dp=0
Environment::ObjectProvider<ICPInductionConvVarSet<Euler2DNEQPivtTv>, 
			    ConvectiveVarSet, ICPNEQModule, 1> 
icpNEQ2DPivtTvConvProvider("ICPNEQ2DPivtTv");

Environment::ObjectProvider<ICPInductionConvVarSet<Euler2DNEQRhoivt>, 
			    ConvectiveVarSet, ICPNEQModule, 1> 
icpNEQ2DRhoivtConvProvider("ICPNEQ2DRhoivt");

Environment::ObjectProvider<ICPInductionConvVarSet<Euler2DNEQRhoivtTv>, 
			    ConvectiveVarSet, ICPNEQModule, 1> 
icpNEQ2DRhoivtTvConvProvider("ICPNEQ2DRhoivtTv");

//////////////////////////////////////////////////////////////////////////////

    } // namespace ICP

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
