#include "PoissonNEQ/PoissonNEQ.hh"
#include "PoissonNEQ/PoissonNEQConvVarSet.hh"
#include "NEQ/Euler2DNEQCons.hh"
#include "NEQ/Euler2DNEQRhoivtTv.hh"
#include "NEQ/Euler3DNEQCons.hh"
#include "NEQ/Euler3DNEQRhoivtTv.hh"
#include "NEQ/Euler3DNEQRhoivt.hh"
#include "Environment/ObjectProvider.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Physics::NavierStokes;
using namespace COOLFluiD::Physics::NEQ;
using namespace COOLFluiD::Physics::PoissonNEQ;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace PoissonNEQ {

//////////////////////////////////////////////////////////////////////////////

/// Cons VarSet does nothing but is needed for giving variables so it 
/// is parametrized with the NEQ varset.
Environment::ObjectProvider<PoissonNEQConvVarSet<Euler2DNEQCons>, 
			    ConvectiveVarSet, PoissonNEQModule, 1> 
poissonNEQ2DConsConvProvider("PoissonNEQ2DCons");
    
Environment::ObjectProvider<PoissonNEQConvVarSet<Euler2DNEQRhoivtTv>, 
			    ConvectiveVarSet, PoissonNEQModule, 1> 
poissonNEQ2DRhoivtTvConvProvider("PoissonNEQ2DRhoivtTv");

Environment::ObjectProvider<PoissonNEQConvVarSet<Euler3DNEQCons>, 
			    ConvectiveVarSet, PoissonNEQModule, 1> 
poissonNEQ3DConsConvProvider("PoissonNEQ3DCons");
    
Environment::ObjectProvider<PoissonNEQConvVarSet<Euler3DNEQRhoivtTv>, 
			    ConvectiveVarSet, PoissonNEQModule, 1> 
poissonNEQ3DRhoivtTvConvProvider("PoissonNEQ3DRhoivtTv");

Environment::ObjectProvider<PoissonNEQConvVarSet<Euler3DNEQRhoivt>, 
			    ConvectiveVarSet, PoissonNEQModule, 1> 
poissonNEQ3DRhoivtConvProvider("PoissonNEQ3DRhoivt"); // VS: New Addition

//////////////////////////////////////////////////////////////////////////////

    } // namespace PoissonNEQ

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
