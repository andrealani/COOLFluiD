#include "NEQKOmega/NEQKOmega.hh"
#include "NEQKOmega/EulerNEQKOmegaConsVarSet.hh"
#include "Environment/ObjectProvider.hh"
#include "NEQ/Euler2DNEQCons.hh"
#include "NEQ/Euler3DNEQCons.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Physics::NavierStokes;
using namespace COOLFluiD::Physics::NEQ;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace NEQKOmega {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<EulerNEQKOmegaConsVarSet<Euler2DNEQCons,1>, 
			    ConvectiveVarSet, NEQKOmegaModule, 1>
euler2DNEQKOmegaConsProvider("Euler2DNEQKOmegaCons");
   
Environment::ObjectProvider<EulerNEQKOmegaConsVarSet<Euler3DNEQCons,1>, 
			    ConvectiveVarSet, NEQKOmegaModule, 1>
euler3DNEQKOmegaConsProvider("Euler3DNEQKOmegaCons");
      
//////////////////////////////////////////////////////////////////////////////

    } // namespace NEQKOmega

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
