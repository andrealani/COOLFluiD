#include "NEQKOmega/NEQKOmega.hh"
#include "KOmega/EulerKOmegaVarSet.hh"
#include "Environment/ObjectProvider.hh"
#include "NEQ/Euler2DNEQRhoivt.hh"
#include "NEQ/Euler3DNEQRhoivt.hh"
#include "NEQ/Euler2DNEQRhoivtTv.hh"
#include "NEQ/Euler3DNEQRhoivtTv.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Physics::NavierStokes;
using namespace COOLFluiD::Physics::NEQ;
using namespace COOLFluiD::Physics::KOmega;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace NEQKOmega {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<EulerKOmegaVarSet<Euler2DNEQRhoivt,2>, 
			    ConvectiveVarSet, NEQKOmegaModule, 1>
euler2DNEQKOmegaRhoivtProvider("Euler2DNEQKOmegaRhoivt");
   
Environment::ObjectProvider<EulerKOmegaVarSet<Euler3DNEQRhoivt,2>, 
			    ConvectiveVarSet, NEQKOmegaModule, 1>
euler3DNEQKOmegaRhoivtProvider("Euler3DNEQKOmegaRhoivt");

Environment::ObjectProvider<EulerKOmegaVarSet<Euler2DNEQRhoivtTv,3>, 
			    ConvectiveVarSet, NEQKOmegaModule, 1>
euler2DNEQKOmegaRhoivtTvProvider("Euler2DNEQKOmegaRhoivtTv");
   
Environment::ObjectProvider<EulerKOmegaVarSet<Euler3DNEQRhoivtTv,3>, 
			    ConvectiveVarSet, NEQKOmegaModule, 1>
euler3DNEQKOmegaRhoivtTvProvider("Euler3DNEQKOmegaRhoivtTv");
      
//////////////////////////////////////////////////////////////////////////////

    } // namespace NEQKOmega

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
