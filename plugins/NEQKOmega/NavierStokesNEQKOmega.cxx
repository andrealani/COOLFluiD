#include "NEQKOmega/NEQKOmega.hh"
#include "NEQKOmega/NavierStokesNEQKOmega.hh"
#include "NEQ/NEQReactionTerm.hh"
#include "Environment/ObjectProvider.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace NEQKOmega {

//////////////////////////////////////////////////////////////////////////////
      
Environment::ObjectProvider<NavierStokesNEQKOmega<DIM_2D, NEQ::NEQReactionTerm>, 
			    PhysicalModelImpl, NEQKOmegaModule, 1>
navierStokes2DNEQOmegaProvider("NavierStokes2DNEQKOmega");
      
Environment::ObjectProvider<NavierStokesNEQKOmega<DIM_3D, NEQ::NEQReactionTerm>,
			    PhysicalModelImpl, NEQKOmegaModule, 1>
navierStokes3DNEQOmegaProvider("NavierStokes3DNEQKOmega");
    
//////////////////////////////////////////////////////////////////////////////

    } // namespace NEQKOmega

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
