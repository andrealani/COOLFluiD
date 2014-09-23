#include "NEQKOmega.hh"
#include "NEQKOmega/NavierStokesNEQKOmegaRhoivt.hh"
#include "KOmega/NavierStokesKOmegaBSLVarSet.hh"
#include "KOmega/NavierStokesKOmegaKEpsilonVarSet.hh"
#include "KOmega/NavierStokesKOmegaSSTVarSet.hh"
#include "NEQ/NavierStokesCNEQVarSet.hh"
#include "NEQ/NavierStokesTCNEQVarSet.hh"
#include "NavierStokes/NavierStokes2DVarSet.hh"
#include "NavierStokes/NavierStokes3DVarSet.hh"
#include "Environment/ObjectProvider.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Physics::NavierStokes;
using namespace COOLFluiD::Physics::KOmega;
using namespace COOLFluiD::Physics::NEQ;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace NEQKOmega {

//////////////////////////////////////////////////////////////////////////////

// chemical NEQ
Environment::ObjectProvider<NavierStokesNEQKOmegaRhoivt
			    <NavierStokesKOmegaVarSet<NavierStokesCNEQVarSet<NavierStokes2DVarSet>, 2> >, 
			    DiffusiveVarSet,
			    NEQKOmegaModule, 2>
ns2DNEQKOmegaRhoivtProvider("NavierStokes2DNEQKOmegaRhoivt");
      
Environment::ObjectProvider<NavierStokesNEQKOmegaRhoivt
			    <NavierStokesKOmegaVarSet<NavierStokesCNEQVarSet<NavierStokes3DVarSet>, 2> >, 
			    DiffusiveVarSet,
			    NEQKOmegaModule, 2>
ns3DNEQNEQKOmegaRhoivtProvider("NavierStokes3DNEQKOmegaRhoivt");

Environment::ObjectProvider<NavierStokesNEQKOmegaRhoivt
			    <NavierStokesKOmegaBSLVarSet<NavierStokesCNEQVarSet<NavierStokes2DVarSet>, 2> >, DiffusiveVarSet,
			    NEQKOmegaModule, 2>
ns2DNEQKOmegaBSLRhoivtProvider("NavierStokes2DNEQKOmegaBSLRhoivt");

Environment::ObjectProvider<NavierStokesNEQKOmegaRhoivt
			    <NavierStokesKOmegaBSLVarSet<NavierStokesCNEQVarSet<NavierStokes3DVarSet>, 2> >, DiffusiveVarSet,
			    NEQKOmegaModule, 2>
ns3DNEQKOmegaBSLRhoivtProvider("NavierStokes3DNEQKOmegaBSLRhoivt");
      
Environment::ObjectProvider<NavierStokesNEQKOmegaRhoivt
			    <NavierStokesKOmegaSSTVarSet<NavierStokesCNEQVarSet<NavierStokes2DVarSet>, 2> >, DiffusiveVarSet,
			    NEQKOmegaModule, 2>
ns2DNEQKOmegaSSTRhoivtProvider("NavierStokes2DNEQKOmegaSSTRhoivt");

Environment::ObjectProvider<NavierStokesNEQKOmegaRhoivt
			    <NavierStokesKOmegaSSTVarSet<NavierStokesCNEQVarSet<NavierStokes3DVarSet>, 2> >, DiffusiveVarSet,
			    NEQKOmegaModule, 2>
ns3DNEQKOmegaSSTRhoivtProvider("NavierStokes3DNEQKOmegaSSTRhoivt");
      
Environment::ObjectProvider<NavierStokesNEQKOmegaRhoivt
			    <NavierStokesKOmegaKEpsilonVarSet<NavierStokesCNEQVarSet<NavierStokes2DVarSet>, 2> >, DiffusiveVarSet,
			    NEQKOmegaModule, 2>
ns2DNEQKOmegaKEpsilonRhoivtProvider("NavierStokes2DNEQKOmegaKEpsilonRhoivt");
      
Environment::ObjectProvider<NavierStokesNEQKOmegaRhoivt
			    <NavierStokesKOmegaKEpsilonVarSet<NavierStokesCNEQVarSet<NavierStokes3DVarSet>, 2> >, DiffusiveVarSet,
			    NEQKOmegaModule, 2>
ns3DNEQKOmegaKEpsilonRhoivtProvider("NavierStokes3DNEQKOmegaKEpsilonRhoivt"); 

// thermo-chemical NEQ
Environment::ObjectProvider<NavierStokesNEQKOmegaRhoivt
			    <NavierStokesKOmegaVarSet<NavierStokesCNEQVarSet<NavierStokes2DVarSet>, 3> >, 
			    DiffusiveVarSet,
			    NEQKOmegaModule, 2>
ns2DNEQKOmegaRhoivtTvProvider("NavierStokes2DNEQKOmegaRhoivtTv");

Environment::ObjectProvider<NavierStokesNEQKOmegaRhoivt
			    <NavierStokesKOmegaVarSet<NavierStokesCNEQVarSet<NavierStokes3DVarSet>, 3> >, 
			    DiffusiveVarSet,
			    NEQKOmegaModule, 2>
ns3DNEQNEQKOmegaRhoivtTvProvider("NavierStokes3DNEQKOmegaRhoivtTv");

Environment::ObjectProvider<NavierStokesNEQKOmegaRhoivt
			    <NavierStokesKOmegaBSLVarSet<NavierStokesCNEQVarSet<NavierStokes2DVarSet>, 3> >, DiffusiveVarSet,
			    NEQKOmegaModule, 2>
ns2DNEQKOmegaBSLRhoivtTvProvider("NavierStokes2DNEQKOmegaBSLRhoivtTv");

Environment::ObjectProvider<NavierStokesNEQKOmegaRhoivt
			    <NavierStokesKOmegaBSLVarSet<NavierStokesCNEQVarSet<NavierStokes3DVarSet>, 3> >, DiffusiveVarSet,
			    NEQKOmegaModule, 2>
ns3DNEQKOmegaBSLRhoivtTvProvider("NavierStokes3DNEQKOmegaBSLRhoivtTv");
      
Environment::ObjectProvider<NavierStokesNEQKOmegaRhoivt
			    <NavierStokesKOmegaSSTVarSet<NavierStokesCNEQVarSet<NavierStokes2DVarSet>, 3> >, DiffusiveVarSet,
			    NEQKOmegaModule, 2>
ns2DNEQKOmegaSSTRhoivtTvProvider("NavierStokes2DNEQKOmegaSSTRhoivtTv");

Environment::ObjectProvider<NavierStokesNEQKOmegaRhoivt
			    <NavierStokesKOmegaSSTVarSet<NavierStokesCNEQVarSet<NavierStokes3DVarSet>, 3> >, DiffusiveVarSet,
			    NEQKOmegaModule, 2>
ns3DNEQKOmegaSSTRhoivtTvProvider("NavierStokes3DNEQKOmegaSSTRhoivtTv");
      
Environment::ObjectProvider<NavierStokesNEQKOmegaRhoivt
			    <NavierStokesKOmegaKEpsilonVarSet<NavierStokesCNEQVarSet<NavierStokes2DVarSet>, 3> >, DiffusiveVarSet,
			    NEQKOmegaModule, 2>
ns2DNEQKOmegaKEpsilonRhoivtTvProvider("NavierStokes2DNEQKOmegaKEpsilonRhoivtTv");
      
Environment::ObjectProvider<NavierStokesNEQKOmegaRhoivt
			    <NavierStokesKOmegaKEpsilonVarSet<NavierStokesCNEQVarSet<NavierStokes3DVarSet>, 3> >, DiffusiveVarSet,
			    NEQKOmegaModule, 2>
ns3DNEQKOmegaKEpsilonRhoivtTvProvider("NavierStokes3DNEQKOmegaKEpsilonRhoivtTv"); 
      
//////////////////////////////////////////////////////////////////////////////

    } // namespace NEQKOmega

  } // namespace Physics

} // namespace COOLFluiD+

//////////////////////////////////////////////////////////////////////////////
