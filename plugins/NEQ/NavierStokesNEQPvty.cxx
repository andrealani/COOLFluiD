#include "NEQ/NEQ.hh"
#include "NEQ/NavierStokesCNEQVarSet.hh"
#include "NEQ/NavierStokesNEQPvty.hh"
#include "NavierStokes/NavierStokes2DVarSet.hh"
#include "NavierStokes/NavierStokes3DVarSet.hh"
#include "Environment/ObjectProvider.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Physics::NavierStokes;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace NEQ {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<NavierStokesNEQPvty<NavierStokesCNEQVarSet<NavierStokes2DVarSet> >,
			    DiffusiveVarSet,
			    NEQModule, 2>
ns2DNEQPvtyProvider("NavierStokes2DNEQPvty");
 

// to be removed
Environment::ObjectProvider<NavierStokesNEQPvty<NavierStokesCNEQVarSet<NavierStokes2DVarSet> >,
			    DiffusiveVarSet,
			    NEQModule, 2>
ns2DNEQPuvtyProvider("NavierStokes2DNEQPuvty");
        
Environment::ObjectProvider<NavierStokesNEQPvty<NavierStokesCNEQVarSet<NavierStokes3DVarSet> >,
			    DiffusiveVarSet,
			    NEQModule, 2>
ns3DNEQPvtyProvider("NavierStokes3DNEQPvty");
      
//////////////////////////////////////////////////////////////////////////////
      
    } // namespace NEQ

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
