#include "NEQ/NEQ.hh"
#include "NEQ/NavierStokesCNEQVarSet.hh"
#include "NEQ/NavierStokesTCNEQVarSet.hh"
#include "NEQ/NavierStokesNEQPivt.hh"
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

Environment::ObjectProvider<NavierStokesNEQPivt<NavierStokesCNEQVarSet<NavierStokes2DVarSet> >,
			    DiffusiveVarSet,
			    NEQModule, 2>
ns2DNEQPivtProvider("NavierStokes2DNEQPivt");
      
Environment::ObjectProvider<NavierStokesNEQPivt<NavierStokesTCNEQVarSet<NavierStokes2DVarSet> >,
			    DiffusiveVarSet,
			    NEQModule, 2>
ns2DNEQPivtTvProvider("NavierStokes2DNEQPivtTv");

Environment::ObjectProvider<NavierStokesNEQPivt<NavierStokesCNEQVarSet<NavierStokes3DVarSet> >,
			    DiffusiveVarSet,
			    NEQModule, 2>
ns3DNEQPivtProvider("NavierStokes3DNEQPivt");
      
Environment::ObjectProvider<NavierStokesNEQPivt<NavierStokesTCNEQVarSet<NavierStokes3DVarSet> >,
			    DiffusiveVarSet,
			    NEQModule, 2>
ns3DNEQPivtTvProvider("NavierStokes3DNEQPivtTv");
      
//////////////////////////////////////////////////////////////////////////////
      
    } // namespace NEQ

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
