#include "NEQ/NEQ.hh"
#include "NEQ/NavierStokesCNEQVarSet.hh"
#include "NEQ/NavierStokesTCNEQVarSet.hh"
#include "NEQ/NavierStokesNEQRhoivt.hh"
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

Environment::ObjectProvider<NavierStokesNEQRhoivt<NavierStokesCNEQVarSet<NavierStokes2DVarSet> >,
			    DiffusiveVarSet,
			    NEQModule, 2>
ns2DNEQRhoivtProvider("NavierStokes2DNEQRhoivt");
    
Environment::ObjectProvider<NavierStokesNEQRhoivt<NavierStokesTCNEQVarSet<NavierStokes2DVarSet> >,
			    DiffusiveVarSet,
			    NEQModule, 2>
ns2DNEQRhoivtTvProvider("NavierStokes2DNEQRhoivtTv");

Environment::ObjectProvider<NavierStokesNEQRhoivt<NavierStokesCNEQVarSet<NavierStokes3DVarSet> >,
			    DiffusiveVarSet,
			    NEQModule, 2>
ns3DNEQRhoivtProvider("NavierStokes3DNEQRhoivt");
      
Environment::ObjectProvider<NavierStokesNEQRhoivt<NavierStokesTCNEQVarSet<NavierStokes3DVarSet> >,
			    DiffusiveVarSet,
			    NEQModule, 2>
ns3DNEQRhoivtTvProvider("NavierStokes3DNEQRhoivtTv");
      
//////////////////////////////////////////////////////////////////////////////
      
    } // namespace NEQ

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
