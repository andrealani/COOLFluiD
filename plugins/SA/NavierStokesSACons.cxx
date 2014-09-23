#include "SA.hh"
#include "SA/NavierStokesSACons.hh"
#include "SA/NavierStokesSAVarSet.hh"
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

    namespace SA {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<NavierStokesSACons<NavierStokesSAVarSet<NavierStokes2DVarSet, 0> >, DiffusiveVarSet,
			    SAModule, 2>
ns2DSAConsProvider("NavierStokes2DSACons");

Environment::ObjectProvider<NavierStokesSACons<NavierStokesSAVarSet<NavierStokes3DVarSet, 0> >, DiffusiveVarSet,
			    SAModule, 2>
ns3DSAConsProvider("NavierStokes3DSACons");
      
//////////////////////////////////////////////////////////////////////////////

    } // namespace SA

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
