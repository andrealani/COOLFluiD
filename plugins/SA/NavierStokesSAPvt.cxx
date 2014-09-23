#include "SA.hh"
#include "SA/NavierStokesSAPvt.hh"
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

Environment::ObjectProvider<NavierStokesSAPvt<NavierStokesSAVarSet<NavierStokes2DVarSet, 0> >, DiffusiveVarSet,
			    SAModule, 2>
ns2DSAPuvtProvider("NavierStokes2DSAPuvt");

Environment::ObjectProvider<NavierStokesSAPvt<NavierStokesSAVarSet<NavierStokes3DVarSet, 0> >, DiffusiveVarSet,
			    SAModule, 2>
ns3DSAPvtProvider("NavierStokes3DSAPvt");
      
//////////////////////////////////////////////////////////////////////////////

    } // namespace SA

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
