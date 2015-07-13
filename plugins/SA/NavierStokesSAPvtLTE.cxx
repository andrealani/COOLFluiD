#include "SALTE.hh"
#include "SA/NavierStokesSAPvtLTE.hh"
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

Environment::ObjectProvider<NavierStokesSAPvtLTE<NavierStokesSAVarSet<NavierStokes2DVarSet, 0> >,
			    DiffusiveVarSet, SALTEModule, 2>
ns2DSAPuvtLTEProvider("NavierStokes2DSAPuvtLTE");

Environment::ObjectProvider<NavierStokesSAPvtLTE<NavierStokesSAVarSet<NavierStokes3DVarSet, 0> >, 
			    DiffusiveVarSet, SALTEModule, 2>
ns3DSAPvtLTEProvider("NavierStokes3DSAPvtLTE");
      
//////////////////////////////////////////////////////////////////////////////

    } // namespace SA

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
