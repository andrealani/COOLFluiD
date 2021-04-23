#include "KOmega.hh"
#include "KOmega/NavierStokesKLogOmegaPvt.hh"
#include "KOmega/NavierStokesKLogOmegaBSLVarSet.hh"
#include "KOmega/NavierStokesKLogOmegaSSTVarSet.hh"
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

    namespace KOmega {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<NavierStokesKLogOmegaPvt<NavierStokesKLogOmegaVarSet<NavierStokes2DVarSet, 0> >, DiffusiveVarSet,
			    KOmegaModule, 2>
ns2DKLogOmegaPuvtProvider("NavierStokes2DKLogOmegaPuvt");

Environment::ObjectProvider<NavierStokesKLogOmegaPvt<NavierStokesKLogOmegaVarSet<NavierStokes3DVarSet, 0> >, DiffusiveVarSet,
			    KOmegaModule, 2>
ns3DKLogOmegaPvtProvider("NavierStokes3DKLogOmegaPvt"); // change this to Pvt

Environment::ObjectProvider<NavierStokesKLogOmegaPvt<NavierStokesKLogOmegaBSLVarSet<NavierStokes2DVarSet, 0> >, DiffusiveVarSet,
			    KOmegaModule, 2>
ns2DKLogOmegaBSLPuvtProvider("NavierStokes2DKLogOmegaBSLPuvt");

Environment::ObjectProvider<NavierStokesKLogOmegaPvt<NavierStokesKLogOmegaBSLVarSet<NavierStokes3DVarSet, 0> >, DiffusiveVarSet,
			    KOmegaModule, 2>
ns3DKLogOmegaBSLPvtProvider("NavierStokes3DKLogOmegaBSLPvt"); // change this to Pvt
      
Environment::ObjectProvider<NavierStokesKLogOmegaPvt<NavierStokesKLogOmegaSSTVarSet<NavierStokes2DVarSet, 0> >, DiffusiveVarSet,
			    KOmegaModule, 2>
ns2DKLogOmegaSSTPuvtProvider("NavierStokes2DKLogOmegaSSTPuvt");

Environment::ObjectProvider<NavierStokesKLogOmegaPvt<NavierStokesKLogOmegaSSTVarSet<NavierStokes3DVarSet, 0> >, DiffusiveVarSet,
			    KOmegaModule, 2>
ns3DKLogOmegaSSTPvtProvider("NavierStokes3DKLogOmegaSSTPvt"); // change this to Pvt
      
//////////////////////////////////////////////////////////////////////////////

    } // namespace KOmega

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
