#include "KOmega.hh"
#include "KOmega/NavierStokesKOmegaPvt.hh"
#include "KOmega/NavierStokesKOmegaBSLVarSet.hh"
#include "KOmega/NavierStokesKOmegaKEpsilonVarSet.hh"
#include "KOmega/NavierStokesKOmegaSSTVarSet.hh"
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

Environment::ObjectProvider<NavierStokesKOmegaPvt<NavierStokesKOmegaVarSet<NavierStokes2DVarSet, 0> >, DiffusiveVarSet,
			    KOmegaModule, 2>
ns2DKOmegaPuvtProvider("NavierStokes2DKOmegaPuvt");

Environment::ObjectProvider<NavierStokesKOmegaPvt<NavierStokesKOmegaVarSet<NavierStokes3DVarSet, 0> >, DiffusiveVarSet,
			    KOmegaModule, 2>
ns3DKOmegaPvtProvider("NavierStokes3DKOmegaPuvt"); // change this to Pvt

Environment::ObjectProvider<NavierStokesKOmegaPvt<NavierStokesKOmegaBSLVarSet<NavierStokes2DVarSet, 0> >, DiffusiveVarSet,
			    KOmegaModule, 2>
ns2DKOmegaBSLPuvtProvider("NavierStokes2DKOmegaBSLPuvt");

Environment::ObjectProvider<NavierStokesKOmegaPvt<NavierStokesKOmegaBSLVarSet<NavierStokes3DVarSet, 0> >, DiffusiveVarSet,
			    KOmegaModule, 2>
ns3DKOmegaBSLPvtProvider("NavierStokes3DKOmegaBSLPuvt"); // change this to Pvt
      
Environment::ObjectProvider<NavierStokesKOmegaPvt<NavierStokesKOmegaSSTVarSet<NavierStokes2DVarSet, 0> >, DiffusiveVarSet,
			    KOmegaModule, 2>
ns2DKOmegaSSTPuvtProvider("NavierStokes2DKOmegaSSTPuvt");

Environment::ObjectProvider<NavierStokesKOmegaPvt<NavierStokesKOmegaSSTVarSet<NavierStokes3DVarSet, 0> >, DiffusiveVarSet,
			    KOmegaModule, 2>
ns3DKOmegaSSTPvtProvider("NavierStokes3DKOmegaSSTPuvt"); // change this to Pvt
      
Environment::ObjectProvider<NavierStokesKOmegaPvt<NavierStokesKOmegaKEpsilonVarSet<NavierStokes2DVarSet, 0> >, DiffusiveVarSet,
			    KOmegaModule, 2>
ns2DKOmegaKEpsilonPuvtProvider("NavierStokes2DKOmegaKEpsilonPuvt");

Environment::ObjectProvider<NavierStokesKOmegaPvt<NavierStokesKOmegaKEpsilonVarSet<NavierStokes3DVarSet, 0> >, DiffusiveVarSet,
			    KOmegaModule, 2>
ns3DKOmegaKEpsilonPvtProvider("NavierStokes3DKOmegaKEpsilonPuvt"); // change this to Pvt
      
//////////////////////////////////////////////////////////////////////////////

    } // namespace KOmega

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
