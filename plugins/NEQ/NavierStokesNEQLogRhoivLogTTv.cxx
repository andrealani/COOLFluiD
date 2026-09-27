#include "NEQ/NEQ.hh"
#include "NEQ/NavierStokesTCNEQVarSet.hh"
#include "NEQ/NavierStokesNEQLogRhoivLogTTv.hh"
#include "NavierStokes/NavierStokes2DVarSet.hh"
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

Environment::ObjectProvider<NavierStokesNEQLogRhoivLogTTv<NavierStokesTCNEQVarSet<NavierStokes2DVarSet> >,
                            DiffusiveVarSet,
                            NEQModule, 2>
ns2DNEQLogRhoivLogTTvProvider("NavierStokes2DNEQLogRhoivLogTTv");

//////////////////////////////////////////////////////////////////////////////

    } // namespace NEQ

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
