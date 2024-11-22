#include "PoissonNEQ/PoissonNEQ.hh"
#include "PoissonNEQ/PoissonNEQDiffVarSet.hh"
#include "NEQ/NavierStokesTCNEQVarSet.hh"
#include "NEQ/NavierStokesCNEQVarSet.hh"
#include "NEQ/NavierStokesNEQRhoivt.hh"
#include "NavierStokes/NavierStokes2DVarSet.hh"
#include "NavierStokes/NavierStokes3DVarSet.hh"
#include "Environment/ObjectProvider.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Physics::NavierStokes;
using namespace COOLFluiD::Physics::NEQ;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace PoissonNEQ {

//////////////////////////////////////////////////////////////////////////////
      
Environment::ObjectProvider<PoissonNEQDiffVarSet
			    <NavierStokesNEQRhoivt<NavierStokesTCNEQVarSet<NavierStokes2DVarSet> > >,
			    DiffusiveVarSet, 
			    PoissonNEQModule, 2> 
poissonNEQ2DRhoivtTvDiffProvider("PoissonNEQ2DRhoivtTv");

Environment::ObjectProvider<PoissonNEQDiffVarSet
			    <NavierStokesNEQRhoivt<NavierStokesTCNEQVarSet<NavierStokes3DVarSet> > >,
			    DiffusiveVarSet, 
			    PoissonNEQModule, 2> 
poissonNEQ3DRhoivtTvDiffProvider("PoissonNEQ3DRhoivtTv");

Environment::ObjectProvider<PoissonNEQDiffVarSet
			    <NavierStokesNEQRhoivt<NavierStokesCNEQVarSet<NavierStokes3DVarSet> > >,
			    DiffusiveVarSet, 
			    PoissonNEQModule, 2> 
poissonNEQ3DRhoivtDiffProvider("PoissonNEQ3DRhoivt"); // VS: new addition 

//////////////////////////////////////////////////////////////////////////////

    } // namespace PoissonNEQ

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////