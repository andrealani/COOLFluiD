#include "ICP/ICPNEQ.hh"
#include "ICP/ICPInductionDiffVarSet.hh"
#include "ICP/ICPReactionTerm.hh"
#include "NEQ/NavierStokesCNEQVarSet.hh"
#include "NEQ/NavierStokesTCNEQVarSet.hh"
#include "NEQ/NavierStokesNEQPivt.hh"
#include "NEQ/NavierStokesNEQRhoivt.hh"
#include "NEQ/NEQReactionTerm.hh"
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

    namespace ICP {

//////////////////////////////////////////////////////////////////////////////
      
Environment::ObjectProvider<ICPInductionDiffVarSet
			    <NavierStokesNEQPivt<NavierStokesCNEQVarSet<NavierStokes2DVarSet> >,
			     ICPReactionTerm<NEQ::NEQReactionTerm> >, 
			    DiffusiveVarSet, 
			    ICPNEQModule, 2> 
icpNEQ2DPivtDiffProvider("ICPNEQ2DPivt");

Environment::ObjectProvider<ICPInductionDiffVarSet
			    <NavierStokesNEQPivt<NavierStokesTCNEQVarSet<NavierStokes2DVarSet> >,
			     ICPReactionTerm<NEQ::NEQReactionTerm> >, 
			    DiffusiveVarSet, 
			    ICPNEQModule, 2> 
icpNEQ2DPivtTvDiffProvider("ICPNEQ2DPivtTv");

Environment::ObjectProvider<ICPInductionDiffVarSet
			    <NavierStokesNEQRhoivt<NavierStokesCNEQVarSet<NavierStokes2DVarSet> >,
			     ICPReactionTerm<NEQ::NEQReactionTerm> >, 
			    DiffusiveVarSet, 
			    ICPNEQModule, 2> 
icpNEQ2DRhoivtDiffProvider("ICPNEQ2DRhoivt");

Environment::ObjectProvider<ICPInductionDiffVarSet
			    <NavierStokesNEQRhoivt<NavierStokesTCNEQVarSet<NavierStokes2DVarSet> >,
			     ICPReactionTerm<NEQ::NEQReactionTerm> >, 
			    DiffusiveVarSet, 
			    ICPNEQModule, 2> 
icpNEQ2DRhoivtTvDiffProvider("ICPNEQ2DRhoivtTv");

//////////////////////////////////////////////////////////////////////////////

    } // namespace ICP

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
