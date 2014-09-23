#include "ICP/ICP.hh"
#include "ICPInductionDiffVarSet.hh"
#include "ICP/ICPReactionTerm.hh"
#include "LTE/NavierStokes2DPuvtLTE.hh"
#include "NavierStokes/EulerTerm.hh"
#include "Environment/ObjectProvider.hh"
#include "Framework/MultiScalarTerm.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Physics::NavierStokes;
using namespace COOLFluiD::Physics::LTE;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace ICP {

//////////////////////////////////////////////////////////////////////////////
      
Environment::ObjectProvider<ICPInductionDiffVarSet
			    <NavierStokes2DPuvtLTE<MultiScalarTerm<EulerTerm> >, 
			     ICPReactionTerm<Framework::BaseTerm> >, 
			    DiffusiveVarSet, 
			    ICPModule, 2> 
icpLTE2DPuvtDiffProvider("ICPLTE2DPuvt");

//////////////////////////////////////////////////////////////////////////////

    } // namespace ICP

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
