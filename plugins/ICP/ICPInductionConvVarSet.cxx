#include "ICP/ICP.hh"
#include "ICPInductionConvVarSet.hh"
#include "NavierStokes/Euler2DCons.hh"
#include "LTE/Euler2DPuvtLTE.hh"
#include "Environment/ObjectProvider.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Physics::NavierStokes;
using namespace COOLFluiD::Physics::LTE;
using namespace COOLFluiD::Physics::ICP;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace ICP {

//////////////////////////////////////////////////////////////////////////////

/// Cons VarSet does nothing but is needed for giving variables so it 
/// is parametrized with the non-LTE varset.
Environment::ObjectProvider<ICPInductionConvVarSet<Euler2DCons>, 
			    ConvectiveVarSet, ICPModule, 1> 
icpLTE2DConsConvProvider("ICPLTE2DCons");
    
Environment::ObjectProvider<ICPInductionConvVarSet<Euler2DPuvtLTE>, 
			       ConvectiveVarSet, ICPModule, 1> 
compicpLTE2DPuvtConvProvider("ICPLTE2DPuvt");

//////////////////////////////////////////////////////////////////////////////

    } // namespace ICP

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
