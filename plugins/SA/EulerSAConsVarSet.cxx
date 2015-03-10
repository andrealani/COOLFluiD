#include "SA/SA.hh"
#include "SA/EulerSAConsVarSet.hh"
#include "Environment/ObjectProvider.hh"
#include "NavierStokes/MultiScalarVarSet.hh"
#include "NavierStokes/Euler2DCons.hh"
#include "NavierStokes/Euler3DCons.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Physics::NavierStokes;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace SA {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<EulerSAConsVarSet<MultiScalarVarSet<Euler2DCons>, 0>, 
			    ConvectiveVarSet, SAModule, 1>
euler2DSAConsProvider("Euler2DSACons");

Environment::ObjectProvider<EulerSAConsVarSet<MultiScalarVarSet<Euler3DCons>, 0>, 
			    ConvectiveVarSet, SAModule, 1>
euler3DSAConsProvider("Euler3DSACons");
      
//////////////////////////////////////////////////////////////////////////////

    } // namespace SA

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
