#include "KOmega/KOmega.hh"
#include "KOmega/EulerKLogOmegaConsVarSet.hh"
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

    namespace KOmega {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<EulerKLogOmegaConsVarSet<MultiScalarVarSet<Euler2DCons>, 0>, 
			    ConvectiveVarSet, KOmegaModule, 1>
euler2DKLogOmegaConsProvider("Euler2DKLogOmegaCons");
    
Environment::ObjectProvider<EulerKLogOmegaConsVarSet<MultiScalarVarSet<Euler3DCons>, 0>, 
			    ConvectiveVarSet, KOmegaModule, 1>
euler3DKLogOmegaConsProvider("Euler3DKLogOmegaCons");
      
//////////////////////////////////////////////////////////////////////////////

    } // namespace KOmega

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
