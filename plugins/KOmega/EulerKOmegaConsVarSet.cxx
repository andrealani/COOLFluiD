#include "KOmega/KOmega.hh"
#include "KOmega/EulerKOmegaConsVarSet.hh"
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

Environment::ObjectProvider<EulerKOmegaConsVarSet<MultiScalarVarSet<Euler2DCons>, 0>, 
			    ConvectiveVarSet, KOmegaModule, 1>
euler2DKOmegaConsProvider("Euler2DKOmegaCons");
    
Environment::ObjectProvider<EulerKOmegaConsVarSet<MultiScalarVarSet<Euler3DCons>, 0>, 
			    ConvectiveVarSet, KOmegaModule, 1>
euler3DKOmegaConsProvider("Euler3DKOmegaCons");
      
//////////////////////////////////////////////////////////////////////////////

    } // namespace KOmega

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
