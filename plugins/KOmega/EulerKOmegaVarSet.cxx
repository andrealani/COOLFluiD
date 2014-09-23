#include "KOmega/KOmega.hh"
#include "KOmega/EulerKOmegaVarSet.hh"
#include "Environment/ObjectProvider.hh"
#include "NavierStokes/MultiScalarVarSet.hh"
#include "NavierStokes/Euler2DPuvt.hh"
#include "NavierStokes/Euler3DPvt.hh"
#include "NavierStokes/Euler2DPrim.hh"
#include "NavierStokes/Euler3DPrim.hh"
#include "NavierStokes/Euler2DRhovt.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Physics::NavierStokes;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace KOmega {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<EulerKOmegaVarSet<MultiScalarVarSet<Euler2DPuvt>,0>, 
			    ConvectiveVarSet, KOmegaModule, 1>
euler2DKOmegaPuvtProvider("Euler2DKOmegaPuvt");
     
Environment::ObjectProvider<EulerKOmegaVarSet<MultiScalarVarSet<Euler3DPvt<Euler3DVarSet> >,0>, 
			    ConvectiveVarSet, KOmegaModule, 1>
euler3DKOmegaPvtProvider("Euler3DKOmegaPvt");
      
Environment::ObjectProvider<EulerKOmegaVarSet<MultiScalarVarSet<Euler2DPrim>,0>, 
			    ConvectiveVarSet, KOmegaModule, 1>
euler2DKOmegaPrimProvider("Euler2DKOmegaPrim");
     
Environment::ObjectProvider<EulerKOmegaVarSet<MultiScalarVarSet<Euler3DPrim>, 0>, 
			    ConvectiveVarSet, KOmegaModule, 1>
euler3DKOmegaPrimProvider("Euler3DKOmegaPrim");
      
Environment::ObjectProvider<EulerKOmegaVarSet<MultiScalarVarSet<Euler2DRhovt>, 0>, 
			    ConvectiveVarSet, KOmegaModule, 1>
euler2DKOmegaRhovtProvider("Euler2DKOmegaRhovt");
      
//////////////////////////////////////////////////////////////////////////////

    } // namespace KOmega

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
