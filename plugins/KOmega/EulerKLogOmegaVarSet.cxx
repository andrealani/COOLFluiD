#include "KOmega/KOmega.hh"
#include "KOmega/EulerKLogOmegaVarSet.hh"
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

Environment::ObjectProvider<EulerKLogOmegaVarSet<MultiScalarVarSet<Euler2DPuvt>,0>, 
			    ConvectiveVarSet, KOmegaModule, 1>
euler2DKLogOmegaPuvtProvider("Euler2DKLogOmegaPuvt");
     
Environment::ObjectProvider<EulerKLogOmegaVarSet<MultiScalarVarSet<Euler3DPvt<Euler3DVarSet> >,0>, 
			    ConvectiveVarSet, KOmegaModule, 1>
euler3DKLogOmegaPvtProvider("Euler3DKLogOmegaPvt");
      
Environment::ObjectProvider<EulerKLogOmegaVarSet<MultiScalarVarSet<Euler2DPrim>,0>, 
			    ConvectiveVarSet, KOmegaModule, 1>
euler2DKLogOmegaPrimProvider("Euler2DKLogOmegaPrim");
     
Environment::ObjectProvider<EulerKLogOmegaVarSet<MultiScalarVarSet<Euler3DPrim>, 0>, 
			    ConvectiveVarSet, KOmegaModule, 1>
euler3DKLogOmegaPrimProvider("Euler3DKLogOmegaPrim");
      
Environment::ObjectProvider<EulerKLogOmegaVarSet<MultiScalarVarSet<Euler2DRhovt>, 0>, 
			    ConvectiveVarSet, KOmegaModule, 1>
euler2DKLogOmegaRhovtProvider("Euler2DKLogOmegaRhovt");
      
//////////////////////////////////////////////////////////////////////////////

    } // namespace KOmega

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
