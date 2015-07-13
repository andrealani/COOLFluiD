#include "SA/SALTE.hh"
#include "SA/EulerSAVarSet.hh"
#include "Environment/ObjectProvider.hh"
#include "NavierStokes/MultiScalarVarSet.hh"
#include "LTE/Euler2DPuvtLTE.hh"
#include "LTE/Euler3DPvtLTE.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Physics::NavierStokes;
using namespace COOLFluiD::Physics::LTE;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace SA {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<EulerSAVarSet<MultiScalarVarSet<Euler2DPuvtLTE>,0>, 
			    ConvectiveVarSet, SALTEModule, 1>
euler2DSAPuvtLTEProvider("Euler2DSAPuvtLTE");
    
Environment::ObjectProvider<EulerSAVarSet<MultiScalarVarSet<Euler3DPvtLTE>,0>, 
			    ConvectiveVarSet, SALTEModule, 1>
euler3DSAPvtLTEProvider("Euler3DSAPvtLTE");

//////////////////////////////////////////////////////////////////////////////

    } // namespace SA

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
