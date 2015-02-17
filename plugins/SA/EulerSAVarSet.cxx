#include "SA/SA.hh"
#include "SA/EulerSAVarSet.hh"
#include "Environment/ObjectProvider.hh"
#include "NavierStokes/MultiScalarVarSet.hh"
#include "NavierStokes/Euler2DPuvt.hh"
#include "NavierStokes/Euler2DPrim.hh"
#include "NavierStokes/Euler2DRhovt.hh"
#include "NavierStokes/Euler3DPvt.hh"
#include "NavierStokes/Euler3DPrim.hh"
#include "NavierStokes/Euler3DVarSet.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Physics::NavierStokes;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace SA {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<EulerSAVarSet<MultiScalarVarSet<Euler2DPuvt>,0>, 
			    ConvectiveVarSet, SAModule, 1>
euler2DSAPuvtProvider("Euler2DSAPuvt");
           
Environment::ObjectProvider<EulerSAVarSet<MultiScalarVarSet<Euler2DPrim>,0>, 
			    ConvectiveVarSet, SAModule, 1>
euler2DSAPrimProvider("Euler2DSAPrim");
     
Environment::ObjectProvider<EulerSAVarSet<MultiScalarVarSet<Euler2DRhovt>, 0>, 
			    ConvectiveVarSet, SAModule, 1>
euler2DSARhovtProvider("Euler2DSARhovt");

Environment::ObjectProvider<EulerSAVarSet<MultiScalarVarSet<Euler3DPvt<Euler3DVarSet> >,0>,
                            ConvectiveVarSet, SAModule, 1>
euler3DSAPvtProvider("Euler3DSAPvt");

Environment::ObjectProvider<EulerSAVarSet<MultiScalarVarSet<Euler3DPrim>,0>,
                            ConvectiveVarSet, SAModule, 1>
euler3DSAPrimProvider("Euler3DSAPrim");

//////////////////////////////////////////////////////////////////////////////

    } // namespace SA

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
