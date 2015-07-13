#include "SA/EulerSAPvtToCons.hh"
#include "NavierStokes/EulerPvtToCons.hh"
#include "Environment/ObjectProvider.hh"
#include "SA.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Physics::NavierStokes;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace SA {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<EulerSAPvtToCons<EulerPvtToCons>, VarSetTransformer, SAModule, 1> 
euler1DSAPvtToConsProvider("Euler1DSAPvtToCons");

Environment::ObjectProvider<EulerSAPvtToCons<EulerPvtToCons>, VarSetTransformer, SAModule, 1> 
euler2DSAPuvtToConsProvider("Euler2DSAPuvtToCons");

Environment::ObjectProvider<EulerSAPvtToCons<EulerPvtToCons>, VarSetTransformer, SAModule, 1> 
euler3DSAPvtToConsProvider("Euler3DSAPvtToCons");
      
//////////////////////////////////////////////////////////////////////////////

    } // namespace SA

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
