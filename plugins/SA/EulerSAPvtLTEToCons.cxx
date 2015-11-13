#include "SA/EulerSAPvtToCons.hh"
#include "LTE/EulerPvtLTEToCons.hh"
#include "Environment/ObjectProvider.hh"
#include "SALTE.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Physics::LTE;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace SA {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<EulerSAPvtToCons<EulerPvtLTEToCons>, VarSetTransformer, SALTEModule, 1> 
euler2DSAPuvtLTEToConsProvider("Euler2DSAPuvtLTEToCons");

Environment::ObjectProvider<EulerSAPvtToCons<EulerPvtLTEToCons>, VarSetTransformer, SALTEModule, 1> 
euler3DSAPvtLTEToConsProvider("Euler3DSAPvtLTEToCons");

//////////////////////////////////////////////////////////////////////////////

    } // namespace SA

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
