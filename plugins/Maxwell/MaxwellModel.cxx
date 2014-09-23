#include "Maxwell/Maxwell.hh"
#include "Maxwell/MaxwellModel.hh"
#include "Environment/ObjectProvider.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace Maxwell {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<MaxwellModel<DIM_2D>, PhysicalModelImpl,MaxwellModule, 1>  
maxwellModel2DProvider("Maxwell2D");

Environment::ObjectProvider<MaxwellModel<DIM_3D>, PhysicalModelImpl,MaxwellModule, 1>
maxwellModel3DProvider("Maxwell3D");

//////////////////////////////////////////////////////////////////////////////

    } // namespace Maxwell

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
