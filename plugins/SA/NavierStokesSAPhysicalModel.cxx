#include "SA/SA.hh"
#include "NavierStokesSAPhysicalModel.hh"
#include "Environment/ObjectProvider.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Physics::NavierStokes; 

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace SA {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<NavierStokesSAPhysicalModel<DIM_2D>, PhysicalModelImpl, SAModule, 1>
navierStokes2DSAProvider("NavierStokes2DSA");
      
Environment::ObjectProvider<NavierStokesSAPhysicalModel<DIM_3D>, PhysicalModelImpl, SAModule, 1>
navierStokes3DSAProvider("NavierStokes3DSA");
      
//////////////////////////////////////////////////////////////////////////////

    } // namespace SA

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

