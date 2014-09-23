#include "LESvki/LESvki.hh"
#include "SmagorinskyPhysicalModel.hh"
#include "Environment/ObjectProvider.hh"

//////////////////////////////////////////////////////////////////////////////
using namespace std;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Physics::NavierStokes;
//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace LESvki {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<SmagorinskyPhysicalModel<DIM_2D>, PhysicalModelImpl, LESvkiModule, 1>
sma2DProvider("Smagorinsky2D");
Environment::ObjectProvider<SmagorinskyPhysicalModel<DIM_3D>, PhysicalModelImpl, LESvkiModule, 1>
sma3DProvider("Smagorinsky3D");

//////////////////////////////////////////////////////////////////////////////

    } // namespace NavierStokes

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

