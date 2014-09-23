#include "LESvki/LESvki.hh"
#include "WALESPhysicalModel.hh"
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

Environment::ObjectProvider<WALESPhysicalModel<DIM_2D>, PhysicalModelImpl, LESvkiModule, 1>
walProvider("WALES2D");

Environment::ObjectProvider<WALESPhysicalModel<DIM_3D>, PhysicalModelImpl, LESvkiModule, 1>
wal3DProvider("WALES3D");
//////////////////////////////////////////////////////////////////////////////

    } // namespace NavierStokes

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

