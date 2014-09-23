#include "LinEuler/LinearizedEuler.hh"
#include "LinEulerPhysicalModel.hh"
#include "Environment/ObjectProvider.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace LinearizedEuler {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<LinEulerPhysicalModel<DIM_2D>, PhysicalModelImpl,LinearizedEulerModule, 1>
lineuler2DModelProvider("LinEuler2D");

Environment::ObjectProvider<LinEulerPhysicalModel<DIM_3D>, PhysicalModelImpl,LinearizedEulerModule, 1>
lineuler3DModelProvider("LinEuler3D");

//////////////////////////////////////////////////////////////////////////////

    } // namespace LinearizedEuler

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

