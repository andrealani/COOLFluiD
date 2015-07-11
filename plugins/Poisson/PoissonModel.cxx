#include "Poisson/Poisson.hh"
#include "Poisson/PoissonModel.hh"
#include "Environment/ObjectProvider.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace Poisson {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<PoissonModel<DIM_2D>, PhysicalModelImpl,PoissonModule, 1>
mfModel2DProvider("Poisson2D");

Environment::ObjectProvider<PoissonModel<DIM_3D>, PhysicalModelImpl,PoissonModule, 1>
mfModel3DProvider("Poisson3D");

//////////////////////////////////////////////////////////////////////////////

    } // namespace Poisson

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
