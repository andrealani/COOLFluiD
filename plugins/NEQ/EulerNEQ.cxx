#include "NEQ/NEQ.hh"
#include "NEQ/EulerNEQ.hh"
#include "Environment/ObjectProvider.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace NEQ {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<EulerNEQ<DIM_1D>, PhysicalModelImpl,NEQModule,1>
euler1DNEQProvider("Euler1DNEQ");

Environment::ObjectProvider<EulerNEQ<DIM_2D>, PhysicalModelImpl,NEQModule,1>
euler2DNEQProvider("Euler2DNEQ");

Environment::ObjectProvider<EulerNEQ<DIM_3D>, PhysicalModelImpl,NEQModule,1>
euler3DNEQProvider("Euler3DNEQ");

//////////////////////////////////////////////////////////////////////////////

    } // namespace NEQ

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

