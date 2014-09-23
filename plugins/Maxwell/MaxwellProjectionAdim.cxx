#include "Maxwell/Maxwell.hh"
#include "MaxwellProjectionAdim.hh"
#include "Environment/ObjectProvider.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace Maxwell {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<MaxwellProjectionAdim<DIM_2D>, PhysicalModelImpl, MaxwellModule,1>
Maxwell2DProjectionAdimModelProvider("Maxwell2DProjectionAdim");

Environment::ObjectProvider<MaxwellProjectionAdim<DIM_3D>, PhysicalModelImpl, MaxwellModule,1>
Maxwell3DProjectionAdimModelProvider("Maxwell3DProjectionAdim");

//////////////////////////////////////////////////////////////////////////////

} // namespace Maxwell

} // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

