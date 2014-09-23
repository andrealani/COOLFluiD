#include "MHD/MHD.hh"
#include "MHDProjectionPolytropic.hh"
#include "Environment/ObjectProvider.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace MHD {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<MHDProjectionPolytropic<DIM_2D>, PhysicalModelImpl, MHDModule,1>
mhd2DProjectionPolytropicModelProvider("MHD2DProjectionPolytropic");

Environment::ObjectProvider<MHDProjectionPolytropic<DIM_3D>, PhysicalModelImpl, MHDModule,1>
mhd3DProjectionPolytropicModelProvider("MHD3DProjectionPolytropic");

//////////////////////////////////////////////////////////////////////////////

} // namespace MHD

} // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

