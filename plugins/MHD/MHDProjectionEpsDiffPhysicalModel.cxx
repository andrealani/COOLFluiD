#include "MHD/MHD.hh"
#include "MHD/MHDProjectionEpsDiffPhysicalModel.hh"
#include "Environment/ObjectProvider.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace MHD {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<MHDProjectionEpsDiffPhysicalModel<DIM_2D>, PhysicalModelImpl, MHDModule, 1>
mhd2DProjectionEpsDiffProvider("MHD2DProjectionEpsDiff");

Environment::ObjectProvider<MHDProjectionEpsDiffPhysicalModel<DIM_3D>, PhysicalModelImpl, MHDModule, 1>
mhd3DProjectionEpsDiffProvider("MHD3DProjectionEpsDiff");

//////////////////////////////////////////////////////////////////////////////

    } // namespace MHD

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

