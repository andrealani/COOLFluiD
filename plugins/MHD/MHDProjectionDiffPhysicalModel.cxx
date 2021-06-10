#include "MHD/MHD.hh"
#include "MHD/MHDProjectionDiffPhysicalModel.hh"
#include "Environment/ObjectProvider.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace MHD {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<MHDProjectionDiffPhysicalModel<DIM_2D>, PhysicalModelImpl, MHDModule, 1>
mhd2DProjectionDiffProvider("MHD2DProjectionDiff");

Environment::ObjectProvider<MHDProjectionDiffPhysicalModel<DIM_3D>, PhysicalModelImpl, MHDModule, 1>
mhd3DProjectionDiffProvider("MHD3DProjectionDiff");

//////////////////////////////////////////////////////////////////////////////

    } // namespace MHD

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

