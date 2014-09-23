#include "MultiFluidMHD/MultiFluidMHD.hh"
#include "MultiFluidMHD/MultiFluidMHDModel.hh"
#include "Environment/ObjectProvider.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace MultiFluidMHD {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<MultiFluidMHDModel<DIM_2D>, PhysicalModelImpl,MultiFluidMHDModule, 1>
mfModel2DProvider("MultiFluidMHD2D");

Environment::ObjectProvider<MultiFluidMHDModel<DIM_3D>, PhysicalModelImpl,MultiFluidMHDModule, 1>
mfModel3DProvider("MultiFluidMHD3D");

//////////////////////////////////////////////////////////////////////////////

    } // namespace MultiFluidMHD

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
