#include "FiniteVolumeMHD/HLLETanakaFlux.hh"
#include "FiniteVolumeMHD/FiniteVolumeMHD.hh"
#include "Framework/MethodStrategyProvider.hh"
#include "MHD/MHD3DProjectionVarSet.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Physics::MHD;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

MethodStrategyProvider<HLLETanakaFlux<MHD3DProjectionVarSet>,
                       CellCenterFVMData,
                       FluxSplitter<CellCenterFVMData>,
                       FiniteVolumeMHDModule>
hlleTanakaMHD3DFluxSplitterProvider("HLLETanakaMHD3D");
      
//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
