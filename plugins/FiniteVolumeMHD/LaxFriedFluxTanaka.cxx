#include "LaxFriedFluxTanaka.hh"
#include "FiniteVolumeMHD/FiniteVolumeMHD.hh"
#include "MHD/MHD3DVarSet.hh"
#include "MHD/MHD3DProjectionVarSet.hh"
#include "MHD/MHD3DProjectionPolytropicVarSet.hh"
#include "Framework/MethodStrategyProvider.hh"

//////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Physics::MHD;

//////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////

MethodStrategyProvider<LaxFriedFluxTanaka<MHD3DVarSet>,
		       CellCenterFVMData,
                       FluxSplitter<CellCenterFVMData>,
		       FiniteVolumeMHDModule>
mHD3DConsLaxFriedTanakaFluxProvider("MHD3DConsLaxFriedTanaka");

MethodStrategyProvider<LaxFriedFluxTanaka<MHD3DProjectionVarSet>,
		       CellCenterFVMData,
                       FluxSplitter<CellCenterFVMData>,
		       FiniteVolumeMHDModule>
mHD3DProjectionConsLaxFriedTanakaFluxProvider("MHD3DProjectionConsLaxFriedTanaka");

MethodStrategyProvider<LaxFriedFluxTanaka<MHD3DProjectionPolytropicVarSet>,
                       CellCenterFVMData,
                       FluxSplitter<CellCenterFVMData>,
                       FiniteVolumeMHDModule>
mHD3DProjectionPolytropicConsLaxFriedTanakaFluxProvider("MHD3DProjectionPolytropicConsLaxFriedTanaka");

//////////////////////////////////////////////////////////////////////

} // namespace FiniteVolume

} // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////
