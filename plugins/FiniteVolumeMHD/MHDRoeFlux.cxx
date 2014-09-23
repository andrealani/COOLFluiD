#include "MHDRoeFlux.hh"
#include "FiniteVolumeMHD/FiniteVolumeMHD.hh"
#include "MHD/MHD2DCons.hh"
#include "MHD/MHD3DCons.hh"
#include "MHD/MHD2DProjectionCons.hh"
#include "MHD/MHD3DProjectionCons.hh"
#include "Framework/MethodStrategyProvider.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Physics::MHD;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

MethodStrategyProvider<MHDRoeFlux<MHD2DCons>,
		       CellCenterFVMData,
                       FluxSplitter<CellCenterFVMData>,
		       FiniteVolumeMHDModule>
mhd2DConsRoeFluxProvider("MHD2DConsRoe");

MethodStrategyProvider<MHDRoeFlux<MHD3DCons>,
		       CellCenterFVMData,
                       FluxSplitter<CellCenterFVMData>,
		       FiniteVolumeMHDModule>
mhd3DConsRoeFluxProvider("MHD3DConsRoe");

MethodStrategyProvider<MHDRoeFlux<MHD2DProjectionCons>,
		       CellCenterFVMData,
                       FluxSplitter<CellCenterFVMData>,
		       FiniteVolumeMHDModule>
mhd2DProjectionConsRoeFluxProvider("MHD2DProjectionConsRoe");

MethodStrategyProvider<MHDRoeFlux<MHD3DProjectionCons>,
		       CellCenterFVMData,
                       FluxSplitter<CellCenterFVMData>,
		       FiniteVolumeMHDModule>
mhd3DProjectionConsRoeFluxProvider("MHD3DProjectionConsRoe");

//////////////////////////////////////////////////////////////////////////////

} // namespace FiniteVolume

} // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
