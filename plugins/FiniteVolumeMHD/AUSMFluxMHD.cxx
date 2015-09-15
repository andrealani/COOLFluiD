#include "FiniteVolumeMHD/AUSMPWFluxMHD.hh"
#include "FiniteVolumeMHD/FiniteVolumeMHD.hh"
#include "MHD/MHD2DProjectionVarSet.hh"
#include "MHD/MHD3DProjectionVarSet.hh"
#include "Maxwell/Maxwell2DProjectionVarSet.hh"
#include "Maxwell/Maxwell3DProjectionVarSet.hh"
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

/// AUSM PW

MethodStrategyProvider<AUSMPWFluxMHD<MHD2DProjectionVarSet>,
		       CellCenterFVMData,
                       FluxSplitter<CellCenterFVMData>,
               FiniteVolumeMHDModule>
ausmPWMHD2dProvider("AUSMPWMHD2D");

MethodStrategyProvider<AUSMPWFluxMHD<MHD3DProjectionVarSet>,
		       CellCenterFVMData,
                       FluxSplitter<CellCenterFVMData>,
               FiniteVolumeMHDModule>
ausmPWMHD3dProvider("AUSMPWMHD3D");


//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
