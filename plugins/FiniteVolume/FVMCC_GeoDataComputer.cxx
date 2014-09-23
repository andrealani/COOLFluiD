#include "Framework/MethodStrategyProvider.hh"
#include "Framework/FVMCCGeoDataComputer.hh"
#include "Framework/NullGeoDataComputer.hh"
#include "FiniteVolume/CellCenterFVMData.hh"
#include "FiniteVolume/FiniteVolume.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

MethodStrategyProvider<NullGeoDataComputer<CellCenterFVMData>,
                       CellCenterFVMData,
                       GeoDataComputer<CellCenterFVMData>,
                       FiniteVolumeModule>
nullGeoDataComputerProv("Null");      

MethodStrategyProvider<FVMCCGeoDataComputer<CellCenterFVMData>,
                       CellCenterFVMData,
                       GeoDataComputer<CellCenterFVMData>,
                       FiniteVolumeModule>
fvmccGeoDataComputerProv("FVMCC");      

//////////////////////////////////////////////////////////////////////////////

} // namespace FiniteVolume

} // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
