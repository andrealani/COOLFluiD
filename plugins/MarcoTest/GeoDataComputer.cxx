#include "Framework/MethodStrategyProvider.hh"
#include "Framework/FVMCCGeoDataComputer.hh"
#include "Framework/NullGeoDataComputer.hh"
#include "MarcoTest/MarcoTest.hh"
#include "MarcoTest/MarcoTestMethodData.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace MarcoTest {

//////////////////////////////////////////////////////////////////////////////

MethodStrategyProvider<NullGeoDataComputer<MarcoTestMethodData>,
                       MarcoTestMethodData,
                       GeoDataComputer<MarcoTestMethodData>,
                       MarcoTestModule>
nullGeoDataComputerMarcoTestProvider("Null");     

MethodStrategyProvider<FVMCCGeoDataComputer<MarcoTestMethodData>,
                       MarcoTestMethodData,
                       GeoDataComputer<MarcoTestMethodData>,
                       MarcoTestModule>
geoDataComputerMarcoTestProvider("FVMCC");      
    
//////////////////////////////////////////////////////////////////////////////

 } // namespace MarcoTest

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
