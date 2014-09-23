#include "NSchemeSysT.hh"
#include "Framework/MethodStrategyProvider.hh"
#include "FluctSplit/FluctSplitSystem.hh"
#include "FluctSplit/FluctuationSplitData.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

MethodStrategyProvider<NSchemeSysT<1>,
		       FluctuationSplitData,
		       Splitter,
		       FluctSplitSystemModule>
nSchemeSys1Provider("SysN1");

MethodStrategyProvider<NSchemeSysT<4>,
		       FluctuationSplitData,
		       Splitter,
                       FluctSplitSystemModule>
nSchemeSys4Provider("SysN4");

MethodStrategyProvider<NSchemeSysT<5>,
		       FluctuationSplitData,
		       Splitter,
                       FluctSplitSystemModule>
nSchemeSys5Provider("SysN5");

MethodStrategyProvider<NSchemeSysT<8>,
		       FluctuationSplitData,
		       Splitter,
                       FluctSplitSystemModule>
nSchemeSys8Provider("SysN8");

//////////////////////////////////////////////////////////////////////////////

     } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
