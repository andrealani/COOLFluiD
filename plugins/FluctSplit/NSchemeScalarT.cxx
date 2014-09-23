#include "NSchemeScalarT.hh"
#include "Framework/MethodStrategyProvider.hh"
#include "FluctSplit/FluctSplitScalar.hh"
#include "FluctSplit/FluctuationSplitData.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

MethodStrategyProvider<NSchemeScalarT<1>,
		       FluctuationSplitData,
		       Splitter,
                       FluctSplitScalarModule>
nSchemeScalar1Provider("ScalarN1");

MethodStrategyProvider<NSchemeScalarT<4>,
		       FluctuationSplitData,
		       Splitter,
                       FluctSplitScalarModule>
nSchemeScalar4Provider("ScalarN4");

MethodStrategyProvider<NSchemeScalarT<5>,
		       FluctuationSplitData,
		       Splitter,
                       FluctSplitScalarModule>
nSchemeScalar5Provider("ScalarN5");

//////////////////////////////////////////////////////////////////////////////

     } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
