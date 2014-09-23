#include "Common/CFLog.hh"
#include "Framework/MethodStrategyProvider.hh"
#include "FluctSplit/FluctSplit.hh"
#include "FluctSplit/FluctuationSplitData.hh"
#include "FluctSplit/NullSplitter.hh"
#include "FluctSplit/SourceTermSplitter.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

MethodStrategyProvider< NullSplitter<Splitter>,
                        FluctuationSplitData,
                        Splitter,
                        FluctSplitModule>
aNullSplitterProvider("Null");

MethodStrategyProvider< NullSplitter<SourceTermSplitter>,
                        FluctuationSplitData,
                        SourceTermSplitter,
                        FluctSplitModule>
aNullSourceTermSplitterScalarProvider("Null");
      
//////////////////////////////////////////////////////////////////////////////

     } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
