#include "SplitStrategyNEQ.hh"
#include "Framework/MethodStrategyProvider.hh"
#include "FluctSplitNEQ/FluctSplitNEQ.hh"
#include "FluctSplit/RD_SplitStrategy.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::FluctSplit;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace FluctSplitNEQ {

//////////////////////////////////////////////////////////////////////////////

// AL: for the moment let's forget this (gets messy because of the update 
//     states that need to be backed up
// MethodStrategyProvider<SplitStrategyNEQ<CRD_SplitStrategy>,
//                        FluctuationSplitData,
//                        FluctuationSplitStrategy,
//                        FluctSplitModule>
// crdNEQFluctSplitStrategyProvider("CRD_NEQ");

MethodStrategyProvider<SplitStrategyNEQ<RD_SplitStrategy>,
                       FluctuationSplitData,
                       FluctuationSplitStrategy,
                       FluctSplitModule>
rdNEQFluctSplitStrategyProvider("RD_NEQ");
    
//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplitNEQ



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
