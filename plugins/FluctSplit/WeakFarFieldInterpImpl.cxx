#include "Framework/MethodCommandProvider.hh"

#include "FluctSplit/FluctSplit.hh"
#include "FluctSplit/WeakFarFieldInterpImpl.hh"
#include "FluctSplit/WeakBC2DImpl.hh"
#include "FluctSplit/WeakBC3DImpl.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<WeakFarFieldInterpImpl<WeakBC2DImpl>,
                      FluctuationSplitData,
                      FluctSplitModule>
                      weakFarFieldInterpImpl2DProvider("WeakFarFieldInterp2DImpl");

MethodCommandProvider<WeakFarFieldInterpImpl<WeakBC3DImpl>,
                      FluctuationSplitData,
                      FluctSplitModule>
                      weakFarFieldInterpImpl3DProvider("WeakFarFieldInterp3DImpl");

//////////////////////////////////////////////////////////////////////////////

} // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
