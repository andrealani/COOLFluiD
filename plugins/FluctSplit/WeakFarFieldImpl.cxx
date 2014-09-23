#include "Framework/MethodCommandProvider.hh"

#include "FluctSplit/FluctSplit.hh"
#include "FluctSplit/WeakFarFieldImpl.hh"
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

MethodCommandProvider<WeakFarFieldImpl<WeakBC2DImpl>,
                      FluctuationSplitData,
                      FluctSplitModule>
                      weakFarFieldImpl2DProvider("WeakFarField2DImpl");

MethodCommandProvider<WeakFarFieldImpl<WeakBC3DImpl>,
                      FluctuationSplitData,
                      FluctSplitModule>
                      weakFarFieldImpl3DProvider("WeakFarField3DImpl");

//////////////////////////////////////////////////////////////////////////////

} // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
