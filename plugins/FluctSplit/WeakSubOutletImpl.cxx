#include "Framework/MethodCommandProvider.hh"

#include "FluctSplit/FluctSplit.hh"
#include "FluctSplit/WeakSubOutletImpl.hh"
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

MethodCommandProvider<WeakSubOutletImpl<WeakBC2DImpl>,
                      FluctuationSplitData,
                      FluctSplitModule>
                      weakSubOutletImpl2DProvider("WeakSubOutlet2DImpl");

MethodCommandProvider<WeakSubOutletImpl<WeakBC3DImpl>,
                      FluctuationSplitData,
                      FluctSplitModule>
                      weakSubOutletImpl3DProvider("WeakSubOutlet3DImpl");

//////////////////////////////////////////////////////////////////////////////

} // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
