#include "Framework/MethodCommandProvider.hh"
#include "FluctSplit/FluctSplit.hh"
#include "FluctSplit/WeakFarField.hh"
#include "FluctSplit/WeakBC2D.hh"
#include "FluctSplit/WeakBC3D.hh"
#include "FluctSplit/WeakBC3DHO.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<WeakFarField<WeakBC2D>,
                      FluctuationSplitData,
                      FluctSplitModule>
                      weakFarField2DProvider("WeakFarField2D");

MethodCommandProvider<WeakFarField<WeakBC3D>,
                      FluctuationSplitData,
                      FluctSplitModule>
                      weakFarField3DProvider("WeakFarField3D");

MethodCommandProvider<WeakFarField<WeakBC3DHO>,
                      FluctuationSplitData,
                      FluctSplitModule>
                      weakFarField3DHOProvider("WeakFarField3DHO");

//////////////////////////////////////////////////////////////////////////////

} // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
