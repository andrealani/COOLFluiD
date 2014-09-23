#include "FluctSplit/FluctSplit.hh"
#include "FluctSplit/SuperInletMHD2DProjection.hh"

#include "Framework/MethodCommandProvider.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<SuperInletMHD2DProjection,
                      FluctuationSplitData,
                      FluctSplitModule>
superInletMHD2DProjectionProvider("SuperInletMHD2DProjection");


//////////////////////////////////////////////////////////////////////////////

SuperInletMHD2DProjection::SuperInletMHD2DProjection(const std::string& name) :
  SuperInlet(name)
{
  _nbEquationsToSkip = 1;
}

//////////////////////////////////////////////////////////////////////////////

SuperInletMHD2DProjection::~SuperInletMHD2DProjection()
{
}

//////////////////////////////////////////////////////////////////////////////


    } // namespace FluctSplit



} // namespace COOLFluiD
