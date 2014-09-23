#include "FluctSplit/FluctSplit.hh"
#include "SuperOutletMHD2DProjection.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/MeshData.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<SuperOutletMHD2DProjection, FluctuationSplitData, FluctSplitModule> superOutletMHD2DProjectionProvider("SuperOutletMHD2DProjection");

//////////////////////////////////////////////////////////////////////////////

void SuperOutletMHD2DProjection::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOption< CFreal >("refPhi","Reference phi value imposed at the supersonic outlet.");
}

//////////////////////////////////////////////////////////////////////////////

SuperOutletMHD2DProjection::SuperOutletMHD2DProjection(const std::string& name) :
  FluctuationSplitCom(name),
  socket_rhs("rhs"),
  socket_isUpdated("isUpdated"),
  socket_states("states")
{
   addConfigOptionsTo(this);
  _refPhi = 0.;
   setParameter("refPhi",&_refPhi);
}

//////////////////////////////////////////////////////////////////////////////

SuperOutletMHD2DProjection::~SuperOutletMHD2DProjection()
{
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
SuperOutletMHD2DProjection::needsSockets()
{

  std::vector<Common::SafePtr<BaseDataSocketSink> > result;

  result.push_back(&socket_rhs);
  result.push_back(&socket_isUpdated);
  result.push_back(&socket_states);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

void SuperOutletMHD2DProjection::executeOnTrs()
{
SafePtr<TopologicalRegionSet> trs = getCurrentTRS();
  CFLogDebugMax( "SuperOutletMHD2DProjection::execute() called for TRS: "
  << trs->getName() << "\n");

  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();

  DataHandle<bool> isUpdated = socket_isUpdated.getDataHandle();

  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();

  DataHandle< CFreal> rhs = socket_rhs.getDataHandle();

  Common::SafePtr< vector<CFuint> > const statesIdx = trs->getStatesInTrs();

  for (CFuint iState = 0; iState < statesIdx->size(); ++iState) {
    const CFuint stateID = (*statesIdx)[iState];
    State& state = *(states[stateID]);
    if (!isUpdated[stateID]) {
        state[8] = _refPhi;
        rhs(stateID, 8, nbEqs) = 0.;
      }
    }
}

//////////////////////////////////////////////////////////////////////////////

void SuperOutletMHD2DProjection::configure ( Config::ConfigArgs& args )
{
  FluctuationSplitCom::configure(args);
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD
