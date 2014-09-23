#include "FluctSplit/FluctSplit.hh"
#include "SuperOutletMHD3DProjection.hh"
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

MethodCommandProvider<SuperOutletMHD3DProjection, FluctuationSplitData, FluctSplitModule> superOutletMHD3DProjectionProvider("SuperOutletMHD3DProjection");

//////////////////////////////////////////////////////////////////////////////

void SuperOutletMHD3DProjection::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOption< CFreal >("refPhi","Reference phi value imposed at the supersonic outlet.");
}

//////////////////////////////////////////////////////////////////////////////

SuperOutletMHD3DProjection::SuperOutletMHD3DProjection(const std::string& name) :
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

SuperOutletMHD3DProjection::~SuperOutletMHD3DProjection()
{
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
SuperOutletMHD3DProjection::needsSockets()
{

  std::vector<Common::SafePtr<BaseDataSocketSink> > result;

  result.push_back(&socket_rhs);
  result.push_back(&socket_isUpdated);
  result.push_back(&socket_states);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

void SuperOutletMHD3DProjection::executeOnTrs()
{
SafePtr<TopologicalRegionSet> trs = getCurrentTRS();
  CFLogDebugMax( "SuperOutletMHD3DProjection::execute() called for TRS: "
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

void SuperOutletMHD3DProjection::configure ( Config::ConfigArgs& args )
{
  FluctuationSplitCom::configure(args);
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD
