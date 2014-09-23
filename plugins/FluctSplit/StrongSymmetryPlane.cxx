#include "FluctSplit/FluctSplit.hh"
#include "StrongSymmetryPlane.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/MeshData.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<StrongSymmetryPlane,
		      FluctuationSplitData,
		      FluctSplitModule>
strongSymmetryPlaneProvider("StrongSymmetryPlane");

//////////////////////////////////////////////////////////////////////////////

void StrongSymmetryPlane::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOption< vector<CFuint> >
     ("annullVarID","Tells which variables should be annulled.");
}

//////////////////////////////////////////////////////////////////////////////

StrongSymmetryPlane::StrongSymmetryPlane(const std::string& name) :
  FluctuationSplitCom(name),
  socket_rhs("rhs"),
  socket_states("states"),
  socket_isUpdated("isUpdated"),
  _useForInitialization(false)
{
  addConfigOptionsTo(this);

  _varIDToAnnull = vector<CFuint>();
  setParameter("annullVarID",&_varIDToAnnull);
}

//////////////////////////////////////////////////////////////////////////////

StrongSymmetryPlane::~StrongSymmetryPlane()
{
}

//////////////////////////////////////////////////////////////////////////////

void StrongSymmetryPlane::setup()
{

  FluctuationSplitCom::setup();

  // fix for restart
  // apply the initial condition
  vector<Common::SafePtr<TopologicalRegionSet> >& trsList = getTrsList();
  const CFuint nbAnnullVars = _varIDToAnnull.size();
  cf_assert(nbAnnullVars > 0);

  for(CFuint iTRS = 0; iTRS < trsList.size(); ++iTRS) {
    SafePtr<TopologicalRegionSet> trs = trsList[iTRS];

    CFLogDebugMax("StrongSymmetryPlane::setup() called for TRS: "
		  << trs->getName() << "\n");

    DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();
    SafePtr< vector<CFuint> > const statesIdx = trs->getStatesInTrs();

    for (CFuint iState = 0; iState < statesIdx->size(); ++iState) {
      const CFuint localStateID = (*statesIdx)[iState];
      State *const currState = states[localStateID];
      for (CFuint iv = 0; iv < nbAnnullVars; ++iv) {
	(*currState)[_varIDToAnnull[iv]] = 0.0;
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
StrongSymmetryPlane::needsSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSink> > result;

  result.push_back(&socket_rhs);
  result.push_back(&socket_states);
  result.push_back(&socket_isUpdated);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

void StrongSymmetryPlane::executeOnTrs()
{
  SafePtr<TopologicalRegionSet> trs = getCurrentTRS();

  CFLogDebugMax
    ( "StrongNoSlipWallAdiabaticNS2DImpl::execute() called for TRS: "
      << trs->getName() << "\n");

  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();
  DataHandle<bool> isUpdated = socket_isUpdated.getDataHandle();
  DataHandle<CFreal> rhs = socket_rhs.getDataHandle();

  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  SafePtr< vector<CFuint> > const statesIdx = trs->getStatesInTrs();
  const CFuint nbAnnullVars = _varIDToAnnull.size();
  _useForInitialization = getMethodData().isInitializationPhase();

  for (CFuint iState = 0; iState < statesIdx->size(); ++iState) {
    const CFuint localStateID = (*statesIdx)[iState];
    State *const currState = states[localStateID];
    if (!_useForInitialization) {
      if (currState->isParUpdatable()) {
	if (!isUpdated[localStateID]) {
	  // annull the RHS components and the selected components
	  for (CFuint iv = 0; iv < nbAnnullVars; ++iv) {
	    rhs(localStateID, _varIDToAnnull[iv], nbEqs) = 0.0;
	    (*currState)[_varIDToAnnull[iv]] = 0.0;
	  }

	  isUpdated[localStateID] = true; // flagging is important!!!!!
	}
      }
    }
    else {
      for (CFuint iv = 0; iv < nbAnnullVars; ++iv) {
	(*currState)[_varIDToAnnull[iv]] = 0.0;
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
