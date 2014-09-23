#include "FluctSplit/FluctSplitNavierStokes.hh"
#include "StrongMirrorAxisymmEuler2DCons.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/MeshData.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<StrongMirrorAxisymmEuler2DCons, FluctuationSplitData, FluctSplitNavierStokesModule> strongMirrorAxisymmEuler2DConsProvider("StrongMirrorAxisymmEuler2DCons");

//////////////////////////////////////////////////////////////////////////////

StrongMirrorAxisymmEuler2DCons::StrongMirrorAxisymmEuler2DCons
  (const std::string& name) :
  FluctuationSplitCom(name),
  socket_rhs("rhs"),
  socket_states("states")
{
}

//////////////////////////////////////////////////////////////////////////////

StrongMirrorAxisymmEuler2DCons::~StrongMirrorAxisymmEuler2DCons()
{
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
StrongMirrorAxisymmEuler2DCons::needsSockets()
{

  std::vector<Common::SafePtr<BaseDataSocketSink> > result;

  result.push_back(&socket_rhs);
  result.push_back(&socket_states);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

void StrongMirrorAxisymmEuler2DCons::setup()
{
  FluctuationSplitCom::setup();
}

//////////////////////////////////////////////////////////////////////////////

void StrongMirrorAxisymmEuler2DCons::executeOnTrs()
{
  DataHandle< CFreal> rhs = socket_rhs.getDataHandle();
  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();

  Common::SafePtr< vector<CFuint> > const statesIdx = getCurrentTRS()->getStatesInTrs();
  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();

  for (CFuint iState = 0; iState < statesIdx->size(); ++iState) {
    const CFuint stateID = (*statesIdx)[iState];

    // set to 0 the residual for the momentum y component
    rhs(stateID, 2, nbEqs) = 0.0;
    // set to 0 the state for the momentum y component
    (*(states[stateID]))[2] = 0.0;
  }
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
