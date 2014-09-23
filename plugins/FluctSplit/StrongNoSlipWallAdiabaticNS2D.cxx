#include "FluctSplit/FluctSplit.hh"
#include "StrongNoSlipWallAdiabaticNS2D.hh"
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

MethodCommandProvider<StrongNoSlipWallAdiabaticNS2D, FluctuationSplitData, FluctSplitModule> strongNoSlipWallAdiabaticNS2DProvider("StrongNoSlipWallAdiabaticNS2D");

//////////////////////////////////////////////////////////////////////////////
void StrongNoSlipWallAdiabaticNS2D::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOption< CFreal >("u0","x-coponent of the velocity of the wall");
   options.addConfigOption< CFreal >("v0","y-coponent of the velocity of the wall");

}


//////////////////////////////////////////////////////////////////////////////

StrongNoSlipWallAdiabaticNS2D::StrongNoSlipWallAdiabaticNS2D(const std::string& name) :
  FluctuationSplitCom(name),
  socket_rhs("rhs"),
  socket_states("states"),
  socket_isUpdated("isUpdated")
{
   addConfigOptionsTo(this);

  m_u0 = 0.0;
  setParameter("u0",&m_u0);

  m_v0 = 0.0;
  setParameter("v0",&m_v0);

}

//////////////////////////////////////////////////////////////////////////////

StrongNoSlipWallAdiabaticNS2D::~StrongNoSlipWallAdiabaticNS2D()
{
}

//////////////////////////////////////////////////////////////////////////////

void StrongNoSlipWallAdiabaticNS2D::setup()
{
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
StrongNoSlipWallAdiabaticNS2D::needsSockets()
{

  std::vector<Common::SafePtr<BaseDataSocketSink> > result;

  result.push_back(&socket_rhs);
  result.push_back(&socket_states);
  result.push_back(&socket_isUpdated);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

void StrongNoSlipWallAdiabaticNS2D::executeOnTrs()
{
  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();
  DataHandle<bool> isUpdated = socket_isUpdated.getDataHandle();
  DataHandle<CFreal> rhs = socket_rhs.getDataHandle();

  Common::SafePtr< vector<CFuint> > const statesIdx = getCurrentTRS()->getStatesInTrs();
  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  const bool initialization = getMethodData().isInitializationPhase();

  for (CFuint iState = 0; iState < statesIdx->size(); ++iState)
  {
    const CFuint stateID = (*statesIdx)[iState];
    State *const state = states[stateID];

    if (!initialization)
    {
      if (!isUpdated[stateID])
      {
        (*state)[1] = m_u0;
        (*state)[2] = m_v0;

        // reset to 0 the momentum residual for BC nodes
        rhs(stateID, 1, nbEqs) = 0.;
        rhs(stateID, 2, nbEqs) = 0.;
      }
      isUpdated[stateID] = true; // flagging is important!!!!!
    }
   else {
     (*state)[1] = m_u0;
     (*state)[2] = m_v0;
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
