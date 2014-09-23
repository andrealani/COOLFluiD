#include "FluctSplit/FluctSplit.hh"
#include "StrongNoSlipWallAdiabaticNSKOmega2D.hh"
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

MethodCommandProvider<StrongNoSlipWallAdiabaticNSKOmega2D, FluctuationSplitData, FluctSplitModule> strongNoSlipWallAdiabaticKOmegaNS2DProvider("StrongNoSlipWallAdiabaticNSKOmega2D");

//////////////////////////////////////////////////////////////////////////////
void StrongNoSlipWallAdiabaticNSKOmega2D::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOption< CFreal >("u0","x-coponent of the velocity of the wall");
   options.addConfigOption< CFreal >("v0","y-coponent of the velocity of the wall");

}


//////////////////////////////////////////////////////////////////////////////

StrongNoSlipWallAdiabaticNSKOmega2D::StrongNoSlipWallAdiabaticNSKOmega2D(const std::string& name) :
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

StrongNoSlipWallAdiabaticNSKOmega2D::~StrongNoSlipWallAdiabaticNSKOmega2D()
{
}

//////////////////////////////////////////////////////////////////////////////

void StrongNoSlipWallAdiabaticNSKOmega2D::setup()
{
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
StrongNoSlipWallAdiabaticNSKOmega2D::needsSockets()
{

  std::vector<Common::SafePtr<BaseDataSocketSink> > result;

  result.push_back(&socket_rhs);
  result.push_back(&socket_states);
  result.push_back(&socket_isUpdated);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

void StrongNoSlipWallAdiabaticNSKOmega2D::executeOnTrs()
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
        (*state)[4] = 0.0;
        (*state)[5] = 0.0;
	

        // reset to 0 the momentum residual for BC nodes
        rhs(stateID, 1, nbEqs) = 0.;
        rhs(stateID, 2, nbEqs) = 0.;
        rhs(stateID, 4, nbEqs) = 0.;
        rhs(stateID, 5, nbEqs) = 0.;
      }
      isUpdated[stateID] = true; // flagging is important!!!!!
    }
   else {
     (*state)[1] = m_u0;
     (*state)[2] = m_v0;
     (*state)[4] = 0.0;
     (*state)[5] = 0.0;    
     
     
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
