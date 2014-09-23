#include "FiniteVolume/FiniteVolume.hh"


#include "InitStateTurb.hh"
#include "Common/CFLog.hh"
#include "Framework/MethodCommandProvider.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<InitStateTurb, CellCenterFVMData, FiniteVolumeModule>
InitStateTurbProvider("InitStateTurb");

//////////////////////////////////////////////////////////////////////////////

void InitStateTurb::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< CFreal >("Intensity","Intensity of turbulent fluctuations in percentage of the mean flow");
  options.addConfigOption< CFuint >("VelocityID","ID of velocity component to apply fluctuations to");
}

//////////////////////////////////////////////////////////////////////////////

InitStateTurb::InitStateTurb(const std::string& name) :
  InitState(name),
  m_intensity(0.),
  m_velocityID(1)
{
   addConfigOptionsTo(this);
   setParameter("Intensity",&m_intensity);
   setParameter("VelocityID",&m_velocityID);
}

//////////////////////////////////////////////////////////////////////////////

InitStateTurb::~InitStateTurb()
{
}

//////////////////////////////////////////////////////////////////////////////

void InitStateTurb::configure ( Config::ConfigArgs& args )
{
  CFAUTOTRACE;
  InitState::configure(args);
}

//////////////////////////////////////////////////////////////////////////////

CFreal InitStateTurb::rand(const CFreal& a, const CFreal& b)
{
  return ((b-a)*((CFreal)std::rand()/RAND_MAX))+a;
}

//////////////////////////////////////////////////////////////////////////////

void InitStateTurb::executeOnTrs()
{  
  SafePtr<TopologicalRegionSet> trs = getCurrentTRS();
  CFLogDebugMax( "InitStateTurb::executeOnTrs() called for TRS: "
  << trs->getName() << "\n");

  SafePtr<std::vector<CFuint> > trsStates = trs->getStatesInTrs();
  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();
  
  std::vector<CFuint>::iterator itd;
  
  if(_inputAdimensionalValues)
  {
    for (itd = trsStates->begin(); itd != trsStates->end(); ++itd) {
      State* const currState = states[(*itd)];
      _vFunction.evaluate(currState->getCoordinates(), *_input);
      *currState = *_inputToUpdateVar->transform(_input);
      (*currState)[m_velocityID] += m_intensity*rand(-(*currState)[m_velocityID],(*currState)[m_velocityID]);
    }
  }
  else
  {
    State dimState;
    for (itd = trsStates->begin(); itd != trsStates->end(); ++itd) {
      State* const currState = states[(*itd)];
      _vFunction.evaluate(currState->getCoordinates(), *_input);
      dimState = *_inputToUpdateVar->transform(_input);
      dimState[m_velocityID] += m_intensity*rand(-dimState[m_velocityID],dimState[m_velocityID]);      
      _varSet->setAdimensionalValues(dimState, *currState);
    }
  }

}

//////////////////////////////////////////////////////////////////////////////

void InitStateTurb::setup()
{
  CFAUTOTRACE;
  InitState::setup();

}

//////////////////////////////////////////////////////////////////////////////

void InitStateTurb::unsetup()
{
  CFAUTOTRACE;
  InitState::unsetup();
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD
