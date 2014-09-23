// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.



#include "Framework/MeshData.hh"
#include "Framework/State.hh"
#include "Framework/PhysicalModel.hh"
#include "Common/CFLog.hh"

#include "Framework/MethodCommandProvider.hh"

#include "BackwardEuler/BackwardEulerMHD.hh"
#include "BackwardEuler/StdSetupOginoMHD.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace BackwardEuler {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<StdSetupOginoMHD, BwdEulerData, BackwardEulerMHDModule> 
stdSetupOginoMHDProvider("StdSetupOginoMHD");

//////////////////////////////////////////////////////////////////////////////

void StdSetupOginoMHD::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< CFreal >("radiusPlanet","Radius of the planet in which dU is not updated.");  
  options.addConfigOption< CFreal >("radiusTransition","Radius of the transition region which is the limit of effect of the planetary state values.");
  options.addConfigOption< vector<CFreal> >("planetSolutionVars","the solution variables inside the planet region.");
}

//////////////////////////////////////////////////////////////////////////////

StdSetupOginoMHD::StdSetupOginoMHD(const std::string& name): StdSetup(name),
  socket_hasToBeUpdated("hasToBeUpdated"),
  _stateCoord(CFNULL)	
{
  addConfigOptionsTo(this);

  _radiusPlanet = 4.0;
  setParameter("radiusPlanet",&_radiusPlanet);

  _radiusTransition = 5.0;
  setParameter("radiusTransition",&_radiusTransition);

  _planetSolutionVars = vector< CFreal >();
  setParameter("planetSolutionVars",&_planetSolutionVars);
}

///////////////////////////////////////////////////////////////////////////////

void StdSetupOginoMHD::setup()
{
  StdSetup::setup();

  cf_assert(_planetSolutionVars.size()==PhysicalModelStack::getActive()->getNbEq());

  _stateCoord.resize(PhysicalModelStack::getActive()->getDim());
}

//////////////////////////////////////////////////////////////////////////////

void StdSetupOginoMHD::execute()
{
  CFAUTOTRACE;

  StdSetup::execute();
  
  DataHandle < Framework::State*, Framework::GLOBAL > states  = socket_states.getDataHandle();
  DataHandle<bool> hasToBeUpdated  = socket_hasToBeUpdated.getDataHandle();
  
  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  const CFuint nbStates = states.size();
  hasToBeUpdated.resize(nbStates);
  
  for (CFuint iState = 0; iState < nbStates; ++iState) {
    hasToBeUpdated[iState] = true;
    _stateCoord = states[iState]->getCoordinates();
    // distance of the state from the origin (center of the planet is assumed to be at the origin)
    CFreal distState;
    if (PhysicalModelStack::getActive()->getDim() == 2)
      distState = sqrt(_stateCoord[0]*_stateCoord[0] +
		       _stateCoord[1]*_stateCoord[1]);
    else
      distState = sqrt(_stateCoord[0]*_stateCoord[0] +
		       _stateCoord[1]*_stateCoord[1] + 
		       _stateCoord[2]*_stateCoord[2]);
    if (distState <= _radiusPlanet) {
      hasToBeUpdated[iState] = false;
      // solution variables inside the planet assigned to fixed user-defined values
      for (CFuint iEq = 0; iEq < nbEqs; ++iEq)
        (*states[iState])[iEq] = _planetSolutionVars[iEq]; 
    }  
  }
}

//////////////////////////////////////////////////////////////////////

void StdSetupOginoMHD::configure ( Config::ConfigArgs& args )
{
  StdSetup::configure(args);
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSource> >
StdSetupOginoMHD::providesSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSource> > result = 
    StdSetup::providesSockets();
  result.push_back(&socket_hasToBeUpdated);
  return result;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace BackwardEuler

  } // namespace Numerics

} // namespace COOLFluiD
