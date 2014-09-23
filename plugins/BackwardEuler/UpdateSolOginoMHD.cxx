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

#include "BackwardEuler/StdSetup.hh"

#include "BackwardEuler/BackwardEulerMHD.hh"
#include "BackwardEuler/UpdateSolOginoMHD.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace BackwardEuler {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<UpdateSolOginoMHD, BwdEulerData, BackwardEulerMHDModule> 
updateSolOginoMHDProvider("UpdateSolOginoMHD");

//////////////////////////////////////////////////////////////////////

void UpdateSolOginoMHD::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< CFreal >("radiusPlanet","Radius of the planet in which dU is not updated.");
  options.addConfigOption< CFreal >("radiusTransition","Radius of the transition region which is the limit of effect of the planetary state values.");
}

//////////////////////////////////////////////////////////////////////////////

UpdateSolOginoMHD::UpdateSolOginoMHD(const std::string& name) : BwdEulerCom(name),
  socket_rhs("rhs"),
  socket_updateCoeff("updateCoeff"),
  socket_hasToBeUpdated("hasToBeUpdated"),
  _stateCoord(CFNULL),
  socket_states("states")
{
  addConfigOptionsTo(this);

  _radiusPlanet = 4.0;
  setParameter("radiusPlanet",&_radiusPlanet);

  _radiusTransition = 5.0;
  setParameter("radiusTransition",&_radiusTransition);
}

///////////////////////////////////////////////////////////////////////////////

void UpdateSolOginoMHD::setup()
{
  BwdEulerCom::setup();

  _stateCoord.resize(PhysicalModelStack::getActive()->getDim());
}

//////////////////////////////////////////////////////////////////////////////

void UpdateSolOginoMHD::execute()
{
  CFAUTOTRACE;

  DataHandle < Framework::State*, Framework::GLOBAL > states  = socket_states.getDataHandle();
  DataHandle<CFreal> rhs = socket_rhs.getDataHandle();

  // rhs is the temporary placeholder for the dU
  DataHandle<CFreal>& dU = rhs;

  DataHandle<CFreal> updateCoeff = socket_updateCoeff.getDataHandle();

  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  const CFuint nbStates = states.size();

  DataHandle<bool> hasToBeUpdated  = socket_hasToBeUpdated.getDataHandle();
  
  CFuint stateInsidePlanetID = 0;

  for (CFuint iState = 0; iState < nbStates; ++iState) {
     if (!hasToBeUpdated[iState]) {
       stateInsidePlanetID = states[iState]->getLocalID();
       break;
     }
  }   

  for (CFuint iState = 0; iState < nbStates; ++iState) {
    // do the update only if the state is parallel updatable
    if (states[iState]->isParUpdatable()) {
      for (CFuint iEq = 0; iEq < nbEqs; ++iEq) {
        (*states[iState])[iEq] += dU(iState, iEq, nbEqs);

        _stateCoord = states[iState]->getCoordinates();
	const CFreal dis = sqrt(_stateCoord[0]*_stateCoord[0] +
			        _stateCoord[1]*_stateCoord[1] +
			        _stateCoord[2]*_stateCoord[2]);
	cout << "Radius " << _radiusTransition << endl;
        if ((dis <= _radiusPlanet) && (dU(iState, iEq, nbEqs) != 0.0))
	  cout << "Problem = " << iState << " " << dis << endl;

	
      }
    }

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

    const CFreal pi = 3.14159265359;

    // interpolation function in the transition region
    const CFreal f = 0.5*(1-cos(pi*(_radiusPlanet-distState)/(_radiusTransition-_radiusPlanet)));
    
    //if (distState < _radiusPlanet)
      //(*states[iState])[0] = 10.0;
      
 
    // solution in the transition region interpolated, effect of the planet taken into account
    if ((hasToBeUpdated[iState]) && (distState <= _radiusTransition)) {
      for (CFuint iEq = 0; iEq < nbEqs; ++iEq) {
        (*states[iState])[iEq] = f*(*states[iState])[iEq] + (1.0-f)*(*states[stateInsidePlanetID])[iEq];
      }
    }
    //if (rState < 4.0)
      //rhs(states[i]->getLocalID(), getMethodData().getVarID(), nbEqs) = 0.0;
    //if ((rState >= 4.0) && (rState <= 5.0))
      //for (CFuint j = 0; j < nbEqs; ++j)
      //rhs(i,j,nbEqs) = f*rhs(i,j,nbEqs);
  }      
      
  // reset to 0 the update coefficient
  updateCoeff = 0.0;
}

//////////////////////////////////////////////////////////////////////

void UpdateSolOginoMHD::configure ( Config::ConfigArgs& args )
{
  BwdEulerCom::configure(args);
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
UpdateSolOginoMHD::needsSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSink> > result;

  result.push_back(&socket_rhs);
  result.push_back(&socket_updateCoeff);
  result.push_back(&socket_hasToBeUpdated);
  result.push_back(&socket_states);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace BackwardEuler

  } // namespace Numerics

} // namespace COOLFluiD
