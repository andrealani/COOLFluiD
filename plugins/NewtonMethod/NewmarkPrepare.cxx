// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "NewtonMethod/NewtonMethod.hh"
#include "NewmarkPrepare.hh"
#include "Framework/PhysicalModel.hh"
#include "Framework/MeshData.hh"
#include "Framework/State.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace NewtonMethod {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<NewmarkPrepare,
                      NewtonIteratorData,
                      NewtonMethodModule>
                      NewmarkPrepareProvider("NewmarkPrepare");

//////////////////////////////////////////////////////////////////////////////

void NewmarkPrepare::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOption< CFreal >("Alpha","Alpha coef");
   options.addConfigOption< CFreal >("Gamma","Gamma coef");
}

//////////////////////////////////////////////////////////////////////////////

NewmarkPrepare::NewmarkPrepare(std::string name) :
  NewtonIteratorCom(name),
  socket_states("states"),
  socket_pastStates("pastStates"),
  socket_pastStatesD("pastStatesD"),
  socket_pastStatesD2("pastStatesD2")
{
   addConfigOptionsTo(this);
  _alpha = 0.5;
   setParameter("Alpha",&_alpha);

 _gamma = 0.5;
   setParameter("Gamma",&_gamma);
}

//////////////////////////////////////////////////////////////////////////////

void NewmarkPrepare::execute()
{

  //only do the prepare step if we are doing normal simulation or
  //on the first step if we are subiterating
  if(SubSystemStatusStack::getActive()->isSubIterationFirstStep()){
    DataHandle < Framework::State*, Framework::GLOBAL > states  = socket_states.getDataHandle();
    DataHandle<State*> pastStates  = socket_pastStates.getDataHandle();
    DataHandle<State*> pastStatesD  = socket_pastStatesD.getDataHandle();
    DataHandle<State*> pastStatesD2 = socket_pastStatesD2.getDataHandle();

    const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
    const CFuint statesSize = states.size();

  //First update the pastStates otherwise at the first iteration,
  //if you restart, you have states=X and pastStates=0, so a non zero speed!!
  //but after, you should update the pastStates afterwards otherwise,
  //states-pastStates will always be zero
  if(SubSystemStatusStack::getActive()->getNbIter() < 2)
  {
    for (CFuint iState = 0; iState < statesSize; ++iState) {
      *(pastStates[iState]) = *(states[iState]);
    }
  }


    // Compute the coefficients
    const CFreal dt = SubSystemStatusStack::getActive()->getDT();
    const CFreal a1 = _alpha * dt;
    const CFreal a2 = (1.- _alpha) * dt;
    const CFreal a3 = 2./(_gamma * dt * dt);
    const CFreal a4 = dt * a3;
    const CFreal a5 = (1./_gamma) - 1.;

    // Compute the first and second derivatives
    State statesD;
    State statesD2;

    for (CFuint iState = 0; iState < statesSize; ++iState) {
      // do the update only if the state is parallel updatable
      if (states[iState]->isParUpdatable()) {
        for (CFuint iEq = 0; iEq < nbEqs; ++iEq) {
          statesD2[iEq] =    a3 * ((*states[iState])[iEq] - (*pastStates[iState])[iEq])
                          - (a4 * (*pastStatesD[iState])[iEq])
                          - (a5 * (*pastStatesD2[iState])[iEq]);
          statesD[iEq]  =         (*pastStatesD[iState])[iEq]
                          + (a2 * (*pastStatesD2[iState])[iEq])
                          + (a1 * statesD2[iEq]);
          }

        *pastStatesD[iState] = statesD;
        *pastStatesD2[iState] = statesD2;
      }
    }

    //Update the pastStates
    for (CFuint iState = 0; iState < statesSize; ++iState) {
      *(pastStates[iState]) = *(states[iState]);
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

vector<SafePtr<BaseDataSocketSink> > NewmarkPrepare::needsSockets()
{
  vector<SafePtr<BaseDataSocketSink> > result;

  result.push_back(&socket_states);
  result.push_back(&socket_pastStates);
  result.push_back(&socket_pastStatesD);
  result.push_back(&socket_pastStatesD2);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace NewtonMethod

  } // namespace Numerics

} // namespace COOLFluiD
