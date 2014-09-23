// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "NewtonMethod/NewtonMethod.hh"
#include "LinearizedPrepare.hh"
#include "Framework/PhysicalModel.hh"
#include "Framework/MeshData.hh"
#include "Framework/State.hh"
#include "Framework/SubSystemStatus.hh"
//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace NewtonMethod {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<LinearizedPrepare, NewtonIteratorData, NewtonMethodModule> LinearizedPrepareProvider("LinearizedPrepare");

//////////////////////////////////////////////////////////////////////////////

void LinearizedPrepare::execute()
{

  if(SubSystemStatusStack::getActive()->isSubIterationFirstStep())
  {

    Common::SafePtr<SubSystemStatus> subSysStatus = SubSystemStatusStack::getActive();
    const CFuint nbIter = subSysStatus->getNbIter();

    // Predictor step
    DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();
    DataHandle<State*> pastStates  = socket_pastStates.getDataHandle();
    DataHandle<State*> linearizedStates = socket_linearizedStates.getDataHandle();
    DataHandle<State*> gstates  = socket_gstates.getDataHandle();
    DataHandle<State*> pastGhostStates  = socket_pastGhostStates.getDataHandle();
    DataHandle<State*> linearizedGhostStates = socket_linearizedGhostStates.getDataHandle();

    // Set Predicted States to current states
    // CFuint nodeID;

    if(nbIter > 1){
      for(CFuint i = 0; i < states.size(); ++i) {
        *(linearizedStates[i]) = 2. * (*(states[i])) - 1. * (*(pastStates[i]));
      }

      for(CFuint i = 0; i < gstates.size(); ++i) {
        *(linearizedGhostStates[i]) = 2. * (*(gstates[i])) - 1. * (*(pastGhostStates[i]));
      }
    }
    else{
      for(CFuint i = 0; i < states.size(); ++i) {
        *(linearizedStates[i]) = (*(states[i]));
      }

      for(CFuint i = 0; i < gstates.size(); ++i) {
        *(linearizedGhostStates[i]) = (*(gstates[i]));
      }

    }

    //Then set :  pastStates = states
    //            pastGhostStates = gstates
    // this needs to be called after computing the predicted states
    // Set Initial States to current states
    //  CFuint nodeID;
    for(CFuint i = 0; i < states.size(); ++i) {
      *(pastStates[i]) = *(states[i]);
    }

    //Do the same for the ghosts
    for(CFuint i = 0; i < gstates.size(); ++i) {
      *(pastGhostStates[i]) = *(gstates[i]);
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

vector<SafePtr<BaseDataSocketSink> > LinearizedPrepare::needsSockets()
{
  vector<SafePtr<BaseDataSocketSink> > result = StdPrepare::needsSockets();

  result.push_back(&socket_linearizedStates);
  result.push_back(&socket_linearizedGhostStates);
  result.push_back(&socket_pastGhostStates);
  result.push_back(&socket_gstates);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace NewtonMethod

  } // namespace Numerics

} // namespace COOLFluiD
