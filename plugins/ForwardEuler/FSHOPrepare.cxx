// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "ForwardEuler/ForwardEuler.hh"
#include "ForwardEuler/ForwardEuler.hh"
#include "FSHOPrepare.hh"
#include "Framework/State.hh"
#include "Framework/PhysicalModel.hh"
#include "Framework/SubSystemStatus.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/MeshData.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



    namespace ForwardEuler {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<FSHOPrepare, FwdEulerData, ForwardEulerLib> fshoPrepareProvider("FSHOPrepare");

//////////////////////////////////////////////////////////////////////////////

void FSHOPrepare::execute()
{
  // Only if unsteady computation
  if (SubSystemStatusStack::getActive()->getDT() > 0.){
    DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();
    DataHandle<State*> pastStates  = socket_pastStates.getDataHandle();
    DataHandle<State*> interStates = socket_interStates.getDataHandle();
    CFuint nbState = states.size();
  //  CF_DEBUG_OBJ(nbState);
    // Set Initial States to current states
    for(CFuint i = 0; i < nbState; ++i) {
      *(pastStates[i])  = *(interStates[i]);
      *(interStates[i]) = *(states[i]);
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

vector<SafePtr<BaseDataSocketSink> > FSHOPrepare::needsSockets()
{
  vector<SafePtr<BaseDataSocketSink> > result;

  CF_DEBUG_POINT;
  result.push_back(&socket_states);
  CF_DEBUG_POINT;
  result.push_back(&socket_pastStates);
  CF_DEBUG_POINT;
  result.push_back(&socket_interStates);
  CF_DEBUG_POINT;

  return result;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace NewtonMethod

} // namespace COOLFluiD
