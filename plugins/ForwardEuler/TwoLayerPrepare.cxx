// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "ForwardEuler/ForwardEuler.hh"
#include "TwoLayerPrepare.hh"
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

MethodCommandProvider<TwoLayerPrepare, FwdEulerData, ForwardEulerLib> TwoLayerPrepareProvider("TwoLayerPrepare");

//////////////////////////////////////////////////////////////////////////////

void TwoLayerPrepare::execute()
{
  // Only if unsteady computation
  if (SubSystemStatusStack::getActive()->getDT() > 0.){
    DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();
    DataHandle<State*> pastStates = socket_pastStates.getDataHandle();
    DataHandle<State*> interStates = socket_interStates.getDataHandle();

    // Set Initial States to current states
    for(CFuint i = 0; i < states.size(); ++i) {
      *(pastStates[i]) = *(states[i]);
      *(interStates[i]) = *(states[i]);
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

vector<SafePtr<BaseDataSocketSink> > TwoLayerPrepare::needsSockets()
{
  vector<SafePtr<BaseDataSocketSink> > result;

  result.push_back(&socket_states);
  result.push_back(&socket_pastStates);
  result.push_back(&socket_interStates);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace NewtonMethod



} // namespace COOLFluiD
