// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "NewtonMethod/NewtonMethod.hh"
#include "NewmarkResetSystem.hh"
#include "Common/SafePtr.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace NewtonMethod {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<NewmarkResetSystem, NewtonIteratorData, NewtonMethodModule> NewmarkResetSystemProvider("NewmarkResetSystem");

//////////////////////////////////////////////////////////////////////////////

void NewmarkResetSystem::execute()
{
  CFAUTOTRACE;

  //1st call the parent
  ResetSystem::execute();

  DataHandle < Framework::State*, Framework::GLOBAL > states  = socket_states.getDataHandle();
  DataHandle<State*> pastStates  = socket_pastStates.getDataHandle();

  const CFuint statesSize = states.size();

  //First update the pastStates otherwise at the first iteration,
  //if you restart, you have states=X and pastStates=0, so a non zero speed!!
  //but after, you should update the pastStates afterwards otherwise,
  //states-pastStates will always be zero
  for (CFuint iState = 0; iState < statesSize; ++iState) {
    *(pastStates[iState]) = *(states[iState]);
  }

}

//////////////////////////////////////////////////////////////////////////////

vector<SafePtr<BaseDataSocketSink> > NewmarkResetSystem::needsSockets()
{
  vector<SafePtr<BaseDataSocketSink> > result = ResetSystem::needsSockets();

  result.push_back(&socket_states);
  result.push_back(&socket_pastStates);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace NewtonMethod

  } // namespace Numerics

} // namespace COOLFluiD
