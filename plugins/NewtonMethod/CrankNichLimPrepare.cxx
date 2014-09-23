// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "NewtonMethod/NewtonMethod.hh"
#include "CrankNichLimPrepare.hh"
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

MethodCommandProvider<CrankNichLimPrepare, NewtonIteratorData, NewtonMethodModule> CrankNichLimPrepareProvider("CrankNichLimPrepare");

//////////////////////////////////////////////////////////////////////////////

void CrankNichLimPrepare::execute()
{

  DataHandle<State*> pastStates  = socket_pastStates.getDataHandle();
  DataHandle<State*> pastPastStates = socket_pastPastStates.getDataHandle();

  // First back up the states at n-1
  if(SubSystemStatusStack::getActive()->isSubIterationFirstStep()){
    for(CFuint i = 0; i < pastStates.size(); ++i) {
      *(pastPastStates[i]) = *(pastStates[i]);
    }
  }

  // Then back up the states at n
  StdPrepare::execute();
}

//////////////////////////////////////////////////////////////////////////////

vector<SafePtr<BaseDataSocketSink> > CrankNichLimPrepare::needsSockets()
{
  vector<SafePtr<BaseDataSocketSink> > result = StdPrepare::needsSockets();

  result.push_back(&socket_pastPastStates);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace NewtonMethod

  } // namespace Numerics

} // namespace COOLFluiD
