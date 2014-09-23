// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "NewtonMethod/NewtonMethod.hh"
#include "NewmarkSetup.hh"
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

MethodCommandProvider<NewmarkSetup, NewtonIteratorData, NewtonMethodModule> newmarkSetupProvider("NewmarkSetup");

//////////////////////////////////////////////////////////////////////////////

NewmarkSetup::NewmarkSetup(std::string name) : 
  StdSetup(name),
  socket_pastStatesD("pastStatesD"),
  socket_pastStatesD2("pastStatesD2")
{
}

//////////////////////////////////////////////////////////////////////////////

NewmarkSetup::~NewmarkSetup()
{
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSource> >
NewmarkSetup::providesSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSource> > result =
    StdSetup::providesSockets();
  
  result.push_back(&socket_pastStatesD);
  result.push_back(&socket_pastStatesD2);
    
  return result;
}

//////////////////////////////////////////////////////////////////////////////

void NewmarkSetup::execute()
{
  CFAUTOTRACE;

  // first call to the parent class
  StdSetup::execute();
  
  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();
  const CFuint nbStates = states.size();
  
  DataHandle<State*> pastStatesD = socket_pastStatesD.getDataHandle();
  pastStatesD.resize(nbStates);
  for (CFuint i = 0; i < nbStates; ++i) {
    pastStatesD[i] = new State();
    *(pastStatesD[i]) = 0.;
  }
  
  DataHandle<State*> pastStatesD2 = socket_pastStatesD2.getDataHandle();
  pastStatesD2.resize(nbStates);
  for (CFuint i = 0; i < nbStates; ++i) {
    pastStatesD2[i] = new State();
    *(pastStatesD2[i]) = 0.;
  }
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace NewtonMethod

  } // namespace Numerics

} // namespace COOLFluiD
