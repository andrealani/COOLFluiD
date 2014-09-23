// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "NewtonMethod/NewtonMethod.hh"
#include "FSHOSetup.hh"
#include "Framework/MeshData.hh"
#include "Framework/PhysicalModel.hh"
#include "Framework/MethodCommandProvider.hh"
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

MethodCommandProvider<FSHOSetup, NewtonIteratorData, NewtonMethodModule> FSHOSetupProvider("FSHOSetup");

//////////////////////////////////////////////////////////////////////////////

FSHOSetup::FSHOSetup(std::string name) : 
  StdSetup(name),
  socket_interStates("interStates")
{
}

//////////////////////////////////////////////////////////////////////////////

FSHOSetup::~FSHOSetup()
{
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSource> >
FSHOSetup::providesSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSource> > result = 
    StdSetup::providesSockets();
  
  result.push_back(&socket_interStates);
  
  return result;
}

//////////////////////////////////////////////////////////////////////////////

void FSHOSetup::execute()
{
  CFAUTOTRACE;
  
  // first call to the parent class
  StdSetup::execute();
  
  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();
  
  const CFuint nbStates = states.size();
  // const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  
  // States at time step u_n
  DataHandle<State*> interStates = socket_interStates.getDataHandle();
  
  interStates.resize(nbStates);
  for (CFuint i = 0; i < nbStates; ++i) {
    interStates[i] = new State();
  }

}

//////////////////////////////////////////////////////////////////////////////

    } // namespace NewtonMethod

  } // namespace Numerics

} // namespace COOLFluiD
