// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "NewtonMethod/NewtonMethod.hh"
#include "UFEMUnSetup.hh"
#include "Framework/MeshData.hh"
#include "Framework/State.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace NewtonMethod {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<UFEMUnSetup, NewtonIteratorData, NewtonMethodModule> UFEMUnSetupProvider("UFEMUnSetup");

//////////////////////////////////////////////////////////////////////////////

UFEMUnSetup::UFEMUnSetup(std::string name) : NewtonIteratorCom(name),
  socket_interStates("interStates"),
  socket_pastStates("pastStates"),
  socket_pastpastStates("pastpastStates")
{
}

//////////////////////////////////////////////////////////////////////////////

UFEMUnSetup::~UFEMUnSetup()
{
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> > UFEMUnSetup::needsSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSink> > result;

  result.push_back(&socket_interStates);
  result.push_back(&socket_pastStates);
  result.push_back(&socket_pastpastStates);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

void UFEMUnSetup::execute()
{
  CFAUTOTRACE;

  DataHandle<State*> interStates = socket_interStates.getDataHandle();

  // deallocate the states in the storage
  for(CFuint i = 0; i < interStates.size(); ++i) {
    delete interStates[i];
  }

  DataHandle<State*> pastStates = socket_pastStates.getDataHandle();

  // deallocate the states in the storage
  for(CFuint i = 0; i < pastStates.size(); ++i) {
    delete pastStates[i];
  }

  DataHandle<State*> pastpastStates = socket_pastpastStates.getDataHandle();

  // deallocate the states in the storage
  for(CFuint i = 0; i < pastpastStates.size(); ++i) {
    delete pastpastStates[i];
  }
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace NewtonMethod

  } // namespace Numerics

} // namespace COOLFluiD
