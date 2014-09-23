// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "NewtonMethod/NewtonMethod.hh"
#include "FSHOUnSetup.hh"
#include "Framework/MeshData.hh"
#include "Framework/State.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace NewtonMethod {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<FSHOUnSetup, NewtonIteratorData, NewtonMethodModule> FSHOUnSetupProvider("FSHOUnSetup");

//////////////////////////////////////////////////////////////////////////////

FSHOUnSetup::FSHOUnSetup(std::string name) : NewtonIteratorCom(name),
socket_pastStates("pastStates"),
socket_interStates("interStates")
{
}

//////////////////////////////////////////////////////////////////////////////

FSHOUnSetup::~FSHOUnSetup()
{
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
FSHOUnSetup::needsSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSink> > result;

  result.push_back(&socket_pastStates);
  result.push_back(&socket_interStates);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

void FSHOUnSetup::execute()
{
  CFAUTOTRACE;

  DataHandle<State*> pastStates = socket_pastStates.getDataHandle();

  // deallocate the states in the storage
  for(CFuint i = 0; i < pastStates.size(); ++i) {
    delete pastStates[i];
  }

  DataHandle<State*> interStates = socket_interStates.getDataHandle();

  // deallocate the states in the storage
  for(CFuint i = 0; i < interStates.size(); ++i) {
    delete interStates[i];
  }
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace NewtonMethod

  } // namespace Numerics

} // namespace COOLFluiD
