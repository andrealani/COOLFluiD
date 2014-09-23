// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "NewtonMethod/NewtonMethod.hh"
#include "NewmarkUnSetup.hh"
#include "Framework/MeshData.hh"
#include "Framework/State.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace NewtonMethod {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<NewmarkUnSetup, NewtonIteratorData, NewtonMethodModule> NewmarkUnSetupProvider("NewmarkUnSetup");

//////////////////////////////////////////////////////////////////////////////

NewmarkUnSetup::NewmarkUnSetup(std::string name) : NewtonIteratorCom(name),
socket_pastStates("pastStates"),
socket_pastStatesD("pastStatesD"),
socket_pastStatesD2("pastStatesD2")
{
}

//////////////////////////////////////////////////////////////////////////////

NewmarkUnSetup::~NewmarkUnSetup()
{
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
NewmarkUnSetup::needsSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSink> > result;

  result.push_back(&socket_pastStates);
  result.push_back(&socket_pastStatesD);
  result.push_back(&socket_pastStatesD2);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

void NewmarkUnSetup::execute()
{
  CFAUTOTRACE;

  DataHandle<State*> pastStates = socket_pastStates.getDataHandle();

  // deallocate the states in the storage
  for(CFuint i = 0; i < pastStates.size(); ++i) {
    delete pastStates[i];
  }

  DataHandle<State*> pastStatesD = socket_pastStatesD.getDataHandle();

  // deallocate the states in the storage
  for(CFuint i = 0; i < pastStatesD.size(); ++i) {
    delete pastStatesD[i];
  }

  DataHandle<State*> pastStatesD2 = socket_pastStatesD2.getDataHandle();

  // deallocate the states in the storage
  for(CFuint i = 0; i < pastStatesD2.size(); ++i) {
    delete pastStatesD2[i];
  }
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace NewtonMethod

  } // namespace Numerics

} // namespace COOLFluiD
