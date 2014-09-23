// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "NewtonMethod/NewtonMethod.hh"
#include "LinearizedUnSetup.hh"
#include "Framework/MeshData.hh"
#include "Framework/State.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace NewtonMethod {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<LinearizedUnSetup, NewtonIteratorData, NewtonMethodModule> LinearizedUnSetupProvider("LinearizedUnSetup");

//////////////////////////////////////////////////////////////////////////////

LinearizedUnSetup::LinearizedUnSetup(std::string name) :
StdUnSetup(name),
socket_linearizedStates("linearizedStates")
{
}

//////////////////////////////////////////////////////////////////////////////

LinearizedUnSetup::~LinearizedUnSetup()
{
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
LinearizedUnSetup::needsSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSink> > result = StdUnSetup::needsSockets();

  result.push_back(&socket_linearizedStates);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

void LinearizedUnSetup::execute()
{
  CFAUTOTRACE;

  StdUnSetup::execute();

  DataHandle<State*> linearStates = socket_linearizedStates.getDataHandle();

  // deallocate the states in the storage
  for(CFuint i = 0; i < linearStates.size(); ++i) {
    delete linearStates[i];
  }

}

//////////////////////////////////////////////////////////////////////////////

    } // namespace NewtonMethod

  } // namespace Numerics

} // namespace COOLFluiD
