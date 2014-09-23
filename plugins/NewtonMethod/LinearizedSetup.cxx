// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "NewtonMethod/NewtonMethod.hh"
#include "LinearizedSetup.hh"
#include "Framework/MeshData.hh"
#include "Framework/PhysicalModel.hh"
#include "Framework/State.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace NewtonMethod {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<LinearizedSetup, NewtonIteratorData, NewtonMethodModule> LinearizedSetupProvider("LinearizedSetup");

//////////////////////////////////////////////////////////////////////////////

LinearizedSetup::LinearizedSetup(std::string name) : StdSetup(name),
  socket_linearizedStates("linearizedStates"),
  socket_pastRhs("pastRhs")
{
}

//////////////////////////////////////////////////////////////////////////////

LinearizedSetup::~LinearizedSetup()
{
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSource> >
LinearizedSetup::providesSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSource> > result = StdSetup::providesSockets();

  result.push_back(&socket_linearizedStates);
  result.push_back(&socket_pastRhs);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

void LinearizedSetup::execute()
{
  CFAUTOTRACE;

  StdSetup::execute();

  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();
  const CFuint nbStates = states.size();

  // States at time step u_n
  DataHandle<State*> linearStates = socket_linearizedStates.getDataHandle();

  linearStates.resize(nbStates);
  for (CFuint i = 0; i < nbStates; ++i) {
    linearStates[i] = new State();
  }

  DataHandle<CFreal> rhs = socket_rhs.getDataHandle();

  DataHandle<CFreal> pastRhs = socket_pastRhs.getDataHandle();
  pastRhs.resize(rhs.size());
  pastRhs = 0.0;

}

//////////////////////////////////////////////////////////////////////////////

    } // namespace NewtonMethod

  } // namespace Numerics

} // namespace COOLFluiD
