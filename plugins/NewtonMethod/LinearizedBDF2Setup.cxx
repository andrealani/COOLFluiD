// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "NewtonMethod/NewtonMethod.hh"
#include "LinearizedBDF2Setup.hh"
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

MethodCommandProvider<LinearizedBDF2Setup, NewtonIteratorData, NewtonMethodModule> LinearizedBDF2SetupProvider("LinearizedBDF2Setup");

//////////////////////////////////////////////////////////////////////////////

LinearizedBDF2Setup::LinearizedBDF2Setup(std::string name) : BDF2Setup(name),
  socket_linearizedStates("linearizedStates"),
  socket_pastRhs("pastRhs")
{
}

//////////////////////////////////////////////////////////////////////////////

LinearizedBDF2Setup::~LinearizedBDF2Setup()
{
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSource> >
LinearizedBDF2Setup::providesSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSource> > result = BDF2Setup::providesSockets();

  result.push_back(&socket_linearizedStates);
  result.push_back(&socket_pastRhs);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

void LinearizedBDF2Setup::execute()
{
  CFAUTOTRACE;

  BDF2Setup::execute();

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
