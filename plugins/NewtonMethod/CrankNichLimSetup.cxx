// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "NewtonMethod/NewtonMethod.hh"
#include "CrankNichLimSetup.hh"
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

MethodCommandProvider<CrankNichLimSetup, NewtonIteratorData, NewtonMethodModule> CrankNichLimSetupProvider("CrankNichLimSetup");

//////////////////////////////////////////////////////////////////////////////

CrankNichLimSetup::CrankNichLimSetup(const std::string& name) :
  CrankNichSetup(name),
  socket_pastPastStates("pastPastStates")
{
}

//////////////////////////////////////////////////////////////////////////////

CrankNichLimSetup::~CrankNichLimSetup()
{
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSource> >
CrankNichLimSetup::providesSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSource> > result = CrankNichSetup::providesSockets();

  result.push_back(&socket_pastPastStates);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

void CrankNichLimSetup::execute()
{
  CFAUTOTRACE;

  CrankNichSetup::execute();


  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();

  const CFuint nbStates = states.size();

  // States at time step u_n-1
  DataHandle<State*> pastPastStates = socket_pastPastStates.getDataHandle();

  pastPastStates.resize(nbStates);
  for (CFuint i = 0; i < nbStates; ++i) {
    pastPastStates[i] = new State();
  }
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace NewtonMethod

  } // namespace Numerics

} // namespace COOLFluiD
