// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "NewtonMethod/NewtonMethod.hh"
#include "StdSetup.hh"
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

MethodCommandProvider<StdSetup, NewtonIteratorData, NewtonMethodModule> stdSetupProvider("StdSetup");

//////////////////////////////////////////////////////////////////////////////

StdSetup::StdSetup(const std::string& name) :
  NewtonIteratorCom(name),
  socket_rhs("rhs"),
  socket_updateCoeff("updateCoeff"),
  socket_pastStates("pastStates"),
  socket_bStatesNeighbors("bStatesNeighbors"),
  socket_states("states")
{
}

//////////////////////////////////////////////////////////////////////////////

StdSetup::~StdSetup()
{
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSource> >
StdSetup::providesSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSource> > result;

  result.push_back(&socket_rhs);
  result.push_back(&socket_updateCoeff);
  result.push_back(&socket_pastStates);
  result.push_back(&socket_bStatesNeighbors);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

void StdSetup::execute()
{
  CFAUTOTRACE;

  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();

  const CFuint nbStates = states.size();
  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();

  // States at time step u_n
  DataHandle<State*> pastStates = socket_pastStates.getDataHandle();

  if(pastStates.size() != nbStates)
  {
    pastStates.resize(nbStates);
    for (CFuint i = 0; i < nbStates; ++i) {
      pastStates[i] = new State();
    }
  }

  DataHandle<CFreal> rhs = socket_rhs.getDataHandle();
  rhs.resize(nbStates*nbEqs);
  rhs = 0.0;

  DataHandle<CFreal> updateCoeff = socket_updateCoeff.getDataHandle();
  const CFuint nbLSS = getMethodData().getLinearSystemSolver().size();
  cf_assert(nbLSS > 0);
  updateCoeff.resize(nbStates*nbLSS);

  DataHandle<std::valarray<State*> > bStatesNeighbors =
    socket_bStatesNeighbors.getDataHandle();
  bStatesNeighbors.resize(nbStates);
}

//////////////////////////////////////////////////////////////////////////////

vector<SafePtr<BaseDataSocketSink> > StdSetup::needsSockets()
{
  vector<SafePtr<BaseDataSocketSink> > result;

  result.push_back(&socket_states);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace NewtonMethod

  } // namespace Numerics

} // namespace COOLFluiD
