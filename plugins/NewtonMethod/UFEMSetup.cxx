// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "NewtonMethod/NewtonMethod.hh"
#include "UFEMSetup.hh"
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

MethodCommandProvider<UFEMSetup, NewtonIteratorData, NewtonMethodModule> UFEMSetupProvider("UFEMSetup");

//////////////////////////////////////////////////////////////////////////////

UFEMSetup::UFEMSetup(std::string name) : NewtonIteratorCom(name),
  socket_interStates("interStates"),
  socket_pastStates("pastStates"),
  socket_pastpastStates("pastpastStates"),
  socket_rhs("rhs"),
  socket_updateCoeff("updateCoeff"),
  socket_bStatesNeighbors("bStatesNeighbors"),
  socket_states("states")
{
}

//////////////////////////////////////////////////////////////////////////////

UFEMSetup::~UFEMSetup()
{
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSource> > UFEMSetup::providesSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSource> > result;

  result.push_back(&socket_interStates);
  result.push_back(&socket_pastStates);
  result.push_back(&socket_pastpastStates);
  result.push_back(&socket_rhs);
  result.push_back(&socket_updateCoeff);
  result.push_back(&socket_bStatesNeighbors);
  
  return result;
}

//////////////////////////////////////////////////////////////////////////////

void UFEMSetup::execute()
{
  CFAUTOTRACE;

  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();

  const CFuint nbStates = states.size();
  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();

  // States at time step u_n+1/2
  DataHandle<State*> interStates = socket_interStates.getDataHandle();
  interStates.resize(nbStates);
  for (CFuint i = 0; i < nbStates; ++i) {
    interStates[i] = new State();
    interStates[i]->setSpaceCoordinates(states[i]->getNodePtr());
  }

  // States at time step u_n-1
  DataHandle<State*> pastStates = socket_pastStates.getDataHandle();
  pastStates.resize(nbStates);
  for (CFuint i = 0; i < nbStates; ++i) {
    pastStates[i] = new State();
  }

  // States at time step u_n-2
  DataHandle<State*> pastpastStates = socket_pastpastStates.getDataHandle();
  pastpastStates.resize(nbStates);
  for (CFuint i = 0; i < nbStates; ++i) {
    pastpastStates[i] = new State();
  }

  // right hand side
  DataHandle<CFreal> rhs = socket_rhs.getDataHandle();
  rhs.resize(nbStates*nbEqs);
  rhs = 0.0;

  // updaetcoeff
  DataHandle<CFreal> updateCoeff = socket_updateCoeff.getDataHandle();
  const CFuint nbLSS = getMethodData().getLinearSystemSolver().size();
  cf_assert(nbLSS > 0);
  updateCoeff.resize(nbStates*nbLSS);

  // bstatesneighbours
  DataHandle<valarray<State*> > bStatesNeighbors = socket_bStatesNeighbors.getDataHandle();
  bStatesNeighbors.resize(nbStates);

}

//////////////////////////////////////////////////////////////////////////////

vector<SafePtr<BaseDataSocketSink> > UFEMSetup::needsSockets()
{
  vector<SafePtr<BaseDataSocketSink> > result;

  result.push_back(&socket_states);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace NewtonMethod

  } // namespace Numerics

} // namespace COOLFluiD
