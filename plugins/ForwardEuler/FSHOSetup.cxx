// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Framework/State.hh"
#include "Framework/PhysicalModel.hh"
#include "Framework/SubSystemStatus.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/MeshData.hh"

#include "ForwardEuler/ForwardEuler.hh"
#include "ForwardEuler/FSHOSetup.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



    namespace ForwardEuler {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<FSHOSetup, FwdEulerData, ForwardEulerLib> fshoSetupProvider("FSHOSetup");

//////////////////////////////////////////////////////////////////////////////

FSHOSetup::FSHOSetup(const std::string& name) :
  FwdEulerCom(name),
  socket_rhs("rhs"),
  socket_updateCoeff("updateCoeff"),
  socket_pastStates("pastStates"),
  socket_interStates("interStates"),
  socket_bStatesNeighbors("bStatesNeighbors"),
  socket_states("states")
{
}

//////////////////////////////////////////////////////////////////////////////

FSHOSetup::~FSHOSetup()
{
}

//////////////////////////////////////////////////////////////////////////////

void FSHOSetup::configure ( Config::ConfigArgs& args )
{
  FwdEulerCom::configure(args);

  //  configureNestedSockets(args);
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSource> >
FSHOSetup::providesSockets()
{
  CF_DEBUG_POINT;
  std::vector<Common::SafePtr<BaseDataSocketSource> > result;

  CF_DEBUG_POINT;
  result.push_back(&socket_rhs);

  CF_DEBUG_POINT;
  result.push_back(&socket_updateCoeff);

  CF_DEBUG_POINT;
  result.push_back(&socket_pastStates);

  CF_DEBUG_POINT;
  result.push_back(&socket_interStates);

  CF_DEBUG_POINT;
  result.push_back(&socket_bStatesNeighbors);

  CF_DEBUG_POINT;

  return result;
}

//////////////////////////////////////////////////////////////////////////////

void FSHOSetup::execute()
{
  CFAUTOTRACE;

  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();

  getMethodData().setNormSize(1);

  const CFuint nbStates = states.size();
  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();

  if (SubSystemStatusStack::getActive()->getDT() > 0.){
    // States at time step u_n
    DataHandle<State*> pastStates = socket_pastStates.getDataHandle();

    if(pastStates.size() != nbStates)
    {
      pastStates.resize(nbStates);
      for (CFuint i = 0; i < nbStates; ++i) {
        pastStates[i] = new State();
      }
    }
    // States at time step u_{n+1/2}
    DataHandle<State*> interStates = socket_interStates.getDataHandle();

    if(interStates.size() != nbStates)
    {
      interStates.resize(nbStates);
      for (CFuint i = 0; i < nbStates; ++i) {
        interStates[i] = new State();
      }
    }
  }

  DataHandle<CFreal> rhs = socket_rhs.getDataHandle();
  rhs.resize(nbStates*nbEqs);
  rhs = 0.0;

  DataHandle<CFreal> updateCoeff = socket_updateCoeff.getDataHandle();
  updateCoeff.resize(nbStates);

  socket_bStatesNeighbors.getDataHandle().resize(nbStates);
}

//////////////////////////////////////////////////////////////////////////////

vector<SafePtr<BaseDataSocketSink> > FSHOSetup::needsSockets()
{
  vector<SafePtr<BaseDataSocketSink> > result;

  result.push_back(&socket_states);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace ForwardEuler



} // namespace COOLFluiD
