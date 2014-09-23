// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Framework/MethodCommandProvider.hh"

#include "EmptySpaceMethod/StdSetup.hh"
#include "EmptySpaceMethod/Empty.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace COOLFluiD::Framework;

namespace COOLFluiD {
  namespace EmptySpaceMethod {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider< StdSetup,EmptySolverData,EmptyModule >
  stdSetupProvider("StdSetup");

//////////////////////////////////////////////////////////////////////////////

StdSetup::StdSetup(const std::string& name) :
  EmptySolverCom(name),
  socket_nstatesProxy("nstatesProxy"),
  socket_states("states")
{
}

//////////////////////////////////////////////////////////////////////////////

StdSetup::~StdSetup()
{
}

//////////////////////////////////////////////////////////////////////////////

std::vector< Common::SafePtr< BaseDataSocketSink > >
  StdSetup::needsSockets()
{
  std::vector< Common::SafePtr< BaseDataSocketSink > > result;
  result.push_back(&socket_states);
  return result;
}

//////////////////////////////////////////////////////////////////////////////

std::vector< Common::SafePtr< BaseDataSocketSource > >
  StdSetup::providesSockets()
{
  std::vector< Common::SafePtr< BaseDataSocketSource > > result;
  result.push_back(&socket_nstatesProxy);
  return result;
}

//////////////////////////////////////////////////////////////////////////////

void StdSetup::execute()
{
  CFAUTOTRACE;

  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();
  DataHandle<ProxyDofIterator< RealVector >* > nstatesProxy =
    socket_nstatesProxy.getDataHandle();

  const CFuint nbStates = states.size();
  nstatesProxy.resize(1);

  // set node to state mapping
  m_nodeIDToStateID.resize(nbStates);
  for (CFuint stateID=0;stateID<nbStates;++stateID) {
    const CFuint nodeID = states[stateID]->getCoordinates().getLocalID();
    cf_assert(nodeID < nbStates);
    m_nodeIDToStateID[nodeID] = stateID;
  }
  nstatesProxy[0] =
    new DofDataHandleIterator< RealVector,State, GLOBAL >(states,&m_nodeIDToStateID);
}

//////////////////////////////////////////////////////////////////////////////

  }  // namespace EmptySpaceMethod
}  // namespace COOLFluiD

