// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Common/NotImplementedException.hh"

#include "Framework/GlobalJacobianSparsity.hh"
#include "Framework/MeshData.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

GlobalJacobianSparsity::GlobalJacobianSparsity() :
  Common::OwnedObject(),
  socket_states("Null"), 
  socket_nodes("Null"),
  socket_bStatesNeighbors("Null")
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

GlobalJacobianSparsity::~GlobalJacobianSparsity()
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

void GlobalJacobianSparsity::setDataSockets
(DataSocketSink<State*, Framework::GLOBAL> statesSocket,
 DataSocketSink<Node*, Framework::GLOBAL> nodesSocket,
 DataSocketSink<std::valarray<State*> > bStatesNeighborsSocket)
{
  socket_states = statesSocket;
  socket_nodes = nodesSocket;
  socket_bStatesNeighbors = bStatesNeighborsSocket;
}

//////////////////////////////////////////////////////////////////////////////

void GlobalJacobianSparsity::computeBoundaryStatesFlag(std::valarray<bool>& isBoundaryState) const
{
  CFAUTOTRACE;

  // loop over all the boundary TRSs to detect all the boundary states
  vector< SafePtr<TopologicalRegionSet> > alltrs =
  MeshDataStack::getActive()->getTrsList();
  vector< Common::SafePtr<TopologicalRegionSet> >::const_iterator itrs;
  for (itrs = alltrs.begin(); itrs != alltrs.end(); ++itrs)
  {
    const SafePtr<TopologicalRegionSet> currTrs = *itrs;
    if (currTrs->hasTag("partition") || currTrs->hasTag("boundary"))
    {
      SafePtr<vector<CFuint> > bStates = currTrs->getStatesInTrs();
      for (CFuint i = 0; i < bStates->size(); ++i)
      {
        isBoundaryState[(*bStates)[i]] = true;
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void GlobalJacobianSparsity::computeMatrixPattern
(DataSocketSink<Framework::State*, Framework::GLOBAL> statesSocket,
 Common::ConnectivityTable<CFuint>& matrixPattern)
{
  throw Common::NotImplementedException(FromHere(), "GlobalJacobianSparsity::computeMatrixPattern()");
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
