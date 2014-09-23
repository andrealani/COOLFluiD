// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Common/CFLog.hh"
#include "Framework/TopologicalRegion.hh"
#include "Framework/PhysicalModel.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

TopologicalRegion::TopologicalRegion() :
  m_startGeoIdx(0),
  m_localNbGeoEnts(0),
  m_globalNbGeoEnts(0),
  m_localIdx(CFNULL),
  m_globalIdx(CFNULL),
  m_geo2states(CFNULL),
  m_geo2nodes(CFNULL)
{
}

//////////////////////////////////////////////////////////////////////////////

TopologicalRegion::~TopologicalRegion()
{
}

//////////////////////////////////////////////////////////////////////////////

void TopologicalRegion::putStatesInTR(vector<CFuint>& statesInTR) const
{
  CFAUTOTRACE;

  vector<CFuint> allStates;
  // preallocate some memory, avegNbPtrs is a just an average number
  // of pointers per geoentity
  const CFuint avegNbPtrs = PhysicalModelStack::getActive()->getDim()+1;
  allStates.reserve(getLocalNbGeoEnts()*avegNbPtrs);

  SafePtr<ConnectivityTable<CFuint> > table = m_geo2states;
  const CFuint nbGeos = m_startGeoIdx + getLocalNbGeoEnts();
  for (CFuint iGeo = m_startGeoIdx; iGeo < nbGeos; ++iGeo) {
    const CFuint nbStatesInGeo = table->nbCols(iGeo);
    for (CFuint is = 0; is < nbStatesInGeo; ++is) {
      const CFuint stateID = (*table)(iGeo,is);
      allStates.push_back(stateID);
    }
  }

  sort(allStates.begin(),allStates.end());

  unique_copy(allStates.begin(),
              allStates.end(),
              back_inserter(statesInTR));
}

//////////////////////////////////////////////////////////////////////////////

void TopologicalRegion::putNodesInTR(vector<CFuint>& nodesInTR) const
{
  CFAUTOTRACE;

  vector<CFuint> allNodes;
  // preallocate some memory, avegNbPtrs is a just an average number
  // of pointers per geoentity
  const CFuint avegNbPtrs = PhysicalModelStack::getActive()->getDim()+1;
  allNodes.reserve(getLocalNbGeoEnts()*avegNbPtrs);

  SafePtr<ConnectivityTable<CFuint> > table = m_geo2nodes;
  const CFuint nbGeos = m_startGeoIdx + getLocalNbGeoEnts();
  for (CFuint iGeo = m_startGeoIdx; iGeo < nbGeos; ++iGeo)
  {
    const CFuint nbNodesInGeo = table->nbCols(iGeo);
    for (CFuint is = 0; is < nbNodesInGeo; ++is)
    {
      const CFuint nodeID = (*table)(iGeo,is);
      allNodes.push_back(nodeID);
    }
  }

  sort(allNodes.begin(),allNodes.end());

  unique_copy(allNodes.begin(),
              allNodes.end(),
              back_inserter(nodesInTR));
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
