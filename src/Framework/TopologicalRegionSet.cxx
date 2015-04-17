// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Framework/TopologicalRegionSet.hh"
#include "Framework/State.hh"
#include "Framework/IndexList.hh"

#include "Common/CFLog.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

TopologicalRegionSet::TopologicalRegionSet(const std::string& nameTRS,
                                           vector<TopologicalRegion*>* listTR)
  : NamedObject(nameTRS),
    TaggedObject(),
    m_listTR(listTR),
    m_geoEntLocalIdx(CFNULL),
    m_geoEntGlobalIdx(CFNULL),
    m_globalNbGeoEnts(0),
    m_geoTypes(CFNULL),
    m_geo2states(CFNULL),
    m_geo2nodes(CFNULL),
    m_statesList(0),
    m_nodesList(0),
    m_hasStatesList(false),
    m_hasNodesList(false)
{
  CFTRACEBEGIN;

  cf_assert(listTR != CFNULL);

  CFLog(NOTICE, "Building TRS: " << getName() << "\n");

  CFTRACEEND;
}

//////////////////////////////////////////////////////////////////////////////

TopologicalRegionSet::~TopologicalRegionSet()
{
  CFTRACEBEGIN;

  // deallocate the list of TRSs
  if(m_listTR != CFNULL) {
    vector<TopologicalRegion*>::iterator it;
    for (it = m_listTR->begin(); it != m_listTR->end(); ++it) {
      deletePtr(*it);
    }
  }
  deletePtr(m_listTR);

  /// @todo check this (use SharedPtr !!!)
  deletePtr(m_geoEntLocalIdx);
  deletePtr(m_geoEntGlobalIdx);
  deletePtr(m_geoTypes);
  
  CFTRACEEND;
}

//////////////////////////////////////////////////////////////////////////////

void TopologicalRegionSet::createStatesList()
{
  CFAUTOTRACE;
  if (!m_hasStatesList)
  {
    CFLogDebugMed("Creating state list in TRS [" + getName() + "]\n");
    cf_assert(m_listTR != CFNULL);

    const CFuint totalNbStates = IndexList<State>::getList().size();
    std::valarray<bool> isInserted(false, totalNbStates);
    vector<CFuint> temp;
    temp.reserve(totalNbStates);

    const CFuint nbGeos = getLocalNbGeoEnts();

    for (CFuint iGeo = 0; iGeo < nbGeos; ++iGeo)
    {
      const CFuint nbStatesInGeo = m_geo2states->nbCols(iGeo);

      for (CFuint is = 0; is < nbStatesInGeo; ++is)
      {
        const CFuint stateID = (*m_geo2states)(iGeo,is);
	
	if (stateID >= totalNbStates) {
	  CFLog(ERROR, "TopologicalRegionSet::createStatesList() => stateID " << 
		stateID << " > " << totalNbStates << "\n");
	  cf_assert(stateID < totalNbStates);
	}
	
        if (!isInserted[stateID])
        {
          temp.push_back(stateID);
          isInserted[stateID] = true;
        }
      }
    }

    const CFuint nbStatesInTrs = temp.size();
    m_statesList.resize(nbStatesInTrs);
    copy(temp.begin(), temp.end(), m_statesList.begin());
    m_hasStatesList = true;
    CFLogDebugMed("Finished state list in TRS [" + getName() + "]\n");
  }
}


//////////////////////////////////////////////////////////////////////////////

void TopologicalRegionSet::createNodesList()
{
  if (!m_hasNodesList) {
    CFLogDebugMed("TopologicalRegionSet()::createNodesList(): BEGIN\n");
    cf_assert(m_listTR != CFNULL);

    const CFuint totalNbNodes = IndexList<Node>::getList().size();
    std::valarray<bool> isInserted(false, totalNbNodes);
    vector<CFuint> temp;
    temp.reserve(totalNbNodes);

    const CFuint nbGeos = getLocalNbGeoEnts();

    for (CFuint iGeo = 0; iGeo < nbGeos; ++iGeo) {
      const CFuint nbNodesInGeo = m_geo2nodes->nbCols(iGeo);
      for (CFuint is = 0; is < nbNodesInGeo; ++is) {
        const CFuint nodeID = (*m_geo2nodes)(iGeo,is);
        cf_assert(nodeID < totalNbNodes);
        if (!isInserted[nodeID])
        {
          temp.push_back(nodeID);
          isInserted[nodeID] = true;
        }
      }
    }

    const CFuint nbNodesInTrs = temp.size();
    m_nodesList.resize(nbNodesInTrs);
    copy(temp.begin(), temp.end(), m_nodesList.begin());
    m_hasNodesList = true;
    CFLogDebugMed( "TopologicalRegionSet()::createNodesList(): END\n");
  }
}

//////////////////////////////////////////////////////////////////////////////

void TopologicalRegionSet::setGeo2NodesConn(Common::SafePtr<Common::ConnectivityTable<CFuint> > geo2nodes)
{
  m_geo2nodes = geo2nodes;
  // set the connectivity acquaintance also in the TR
  cf_assert(m_listTR != CFNULL);
  for (CFuint i = 0; i < m_listTR->size(); ++i) {
    (*m_listTR)[i]->setGeo2NodesConn(m_geo2nodes);
  }
}

//////////////////////////////////////////////////////////////////////////////

void TopologicalRegionSet::setGeo2StatesConn(Common::SafePtr<Common::ConnectivityTable<CFuint> > geo2states)
{
  m_geo2states = geo2states;
  // set the connectivity acquaintance also in the TR
  cf_assert(m_listTR != CFNULL);
  for (CFuint i = 0; i < m_listTR->size(); ++i) {
    (*m_listTR)[i]->setGeo2StatesConn(m_geo2states);
  }
}

//////////////////////////////////////////////////////////////////////////////


  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

