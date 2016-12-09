// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Framework/FaceTrsGeoBuilder.hh"
#include "Framework/MeshData.hh"
#include "Framework/GeometricEntityRegister.hh"
#include "Framework/BaseGeometricEntityProvider.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

FaceTrsGeoBuilder::FaceTrsGeoBuilder() :
  m_data(),
  socket_gstates("Null"),
  socket_nodes("Null"),
  socket_states("Null"),
  m_poolData(),
  m_builtGeo(CFNULL),
  m_isSetup(false),
  m_isSocketsSet(false)
{
}

//////////////////////////////////////////////////////////////////////////////

FaceTrsGeoBuilder::~FaceTrsGeoBuilder()
{
  if ( m_isSetup ) unsetup();
}

//////////////////////////////////////////////////////////////////////////////

void FaceTrsGeoBuilder::setDataSockets
(DataSocketSink<Framework::State*, Framework::GLOBAL> statesSocket,
   DataSocketSink<State*> gstatesSocket,
 DataSocketSink<Framework::Node*, Framework::GLOBAL> nodesSocket)
{

  socket_states = statesSocket;
  socket_gstates = gstatesSocket;
  socket_nodes = nodesSocket;

  m_isSocketsSet = true;
}

//////////////////////////////////////////////////////////////////////////////

void FaceTrsGeoBuilder::setup()
{
  cf_assert(!m_isSetup);

  // build the pool data
  GeometricEntityRegister& geoTypesReg =
    GeometricEntityRegister::getInstance();

  const CFuint nbGeoTypes = geoTypesReg.getSize();
  m_poolData.resize(nbGeoTypes);

  for (CFuint i = 0; i < nbGeoTypes; ++i) {
    BaseGeometricEntityProvider *const prov =  geoTypesReg.getProvider(i);
    GeometricEntity *const geo = prov->create();
    geo->resizeNodes(prov->getGeometryShapeFunctionNbNodes());
    // faces for cell centered FVM have always only left and right state
    geo->resizeStates(2);
    m_poolData[i] = geo;
  }

  m_isSetup = true;
}

//////////////////////////////////////////////////////////////////////////////

void FaceTrsGeoBuilder::unsetup()
{
  cf_assert ( m_isSetup );
  
  vector<GeometricEntity*>::iterator it;
  for (it = m_poolData.begin(); it != m_poolData.end(); ++it) 
    deletePtr(*it);

  m_poolData.resize(0);
        
  m_isSetup = false;
}

//////////////////////////////////////////////////////////////////////////////

GeometricEntity* FaceTrsGeoBuilder::buildGE()
{
  cf_assert(m_isSetup);
  cf_assert(m_isSocketsSet);

  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();
  DataHandle<State*> gstates = socket_gstates.getDataHandle();
  DataHandle < Framework::Node*, Framework::GLOBAL > nodes = socket_nodes.getDataHandle();

  // cache the data
  const CFuint idx = m_data.idx;
  const bool isBFace = m_data.isBFace;
  cf_assert(m_data.trs.isNotNull());
  const TopologicalRegionSet& trs = *m_data.trs;

  // local ID (in all the mesh) of the geometric entity to get
  const CFuint geoType = trs.getGeoType(idx);
  GeometricEntity *const geo = m_poolData[geoType];

  // set the local ID
  geo->setID(trs.getLocalGeoID(idx));
  const CFuint stateID0 = trs.getStateID(idx, 0);
  geo->setState(0, states[stateID0]);

  const CFuint stateID1 = trs.getStateID(idx, 1);
  CFLog(DEBUG_MAX, "FaceTrsGeoBuilder::buildGE() => isBFace[" << isBFace 
	<< "] => stateID(0,1) = (" << stateID0 << ", " << stateID1 << "), " 
	<< states.size() << " - " << gstates.size() << "\n");
  State *const state1 = (!isBFace) ? states[stateID1] : gstates[stateID1];
  geo->setState(1, state1);
  
  const CFuint nbGeoNodes = trs.getNbNodesInGeo(idx);
  for (CFuint in = 0; in < nbGeoNodes; ++in) {
    const CFuint nodeID = trs.getNodeID(idx, in);
    geo->setNode(in, nodes[nodeID]);
  }

  // keep track of the created GeometricEntity
  cf_assert(m_builtGeo == CFNULL);
  m_builtGeo = geo;
  return geo;
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
