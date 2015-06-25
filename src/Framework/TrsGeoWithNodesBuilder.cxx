// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Framework/TrsGeoWithNodesBuilder.hh"
#include "Framework/MeshData.hh"
#include "Framework/GeometricEntityRegister.hh"
#include "Framework/BaseGeometricEntityProvider.hh"
#include "Framework/NamespaceSwitcher.hh"
#include "Framework/SubSystemStatus.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

TrsGeoWithNodesBuilder::TrsGeoWithNodesBuilder() :
  m_data(),
  m_nodes(CFNULL),
  m_poolData(),
  m_builtGeo(CFNULL),
  m_isSetup(false)
{
}

//////////////////////////////////////////////////////////////////////////////

TrsGeoWithNodesBuilder::~TrsGeoWithNodesBuilder()
{
  if ( m_isSetup ) unsetup();
}

//////////////////////////////////////////////////////////////////////////////
    
void TrsGeoWithNodesBuilder::unsetup()
{
  cf_assert ( m_isSetup );
  
  vector<GeometricEntity*>::iterator it;
  for (it = m_poolData.begin(); it != m_poolData.end(); ++it) {
    deletePtr(*it);
  }
  
  m_poolData.resize(0);
  
  m_isSetup = false;
}

//////////////////////////////////////////////////////////////////////////////

void TrsGeoWithNodesBuilder::setupInNamespace(const std::string& namespaceName)
{

  Common::SafePtr<Namespace> nsp = NamespaceSwitcher::getInstance
    (SubSystemStatusStack::getCurrentName()).getNamespace(namespaceName);
  m_nodes  = MeshDataStack::getInstance().getEntryByNamespace(nsp)->
    getNodeDataSocketSink().getDataHandle();

  // build the pool data
  GeometricEntityRegister& geoTypesReg =
    GeometricEntityRegister::getInstance();

  const CFuint nbGeoTypes = geoTypesReg.getSize();
  m_poolData.resize(nbGeoTypes);

  for (CFuint i = 0; i < nbGeoTypes; ++i) {
    BaseGeometricEntityProvider *const prov =  geoTypesReg.getProvider(i);
    GeometricEntity *const geo = prov->create();
    geo->resizeNodes(prov->getGeometryShapeFunctionNbNodes());
    geo->resizeStates(prov->getSolutionShapeFunctionNbNodes());
    m_poolData[i] = geo;
  }
}

//////////////////////////////////////////////////////////////////////////////

void TrsGeoWithNodesBuilder::setup()
{
  cf_assert(!m_isSetup);

  m_nodes  = MeshDataStack::getActive()->getNodeDataSocketSink().getDataHandle();

  // build the pool data
  GeometricEntityRegister& geoTypesReg =
    GeometricEntityRegister::getInstance();

  const CFuint nbGeoTypes = geoTypesReg.getSize();
  m_poolData.resize(nbGeoTypes);

  for (CFuint i = 0; i < nbGeoTypes; ++i) {
    BaseGeometricEntityProvider *const prov =  geoTypesReg.getProvider(i);
    GeometricEntity *const geo = prov->create();
    geo->resizeNodes(prov->getGeometryShapeFunctionNbNodes());
    geo->resizeStates(prov->getSolutionShapeFunctionNbNodes());
    m_poolData[i] = geo;
  }

  m_isSetup = true;
}

//////////////////////////////////////////////////////////////////////////////

GeometricEntity* TrsGeoWithNodesBuilder::buildGE()
{
  // cache the data
  const CFuint idx = m_data.idx;
  const TopologicalRegionSet& trs = *m_data.trs;

  // local ID (in all the mesh) of the geometric entity to get
  const CFuint geoType    = trs.getGeoType(idx);
  GeometricEntity *const geo = m_poolData[geoType];

  // set the local ID
  geo->setID(trs.getLocalGeoID(idx));

  // set the nodes
  const CFuint nbGeoNodes = trs.getNbNodesInGeo(idx);
  for (CFuint in = 0; in < nbGeoNodes; ++in) {
    const CFuint nodeID = trs.getNodeID(idx, in);
    geo->setNode(in, m_nodes[nodeID]);
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
