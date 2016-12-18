// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Framework/StdTrsGeoBuilder.hh"
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

StdTrsGeoBuilder::StdTrsGeoBuilder() :
  _data(),
  _states(CFNULL),
  _nodes(CFNULL),
  _poolData(),
  _builtGeo(CFNULL),
  _isSetup(false),
    pass(0), 
    unpass(0)
{
}

//////////////////////////////////////////////////////////////////////////////

StdTrsGeoBuilder::~StdTrsGeoBuilder()
{
  if (_isSetup) unsetup();
}

//////////////////////////////////////////////////////////////////////////////

void StdTrsGeoBuilder::setupInNamespace(const std::string& namespaceName)
{
  cf_assert ( !_isSetup );
  
  Common::SafePtr<Namespace> nsp = NamespaceSwitcher::getInstance
    (SubSystemStatusStack::getCurrentName()).getNamespace(namespaceName);
  
  _states = MeshDataStack::getInstance().getEntryByNamespace(nsp)->
    getStateDataSocketSink().getDataHandle();
  
  _nodes = MeshDataStack::getInstance().getEntryByNamespace(nsp)->
    getNodeDataSocketSink().getDataHandle();
  
  // build the pool data
  GeometricEntityRegister& geoTypesReg =
    GeometricEntityRegister::getInstance();

  const CFuint nbGeoTypes = geoTypesReg.getSize();
  _poolData.resize(nbGeoTypes);

  for (CFuint i = 0; i < nbGeoTypes; ++i) {
    BaseGeometricEntityProvider *const prov =  geoTypesReg.getProvider(i);
    GeometricEntity *const geo = prov->create();
    geo->resizeNodes(prov->getGeometryShapeFunctionNbNodes());
    geo->resizeStates(prov->getSolutionShapeFunctionNbNodes());
    _poolData[i] = geo;
  }

  _isSetup = true;
}


//////////////////////////////////////////////////////////////////////////////

void StdTrsGeoBuilder::setup()
{
  cf_assert ( !_isSetup );
  
  _states = MeshDataStack::getActive()->getStateDataSocketSink().getDataHandle();
  _nodes  = MeshDataStack::getActive()->getNodeDataSocketSink().getDataHandle();
  
  // build the pool data
  GeometricEntityRegister& geoTypesReg = GeometricEntityRegister::getInstance();
  
  const CFuint nbGeoTypes = geoTypesReg.getSize();
  _poolData.resize(nbGeoTypes);
  
  for (CFuint i = 0; i < nbGeoTypes; ++i) 
  {
    BaseGeometricEntityProvider *const prov =  geoTypesReg.getProvider(i);
    CFLog(VERBOSE, "StdTrsGeoBuilder::setup() => registering provider " 
	  << prov->getProviderName() << "\n");
    GeometricEntity *const geo = prov->create();
    geo->resizeNodes(prov->getGeometryShapeFunctionNbNodes());
    geo->resizeStates(prov->getSolutionShapeFunctionNbNodes());
    _poolData[i] = geo;
  }
  
  _isSetup = true;
}

//////////////////////////////////////////////////////////////////////////////

void StdTrsGeoBuilder::unsetup()
{
  cf_assert ( _isSetup );
  
  // destroy the pool data
  std::vector<GeometricEntity*>::iterator it;
  for (it = _poolData.begin(); it != _poolData.end(); ++it) 
  {
    deletePtr(*it);
  }
  _poolData.resize(0);

  _isSetup = false;
}

//////////////////////////////////////////////////////////////////////////////

GeometricEntity* StdTrsGeoBuilder::buildGE()
{
  cf_assert ( _isSetup );
  
  // cache the data
  const CFuint idx = _data.idx;
  const TopologicalRegionSet& trs = *_data.trs;

  // local ID (in all the mesh) of the geometric entity to get
  const CFuint geoType = trs.getGeoType(idx);
  GeometricEntity *const geo = _poolData[geoType];

  // set the local ID
  geo->setID(trs.getLocalGeoID(idx));

  const CFuint nbGeoStates = trs.getNbStatesInGeo(idx);
  for (CFuint is = 0; is < nbGeoStates; ++is) {
    const CFuint stateID = trs.getStateID(idx, is);
    geo->setState(is, _states[stateID]);
  }

  const CFuint nbGeoNodes = trs.getNbNodesInGeo(idx);
  for (CFuint in = 0; in < nbGeoNodes; ++in) {
    const CFuint nodeID = trs.getNodeID(idx, in);
    geo->setNode(in, _nodes[nodeID]);
  }

  // keep track of the created GeometricEntity
  cf_assert(_builtGeo == CFNULL);
  _builtGeo = geo;
  return geo;
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
