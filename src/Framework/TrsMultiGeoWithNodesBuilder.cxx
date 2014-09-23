// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "TrsMultiGeoWithNodesBuilder.hh"
#include "MeshData.hh"
#include "GeometricEntityRegister.hh"
#include "BaseGeometricEntityProvider.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

TrsMultiGeoWithNodesBuilder::TrsMultiGeoWithNodesBuilder() :
  _data(),
  _nodes(CFNULL),
  _poolData(),
  _countGeo(),
  _builtGeos(),
  _isSetup(false)
{
}

//////////////////////////////////////////////////////////////////////////////

TrsMultiGeoWithNodesBuilder::~TrsMultiGeoWithNodesBuilder()
{
  vector<vector<GeometricEntity*> >::iterator it;
  for (it = _poolData.begin(); it != _poolData.end(); ++it) {
    vector<GeometricEntity*>::iterator itg;
    for (itg = (*it).begin(); itg != (*it).end(); ++itg) {
      deletePtr(*itg);
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void TrsMultiGeoWithNodesBuilder::setup()
{
  cf_assert(!_isSetup);
  
  _nodes  = MeshDataStack::getActive()->getNodeDataSocketSink().getDataHandle();
  
  /// @todo resize and set _poolData ...
  _countGeo.resize(_poolData.size());

  // build the pool data
  GeometricEntityRegister& geoTypesReg =
    GeometricEntityRegister::getInstance();

  const CFuint nbGeoTypes = geoTypesReg.getSize();
  _poolData.resize(nbGeoTypes);

  const CFuint MAX_NB_GEOS_PER_TYPE = 10;
  for (CFuint i = 0; i < nbGeoTypes; ++i) {
    _poolData[i].resize(MAX_NB_GEOS_PER_TYPE);

    for (CFuint j = 0; j < MAX_NB_GEOS_PER_TYPE; ++j) {
      BaseGeometricEntityProvider *const prov =  geoTypesReg.getProvider(i);
      GeometricEntity *const geo = prov->create();
      geo->resizeNodes(prov->getGeometryShapeFunctionNbNodes());
      geo->resizeStates(prov->getSolutionShapeFunctionNbNodes());
      _poolData[i][j] = geo;
    }
  }

  // the maximum number of built geos must NEVER exceed the
  // maximum number of stored GeometricEntity's
  CFuint maxSize = 0;
  vector<vector<GeometricEntity*> >::iterator it;
  for (it = _poolData.begin(); it != _poolData.end(); ++it) {
    maxSize += (*it).size();
  }
  _builtGeos.reserve(maxSize);

  _isSetup = true;
}

//////////////////////////////////////////////////////////////////////////////

GeometricEntity* TrsMultiGeoWithNodesBuilder::buildGE()
{
  // cache the data
  const CFuint idx = _data.idx;
  const TopologicalRegionSet& trs = *_data.trs;

  // local ID (in all the mesh) of the geometric entity to get
  const CFuint geoType    = trs.getGeoType(idx);
  GeometricEntity *const geo = _poolData[geoType][_countGeo[geoType]++];

  // set the local ID
  geo->setID(trs.getLocalGeoID(idx));

  const CFuint nbGeoNodes = trs.getNbNodesInGeo(idx);
  for (CFuint in = 0; in < nbGeoNodes; ++in) {
    const CFuint nodeID = trs.getNodeID(idx, in);
    geo->setNode(in, _nodes[nodeID]);
  }

  _builtGeos.push_back(geo);
  return geo;
}

    //////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
