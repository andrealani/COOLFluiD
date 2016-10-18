// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Common/NotImplementedException.hh"

#include "Framework/CellTrsGeoBuilder.hh"
#include "Framework/MeshData.hh"
#include "Framework/GeometricEntityRegister.hh"
#include "Framework/BaseGeometricEntityProvider.hh"
#include "Framework/LocalConnectionData.hh"
#include "Framework/MapGeoToTrsAndIdx.hh"
#include "Framework/NamespaceSwitcher.hh"
#include "Framework/SubSystemStatus.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

CellTrsGeoBuilder::CellTrsGeoBuilder() :
  _data(),
  socket_gstates("Null"),
  socket_nodes("Null"),
  socket_states("Null"),
  _mapGeoToTrs(CFNULL),
  _cellFaces(CFNULL),
  _poolData(),
  _countGeo(),
  _builtGeos(),
  _isSetup(false),
  _isSocketsSet(false)
{
}

//////////////////////////////////////////////////////////////////////////////

CellTrsGeoBuilder::~CellTrsGeoBuilder()
{
  if (_isSetup) unsetup();
}
  
//////////////////////////////////////////////////////////////////////////////
    
void CellTrsGeoBuilder::unsetup()
{
  cf_assert ( _isSetup );
  
  vector<vector<GeometricEntity*> >::iterator it;
  for (it = _poolData.begin(); it != _poolData.end(); ++it) {
    vector<GeometricEntity*>::iterator itg;
    for (itg = (*it).begin(); itg != (*it).end(); ++itg) {
      deletePtr(*itg);
    }
  }
  
  _poolData.resize(0);
  
  _isSetup = false;
}
    
//////////////////////////////////////////////////////////////////////////////

void CellTrsGeoBuilder::setDataSockets
(DataSocketSink<Framework::State*,Framework::GLOBAL> statesSocket,
   DataSocketSink<State*> gstatesSocket,
   DataSocketSink<Framework::Node*,Framework::GLOBAL> nodesSocket)
{
  
  socket_states = statesSocket;
  socket_gstates = gstatesSocket;
  socket_nodes = nodesSocket;

  _isSocketsSet = true;
}

//////////////////////////////////////////////////////////////////////////////

void CellTrsGeoBuilder::setup()
{
  cf_assert(!_isSetup);

  _mapGeoToTrs = MeshDataStack::getActive()->getMapGeoToTrs("MapFacesToTrs");
  _cellFaces = MeshDataStack::getActive()->getConnectivity("cellFaces");

  // build the pool data
  GeometricEntityRegister& geoTypesReg = GeometricEntityRegister::getInstance();
  
  const CFuint nbGeoTypes = geoTypesReg.getSize();
  _poolData.resize(nbGeoTypes);
  
  const CFuint MAX_NB_GEOS_PER_TYPE = 30;
  for (CFuint i = 0; i < nbGeoTypes; ++i)
  {
    _poolData[i].resize(MAX_NB_GEOS_PER_TYPE);

    BaseGeometricEntityProvider *const prov =  geoTypesReg.getProvider(i);

    for (CFuint j = 0; j < MAX_NB_GEOS_PER_TYPE; ++j)
    {

      GeometricEntity *const geo = prov->create();
      geo->resizeNodes(prov->getGeometryShapeFunctionNbNodes());

      switch (prov->getGeomType())
      {

      case CFGeoEnt::CELL:
      {
        geo->resizeStates(prov->getSolutionShapeFunctionNbNodes());
        const CFGeoShape::Type shape = prov->getShape();
        const CFuint nbFacesInCell = LocalConnectionData::getInstance().getNbFacesInShape(shape);
        geo->resizeNeighborGeos(nbFacesInCell);
      } break;

      case CFGeoEnt::FACE:
      {
        // here if it is a face allocate with 2 states else 1 !!
        // faces for cell centered FVM have always only left and right state
        geo->resizeStates(2);
	geo->resizeNeighborGeos(2);
      } break;
      
      case CFGeoEnt::EDGE:
      {
        throw Common::NotImplementedException(FromHere(),"Edges are currently not supported in CellTrsGeoBuilder");
      } break;

      case CFGeoEnt::INVALID:
      {
        throw Common::NotImplementedException(FromHere(),"Found geometry entity with no geometric type? Somethign has gone really bad.");
      } break;

      } // end switch CFGeoEnt::Type

      _poolData[i][j] = geo;
    }
  }

  _countGeo.resize(_poolData.size());
  _countGeo = 0;

  // the maximum number of built geos must NEVER exceed the
  // maximum number of stored GeometricEntity's
  CFuint maxSize = 0;
  vector<vector<GeometricEntity*> >::iterator it;
  for (it = _poolData.begin(); it != _poolData.end(); ++it) {
    maxSize += (*it).size();
  }
  
  _builtGeos.reserve(maxSize);
  cf_assert(maxSize > 0);
  _isSetup = true;
}

//////////////////////////////////////////////////////////////////////////////

void CellTrsGeoBuilder::setupInNamespace(const std::string& namespaceName)
{

  Common::SafePtr<Namespace> nsp = NamespaceSwitcher::getInstance
    (SubSystemStatusStack::getCurrentName()).getNamespace(namespaceName);
  
  _mapGeoToTrs = MeshDataStack::getInstance().getEntryByNamespace(nsp)->getMapGeoToTrs("MapFacesToTrs");
  _cellFaces = MeshDataStack::getInstance().getEntryByNamespace(nsp)->getConnectivity("cellFaces");

  // build the pool data
  GeometricEntityRegister& geoTypesReg =
    GeometricEntityRegister::getInstance();

  const CFuint nbGeoTypes = geoTypesReg.getSize();
  _poolData.resize(nbGeoTypes);

  const CFuint MAX_NB_GEOS_PER_TYPE = 30;
  for (CFuint i = 0; i < nbGeoTypes; ++i)
  {
    _poolData[i].resize(MAX_NB_GEOS_PER_TYPE);
    BaseGeometricEntityProvider *const prov =  geoTypesReg.getProvider(i);

    for (CFuint j = 0; j < MAX_NB_GEOS_PER_TYPE; ++j)
    {

      GeometricEntity *const geo = prov->create();
      geo->resizeNodes(prov->getGeometryShapeFunctionNbNodes());

      switch (prov->getGeomType())
      {

      case CFGeoEnt::CELL:
      {
        geo->resizeStates(prov->getSolutionShapeFunctionNbNodes());
        const CFGeoShape::Type shape = prov->getShape();
        const CFuint nbFacesInCell = LocalConnectionData::getInstance().getNbFacesInShape(shape);
        geo->resizeNeighborGeos(nbFacesInCell);
      } break;

      case CFGeoEnt::FACE:
      {
        // here if it is a face allocate with 2 states else 1 !!
        // faces for cell centered FVM have always only left and right state
        geo->resizeStates(2);
	geo->resizeNeighborGeos(2);
      } break;

      case CFGeoEnt::EDGE:
      {
        throw Common::NotImplementedException (FromHere(),"Edges are currently not supported in CellTrsGeoBuilder");
      } break;

      case CFGeoEnt::INVALID:
      {
        throw Common::NotImplementedException (FromHere(),"Found geometry entity with no geometric type? Somethign has gone really bad.");
      } break;

      } // end switch CFGeoEnt::Type

      _poolData[i][j] = geo;
    }
  }

  _countGeo.resize(_poolData.size());
  _countGeo = 0;

  // the maximum number of built geos must NEVER exceed the
  // maximum number of stored GeometricEntity's
  CFuint maxSize = 0;
  vector<vector<GeometricEntity*> >::iterator it;
  for (it = _poolData.begin(); it != _poolData.end(); ++it) {
    maxSize += (*it).size();
  }
  
  _builtGeos.reserve(maxSize);
  cf_assert(maxSize > 0);
  _isSetup = true;
}

//////////////////////////////////////////////////////////////////////////////

GeometricEntity* CellTrsGeoBuilder::buildGE()
{
  cf_assert(_builtGeos.size() == 0);
  cf_assert(_isSocketsSet);
  cf_assert(_isSetup);
  
  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();
  DataHandle<State*> gstates = socket_gstates.getDataHandle();
  DataHandle < Framework::Node*, Framework::GLOBAL > nodes = socket_nodes.getDataHandle();

  // cache the data
  const CFuint idx = _data.idx;
  cf_assert(_data.trs.isNotNull());
  const TopologicalRegionSet& trs = *_data.trs;
  
  // first create the cell
  // local ID (in all the mesh) of the geometric entity to get
  const CFuint geoType = trs.getGeoType(idx);
  GeometricEntity *const cell = _poolData[geoType][_countGeo[geoType]++];
  
  // set the local ID
  cell->setID(trs.getLocalGeoID(idx));
  
  cf_assert(trs.getNbStatesInGeo(idx) == 1);
  const CFuint stateID0 = trs.getStateID(idx, 0);
  cell->setState(0, states[stateID0]);
  
  const CFuint nbGeoNodes = trs.getNbNodesInGeo(idx);
  for (CFuint in = 0; in < nbGeoNodes; ++in) {
    const CFuint nodeID = trs.getNodeID(idx, in);
    cell->setNode(in, nodes[nodeID]);
  }
  
  // keep track of the created GeometricEntity
  _builtGeos.push_back(cell);

  // create all the neighbor faces for the current cell
  const CFuint nbFacesInCell = _cellFaces->nbCols(idx);
  cf_assert(nbFacesInCell == cell->nbNeighborGeos());
  
  CFLog(DEBUG_MAX, "CellTrsGeoBuilder::buildGE() => nbFacesInCell = " << nbFacesInCell << "\n");
  
  for (CFuint iFace = 0; iFace < nbFacesInCell; ++iFace)
  {
    CFLog(DEBUG_MAX, "CellTrsGeoBuilder::buildGE() => iFace = " << iFace << "\n");
    const CFuint faceID = (*_cellFaces)(idx, iFace);
    CFLog(DEBUG_MAX, "CellTrsGeoBuilder::buildGE() => faceID = " << faceID << "\n");
    const TopologicalRegionSet& faceTrs = *_mapGeoToTrs->getTrs(faceID);
    const CFuint faceIdx = _mapGeoToTrs->getIdxInTrs(faceID);
    CFLog(DEBUG_MAX, "CellTrsGeoBuilder::buildGE() => faceIdx = " << faceIdx << "\n");
    const bool isBFace = _mapGeoToTrs->isBGeo(faceID);
    // local ID (in all the mesh) of the geometric entity to get
    const CFuint faceType = faceTrs.getGeoType(faceIdx);
    CFLog(DEBUG_MAX, "CellTrsGeoBuilder::buildGE() => faceType = " << faceType << "\n");
    GeometricEntity *const face = _poolData[faceType][_countGeo[faceType]++];
    // set the local ID
    face->setID(faceTrs.getLocalGeoID(faceIdx));
    
    // place the cell-centered states in those faces
    const CFuint sID0 = faceTrs.getStateID(faceIdx, 0);
    CFLog(DEBUG_MAX, "CellTrsGeoBuilder::buildGE() => sID0 = " << sID0 << "\n");
    face->setState(0, states[sID0]);
    
    const CFuint sID1 = faceTrs.getStateID(faceIdx, 1);
    CFLog(DEBUG_MAX, "CellTrsGeoBuilder::buildGE() => sID1 = " << sID1 << "\n");
    CFLog(DEBUG_MAX, "CellTrsGeoBuilder::buildGE() => isBFace = " << isBFace << "\n");
    CFLog(DEBUG_MAX, "CellTrsGeoBuilder::buildGE() => states.size() = " << states.size() << "\n");
    CFLog(DEBUG_MAX, "CellTrsGeoBuilder::buildGE() => gstates.size() = " << gstates.size() << "\n");
    cf_assert((!isBFace && sID1 < states.size()) || (isBFace && sID1 < gstates.size()));
    State *const state1 = (!isBFace) ? states[sID1] : gstates[sID1];
    face->setState(1, state1);
    
    // place the nodes in those faces
    const CFuint nbFaceNodes = faceTrs.getNbNodesInGeo(faceIdx);
    CFLog(DEBUG_MAX, "CellTrsGeoBuilder::buildGE() => nbFaceNodes = " << nbFaceNodes << "\n");
    for (CFuint in = 0; in < nbFaceNodes; ++in) {
      const CFuint nodeID = faceTrs.getNodeID(faceIdx, in);
      CFLog(DEBUG_MAX, "CellTrsGeoBuilder::buildGE() => nodeID = " << nodeID << "\n");
      face->setNode(in, nodes[nodeID]);
    }

    cell->setNeighborGeo(iFace, face);
    // keep track of the created GeometricEntity
    _builtGeos.push_back(face);

    // cout << "trs.getName() = " << trs.getName() << endl;
    // cout << "trs.getLocalNbGeoEnts() = " << trs.getLocalNbGeoEnts() << endl;
    // cout << "trs.getNbStatesInGeo() = " << trs.getNbStatesInGeo(faceIdx) << endl;
    // cout << "trs.getNbNodesInGeo() = " << trs.getNbNodesInGeo(faceIdx) << endl;
    // cout << "faceIdx = " << faceIdx << endl;
  }

  return cell;
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
