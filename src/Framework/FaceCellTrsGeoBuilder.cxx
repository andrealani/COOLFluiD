// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Framework/FaceCellTrsGeoBuilder.hh"
#include "Framework/MeshData.hh"
#include "Framework/GeometricEntityRegister.hh"
#include "Framework/BaseGeometricEntityProvider.hh"
#include "Framework/LocalConnectionData.hh"
#include "Framework/MapGeoToTrsAndIdx.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

FaceCellTrsGeoBuilder::FaceCellTrsGeoBuilder() : 
  CellTrsGeoBuilder(), 
  socket_cellFlag("Null"),
  m_fcdata()
{
}

//////////////////////////////////////////////////////////////////////////////

FaceCellTrsGeoBuilder::~FaceCellTrsGeoBuilder()
{
}

//////////////////////////////////////////////////////////////////////////////

GeometricEntity* FaceCellTrsGeoBuilder::buildGE()
{
  cf_assert(_builtGeos.size() == 0);
  cf_assert(_isSocketsSet);
  cf_assert(_isSetup);
  
  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();
  DataHandle<State*> gstates = socket_gstates.getDataHandle();
  DataHandle < Framework::Node*, Framework::GLOBAL > nodes = socket_nodes.getDataHandle();
  DataHandle < bool > cellFlag = socket_cellFlag.getDataHandle();
  
  // cache the data
  const CFuint faceidx = m_fcdata.idx;
  const bool isBFace = m_fcdata.isBFace;
  cf_assert(m_fcdata.faces.isNotNull());
  const TopologicalRegionSet& faces = *m_fcdata.faces;
  
  // local ID (in all the mesh) of the geometric entity to get
  const CFuint geoType = faces.getGeoType(faceidx);
  GeometricEntity *const currFace = _poolData[geoType][_countGeo[geoType]++];
  
  // set the local ID
  currFace->setID(faces.getLocalGeoID(faceidx));
  
  // set nodes in current face
  const CFuint nbGeoNodes = faces.getNbNodesInGeo(faceidx);
  for (CFuint in = 0; in < nbGeoNodes; ++in) {
    const CFuint nodeID = faces.getNodeID(faceidx, in);
    currFace->setNode(in, nodes[nodeID]);
  }
  
  // set cells on both sides of  the face
  for (CFuint iCell = 0; iCell < 2; ++iCell) {
    const CFuint stateID = faces.getStateID(faceidx, iCell);
    State* state = CFNULL;
    if (iCell == 0) {
      state = states[stateID];
    }
    else {
      state = (!isBFace) ? states[stateID] : gstates[stateID];
    }
    
    currFace->setState(iCell, state); 
    
    if (!state->isGhost()) {
      if (m_fcdata.allCells || (!cellFlag[stateID])) {
	const TopologicalRegionSet& cells = *m_fcdata.cells;
	
	// first create the cell
	// local ID (in all the mesh) of the geometric entity to get
	const CFuint cellType = cells.getGeoType(stateID); 
	GeometricEntity *const cell = _poolData[cellType][_countGeo[cellType]++];
	
	// set the local ID // here we assume cellID = stateID !!!!
	cell->setID(stateID);
	cf_assert(cells.getNbStatesInGeo(stateID) == 1);
	cell->setState(0, states[stateID]);
	
	const CFuint nbGeoNodes = cells.getNbNodesInGeo(stateID);
	for (CFuint in = 0; in < nbGeoNodes; ++in) {
	  const CFuint nodeID = cells.getNodeID(stateID, in);
	  cell->setNode(in, nodes[nodeID]);
	}
	
	// keep track of the created GeometricEntity
	_builtGeos.push_back(cell);
	
	// create all the neighbor faces for the current cell
	const CFuint nbFacesInCell = _cellFaces->nbCols(stateID);
	cf_assert(nbFacesInCell == cell->nbNeighborGeos());
	
	for (CFuint iFace = 0; iFace < nbFacesInCell; ++iFace)
	{      
	  const CFuint faceID = (*_cellFaces)(stateID, iFace);
	  const TopologicalRegionSet& faceTrs = *_mapGeoToTrs->getTrs(faceID);
	  const CFuint faceIdx = _mapGeoToTrs->getIdxInTrs(faceID);
	  const bool isBFace = _mapGeoToTrs->isBGeo(faceID);
	  
	  // local ID (in all the mesh) of the geometric entity to get
	  const CFuint faceType = faceTrs.getGeoType(faceIdx);
	  GeometricEntity *const face = _poolData[faceType][_countGeo[faceType]++];
	  
	  // set the local ID
	  face->setID(faceTrs.getLocalGeoID(faceIdx));
	  
	  // place the cell-centered states in those faces
	  const CFuint sID0 = faceTrs.getStateID(faceIdx, 0);
	  face->setState(0, states[sID0]);
	  
	  const CFuint sID1 = faceTrs.getStateID(faceIdx, 1);
	  State *const state1 = (!isBFace) ? states[sID1] : gstates[sID1];
	  face->setState(1, state1);
	  
	  // place the nodes in those faces
	  const CFuint nbFaceNodes = faceTrs.getNbNodesInGeo(faceIdx);
	  for (CFuint in = 0; in < nbFaceNodes; ++in) {
	    const CFuint nodeID = faceTrs.getNodeID(faceIdx, in);
	    face->setNode(in, nodes[nodeID]);
	  }
	  
	  cell->setNeighborGeo(iFace, face);
	  // keep track of the created GeometricEntity
	  _builtGeos.push_back(face);
	}
	
	currFace->setNeighborGeo(iCell, cell);  
      }
    }
  }
  
  _builtGeos.push_back(currFace);
  return currFace;
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
