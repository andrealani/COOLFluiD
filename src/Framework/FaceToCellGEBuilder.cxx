// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Common/NotImplementedException.hh"

#include "Framework/MeshData.hh"
#include "Framework/GeometricEntityRegister.hh"
#include "Framework/BaseGeometricEntityProvider.hh"

#include "Framework/CFSide.hh"
#include "Framework/FaceToCellGEBuilder.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

FaceToCellGEBuilder::FaceToCellGEBuilder() :
  m_data(),
  m_states(CFNULL),
  m_nodes(CFNULL),
  m_PoolCells(),
  m_PoolFaces(),
  m_BuiltGeoCells(),
  m_BuiltGeoFaces(),
  m_isSetup(false)
{
  m_BuiltGeoCells.resize(2);
  m_BuiltGeoFaces.resize(1);
}

//////////////////////////////////////////////////////////////////////////////

FaceToCellGEBuilder::~FaceToCellGEBuilder()
{
  // clean the pool of faces
  std::map<CFuint,Framework::GeometricEntity*>::iterator itf;
  for (itf = m_PoolFaces.begin(); itf != m_PoolFaces.end(); ++itf)
  {
    deletePtr(itf->second);
  }

  // clean the pool of cells
  std::map<CFuint,
           std::pair<Framework::GeometricEntity*,Framework::GeometricEntity*>
          >::iterator itc;
  for (itc = m_PoolCells.begin(); itc != m_PoolCells.end(); ++itc)
  {
    deletePtr(itc->second.first);
    deletePtr(itc->second.second);
  }
}

//////////////////////////////////////////////////////////////////////////////

void FaceToCellGEBuilder::setup()
{
  cf_assert(!m_isSetup);

  // get cells-states connectivity
  m_CellToStates = MeshDataStack::getActive()->getConnectivity("cellStates_InnerCells");

  // get cells-nodes connectivity
  m_CellToNodes = MeshDataStack::getActive()->getConnectivity("cellNodes_InnerCells");

  // get datahandles to states and nodes
  m_states = MeshDataStack::getActive()->getStateDataSocketSink().getDataHandle();
  m_nodes  = MeshDataStack::getActive()->getNodeDataSocketSink().getDataHandle();

  // get the types of GeometricEntity
  GeometricEntityRegister& geoTypesReg = GeometricEntityRegister::getInstance();

  // loop over the GeometricEntity types
  // and create a GeometricEntity for each type
  // placing it in the pool
  for (CFuint i = 0; i < geoTypesReg.getSize(); ++i)
  {

    // create the GeometricEntity
    BaseGeometricEntityProvider *const prov =  geoTypesReg.getProvider(i);

    switch (prov->getGeomType())
    {

    case CFGeoEnt::CELL:
    {
      // create two cells of each type because it can happen that
      // the neighbour cells of a face are both of the same type
      GeometricEntity *const geoRight = prov->create();
      GeometricEntity *const geoLeft  = prov->create();

      geoRight->resizeNodes(prov->getGeometryShapeFunctionNbNodes());
      geoRight->resizeStates(prov->getSolutionShapeFunctionNbNodes());

      geoLeft->resizeNodes(prov->getGeometryShapeFunctionNbNodes());
      geoLeft->resizeStates(prov->getSolutionShapeFunctionNbNodes());
      
      std::pair<GeometricEntity*,GeometricEntity*> tp(geoLeft,geoRight);
      m_PoolCells[i] = tp;
      
    } break;

    case CFGeoEnt::FACE:
    {
      // resize the neighbors to have 2 neighbor cells
      // one in each side of the face
      GeometricEntity *const face = prov->create();
      face->resizeNodes(prov->getGeometryShapeFunctionNbNodes());
      face->resizeNeighborGeos(2);
      m_PoolFaces[i] = face;

    } break;

    case CFGeoEnt::EDGE:
    {
      throw Common::NotImplementedException (FromHere(),"Edges are currently not supported in FaceToCellGEBuilder");
    } break;

    case CFGeoEnt::INVALID:
    {
      throw Common::NotImplementedException (FromHere(),"Found geometry entity with no geometric type? Somethign has gone really bad.");
    } break;

    } // end switch CFGeoEnt::Type
  }

  m_isSetup = true;
}

//////////////////////////////////////////////////////////////////////////////

GeometricEntity* FaceToCellGEBuilder::buildGE()
{

  // cache the data
  const CFuint idx = m_data.idx;
  const TopologicalRegionSet& facesTRS = *m_data.facesTRS;
  const bool isNotBoundary = !m_data.isBoundary;

  m_FaceToCells = MeshDataStack::getActive()->getConnectivity(facesTRS.getName()+"-Faces2Cells");

  // Get the ID of the geometric entity type
  // should be a face
  const CFuint geoType = facesTRS.getGeoType(idx);
  m_tmpFace = m_PoolFaces[geoType];

  // set the local ID
  m_tmpFace->setID(facesTRS.getLocalGeoID(idx));

  // place the nodes in the Face
  const CFuint nbGeoNodes = facesTRS.getNbNodesInGeo(idx);
  for (CFuint in = 0; in < nbGeoNodes; ++in)
  {
    const CFuint nodeID = facesTRS.getNodeID(idx, in);
    m_tmpFace->setNode(in, m_nodes[nodeID]);
  }

  // there are two cells for each face
  cf_assert( (!isNotBoundary && m_FaceToCells->nbCols(idx) == 1) ||
          ( isNotBoundary && m_FaceToCells->nbCols(idx) == 2));

  // assemble left cell
  assembleNeighbourCell(LEFT);

  // if not a boundary face, assemble right cell
  if (isNotBoundary)
  {
    assembleNeighbourCell(RIGHT);
  }

  // keep track of the created GeometricEntities
  cf_assert(m_BuiltGeoFaces[0] == CFNULL);
  m_BuiltGeoFaces[0] = m_tmpFace;

  return m_tmpFace;
}

//////////////////////////////////////////////////////////////////////////////

void FaceToCellGEBuilder::assembleNeighbourCell(const CFuint side)
{
  const TopologicalRegionSet& cellsTRS = *m_data.cellsTRS;

  // cell ID and type
  const CFuint cellID = (*m_FaceToCells)(m_data.idx,side);
  const CFuint cellGeoType = cellsTRS.getGeoType(cellID);

  // get cell from the pool
  GeometricEntity *const cell  = side == LEFT ?
    m_PoolCells[cellGeoType ].first :
    m_PoolCells[cellGeoType ].second;

  // set cell ID
  cell->setID(cellID);

  // set the states in the cell
  const CFuint nbrStates = m_CellToStates->nbCols(cellID);
  for (CFuint is = 0; is < nbrStates; ++is)
  {
    const CFuint stateID = (*m_CellToStates)(cellID,is);
    cell->setState(is, m_states[stateID]);
  }

  // set the nodes in the cell
  const CFuint nbrNodes = m_CellToNodes->nbCols(cellID);
  for (CFuint in = 0; in < nbrNodes; ++in)
  {
    const CFuint nodeID = (*m_CellToNodes)(cellID,in);
    cell->setNode(in, m_nodes[nodeID]);
  }

  // set the neighbour cells in this GeometricEntity
  m_tmpFace->setNeighborGeo(side,cell);

  cf_assert(m_BuiltGeoCells[side] == CFNULL);
  m_BuiltGeoCells[side] = cell;
}

//////////////////////////////////////////////////////////////////////////////

void FaceToCellGEBuilder::releaseGE()
{
  fill(m_BuiltGeoCells.begin(),m_BuiltGeoCells.end(),static_cast<GeometricEntity*>(CFNULL));
  fill(m_BuiltGeoFaces.begin(),m_BuiltGeoFaces.end(),static_cast<GeometricEntity*>(CFNULL));
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
