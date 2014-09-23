#include "Framework/BaseGeometricEntityProvider.hh"
#include "Framework/GeometricEntityRegister.hh"
#include "Framework/LocalConnectionData.hh"
#include "Framework/MeshData.hh"

#include "Framework/CFSide.hh"
#include "SpectralFV/CellToFaceGEBuilder.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace SpectralFV {

//////////////////////////////////////////////////////////////////////////////

CellToFaceGEBuilder::CellToFaceGEBuilder() :
  m_data(),
  m_cellToStates(CFNULL),
  m_cellToNodes(CFNULL),
  m_cellToFaces(CFNULL),
  m_intFaceToCells(CFNULL),
  m_facesToIntFaceIdxAndOrientOrBCType(CFNULL),
  m_states(CFNULL),
  m_nodes(CFNULL),
  m_poolCells(),
  m_poolFaces(),
  m_cellFaceNodes(),
  m_cellFaceGeoTypes(),
  m_builtGeoCells(),
  m_builtGeoFaces(),
  m_isFaceOnBoundary(),
  m_neighbrCellSide(),
  m_currentCellSide(),
  m_faceOrient(),
  m_faceBCIdx()
{
}

//////////////////////////////////////////////////////////////////////////////

CellToFaceGEBuilder::~CellToFaceGEBuilder()
{
  // clean the pool of faces
  map< CFuint , vector< Framework::GeometricEntity* > >::iterator itFaceType;
  for (itFaceType = m_poolFaces.begin(); itFaceType != m_poolFaces.end(); ++itFaceType)
  {
    vector< Framework::GeometricEntity* >::iterator itFace;
    for (itFace = itFaceType->second.begin(); itFace != itFaceType->second.end(); ++itFace)
    {
      deletePtr(*itFace);
    }
  }

  // clean the pool of cells
  map< CFuint , vector< Framework::GeometricEntity* > >::iterator itCellType;
  for (itCellType = m_poolCells.begin(); itCellType != m_poolCells.end(); ++itCellType)
  {
    vector< Framework::GeometricEntity* >::iterator itCell;
    for (itCell = itCellType->second.begin(); itCell != itCellType->second.end(); ++itCell)
    {
      deletePtr(*itCell);
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void CellToFaceGEBuilder::setup()
{
  // get cells-states connectivity
  m_cellToStates = MeshDataStack::getActive()->getConnectivity("cellStates_InnerCells");

  // get cells-nodes connectivity
  m_cellToNodes = MeshDataStack::getActive()->getConnectivity("cellNodes_InnerCells");

  // get cells-faces connectivity
  m_cellToFaces = MeshDataStack::getActive()->getConnectivity("InnerCells-Cells2Faces");

  // get internal faces-cells connectivity
  m_intFaceToCells = MeshDataStack::getActive()->getConnectivity("InnerFaces-Faces2Cells");

  // get the faces-internal face idx and orientation or external face boundary condition type connectivity
  m_facesToIntFaceIdxAndOrientOrBCType = MeshDataStack::getActive()->getConnectivity("faceToInFaceIdxOrientOrBCIdx");

  // get datahandles to states and nodes
  m_states = MeshDataStack::getActive()->getStateDataSocketSink().getDataHandle();
  m_nodes  = MeshDataStack::getActive()->getNodeDataSocketSink().getDataHandle();

  // get the types of GeometricEntity
  GeometricEntityRegister& geoTypesReg = GeometricEntityRegister::getInstance();

  // count maximum number of faces on a cell and set the nodes and geoType in each face
  CFuint maxNbrCellFaces = 0;
  for (CFuint i = 0; i < geoTypesReg.getSize(); ++i)
  {
    // create the GeometricEntity
    BaseGeometricEntityProvider *const prov =  geoTypesReg.getProvider(i);

    if (prov->getGeomType() == CFGeoEnt::CELL)
    {
      // cell shape
      const CFGeoShape::Type cellShape = prov->getShape();

      // geometric polynomial order
      const CFPolyOrder::Type geoPolyOrder = prov->getGeometryShapeFunctionOrder();

      // number of faces
      const CFuint nbrCellFaces = LocalConnectionData::getInstance().getNbFacesInShape(cellShape);
      maxNbrCellFaces = nbrCellFaces > maxNbrCellFaces ?
                        nbrCellFaces                   :
                        maxNbrCellFaces;

      // set nodes
      m_cellFaceNodes[i] = LocalConnectionData::getInstance().getFaceDofLocal(cellShape,geoPolyOrder,NODE,CFPolyForm::LAGRANGE);
      cf_assert(nbrCellFaces == m_cellFaceNodes[i]->nbRows());

      // set face geotypes
      m_cellFaceGeoTypes[i].resize(nbrCellFaces);
      for (CFuint iFace = 0; iFace < nbrCellFaces; ++iFace)
      {
        // get face shape
        const CFGeoShape::Type faceShape = LocalConnectionData::getInstance().getFaceShape(cellShape,iFace);

        // search face geoType
        for (CFuint j = 0; j < geoTypesReg.getSize(); ++j)
        {
          BaseGeometricEntityProvider *const prov2 =  geoTypesReg.getProvider(j);
          if (prov2->getGeomType() == CFGeoEnt::FACE)
          {
            if (faceShape == prov2->getShape() && geoPolyOrder == prov2->getGeometryShapeFunctionOrder())
            {
              m_cellFaceGeoTypes[i][iFace] = j;
            }
          }
        }
      }
    } // end switch CFGeoEnt::Type
  }

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
        // create enough cells of each type
        m_poolCells[i].resize(maxNbrCellFaces+1);
        for (CFuint iCell = 0; iCell < maxNbrCellFaces+1; ++iCell)
        {
          GeometricEntity *const cell = prov->create();

          cell->resizeNodes(prov->getGeometryShapeFunctionNbNodes());
          cell->resizeStates(prov->getSolutionShapeFunctionNbNodes());
          cell->resizeNeighborGeos(LocalConnectionData::getInstance().getNbFacesInShape(prov->getShape()));

          m_poolCells[i][iCell] = cell;
        }
      } break;
      case CFGeoEnt::FACE:
      {
        // create enough faces of each type
        m_poolFaces[i].resize(maxNbrCellFaces);
        for (CFuint iFace = 0; iFace < maxNbrCellFaces; ++iFace)
        {
          GeometricEntity *const face = prov->create();

          face->resizeNodes(prov->getGeometryShapeFunctionNbNodes());
          face->resizeNeighborGeos(2);

          m_poolFaces[i][iFace] = face;
        }
      } break;
      case CFGeoEnt::EDGE:
      {
        throw Common::NotImplementedException (FromHere(),"Edges are currently not supported in CellToFaceGEBuilder");
      } break;
      default:
      {
      }
    } // end switch CFGeoEnt::Type
  }

  // resize m_builtGeoCells and m_BuiltGeoFaces
  m_builtGeoCells.resize(maxNbrCellFaces+1);
  m_builtGeoFaces.resize(maxNbrCellFaces);
  m_isFaceOnBoundary.resize(maxNbrCellFaces);
  m_neighbrCellSide.resize(maxNbrCellFaces);
  m_currentCellSide.resize(maxNbrCellFaces);
  m_faceOrient.resize(maxNbrCellFaces);
  m_faceBCIdx.resize(maxNbrCellFaces);
}

//////////////////////////////////////////////////////////////////////////////

GeometricEntity* CellToFaceGEBuilder::buildGE()
{
  // cache the data
  const CFuint idx = m_data.idx;
  const TopologicalRegionSet& trs = *m_data.trs;

  // Get the ID of the geometric entity type (the index in geoTypesReg, see setup())
  // should be a cell
  const CFuint geoType = trs.getGeoType(idx);

  // get the number of faces
  const CFuint nbrFaces = m_cellFaceNodes[geoType]->nbRows();

  // get cell from the pool
  m_tmpCell = m_poolCells[geoType][nbrFaces];

  // set the local ID
  const CFuint cellID = trs.getLocalGeoID(idx);
  m_tmpCell->setID(cellID);

  // set the states in the Cell
  const CFuint nbrStates = m_cellToStates->nbCols(cellID);
  for (CFuint iState = 0; iState < nbrStates; ++iState)
  {
    const CFuint stateID = (*m_cellToStates)(cellID,iState);
    m_tmpCell->setState(iState, m_states[stateID]);
  }

  // set the nodes in the Cell
  const CFuint nbrNodes = m_cellToNodes->nbCols(cellID);
  for (CFuint iNode = 0; iNode < nbrNodes; ++iNode)
  {
    const CFuint nodeID = (*m_cellToNodes)(cellID,iNode);
    m_tmpCell->setNode(iNode, m_nodes[nodeID]);
  }

  // create neighbouring faces
  for (CFuint iFace = 0; iFace < nbrFaces; ++iFace)
  {
    assembleNeighbourFace(geoType,iFace);
  }

  // keep track of the created GeometricEntities
  cf_assert(m_builtGeoCells[nbrFaces] == CFNULL);
  m_builtGeoCells[nbrFaces] = m_tmpCell;

  return m_tmpCell;
}

//////////////////////////////////////////////////////////////////////////////

void CellToFaceGEBuilder::assembleNeighbourFace(const CFuint cellGeoType, const CFuint faceLocalID)
{
  // face geo type
  const CFuint faceGeoType = m_cellFaceGeoTypes[cellGeoType][faceLocalID];

  // get face from the pool
  GeometricEntity *const face = m_poolFaces[faceGeoType][faceLocalID];

  // set the local ID
  const CFuint faceID = (*m_cellToFaces)(m_data.idx,faceLocalID);
  face->setID(faceID);

  // set boundary face boolean
  const bool isFaceOnBoundary = m_facesToIntFaceIdxAndOrientOrBCType->nbCols(faceID) == 1;
  m_isFaceOnBoundary[faceLocalID] = isFaceOnBoundary;

  // boolean telling whether to invert the face node order
  bool invertFaceNodeOrder = false;

  // assemble neighbour cell or set the boundary face data
  if (isFaceOnBoundary)
  {
    // set the BC index
    m_faceBCIdx[faceLocalID] = (*m_facesToIntFaceIdxAndOrientOrBCType)(faceID,0);

    // set a large value in m_faceOrient, orientation is not needed for a boundary face
    m_faceOrient[faceLocalID] = 10000;

    // set the neighbouring cell in the face
    face->setNeighborGeo(0,m_tmpCell);
  }
  else
  {
    // get internal face index
    const CFuint intFaceIdx = (*m_facesToIntFaceIdxAndOrientOrBCType)(faceID,0);

    // get the face orientation
    m_faceOrient[faceLocalID] = (*m_facesToIntFaceIdxAndOrientOrBCType)(faceID,1);

    // determine wether to invert the face node order
    invertFaceNodeOrder = (*m_intFaceToCells)(intFaceIdx,LEFT) != m_data.idx;

    // determine cell sides with respect to face
    const CFuint neighbrCellSide = invertFaceNodeOrder ? LEFT : RIGHT;
    m_neighbrCellSide[faceLocalID] = neighbrCellSide;
    const CFuint currentCellSide = invertFaceNodeOrder ? RIGHT : LEFT;
    m_currentCellSide[faceLocalID] = currentCellSide;

    // get neighbouring cell
    const CFuint neighbourCellIdx = (*m_intFaceToCells)(intFaceIdx,neighbrCellSide);

    // neigbour cell geoType
    const CFuint neighbourCellGeoType = (m_data.trs)->getGeoType(neighbourCellIdx);

    // get cell from the pool
    GeometricEntity *const cell = m_poolCells[neighbourCellGeoType][faceLocalID];

    // set the local ID
    cell->setID(neighbourCellIdx);

    // set the neigbouring cell states
    const CFuint nbrStates = m_cellToStates->nbCols(neighbourCellIdx);
    for (CFuint iState = 0; iState < nbrStates; ++iState)
    {
      const CFuint stateID = (*m_cellToStates)(neighbourCellIdx,iState);
      cell->setState(iState, m_states[stateID]);
    }

    // set the neigbouring cell nodes
    const CFuint nbrNodes = m_cellToNodes->nbCols(neighbourCellIdx);
    for (CFuint iNode = 0; iNode < nbrNodes; ++iNode)
    {
      const CFuint nodeID = (*m_cellToNodes)(neighbourCellIdx,iNode);
      cell->setNode(iNode, m_nodes[nodeID]);
    }

    // set the neighbouring cells in the face
    face->setNeighborGeo(currentCellSide,m_tmpCell);
    face->setNeighborGeo(neighbrCellSide,cell     );

    cf_assert(m_builtGeoCells[faceLocalID] == CFNULL);
    m_builtGeoCells[faceLocalID] = face;

    // set a large value in m_bndFaceBCType, there is no BC in this face
    m_faceBCIdx[faceLocalID] = 10000;
  }

  if (invertFaceNodeOrder)
  {
    // set the nodes in the face
    Table< CFuint >& cellFaceNodes = *m_cellFaceNodes[cellGeoType];
    const CFuint nbrNodes = cellFaceNodes.nbCols(faceLocalID);
    for (CFuint iNode = 0; iNode < nbrNodes; ++iNode)
    {
      const CFuint nodeIdx = cellFaceNodes(faceLocalID,nbrNodes-1-iNode);
      face->setNode(iNode,m_tmpCell->getNode(nodeIdx));
    }
  }
  else
  {
    // set the nodes in the face
    Table< CFuint >& cellFaceNodes = *m_cellFaceNodes[cellGeoType];
    const CFuint nbrNodes = cellFaceNodes.nbCols(faceLocalID);
    for (CFuint iNode = 0; iNode < nbrNodes; ++iNode)
    {
      const CFuint nodeIdx = cellFaceNodes(faceLocalID,iNode);
      face->setNode(iNode,m_tmpCell->getNode(nodeIdx));
    }
  }

  // set the neighbour faces in the main cell
  m_tmpCell->setNeighborGeo(faceLocalID,face);

  cf_assert(m_builtGeoFaces[faceLocalID] == CFNULL);
  m_builtGeoFaces[faceLocalID] = face;
}

//////////////////////////////////////////////////////////////////////////////

void CellToFaceGEBuilder::releaseGE()
{
  fill(m_builtGeoCells.begin(),m_builtGeoCells.end(),static_cast<GeometricEntity*>(CFNULL));
  fill(m_builtGeoFaces.begin(),m_builtGeoFaces.end(),static_cast<GeometricEntity*>(CFNULL));
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace SpectralFV

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
