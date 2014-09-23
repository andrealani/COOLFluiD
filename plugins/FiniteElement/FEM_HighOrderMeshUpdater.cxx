#include <ctime>
#include <set>


#include "Environment/ObjectProvider.hh"
#include "Common/CFLog.hh"
#include "Common/BadValueException.hh"
#include "Framework/MeshData.hh"
#include "Framework/LocalConnectionData.hh"
#include "Framework/GeometricEntityProvider.hh"
#include "Framework/Face.hh"
#include "Framework/Cell.hh"
#include "Framework/PhysicalModel.hh"
#include "Framework/BadFormatException.hh"
#include "Framework/MapGeoToTrsAndIdx.hh"
#include "FiniteElement/FiniteElement.hh"
#include "FiniteElement/FEM_HighOrderMeshUpdater.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteElement {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<FEM_HighOrderMeshUpdater,
                            MeshDataBuilder,
                            FiniteElementModule,
                            1>
FEM_HighOrderMeshUpdaterProvider("FiniteElementHO");

//////////////////////////////////////////////////////////////////////////////
void FEM_HighOrderMeshUpdater::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOption< CFuint >("order","Order of the new mesh");
}
//////////////////////////////////////////////////////////////////////////////

FEM_HighOrderMeshUpdater::FEM_HighOrderMeshUpdater(const std::string& name) :
  FEM_MeshDataBuilder(name)
{
 addConfigOptionsTo(this);
  m_newPolyOrder_int = 2;
 setParameter("order",&m_newPolyOrder_int);

  m_newPolyOrder = static_cast<COOLFluiD::CFPolyOrder::Type>(m_newPolyOrder_int);

}

//////////////////////////////////////////////////////////////////////////////

FEM_HighOrderMeshUpdater::~FEM_HighOrderMeshUpdater()
{
}

//////////////////////////////////////////////////////////////////////////////

void FEM_HighOrderMeshUpdater::releaseMemory()
{
  FEM_MeshDataBuilder::releaseMemory();
}

//////////////////////////////////////////////////////////////////////////////

void FEM_HighOrderMeshUpdater::computeGeoTypeInfo()
{
  CFAUTOTRACE;

  //cf_assert(m_newPolyOrder.isNotNull());

  // get the element type data
  SafePtr< std::vector< ElementTypeData > > elementType = getCFmeshData().getElementTypeData();

  std::vector< ElementTypeData >::iterator type_itr = elementType->begin();
  for (; type_itr != elementType->end(); ++type_itr)
  {
    // get the number of control volumes (states) in this element type
    const CFuint nbStatesPerElem = getNbStatesInNewType(type_itr->getGeoShape(),m_newPolyOrder);


    type_itr->setNbStates(nbStatesPerElem);
    type_itr->setSolOrder(m_newPolyOrder);
  }

  /// @todo since only meshes with the same order are currently supported
  ///       we have to change the order in the CFmeshData and
  ///       not only on the ElementTypeData
  getCFmeshData().setSolutionPolyOrder(m_newPolyOrder);

  // continue with the standard algorithm
  FEM_MeshDataBuilder::computeGeoTypeInfo();
}


//////////////////////////////////////////////////////////////////////////////

CFuint FEM_HighOrderMeshUpdater::getNbStatesInNewType(const CFGeoShape::Type geoShape, const CFPolyOrder::Type polyOrder)
{
  cf_assert(polyOrder == 2);

  switch(polyOrder)
  {
  case CFPolyOrder::ORDER2:
      switch (geoShape)
     {
      case CFGeoShape::LINE:
        return 3;
      case CFGeoShape::TRIAG:
        return 6;
      case CFGeoShape::QUAD:
        return 9;
      case CFGeoShape::TETRA:
        return 10;
      case CFGeoShape::PRISM:
        throw Common::NotImplementedException (FromHere(),"Element type PRISM not implemented...");
      case CFGeoShape::PYRAM:
        throw Common::NotImplementedException (FromHere(),"Element type PYRAM not implemented...");
      case CFGeoShape::HEXA:
        return 20;
      default:
        throw Common::NotImplementedException (FromHere(),"Element type not found...");
    }
  break;
  case CFPolyOrder::ORDER3:
   switch (geoShape)
     {
      case CFGeoShape::LINE:
        return 4;
      case CFGeoShape::TRIAG:
        return 9;
      case CFGeoShape::QUAD:
        throw Common::NotImplementedException (FromHere(),"Element type QUAD not implemented...");
      case CFGeoShape::TETRA:
        throw Common::NotImplementedException (FromHere(),"Element type TETRA not implemented...");
      case CFGeoShape::PRISM:
        throw Common::NotImplementedException (FromHere(),"Element type PRISM not implemented...");
      case CFGeoShape::PYRAM:
        throw Common::NotImplementedException (FromHere(),"Element type PYRAM not implemented...");
      case CFGeoShape::HEXA:
        throw Common::NotImplementedException (FromHere(),"Element type HEXA not implemented...");
      default:
        throw Common::NotImplementedException (FromHere(),"Element type not found...");
    }
  break;
  default: { throw Common::NotImplementedException (FromHere(),"PolyOrder not handled" + StringOps::to_str( (CFuint) polyOrder )); break; }
  }

  return 0; // just to avoid warnings
}

//////////////////////////////////////////////////////////////////////////////

CFuint FEM_HighOrderMeshUpdater::getNbInteriorStatesInNewType(const CFGeoShape::Type geoShape, const CFPolyOrder::Type polyOrder)
{

switch(polyOrder)
  {
  case CFPolyOrder::ORDER2:
  switch (geoShape)
  {
    case CFGeoShape::LINE:
      return 1;
    case CFGeoShape::TRIAG:
      return 0;
    case CFGeoShape::QUAD:
      return 1;
    case CFGeoShape::TETRA:
      return 0;
    case CFGeoShape::PRISM:
      throw Common::NotImplementedException (FromHere(),"Element type PRISM not implemented...");
    case CFGeoShape::PYRAM:
      throw Common::NotImplementedException (FromHere(),"Element type PYRAM not implemented...");
    case CFGeoShape::HEXA:
      return 0;
    default:
      throw Common::NotImplementedException (FromHere(),"Element type not found...");
  }
  case CFPolyOrder::ORDER3:
  switch (geoShape)
  {
    case CFGeoShape::LINE:
      return 2;
    case CFGeoShape::TRIAG:
      return 1;
     case CFGeoShape::QUAD:
        throw Common::NotImplementedException (FromHere(),"Element type QUAD not implemented...");
      case CFGeoShape::TETRA:
        throw Common::NotImplementedException (FromHere(),"Element type TETRA not implemented...");
      case CFGeoShape::PRISM:
        throw Common::NotImplementedException (FromHere(),"Element type PRISM not implemented...");
      case CFGeoShape::PYRAM:
        throw Common::NotImplementedException (FromHere(),"Element type PYRAM not implemented...");
      case CFGeoShape::HEXA:
        throw Common::NotImplementedException (FromHere(),"Element type HEXA not implemented...");
      default:
        throw Common::NotImplementedException (FromHere(),"Element type not found...");
  }
  break;
  default: { throw Common::NotImplementedException (FromHere(),"PolyOrder not handled" + StringOps::to_str( (CFuint) polyOrder )); break; }
}

  return 0; // just to avoid warnings
}

//////////////////////////////////////////////////////////////////////////////

CFuint FEM_HighOrderMeshUpdater::getNewNbStatesInBoundaryGeo(const CFuint oldNbStates, const CFPolyOrder::Type polyOrder)
{

  const CFuint dim = getCFmeshData().getDimension();

  switch(polyOrder)
  {
  case CFPolyOrder::ORDER2:
    if(dim == DIM_2D)
      {
        cf_assert(oldNbStates == 2);
        return 3;
      }
      else
      {
        cf_assert(dim == DIM_3D);
        switch (oldNbStates)
      {
          case 3: // we have a triangular face
            return 6;
          case 4: // we have a quad face
            return 8;
          default:
            throw Common::NotImplementedException (FromHere(),"Face type not found...");
        }
      }
  break;
  case CFPolyOrder::ORDER3:
    if(dim == DIM_2D)
      {
        cf_assert(oldNbStates == 2);
        return 4;
      }
      else
      {
        cf_assert(dim == DIM_3D);
        switch (oldNbStates)
      {
          case 3: // we have a triangular face
            return 9;
          case 4: // we have a quad face
          throw Common::NotImplementedException (FromHere(),"Element type PRISM not implemented...");
          default:
            throw Common::NotImplementedException (FromHere(),"Face type not found...");
        }
      }
  break;
  default: { throw Common::NotImplementedException (FromHere(),"PolyOrder not handled" + StringOps::to_str( (CFuint) polyOrder )); break; }
}

  return 0; // just to avoid warnings
}

//////////////////////////////////////////////////////////////////////////////

void FEM_HighOrderMeshUpdater::createTopologicalRegionSets()
{
  CFAUTOTRACE;

  CFout << "Converting the P1P1 mesh to a P1Pk\n";
  // first transform the cell-states connectivity
  // from P1P1 to P1P2
  upgradeStateConnectivity();

  // recreate the states
  recreateStates();

  // modify the TRS states
  updateTRSData();

  // continue with the standard algorithm
  FEM_MeshDataBuilder::createTopologicalRegionSets();
}

//////////////////////////////////////////////////////////////////////////////

void FEM_HighOrderMeshUpdater::upgradeStateConnectivity()
{
  CFAUTOTRACE;

CFout << "Updating existing states connectivity\n";
  // get the element type data
  SafePtr< std::vector< ElementTypeData > > elementType = getCFmeshData().getElementTypeData();
  std::vector< ElementTypeData >::iterator type_itr;

  // get the nodal connectivity
  Common::SafePtr<MeshData::ConnTable> cellNodes = MeshDataStack::getActive()->getConnectivity("cellNodes_InnerCells");

  // get the connectivity that we will change
  Common::SafePtr<MeshData::ConnTable> cellStates = MeshDataStack::getActive()->getConnectivity("cellStates_InnerCells");

  const CFuint nbElems = getNbElements();

  // we will create a mesh with equal order on all the elements
  // and we reset the connectivity accordingly, using a std::valarray
  // first create the std::valarray
  std::valarray<CFuint> newPattern(nbElems);
  for (type_itr = elementType->begin(); type_itr != elementType->end(); ++type_itr)
  {
    // get the number of control volumes (states) in this element type
    const CFuint newNbStatesPerElem = type_itr->getNbStates();

    // get the number of elements
    const CFuint nbElemPerType = type_itr->getNbElems();

    // get start index of this element type in global element list
    CFuint globalIdx = type_itr->getStartIdx();

    // loop over elements of this type
    for (CFuint iElem = 0; iElem < nbElemPerType; ++iElem, ++globalIdx)
    {
      newPattern[globalIdx] = newNbStatesPerElem;
    }
  }
  // then resize the connectivity
  cellStates->resize(newPattern);

  // loop on all the elements and reassign the IDs of the nodal states
  CFuint maxStateID = 0;
  for (type_itr = elementType->begin(); type_itr != elementType->end(); ++type_itr)
  {
    // get the number of nodes in this element type
    const CFuint nbNodesPerElem = type_itr->getNbNodes();

    // get the number of elements
    const CFuint nbElemPerType = type_itr->getNbElems();

    // get start index of this element type in global element list
    CFuint globalIdx = type_itr->getStartIdx();

    // loop over elements of this type and set the P1P1 states
    for (CFuint iElem = 0; iElem < nbElemPerType; ++iElem, ++globalIdx)
    {
// CFout << "Cell: " << globalIdx << "\n";
      for (CFuint jState = 0; jState < nbNodesPerElem; ++jState)
      {
        (*cellStates)(globalIdx,jState) = (*cellNodes)(globalIdx,jState);
// CFout << "Recopying the state: " << (*cellStates)(globalIdx, jState) <<"\n";
        maxStateID = std::max((*cellStates)(globalIdx,jState),maxStateID);
// CFout << "Checking nodal state: " << (*cellStates)(globalIdx,jState) << "\n";
      }
    }
  }

  maxStateID++;
  setExtraStates(maxStateID);

  _totalNewNbStates = maxStateID;
//  cf_assert(maxStateID == newPattern.sum());
}

//////////////////////////////////////////////////////////////////////////////

void FEM_HighOrderMeshUpdater::setExtraStates(CFuint& nextStateID)
{
  // get the element type data
  SafePtr< std::vector< ElementTypeData > > elementType = getCFmeshData().getElementTypeData();
  std::vector< ElementTypeData >::iterator type_itr;
  const CFuint nbElemTypes = elementType->size();
  const CFuint dim = getCFmeshData().getDimension();

  // get the connectivity that we will change
  Common::SafePtr<MeshData::ConnTable> cellStates = MeshDataStack::getActive()->getConnectivity("cellStates_InnerCells");

  // get the nodal connectivity
  Common::SafePtr<MeshData::ConnTable> cellNodes = MeshDataStack::getActive()->getConnectivity("cellNodes_InnerCells");

  const CFuint nbElements = getNbElements();

  //marking if it is updated and the element type
  std::vector< CFuint > isUpdated(nbElements);
  for(CFuint iElem=0; iElem < nbElements; ++iElem)
  {
    isUpdated[iElem] = std::numeric_limits<CFuint>::max();
  }

CFout << "Creating neighbour connectivity\n";
  //Create neighbour connectivity
  const CFuint nbCells = cellNodes->nbRows();
  CFuint totalMapSize = 0;
  for(CFuint iCell = 0; iCell < nbCells; ++iCell)
  {
    const CFuint nbNodesInCell = cellNodes->nbCols(iCell);
    totalMapSize += nbNodesInCell;
  }

  //Create a CFMultiMap that will store for each node, the related cells
  Common::CFMultiMap<CFuint,CFuint> mapNode2Cells;

  //Reserve the memory of totalSize
  mapNode2Cells.reserve(totalMapSize);

  //Loop over the cells
  for(CFuint iCell = 0; iCell < nbCells; ++iCell)
  {
    const CFuint nbNodesInCell = cellNodes->nbCols(iCell);
    //Loop over the nodes of the cell and insert
    for(CFuint iNode = 0; iNode < nbNodesInCell; ++iNode)
    {
      mapNode2Cells.insert((*cellNodes)(iCell,iNode), iCell);
    }
  }
  mapNode2Cells.sortKeys();

  ///Store the element local connectivity in an easy way
  // local connectivity face-node for each element type
  // local connectivity face-state for each element type
  vector<Table<CFuint>*> faceNodeElement(nbElemTypes);
  vector<Table<CFuint>*> faceStateElement(nbElemTypes);
  for (CFuint iType = 0; iType < nbElemTypes; ++iType) {
    faceNodeElement[iType] = LocalConnectionData::getInstance().getFaceDofLocal
      ((*elementType)[iType].getGeoShape(),
       getGeometricPolyOrder(),
       NODE,
       CFPolyForm::LAGRANGE);
    faceStateElement[iType] = LocalConnectionData::getInstance().getFaceDofLocal
      ((*elementType)[iType].getGeoShape(),
       getSolutionPolyOrder(),
       STATE,
       CFPolyForm::LAGRANGE);
  }

  // local connectivity edge-node for each element type
  // local connectivity edge-state for each element type
  vector<Table<CFuint>*> edgeNodeElement(nbElemTypes);
  vector<Table<CFuint>*> edgeStateElement(nbElemTypes);
  if(dim == DIM_3D)
  {
    for (CFuint iType = 0; iType < nbElemTypes; ++iType) {
      edgeNodeElement[iType] = LocalConnectionData::getInstance().getEdgeDofLocal
        ((*elementType)[iType].getGeoShape(),
         getGeometricPolyOrder(),
         NODE,
         CFPolyForm::LAGRANGE);
      edgeStateElement[iType] = LocalConnectionData::getInstance().getEdgeDofLocal
        ((*elementType)[iType].getGeoShape(),
        getSolutionPolyOrder(),
        STATE,
        CFPolyForm::LAGRANGE);
    }
  }


if(dim == DIM_3D)
{
CFout << "We are in 3D - Adding extra states on the edges\n";
  //Loop over cells
  for (CFuint iType = 0; iType < nbElemTypes; ++iType)
  {
    const CFuint nbNodesPerElem = (*elementType)[iType].getNbNodes();

    // get the number of elements
    const CFuint nbElemPerType = (*elementType)[iType].getNbElems();

    // get start index of this element type in global element list
    CFuint globalIdx = (*elementType)[iType].getStartIdx();

    // loop over elements of this type and set the P1P1 states
    for (CFuint iElem = 0; iElem < nbElemPerType; ++iElem, ++globalIdx)
    {
      const CFuint elemID = globalIdx;

//CFout << "Cell: " << elemID <<"\n";
//CFout << "nbNeighCells: " << nbNeighCells <<"\n";
      //----------------------------------------------
      //Loop over the element edges to insert the edge extra states
      //----------------------------------------------
      const CFuint nbEdges = LocalConnectionData::getInstance().getNbEdgesInShape
        ((*elementType)[iType].getGeoShape());
// CFout << "nbEdges: " << nbEdges <<"\n";
      for(CFuint iEdge = 0; iEdge < nbEdges; ++iEdge)
      {
        const CFuint nbNodesPerEdge = edgeNodeElement[iType]->nbCols(iEdge);
        const CFuint nbStatesPerEdge = edgeStateElement[iType]->nbCols(iEdge);

        CFuint matchingCellID = std::numeric_limits<CFuint>::max();

        //Loop over the neighbours of the nodes of the edge
        //To check if one or more of the neighborCell of the edge has been updated
        const CFuint firstNodeID = (*cellNodes)(elemID, (*edgeNodeElement[iType])(iEdge, 0));
        typedef CFMultiMap<CFuint, CFuint>::MapIterator MapIterator;
        
	bool fo = false;
	pair<MapIterator, MapIterator> neighCellsOfNode0 = mapNode2Cells.find(firstNodeID, fo);
	if (!fo) cout << "FEM_HighOrderMeshUpdater::setExtraStates() => node " << firstNodeID  << " not found!\n";
	
        for (MapIterator neighCellItr1 = neighCellsOfNode0.first;
               neighCellItr1 != neighCellsOfNode0.second;
               ++neighCellItr1)
        {
          CFuint nbCommonNodes = 1;
          const CFuint neighCellID = neighCellItr1->second;
          if(neighCellID != elemID)
          {
            for (CFuint iNode = 1; iNode < nbNodesPerEdge; ++iNode)
            {
              const CFuint localNodeID = (*edgeNodeElement[iType])(iEdge, iNode);
              const CFuint nodeID = (*cellNodes)(elemID, localNodeID);
	      
	      bool fo1 = false;
              pair<MapIterator, MapIterator> neighCellsOfNode = mapNode2Cells.find(nodeID, fo1);
	      if (!fo1) cout << "FEM_HighOrderMeshUpdater::setExtraStates() => node " << nodeID << " not found!\n"; 
		
	      for (MapIterator neighCellItr2 = neighCellsOfNode.first;
                     neighCellItr2 != neighCellsOfNode.second;
                     ++neighCellItr2)
              {
                const CFuint otherNeighCellID = neighCellItr2->second;
                if(otherNeighCellID == neighCellID) ++nbCommonNodes;
              }
            }
            if((nbCommonNodes == nbNodesPerEdge) && (isUpdated[neighCellID]!=  std::numeric_limits<CFuint>::max()))
              matchingCellID = neighCellID;
          }
        }

        // Check if one of the neighborCells has already been updated
        bool highOrderNeighbor = false;
        if(matchingCellID != std::numeric_limits<CFuint>::max()) highOrderNeighbor = true;

        if(!highOrderNeighbor)
        {
          // Create new states on the edge
          for (CFuint iState = 0; iState < nbStatesPerEdge; ++iState)
          {
            const CFuint localStateID = (*edgeStateElement[iType])(iEdge, iState);
// unused // const CFuint nbNodesInElement = cellNodes->nbCols(elemID);
            bool nodalState(false);
            ///@todo find a safer way to check that it is a nodalState or not...
            if(localStateID < nbNodesPerElem) nodalState = true;

            if(!nodalState)
            {
              (*cellStates)(elemID, localStateID) = nextStateID;
// CFout << "Creating the state: " << nextStateID <<" in Cell: " << elemID<<"\n";
              nextStateID++;
            }
          }
        }
        else{
          // Add the new states created by cell on the edge

          //find which edge on the neighbor is matching
          CFuint matchingEdgeID = std::numeric_limits<CFuint>::max();
          const CFuint neighborType = isUpdated[matchingCellID];
          const CFuint nbNeighborEdges = edgeNodeElement[neighborType]->nbRows();
          for(CFuint neighEdge = 0; neighEdge < nbNeighborEdges; ++neighEdge)
          {
            const CFuint nbNodesInNeighborEdge = edgeNodeElement[neighborType]->nbCols(neighEdge);
            CFuint nbMatchingNodes = 0;
            for(CFuint iNodeNeighEdge = 0; iNodeNeighEdge < nbNodesInNeighborEdge; ++iNodeNeighEdge)
            {
              const CFuint neighNodeID = (*edgeNodeElement[neighborType])(neighEdge, iNodeNeighEdge);
              const CFuint neighNode = (*cellNodes)(matchingCellID, neighNodeID);
              for(CFuint iNodeEdge = 0; iNodeEdge < nbNodesPerEdge; ++iNodeEdge)
              {
                const CFuint currentNodeID = (*edgeNodeElement[iType])(iEdge, iNodeEdge);
                const CFuint currentNode = (*cellNodes)(elemID, currentNodeID);
                if(currentNode == neighNode) nbMatchingNodes++;
              }
            }
            if(nbMatchingNodes == nbNodesPerEdge) matchingEdgeID = neighEdge;
          }

          cf_assert(matchingEdgeID != std::numeric_limits<CFuint>::max());

          //Insert the states of the neighbour edge on the current edge
          //excepted if they are nodal states

          const CFuint nbStatesInNeighborEdge = edgeStateElement[neighborType]->nbCols(matchingEdgeID);
// unused // const CFuint nbNodesInNeighborEdge = edgeNodeElement[neighborType]->nbCols(matchingEdgeID);
          CFuint idx = 1;
//             CFout << "NeighbourCell: " << matchingCellID <<"\n";
//             CFout << "NeighbourEdge: " << matchingEdgeID <<"\n";
          for (CFuint iNeighState = 0; iNeighState < nbStatesInNeighborEdge; ++iNeighState)
          {
            bool nodalState(false);

            const CFuint neighLocalStateID = (*edgeStateElement[iType])(matchingEdgeID, iNeighState);
//             CFout << "EdgeState: " << neighLocalStateID <<"\n";
            const CFuint neighStateID = (*cellStates)(matchingCellID, neighLocalStateID);
/*            CFout << "Should we add the state: " << (*cellStates)(matchingCellID, neighLocalStateID) <<"\n";
            CFout << "Which is localID: " << neighLocalStateID <<" in Cell: " <<matchingCellID <<"\n";*/
            if(neighLocalStateID < nbNodesPerElem) nodalState = true;

            if(!nodalState)
            {
              CFuint newStateLocalID = (*edgeStateElement[iType])(iEdge, nbStatesPerEdge - idx);
              (*cellStates)(elemID, newStateLocalID) = neighStateID;
// CFout << "Adding the state: " << (*cellStates)(elemID, newStateLocalID) <<" at location: " << newStateLocalID << "\n";
              ++idx;
            }
          }
        }
      } // end loop over edges
      isUpdated[elemID] = iType;
    } // end loop over elements
  } // end loop over elementType
} // end of if(DIM_3D)

CFout << "Adding extra states on the faces\n";
  //Reset the isUpdated counter
  for(CFuint iElem=0; iElem < nbElements; ++iElem)
  {
    isUpdated[iElem] = std::numeric_limits<CFuint>::max();
  }

  //Loop over cells
  std::deque<CFuint> nextStateIDReserve(0);
  std::deque<CFuint> lastCreatedStateCell(0);
  std::deque<CFuint> lastCreatedStateLocalID(0);
  for (CFuint iType = 0; iType < nbElemTypes; ++iType)
  {
    const CFuint nbStatesPerElem = (*elementType)[iType].getNbStates();
    const CFuint nbElemPerType = (*elementType)[iType].getNbElems();

    // get start index of this element type in global element list
    CFuint globalIdx = (*elementType)[iType].getStartIdx();

    // If 3D, the edge states have already been added, so discard them
    vector<bool> isOnEdge(nbStatesPerElem, false);
    if(dim == DIM_3D)
    {
      const CFuint nbEdgesPerElement = edgeStateElement[iType]->nbRows();
      for(CFuint iEdge=0; iEdge<nbEdgesPerElement; ++iEdge)
      {
        const CFuint nbStatesPerEdge = edgeStateElement[iType]->nbCols(iEdge);
        for(CFuint iState=0; iState<nbStatesPerEdge; ++iState)
        {
          const CFuint localStateID = (*edgeStateElement[iType])(iEdge, iState);
          isOnEdge[localStateID] = true;
        }
      }
    }

    // loop over elements of this type and set the P1P1 states
    for (CFuint iElem = 0; iElem < nbElemPerType; ++iElem, ++globalIdx)
    {
      const CFuint elemID = globalIdx;
      const CFuint nbFaces = LocalConnectionData::getInstance().getNbFacesInShape
        ((*elementType)[iType].getGeoShape());
      vector <bool> isAlreadyCreatedState(nbStatesPerElem, false);

      //----------------------------------------------
      //Loop over the element faces to insert the face extra states
      //----------------------------------------------
      for(CFuint iFace = 0; iFace < nbFaces; ++iFace)
      {
//  CFout << "Face: " << iFace <<"\n";

        const CFuint nbNodesPerFace = faceNodeElement[iType]->nbCols(iFace);
        const CFuint nbStatesPerFace = faceStateElement[iType]->nbCols(iFace);

        CFuint matchingCellID = std::numeric_limits<CFuint>::max();


        // Get the neighborCell of the face
        const CFuint firstNodeID = (*cellNodes)(elemID, (*faceNodeElement[iType])(iFace, 0));
        typedef CFMultiMap<CFuint, CFuint>::MapIterator MapIterator;
	
	bool fo2 = false;
        pair<MapIterator, MapIterator> neighCellsOfNode0 = mapNode2Cells.find(firstNodeID, fo2);
	if (!fo2) cout << "FEM_HighOrderMeshUpdater::setExtraStates() => node " << firstNodeID << " not found!\n"; 
	
        for (MapIterator neighCellItr1 = neighCellsOfNode0.first;
               neighCellItr1 != neighCellsOfNode0.second;
               ++neighCellItr1)
        {
          CFuint nbCommonNodes = 1;
          const CFuint neighCellID = neighCellItr1->second;
          if(neighCellID != elemID)
          {
            for (CFuint iNode = 1; iNode < nbNodesPerFace; ++iNode)
            {
              const CFuint localNodeID = (*faceNodeElement[iType])(iFace, iNode);
              const CFuint nodeID = (*cellNodes)(elemID, localNodeID);
	      
	      bool fo3 = false;
              pair<MapIterator, MapIterator> neighCellsOfNode = mapNode2Cells.find(nodeID, fo3);
	      if (!fo3) cout << "FEM_HighOrderMeshUpdater::setExtraStates() => node " << nodeID << " not found!\n"; 
	      
	      for (MapIterator neighCellItr2 = neighCellsOfNode.first;
                     neighCellItr2 != neighCellsOfNode.second;
                     ++neighCellItr2)
              {
                const CFuint otherNeighCellID = neighCellItr2->second;
                if(otherNeighCellID == neighCellID) ++nbCommonNodes;
              }
            }
            if(nbCommonNodes == nbNodesPerFace) matchingCellID = neighCellID;
          }
        }

        // Check if the neighborCell has already been updated
        bool highOrderP2 = false;
        if((matchingCellID != std::numeric_limits<CFuint>::max()) &&
           (isUpdated[matchingCellID] !=  std::numeric_limits<CFuint>::max())) highOrderP2 = true;

        // If the neighbour is still P1, we should create states on the common face
        if(!highOrderP2)
        {
          for (CFuint iState = 0; iState < nbStatesPerFace; ++iState)
          {
            const CFuint localStateID = (*faceStateElement[iType])(iFace, iState);
            const CFuint nbNodesInElement = cellNodes->nbCols(elemID);
            bool nodalState(false);
            /// @todo find a safer way to check that it is a nodalState or not...
            if(localStateID < nbNodesInElement) nodalState = true;
            if(isOnEdge[localStateID]) nodalState = true;

            if((isAlreadyCreatedState[localStateID] == false) && (!nodalState))
            {
              if(nextStateIDReserve.size() > 0){
                (*cellStates)(elemID, localStateID) = nextStateIDReserve[0];
                nextStateIDReserve.pop_front();
              }
              else{
                (*cellStates)(elemID, localStateID) = nextStateID;
                nextStateID++;
                lastCreatedStateCell.push_back(elemID);
                lastCreatedStateLocalID.push_back(localStateID);

                if(lastCreatedStateCell.size() > 10) lastCreatedStateCell.pop_front();
                if(lastCreatedStateLocalID.size() > 10) lastCreatedStateLocalID.pop_front();
              }
              isAlreadyCreatedState[localStateID] = true;
            }
          }
        }
        // if the neighbour is already P2, we should share the states on the common face
        else
        {
          // check which face of the neighbor is matching with the current face
          CFuint matchingNeighborFaceID = std::numeric_limits<CFuint>::max();

          const CFuint neighborType = isUpdated[matchingCellID];
          const CFuint nbNeighborFaces = faceNodeElement[neighborType]->nbRows();
          for(CFuint neighFace = 0; neighFace < nbNeighborFaces; ++neighFace)
          {
            const CFuint nbNodesInNeighborFace = faceNodeElement[neighborType]->nbCols(neighFace);
            CFuint nbMatchingNodes = 0;
            for(CFuint iNodeNeighFace = 0; iNodeNeighFace < nbNodesInNeighborFace; ++iNodeNeighFace)
            {
              const CFuint neighNodeID = (*faceNodeElement[neighborType])(neighFace, iNodeNeighFace);
              const CFuint neighNode = (*cellNodes)(matchingCellID, neighNodeID);
              for(CFuint iNodeFace = 0; iNodeFace < nbNodesPerFace; ++iNodeFace)
              {
                const CFuint currentNodeID = (*faceNodeElement[iType])(iFace, iNodeFace);
                const CFuint currentNode = (*cellNodes)(elemID, currentNodeID);
                if(currentNode == neighNode) nbMatchingNodes++;
              }
            }
            if(nbMatchingNodes == nbNodesPerFace) matchingNeighborFaceID = neighFace;
          }

          cf_assert(matchingNeighborFaceID != std::numeric_limits<CFuint>::max());

          //we have the matching face of the matching neighbor
          // for each state, we have to check if it is a nodal state
          //if not we can add it
          const CFuint nbStatesInNeighborFace = faceStateElement[neighborType]->nbCols(matchingNeighborFaceID);
          CFuint idx = 1;

          for (CFuint iNeighState = 0; iNeighState < nbStatesInNeighborFace; ++iNeighState)
          {
            bool alreadyThere(false);

            const CFuint neighLocalStateID = (*faceStateElement[iType])(matchingNeighborFaceID, iNeighState);
            const CFuint neighStateID = (*cellStates)(matchingCellID, neighLocalStateID);
//             CFout << "Should we add the state: " << (*cellStates)(matchingCellID, neighLocalStateID) <<"\n";

            for (CFuint iState = 0; iState < nbStatesPerFace; ++iState)
            {
              const CFuint currentLocalStateID = (*faceStateElement[iType])(iFace, iState);
              const CFuint currentStateID = (*cellStates)(elemID, currentLocalStateID);
              if( neighStateID == currentStateID) alreadyThere = true;
              if(isOnEdge[currentLocalStateID]) alreadyThere = true;
            }

            //the state localNeighStateID should be added to the face iFace
            // but to which state in the face??
            if(!alreadyThere)
            {
//                CFout << "YES it is not yet there...\n";
//  CFout << "idx: " << nbStatesPerFace - idx <<" \n";
              CFuint newStateLocalID = (*faceStateElement[iType])(iFace, nbStatesPerFace - idx);

              ///@todo here we are not sure about the order of the states to be added
              ///->to be modified or make sure that it is ok!!!!
//   CFout << "Trying to insert it at location: " << newStateLocalID << " in Cell: " <<elemID << " \n";
//   CFout << "CellStates connectivity of cell: " << elemID << " has size: " << cellStates->nbCols(elemID) <<"\n";
              if(isAlreadyCreatedState[newStateLocalID] == true)
              {
                nextStateIDReserve.push_back((*cellStates)(elemID, newStateLocalID));
              }

              (*cellStates)(elemID, newStateLocalID) = neighStateID;
// CFout << "Adding the state: " << (*cellStates)(elemID, newStateLocalID) <<" at location: " << newStateLocalID << "\n";
              ++idx;
            }
          }
        }

        //if there is no neighbour sharing the face iFace,
        //then it is a boundary face, so save its connectivity for later
        if(matchingCellID == std::numeric_limits<CFuint>::max())
        {
// CFout << "This is a boundary face!!!\n";
          BoundaryFaceConnect bFaceConnect;
          vector<CFuint*> faceNodeConn(nbNodesPerFace);
          vector<CFuint*> faceStateConn(nbStatesPerFace);

          for (CFuint iNode = 0; iNode < nbNodesPerFace; ++iNode)
          {
            const CFuint localNodeID = (*faceNodeElement[iType])(iFace, iNode);
            faceNodeConn[iNode] = &(*cellNodes)(elemID, localNodeID);
          }
          for (CFuint iState = 0; iState < nbStatesPerFace; ++iState)
          {
            const CFuint localStateID = (*faceStateElement[iType])(iFace, iState);
            faceStateConn[iState] = &(*cellStates)(elemID, localStateID);
          }
          bFaceConnect.first = faceNodeConn;
          bFaceConnect.second = faceStateConn;

          _boundaryFaces.push_back(bFaceConnect);
        }
      }

        //
        if((elemID == nbCells-1) && (nextStateIDReserve.size() > 0))
        {
            CFout << "There are " << nextStateIDReserve.size() <<"remaining unused stateID\n";

cf_assert(lastCreatedStateCell.size() > nextStateIDReserve.size());
cf_assert(lastCreatedStateLocalID.size() > nextStateIDReserve.size());
for(CFuint iState = 0 ; iState < nextStateIDReserve.size(); ++iState)
{
const CFuint cellID2Replace = lastCreatedStateCell.size()-1-iState;
const CFuint localStateID2Replace = lastCreatedStateLocalID.size()-1-iState;

            CFout << "Last Created state is in cell: " << lastCreatedStateCell[cellID2Replace]<< " at localID: " << lastCreatedStateLocalID[localStateID2Replace]<<"\n";
            CFout << "Replacing: "<<(*cellStates)(lastCreatedStateCell[cellID2Replace],lastCreatedStateLocalID[localStateID2Replace])<<"\n";

(*cellStates)(lastCreatedStateCell[cellID2Replace],lastCreatedStateLocalID[localStateID2Replace]) = nextStateIDReserve[iState];

            CFout << "by: "<<(*cellStates)(lastCreatedStateCell[cellID2Replace],lastCreatedStateLocalID[localStateID2Replace])<<"\n";
  nextStateID -= 1;
}
CFout << "Total nb states: "<< nextStateID<<"\n";
nextStateIDReserve.resize(0);
        }

      const CFuint nbInteriorStates = getNbInteriorStatesInNewType((*elementType)[iType].getGeoShape(),m_newPolyOrder);
      for(CFuint iInteriorState=0; iInteriorState < nbInteriorStates ; ++iInteriorState)
      {
        (*cellStates)(elemID,nbStatesPerElem-nbInteriorStates) = nextStateID;
// CFout << "Creating the central state: " << nextStateID <<"\n";
        nextStateID++;
      }

      isUpdated[elemID] = iType;
    }
  }

}

//////////////////////////////////////////////////////////////////////////////

void FEM_HighOrderMeshUpdater::updateTRSData()
{
  CFAUTOTRACE;

  CFout << "Converting the TRS Data to P2\n";

  SafePtr< vector<CFuint> > nbTRs = getCFmeshData().getNbTRs();
  SafePtr< vector<vector<CFuint> > > nbGeomEntsPerTR = getCFmeshData().getNbGeomEntsPerTR();

  for(CFuint iTRS = 0; iTRS < getCFmeshData().getNbTRSs(); ++iTRS )
  {
    TRGeoConn& geoConn = getCFmeshData().getTRGeoConn(iTRS);
    for(CFuint iTR = 0; iTR < (*nbTRs)[iTRS]; ++iTR )
    {
      GeoConn& geoC = geoConn[iTR];
      const CFuint nbGeos = geoC.size();
      for (CFuint iGeo = 0; iGeo < nbGeos; ++iGeo)
      {
        std::valarray<CFuint>& nodeCon = geoC[iGeo].first;
        std::valarray<CFuint>& stateCon = geoC[iGeo].second;
        const CFuint oldNbStatesInGeo = stateCon.size();
        const CFuint nbStatesInGeo = getNewNbStatesInBoundaryGeo(oldNbStatesInGeo, m_newPolyOrder);

//         std::valarray<CFuint> oldStateCon = geoC[iGeo].second;
// CFout << "iGeo: " << iGeo <<"\n";
        stateCon.resize(nbStatesInGeo);
// CFout << "stateCon resized to : " << nbStatesInGeo <<"\n";
//         for(CFuint stateID = 0; stateID < oldNbStatesInGeo; ++stateID)
//         {
//           CFout << "states: " << stateCon[stateID] <<"\n";
//           stateCon[stateID] = oldStateCon[stateID];
//         }
        CFuint matchingFaceID = std::numeric_limits<CFuint>::max();
        const CFuint nbBoundaryFaces = _boundaryFaces.size();
// CFout << "Trying to get the corresponding boundaryFace among "<< nbBoundaryFaces<< " faces \n";
        //Look for the matching neighbor face by matching the nodes
        for(CFuint iFace = 0; iFace < nbBoundaryFaces; ++iFace)
        {
// CFout << "Trying to match face: "<< iFace<< "\n";
          vector<CFuint*>& faceNodes = _boundaryFaces[iFace].first;
          CFuint matchingNodes = 0;
          for(CFuint iNode = 0; iNode < nodeCon.size(); iNode++)
          {
// CFout << "Node: "<< nodeCon[iNode] << "\n";

            for(CFuint jNode = 0; jNode < faceNodes.size(); jNode++)
            {
// CFout << "BFace Node: "<< *(faceNodes[jNode]) << "\n";
              if(*(faceNodes[jNode]) == nodeCon[iNode]) matchingNodes++;
            }
          }

          if(matchingNodes == nodeCon.size())
          {
            matchingFaceID = iFace;
          }
        }

        cf_assert(matchingFaceID != std::numeric_limits<CFuint>::max());

        vector<CFuint*>& faceStates = _boundaryFaces[matchingFaceID].second;
        for(CFuint iState = 0; iState < faceStates.size(); ++iState)
        {
          stateCon[iState] = *(faceStates[iState]);
//           CFout << "states: " << stateCon[iState] <<"\n";
        }

        vector<CFuint*>& faceNodes = _boundaryFaces[matchingFaceID].first;
        for(CFuint iNode = 0; iNode < faceNodes.size(); ++iNode)
        {
           nodeCon[iNode] = *(faceNodes[iNode]);
//           CFout << "states: " << stateCon[iState] <<"\n";
        }
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void FEM_HighOrderMeshUpdater::recreateStates()
{
  CFAUTOTRACE;

  Common::SafePtr<MeshData::ConnTable> cellStates = MeshDataStack::getActive()->getConnectivity("cellStates_InnerCells");

  DataHandle < Framework::State*, Framework::GLOBAL > states = getCFmeshData().getStatesHandle();

  const CFuint newNbStates = _totalNewNbStates;
  const CFuint oldNbStates = states.size();
  CFout << "Updating the nb of states from: " << oldNbStates << "  to  " << newNbStates <<"\n";

  std::vector<State*> oldStates(oldNbStates);
  // backup the existing states
  for (CFuint i = 0; i < oldNbStates; ++i)
  {
    oldStates[i] = states[i];
  }

  // Resize the datahandle for the states
  states.resize(newNbStates);

  // backup the existing states
  for (CFuint i = 0; i < oldNbStates; ++i)
  {
    states[i] = oldStates[i];
  }

  // allocate the new states
  const CFuint nbeq = PhysicalModelStack::getActive()->getNbEq();
  RealVector stateData (nbeq);
  for (CFuint iState = oldNbStates; iState < states.size(); ++iState)
  {
    getCFmeshData().createState(iState,stateData);
    ///@todo modify this for the parallel
    states[iState]->setGlobalID(iState);
  }
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteElement

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
