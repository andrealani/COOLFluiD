#include <numeric>

#include "Common/CFLog.hh"
#include "Environment/ObjectProvider.hh"

#include "Common/BadValueException.hh"
#include "Framework/MeshData.hh"
#include "Framework/PhysicalModel.hh"
#include "Framework/BadFormatException.hh"

#include "FluxReconstructionMethod/FluxReconstruction.hh"
#include "FluxReconstructionMethod/MeshUpgradeBuilder.hh"
#include "FluxReconstructionMethod/FluxReconstructionElementData.hh"
#include "FluxReconstructionMethod/QuadFluxReconstructionElementData.hh"
#include "FluxReconstructionMethod/HexaFluxReconstructionElementData.hh"
#include "FluxReconstructionMethod/TriagFluxReconstructionElementData.hh"
#include "FluxReconstructionMethod/TetraFluxReconstructionElementData.hh"

#include <algorithm>
#include <iostream>

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////


Environment::ObjectProvider<MeshUpgradeBuilder,
               MeshDataBuilder,
               FluxReconstructionModule,
               1>
meshUpgradeBuilderProvider("MeshUpgrade");

//////////////////////////////////////////////////////////////////////////////

void MeshUpgradeBuilder::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< std::string >("PolynomialOrder","Flux Reconstruction polynomial order.");
  options.addConfigOption< std::string >("GeoPolynomialOrder","Geometrical polynomial order.");
  options.addConfigOption< CFuint >("DivideElements","Divide elements on equal parts to form new cells. This number is equal to te number of element adjescend to an old element face, so for 1 nothing happens.");
  options.addConfigOption< bool >("UpgradeInit","Flag telling whether the initialization needs to be upgraded too.");
}

//////////////////////////////////////////////////////////////////////////////

MeshUpgradeBuilder::MeshUpgradeBuilder(const std::string& name) :
  FluxReconstructionBuilder(name),
  m_solPolyOrder(),
  m_elementDivision(),
  m_geoPolyOrder(),
  m_prevGeoPolyOrder(),
  m_bndFacesNodes(),
  m_updatables(),
  m_globalIDs(),
  m_elemIDOfState(),
  m_elemLocalIDOfState(),
  m_elemFirstStateLocalID(),
  m_newToOldNodeID(),
  m_newStatesVal()
{
  addConfigOptionsTo(this);

  m_solPolyOrderStr = "P1";
  setParameter( "PolynomialOrder", &m_solPolyOrderStr);

  m_geoPolyOrderStr = "P1";
  setParameter( "GeoPolynomialOrder", &m_geoPolyOrderStr);
  
  m_elementDivision = 1;
  setParameter( "DivideElements", &m_elementDivision);
  
  m_upgradeInit = false;
  setParameter( "UpgradeInit", &m_upgradeInit);
}

//////////////////////////////////////////////////////////////////////////////

MeshUpgradeBuilder::~MeshUpgradeBuilder()
{
}

//////////////////////////////////////////////////////////////////////////////

void MeshUpgradeBuilder::releaseMemory()
{
  FluxReconstructionBuilder::releaseMemory();
}

//////////////////////////////////////////////////////////////////////////////

void MeshUpgradeBuilder::configure ( Config::ConfigArgs& args )
{
  FluxReconstructionBuilder::configure(args);

  m_solPolyOrder = CFPolyOrder::Convert::to_enum( m_solPolyOrderStr );

  m_geoPolyOrder = CFPolyOrder::Convert::to_enum( m_geoPolyOrderStr );
}

//////////////////////////////////////////////////////////////////////////////

CFuint MeshUpgradeBuilder::getNbrOfStatesInSDCellType(CFGeoShape::Type geoShape, CFPolyOrder::Type polyOrder)
{
  switch (geoShape)
  {
    case CFGeoShape::LINE:
      return polyOrder + 1;
    case CFGeoShape::TRIAG:
      return (polyOrder + 1)*(polyOrder + 2)/2;
    case CFGeoShape::QUAD:
      return (polyOrder + 1)*(polyOrder + 1);
    case CFGeoShape::TETRA:
      return (polyOrder + 1)*(polyOrder + 2)*(polyOrder + 3)/6;
    case CFGeoShape::HEXA:
      return (polyOrder + 1)*(polyOrder + 1)*(polyOrder + 1);
    default:
      throw Common::NotImplementedException (FromHere(),"Element type not implemented...");
  }
}

//////////////////////////////////////////////////////////////////////////////

CFuint MeshUpgradeBuilder::getNbrOfNodesInCellType(CFGeoShape::Type geoShape, CFPolyOrder::Type polyOrder)
{
  cout<<"getNbrOfNodesInCellType " << polyOrder << endl;
  return getNbrOfStatesInSDCellType(geoShape,polyOrder);
}

//////////////////////////////////////////////////////////////////////////////

CFuint MeshUpgradeBuilder::getNbrOfInternalNodesInCellType(CFGeoShape::Type geoShape, CFPolyOrder::Type polyOrder)
{
  if (polyOrder == CFPolyOrder::ORDER0)
  {
    throw BadValueException(FromHere(),"Cell with geometrical polynomial order of 0");
  }

  switch (geoShape)
  {
    case CFGeoShape::LINE:
      return polyOrder - 1;
    case CFGeoShape::TRIAG:
    {
      switch (polyOrder)
      {
        case CFPolyOrder::ORDER1:
        {
          return 0;
        }
        case CFPolyOrder::ORDER2:
        {
          return 0;
        }
        default:
        {
          throw Common::NotImplementedException (FromHere(),"Element type not implemented...");
        }
      }
    }
    case CFGeoShape::QUAD:
    {
      switch (polyOrder)
      {
        case CFPolyOrder::ORDER1:
        {
          return 0;
        }
        case CFPolyOrder::ORDER2:
        {
          return 1;
        }
        default:
        {
          throw Common::NotImplementedException (FromHere(),"Element type not implemented...");
        }
      }
    }
    case CFGeoShape::TETRA:
    {
      switch (polyOrder)
      {
        case CFPolyOrder::ORDER1:
        {
          return 0;
        }
        case CFPolyOrder::ORDER2:
        {
          return 0;
        }
        default:
        {
          throw Common::NotImplementedException (FromHere(),"Element type not implemented...");
        }
      }
    }
    case CFGeoShape::HEXA:
    {
      switch (polyOrder)
      {
        case CFPolyOrder::ORDER1:
        {
          return 0;
        }
        case CFPolyOrder::ORDER2:
        {
          return 1;
        }
        default:
        {
          throw Common::NotImplementedException (FromHere(),"Element type not implemented...");
        }
      }
    }
    default:
    {
      throw Common::NotImplementedException (FromHere(),"Element type not implemented...");
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

CFuint MeshUpgradeBuilder::getNbrOfNodesInFaceType(CFGeoShape::Type geoShape, CFPolyOrder::Type polyOrder)
{
  switch (geoShape)
  {
    case CFGeoShape::LINE:
      return polyOrder + 1;
    case CFGeoShape::TRIAG:
      return (polyOrder + 1)*(polyOrder + 2)/2;
    case CFGeoShape::QUAD:
      return (polyOrder + 1)*(polyOrder + 1);
    default:
      throw Common::NotImplementedException (FromHere(),"Face type not implemented...");
  }
}

//////////////////////////////////////////////////////////////////////////////

CFuint MeshUpgradeBuilder::getNbrOfInternalNodesInFaceType(CFGeoShape::Type geoShape, CFPolyOrder::Type polyOrder)
{
  if (polyOrder == CFPolyOrder::ORDER0)
  {
    throw BadValueException(FromHere(),"Face with geometrical polynomial order of 0");
  }
  switch (geoShape)
  {
    case CFGeoShape::LINE:
      return polyOrder - 1;
    case CFGeoShape::TRIAG:
    {
      switch (polyOrder)
      {
        case CFPolyOrder::ORDER1:
        {
          return 0;
        }
        case CFPolyOrder::ORDER2:
        {
          return 1;
        }
        default:
        {
          throw Common::NotImplementedException (FromHere(),"Element type not implemented...");
        }
      }
    }
    case CFGeoShape::QUAD:
    {
      switch (polyOrder)
      {
        case CFPolyOrder::ORDER1:
        {
          return 0;
        }
        case CFPolyOrder::ORDER2:
        {
          return 1;
        }
        default:
        {
          throw Common::NotImplementedException (FromHere(),"Face type not implemented...");
        }
      }
    }
    default:
    {
      throw Common::NotImplementedException (FromHere(),"Face type not implemented...");
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void MeshUpgradeBuilder::computeGeoTypeInfo()
{
  CFAUTOTRACE;

  // check requested geometrical polynomial order
  if (m_geoPolyOrder > CFPolyOrder::ORDER2)
  {
    throw Common::NotImplementedException (FromHere(),"Geometrical polynomial order higher than P2 not supported...");
  }

  // get previous geometrical polynomial order
  m_prevGeoPolyOrder = getCFmeshData().getGeometricPolyOrder();

  // get the element type data
  SafePtr< vector< ElementTypeData > > elementType = getCFmeshData().getElementTypeData();

  vector< ElementTypeData >::iterator type_itr = elementType->begin();
  for (; type_itr != elementType->end(); ++type_itr)
  {
    // get the number of solution points (states) in this element type
    const CFuint nbStatesPerElem = getNbrOfStatesInSDCellType(type_itr->getGeoShape(),m_solPolyOrder);
    cout<<"nbStatesPerElem = getNbrOfStatesInSDCellType  "<< nbStatesPerElem<< endl;
    // set the solution polynomial order and number of states
    type_itr->setSolOrder(m_solPolyOrder);
    type_itr->setNbStates(nbStatesPerElem);

    // get the number of nodes in this element type (is computed by the same expression as the number of states)
    const CFuint nbNodesPerElem = getNbrOfStatesInSDCellType(type_itr->getGeoShape(),m_geoPolyOrder);

    // set the geometric polynomial order and the number of nodes
    if (static_cast<CFuint>(m_geoPolyOrder) > type_itr->getGeoOrder())
    {
      type_itr->setGeoOrder(m_geoPolyOrder);
      type_itr->setNbNodes (nbNodesPerElem);
    }
  }

  /// @todo since only meshes with the same order are currently supported
  ///       we have to change the order in the CFmeshData and
  ///       not only on the ElementTypeData
  getCFmeshData().setSolutionPolyOrder (m_solPolyOrder);
  if (m_geoPolyOrder > m_prevGeoPolyOrder)
  {
    getCFmeshData().setGeometricPolyOrder(m_geoPolyOrder);
  }

  // continue with the standard algorithm
  FluxReconstructionBuilder::computeGeoTypeInfo();
}

//////////////////////////////////////////////////////////////////////////////

void MeshUpgradeBuilder::createTopologicalRegionSets()
{
  CFAUTOTRACE;
  
  if (m_elementDivision != 1)
  {
    CFLog(NOTICE, "MeshUpgradeBuilder: dividing elements in " << pow(m_elementDivision,getCFmeshData().getDimension()) << " equal parts to form new elements\n");
    if (m_prevGeoPolyOrder != CFPolyOrder::ORDER1)
    {
      CFLog(NOTICE, "For now, by dividing the elements, the geometric poly order will be reduced to P1!\n");
    }
    
    divideElements();
  }

  CFLog(NOTICE,"MeshUpgradeBuilder: upgrading mesh to solution polynomial order " << m_solPolyOrderStr << "\n");

  
  // first transform the cell-states connectivity
  // from cell centered to a FluxReconstruction
  upgradeStateConnectivity();

  // recreate the states
  recreateStates();

  // if necessary, add high order nodes
  if (m_geoPolyOrder > m_prevGeoPolyOrder)
  {
    CFLog(NOTICE,"MeshUpgradeBuilder: upgrading mesh to geometric polynomial order " << m_geoPolyOrderStr << "\n");

    upgradeNodeConnectivity();
    recreateNodes();
  }

  // continue with the standard algorithm
  FluxReconstructionBuilder::createTopologicalRegionSets();
}

//////////////////////////////////////////////////////////////////////////////

void MeshUpgradeBuilder::divideElements()
{
  // get the element type data
  SafePtr< vector<ElementTypeData> > elementType =
    getCFmeshData().getElementTypeData();
    
  // get the old state and node connectivities
  SafePtr<MeshData::ConnTable> cellStates = MeshDataStack::getActive()->getConnectivity("cellStates_InnerCells");
  SafePtr<MeshData::ConnTable> cellNodes = MeshDataStack::getActive()->getConnectivity("cellNodes_InnerCells");
  
  // get the nodes
  DataHandle < Framework::Node*, Framework::GLOBAL > nodes = getCFmeshData().getNodesHandle();
  
  // get the states
  DataHandle < Framework::State*, Framework::GLOBAL > states = getCFmeshData().getStatesHandle();
     
  const CFuint nbElements = getCFmeshData().getNbElements();
  const CFuint nbNodes = getCFmeshData().getTotalNbNodes();

  const CFuint nbElementTypes = getCFmeshData().getNbElementTypes();
  cf_assert(nbElementTypes == elementType->size());
  
  // the number of new cells that will be created in one old cell
  const CFuint nbNewCellsPerOldCell = pow(m_elementDivision,getCFmeshData().getDimension());
  
  // total nb of new elements
  const CFuint newNbElements = nbElements*nbNewCellsPerOldCell;
  
  // temporary storage for original nodes and state connectivities
  MeshData::ConnTable initialCellNodes = (*cellNodes);
  const MeshData::ConnTable initialCellStates = (*cellStates);
  
  // the new node coordinates
  vector< vector< RealVector > > newNodeCoords;
  
  // vector holding whether the old cell was updatable
  vector< bool > updatables;
  
  // this pattern will be used to resize the connectivities
  m_pattern.resize(newNbElements);
  
  // iterator over the elementTypes, for now there should be only one
  vector< ElementTypeData >::iterator type_itr;
  
  // first free idx for a cell
  CFuint firstFreeIdx = 0;
  // needed for a new node to cell map
  CFuint nodeToCellsMapSize = 0;
  // nb of new states
  CFuint nbNewStates = 0;
  
  // nodes to cells multimap for the old nodes
  CFMultiMap<CFuint,CFuint> mapNode2CellsOld;
  mapNode2CellsOld.reserve(nbNodes);

  // set the correct number of nodes per element in m_pattern
  for (type_itr = elementType->begin(); type_itr != elementType->end(); ++type_itr)
  {
    // get the current geo shape
    const CFGeoShape::Type elemGeoShape = type_itr->getGeoShape();
    
    // get the number of elements
    const CFuint nbrCellsPerType = type_itr->getNbElems();
    
    // loop over elements of this type
    CFuint globalIdx = type_itr->getStartIdx();
    for (CFuint iCell = 0; iCell < nbrCellsPerType; ++iCell, ++globalIdx)
    {
      // get the number of nodes in this cell type originally
      const CFuint nbNodesPerCell = getNbrOfNodesInCellType(elemGeoShape,m_prevGeoPolyOrder);
      
      // add new nodes to the new node to cell map
      nodeToCellsMapSize += nbNodesPerCell*nbNewCellsPerOldCell;

      // compute the coordinates of the new nodes
      newNodeCoords.push_back(getNewNodesCoords(elemGeoShape,m_elementDivision,globalIdx,initialCellNodes));

      // fill in updatables
      updatables.push_back(states[(*cellStates)(globalIdx,0)]->isParUpdatable());
      
      // create the old node to cell map
      for (CFuint iNode = 0; iNode < nbNodesPerCell; ++iNode)
      {
        mapNode2CellsOld.insert((*cellNodes)(globalIdx,iNode),globalIdx);
      }

      // set the nb of nodes in m_pattern for the old cell and add the new states to nbNewStates
      switch(elemGeoShape) 
      {

      case CFGeoShape::TRIAG:
	for (CFuint iNewCell = 0; iNewCell < nbNewCellsPerOldCell; ++iNewCell,++firstFreeIdx)
	{
	  m_pattern[firstFreeIdx] = 3;
	  nbNewStates += 3;
	}
        break;

      case CFGeoShape::QUAD:
        for (CFuint iNewCell = 0; iNewCell < nbNewCellsPerOldCell; ++iNewCell,++firstFreeIdx)
	{
	  m_pattern[firstFreeIdx] = 4;
	  nbNewStates += 4;
	}
        break;
	
      case CFGeoShape::HEXA:
        for (CFuint iNewCell = 0; iNewCell < nbNewCellsPerOldCell; ++iNewCell,++firstFreeIdx)
	{
	  m_pattern[firstFreeIdx] = 8;
	  nbNewStates += 8;
	}
        break;

      case CFGeoShape::TETRA:
        for (CFuint iNewCell = 0; iNewCell < nbNewCellsPerOldCell; ++iNewCell,++firstFreeIdx)
	{
	  m_pattern[firstFreeIdx] = 4;
	  nbNewStates += 4;
	}
        break;

      default:
        std::string shape =
          CFGeoShape::Convert::to_str(elemGeoShape);
        std::string msg = std::string("Element type not implemented: ") + shape;
        throw Common::NotImplementedException (FromHere(),msg);
      }
    }
  }
  cf_assert(firstFreeIdx == m_pattern.size());
  
  // sort the nodes to cells map
  mapNode2CellsOld.sortKeys();
  
  // set the connectivities to their new sizes
  cellNodes->resize(m_pattern);
  cellStates->resize(m_pattern);

  // resize the states data handle to put dummy states in it 
  states.resize(nbNewStates);
  // resize the node data handle for now to the new nb of states which is larger than the new nb of nodes, same for m_newToOldNodeID
  nodes.resize(nbNewStates);
  m_newToOldNodeID.resize(nbNewStates);
  
  // reset firstFreeIdx
  firstFreeIdx = 0;
  CFuint nodeIdx = 0;
  CFuint stateIdx = 0;
  const CFuint nodes1D = m_elementDivision + 1;
  
  // variable marking whether a cell has been updated
  // the element type is also stored here, though this is not really necessary since there is only one type of element
  vector<bool> isUpdated(nbElements);
  for (CFuint iElem = 0; iElem < nbElements; ++iElem)
  {
    isUpdated[iElem] = false;
  }
  
  // new nodes to cells multimap
  CFMultiMap<CFuint,CFuint> mapNode2Cells;
  typedef CFMultiMap<CFuint, CFuint>::MapIterator MapIterator;
  mapNode2Cells.reserve(nodeToCellsMapSize);
  
  const CFuint nbeq = PhysicalModelStack::getActive()->getNbEq();
  
  SafePtr<vector<ElementTypeData> > elementType2 = MeshDataStack::getActive()->getElementTypeData();
  
  vector< ElementTypeData >::iterator type_itr2 = elementType2->begin();
  
  // loop over element types to create new connectivities and nodes
  for (type_itr = elementType->begin(); type_itr != elementType->end(); ++type_itr, ++type_itr2)
  {
    // get the current geo shape
    const CFGeoShape::Type elemGeoShape = type_itr->getGeoShape();
    
    // get the number of elements
    const CFuint nbrCellsPerType = type_itr->getNbElems();
    
    // loop over elements of this type
    CFuint globalIdx = type_itr->getStartIdx();
    for (CFuint iCell = 0; iCell < nbrCellsPerType; ++iCell, ++globalIdx)
    {
      // get the number of nodes in this cell type originally
      const CFuint nbNodesPerCell = getNbrOfNodesInCellType(elemGeoShape,m_prevGeoPolyOrder);
      
      // first free idx for the new element
      const CFuint idxFirstNewElement = firstFreeIdx;
      
      // iterator over the neighbors of a given old node
      vector< pair<MapIterator, MapIterator> > neighbors;
      
      // neighbors of a face
      vector< CFint > faceNeighbors;
      
      // flag if neighbour cell has been found
      bool foundNghb = false;
      
      // get face to node connectivity table
      Table<CFuint>* face2Node = LocalConnectionData::getInstance().getFaceDofLocal(elemGeoShape,m_prevGeoPolyOrder,NODE,CFPolyForm::LAGRANGE);
      
      // neighbours of a given old node 
      vector< CFint > nodeNeighbors;

      switch(elemGeoShape) 
      {
      case CFGeoShape::TRIAG:

        break;

      case CFGeoShape::QUAD:
	
	// loop over the old nodes, to find their neighbour elements
	for (CFuint iNode = 0; iNode < 4; ++iNode)
	{
	  const CFuint id = initialCellNodes(globalIdx,iNode);

	  bool fo = false;
	  neighbors.push_back(mapNode2CellsOld.find(id,fo));
	  cf_assert(fo);
	}
	
	// for each old face, find the old neighbour elements
	for (CFuint iFace = 0; iFace < 4; ++iFace)
	{
	  for (MapIterator nghbrCellItr = neighbors[(*face2Node)(iFace,0)].first;
               nghbrCellItr != neighbors[(*face2Node)(iFace,0)].second;
               ++nghbrCellItr)
          {
            const CFuint currCellIdx = nghbrCellItr->second;
            if (currCellIdx != globalIdx)
	    {
	      for (MapIterator nghbrCellItr2 = neighbors[(*face2Node)(iFace,1)].first;
                   nghbrCellItr2 != neighbors[(*face2Node)(iFace,1)].second;
                   ++nghbrCellItr2)
	      {
		if (currCellIdx == nghbrCellItr2->second)
		{
		  cf_assert(!foundNghb);
		  // current cell is a neighbour of both nodes, so also of the face
		  faceNeighbors.push_back(currCellIdx);
		  foundNghb = true;
		}
	      }  
	    }
          }
          // if no neighbor is found, it is a bnd face, the negative value is used to detect this later
          if (!foundNghb)
	  {
	    faceNeighbors.push_back(-1);
	  }
	  foundNghb = false;
	}
	cf_assert(faceNeighbors.size() == 4);
	
	// find the neighbouring cells of the nodes in the same way as for the faces
	for (CFuint iNode = 0; iNode < 4; ++iNode)
	{
	  for (MapIterator nghbrCellItr = neighbors[iNode].first;
               nghbrCellItr != neighbors[iNode].second;
               ++nghbrCellItr)
          {
            const CFuint currCellIdx = nghbrCellItr->second;
            if (currCellIdx != globalIdx)
	    {
	      bool isInNghbrs = false;
	      for (CFuint iNghb = 0; iNghb < 4; ++iNghb)
	      {
		if (currCellIdx == static_cast< CFuint >(faceNeighbors[iNghb]))
		{
		  isInNghbrs = true;
		}
	      }  
	      if (!isInNghbrs)
	      {
		cf_assert(!foundNghb);
	        nodeNeighbors.push_back(currCellIdx);
		foundNghb = true;
	      }
	    }
          }
          // if no neighbor is found, it is a corner node
          if (!foundNghb)
	  {
	    nodeNeighbors.push_back(-1);
	  }
	  foundNghb = false;
	}
	cf_assert(nodeNeighbors.size() == 4);
	
	// loop over the new nodes and depending on where the new node is situated (corner, face, internal), create the node if it hasn't been created yet and 
	// fill in the new node connectivity
	for (CFuint iEta = 0; iEta < nodes1D; ++iEta)
	{
	  for (CFuint iKsi = 0; iKsi < nodes1D; ++iKsi)
	  {
	    // internal nodes, these must always be created
	    if (iKsi != 0 && iKsi != nodes1D-1 && iEta != 0 && iEta != nodes1D-1)
	    {
	      // insert in node to cell to connectivity
              mapNode2Cells.insert(nodeIdx, idxFirstNewElement+iKsi-1+(iEta-1)*m_elementDivision);
	      (*cellNodes)(idxFirstNewElement+iKsi-1+(iEta-1)*m_elementDivision,2) = nodeIdx;
	      mapNode2Cells.insert(nodeIdx, idxFirstNewElement+iKsi+(iEta-1)*m_elementDivision);
	      (*cellNodes)(idxFirstNewElement+iKsi+(iEta-1)*m_elementDivision,3) = nodeIdx;
	      mapNode2Cells.insert(nodeIdx, idxFirstNewElement+iKsi+(iEta)*m_elementDivision);
	      (*cellNodes)(idxFirstNewElement+iKsi+(iEta)*m_elementDivision,0) = nodeIdx;
	      mapNode2Cells.insert(nodeIdx, idxFirstNewElement+iKsi-1+(iEta)*m_elementDivision);
	      (*cellNodes)(idxFirstNewElement+iKsi-1+(iEta)*m_elementDivision,1) = nodeIdx;
	      
	      // create node 
	      getCFmeshData().createNode(nodeIdx,newNodeCoords[globalIdx][iKsi+nodes1D*iEta]);
	      nodes[nodeIdx]->setGlobalID(nodeIdx);
              nodes[nodeIdx]->setIsOnMesh(true);
              nodes[nodeIdx]->setIsOwnedByState(false);
	      ++nodeIdx;
	    }
	    // nodes on the 4 faces, but not on a corner
	    else if (iKsi == 0 && iEta != 0 && iEta != nodes1D-1)
	    {
	      // check if it is a bnd face
	      if (faceNeighbors[3] != -1)
	      {
		// check if the node has already been created
		if (!(isUpdated[faceNeighbors[3]]))
		{
		  // insert in nodes to cells map
	          mapNode2Cells.insert(nodeIdx, idxFirstNewElement+iKsi+(iEta-1)*m_elementDivision);
		  (*cellNodes)(idxFirstNewElement+iKsi+(iEta-1)*m_elementDivision,3) = nodeIdx;
	          mapNode2Cells.insert(nodeIdx, idxFirstNewElement+iKsi+(iEta)*m_elementDivision);
		  (*cellNodes)(idxFirstNewElement+iKsi+(iEta)*m_elementDivision,0) = nodeIdx;
		  
		  // create node 
	          getCFmeshData().createNode(nodeIdx,newNodeCoords[globalIdx][iKsi+nodes1D*iEta]);
	          nodes[nodeIdx]->setGlobalID(nodeIdx);
                  nodes[nodeIdx]->setIsOnMesh(true);
                  nodes[nodeIdx]->setIsOwnedByState(false);
	          ++nodeIdx;
		}
		else
		{
		  const CFuint oldNghbCellID = nbNewCellsPerOldCell*faceNeighbors[3];
		  const CFuint existingNodeIdx1 = (*cellNodes)(oldNghbCellID+nodes1D-2+(iEta-1)*m_elementDivision,2);
		  const CFuint existingNodeIdx2 = (*cellNodes)(oldNghbCellID+nodes1D-2+(iEta)*m_elementDivision,1);
		  (*cellNodes)(idxFirstNewElement+iKsi+(iEta-1)*m_elementDivision,3) = existingNodeIdx1;
		  (*cellNodes)(idxFirstNewElement+iKsi+(iEta)*m_elementDivision,0) = existingNodeIdx2;
		}
	      }
	      else 
	      {
		// insert in nodes to cells map
	        mapNode2Cells.insert(nodeIdx, idxFirstNewElement+iKsi+(iEta-1)*m_elementDivision);
		(*cellNodes)(idxFirstNewElement+iKsi+(iEta-1)*m_elementDivision,3) = nodeIdx;
	        mapNode2Cells.insert(nodeIdx, idxFirstNewElement+iKsi+(iEta)*m_elementDivision);
		(*cellNodes)(idxFirstNewElement+iKsi+(iEta)*m_elementDivision,0) = nodeIdx;
		
		m_newToOldNodeID[nodeIdx].push_back(initialCellNodes(globalIdx,0));
		m_newToOldNodeID[nodeIdx].push_back(initialCellNodes(globalIdx,3));
		  
		// create node 
	        getCFmeshData().createNode(nodeIdx,newNodeCoords[globalIdx][iKsi+nodes1D*iEta]);
	        nodes[nodeIdx]->setGlobalID(nodeIdx);
                nodes[nodeIdx]->setIsOnMesh(true);
                nodes[nodeIdx]->setIsOwnedByState(false);
	        ++nodeIdx;
	      } 
	    }
	    else if (iKsi == nodes1D-1 && iEta != 0 && iEta != nodes1D-1)
	    {
	      if (faceNeighbors[1] != -1)
	      {
		if (!(isUpdated[faceNeighbors[1]]))
		{
		  // insert in nodes to cells map
	          mapNode2Cells.insert(nodeIdx, idxFirstNewElement+iKsi-1+(iEta-1)*m_elementDivision);
		  (*cellNodes)(idxFirstNewElement+iKsi-1+(iEta-1)*m_elementDivision,2) = nodeIdx;
	          mapNode2Cells.insert(nodeIdx, idxFirstNewElement+iKsi-1+(iEta)*m_elementDivision);
		  (*cellNodes)(idxFirstNewElement+iKsi-1+(iEta)*m_elementDivision,1) = nodeIdx;
		  
		  // create node 
	          getCFmeshData().createNode(nodeIdx,newNodeCoords[globalIdx][iKsi+nodes1D*iEta]);
	          nodes[nodeIdx]->setGlobalID(nodeIdx);
                  nodes[nodeIdx]->setIsOnMesh(true);
                  nodes[nodeIdx]->setIsOwnedByState(false);
	          ++nodeIdx;
		}
		else
		{
		  const CFuint oldNghbCellID = nbNewCellsPerOldCell*faceNeighbors[1];
		  const CFuint existingNodeIdx1 = (*cellNodes)(oldNghbCellID+0+(iEta-1)*m_elementDivision,3);
		  const CFuint existingNodeIdx2 = (*cellNodes)(oldNghbCellID+0+(iEta)*m_elementDivision,0);
		  (*cellNodes)(idxFirstNewElement+iKsi-1+(iEta-1)*m_elementDivision,2) = existingNodeIdx1;
		  (*cellNodes)(idxFirstNewElement+iKsi-1+(iEta)*m_elementDivision,1) = existingNodeIdx2;
		}
	      }
	      else 
	      {
		// insert in nodes to cells map
	        mapNode2Cells.insert(nodeIdx, idxFirstNewElement+iKsi-1+(iEta-1)*m_elementDivision);
		(*cellNodes)(idxFirstNewElement+iKsi-1+(iEta-1)*m_elementDivision,2) = nodeIdx;
	        mapNode2Cells.insert(nodeIdx, idxFirstNewElement+iKsi-1+(iEta)*m_elementDivision);
		(*cellNodes)(idxFirstNewElement+iKsi-1+(iEta)*m_elementDivision,1) = nodeIdx;
		
		m_newToOldNodeID[nodeIdx].push_back(initialCellNodes(globalIdx,1));
		m_newToOldNodeID[nodeIdx].push_back(initialCellNodes(globalIdx,2));
		  
		// create node 
	        getCFmeshData().createNode(nodeIdx,newNodeCoords[globalIdx][iKsi+nodes1D*iEta]);
	        nodes[nodeIdx]->setGlobalID(nodeIdx);
                nodes[nodeIdx]->setIsOnMesh(true);
                nodes[nodeIdx]->setIsOwnedByState(false);
	        ++nodeIdx;
	      } 
	    }
	    else if (iEta == 0 && iKsi != 0 && iKsi != nodes1D-1)
	    {
	      if (faceNeighbors[0] != -1)
	      {
		if (!(isUpdated[faceNeighbors[0]]))
		{
		  // insert in nodes to cells map
	          mapNode2Cells.insert(nodeIdx, idxFirstNewElement+iKsi-1+(iEta)*m_elementDivision);
		  (*cellNodes)(idxFirstNewElement+iKsi-1+(iEta)*m_elementDivision,1) = nodeIdx;
	          mapNode2Cells.insert(nodeIdx, idxFirstNewElement+iKsi+(iEta)*m_elementDivision);
		  (*cellNodes)(idxFirstNewElement+iKsi+(iEta)*m_elementDivision,0) = nodeIdx;
		  
		  // create node 
	          getCFmeshData().createNode(nodeIdx,newNodeCoords[globalIdx][iKsi+nodes1D*iEta]);
	          nodes[nodeIdx]->setGlobalID(nodeIdx);
                  nodes[nodeIdx]->setIsOnMesh(true);
                  nodes[nodeIdx]->setIsOwnedByState(false);
	          ++nodeIdx;
		}
		else
		{
		  const CFuint oldNghbCellID = nbNewCellsPerOldCell*faceNeighbors[0];
		  const CFuint existingNodeIdx1 = (*cellNodes)(oldNghbCellID+iKsi-1+(nodes1D-2)*m_elementDivision,2);
		  const CFuint existingNodeIdx2 = (*cellNodes)(oldNghbCellID+iKsi+(nodes1D-2)*m_elementDivision,3);
		  (*cellNodes)(idxFirstNewElement+iKsi-1+(iEta)*m_elementDivision,1) = existingNodeIdx1;
		  (*cellNodes)(idxFirstNewElement+iKsi+(iEta)*m_elementDivision,0) = existingNodeIdx2;
		}
	      }
	      else 
	      {
		// insert in nodes to cells map
	        mapNode2Cells.insert(nodeIdx, idxFirstNewElement+iKsi-1+(iEta)*m_elementDivision);
		(*cellNodes)(idxFirstNewElement+iKsi-1+(iEta)*m_elementDivision,1) = nodeIdx;
	        mapNode2Cells.insert(nodeIdx, idxFirstNewElement+iKsi+(iEta)*m_elementDivision);
		(*cellNodes)(idxFirstNewElement+iKsi+(iEta)*m_elementDivision,0) = nodeIdx;
		
		m_newToOldNodeID[nodeIdx].push_back(initialCellNodes(globalIdx,0));
		m_newToOldNodeID[nodeIdx].push_back(initialCellNodes(globalIdx,1));
		  
		// create node 
	        getCFmeshData().createNode(nodeIdx,newNodeCoords[globalIdx][iKsi+nodes1D*iEta]);
	        nodes[nodeIdx]->setGlobalID(nodeIdx);
                nodes[nodeIdx]->setIsOnMesh(true);
                nodes[nodeIdx]->setIsOwnedByState(false);
	        ++nodeIdx;
	      } 
	    }
	    else if (iEta == nodes1D-1 && iKsi != 0 && iKsi != nodes1D-1)
	    {
	      if (faceNeighbors[2] != -1)
	      {
		if (!(isUpdated[faceNeighbors[2]]))
		{
		  // insert in nodes to cells map
	          mapNode2Cells.insert(nodeIdx, idxFirstNewElement+iKsi-1+(iEta-1)*m_elementDivision);
		  (*cellNodes)(idxFirstNewElement+iKsi-1+(iEta-1)*m_elementDivision,2) = nodeIdx;
	          mapNode2Cells.insert(nodeIdx, idxFirstNewElement+iKsi+(iEta-1)*m_elementDivision);
		  (*cellNodes)(idxFirstNewElement+iKsi+(iEta-1)*m_elementDivision,3) = nodeIdx;
		  
		  // create node 
	          getCFmeshData().createNode(nodeIdx,newNodeCoords[globalIdx][iKsi+nodes1D*iEta]);
	          nodes[nodeIdx]->setGlobalID(nodeIdx);
                  nodes[nodeIdx]->setIsOnMesh(true);
                  nodes[nodeIdx]->setIsOwnedByState(false);
	          ++nodeIdx;
		}
		else
		{
		  const CFuint oldNghbCellID = nbNewCellsPerOldCell*faceNeighbors[2];
		  const CFuint existingNodeIdx1 = (*cellNodes)(oldNghbCellID+iKsi-1+(0)*m_elementDivision,1);
		  const CFuint existingNodeIdx2 = (*cellNodes)(oldNghbCellID+iKsi+(0)*m_elementDivision,0);
		  (*cellNodes)(idxFirstNewElement+iKsi-1+(iEta-1)*m_elementDivision,2) = existingNodeIdx1;
		  (*cellNodes)(idxFirstNewElement+iKsi+(iEta-1)*m_elementDivision,3) = existingNodeIdx2;
		}
	      }
	      else 
	      {
		// insert in nodes to cells map
	        mapNode2Cells.insert(nodeIdx, idxFirstNewElement+iKsi-1+(iEta-1)*m_elementDivision);
		(*cellNodes)(idxFirstNewElement+iKsi-1+(iEta-1)*m_elementDivision,2) = nodeIdx;
	        mapNode2Cells.insert(nodeIdx, idxFirstNewElement+iKsi+(iEta-1)*m_elementDivision);
		(*cellNodes)(idxFirstNewElement+iKsi+(iEta-1)*m_elementDivision,3) = nodeIdx;
		
		m_newToOldNodeID[nodeIdx].push_back(initialCellNodes(globalIdx,2));
		m_newToOldNodeID[nodeIdx].push_back(initialCellNodes(globalIdx,3));
		  
		// create node 
	        getCFmeshData().createNode(nodeIdx,newNodeCoords[globalIdx][iKsi+nodes1D*iEta]);
	        nodes[nodeIdx]->setGlobalID(nodeIdx);
                nodes[nodeIdx]->setIsOnMesh(true);
                nodes[nodeIdx]->setIsOwnedByState(false);
	        ++nodeIdx;
	      } 
	    }
	    // the four corner nodes
	    else if (iKsi == 0 && iEta == 0)
	    {
	      bool nodeFound = false;
	      // insert in nodes to cells map
	      if (nodeNeighbors[0] != -1)
	      {
		if (isUpdated[nodeNeighbors[0]])
		{
		  const CFuint oldNghbCellID = nbNewCellsPerOldCell*nodeNeighbors[0];
		  const CFuint existingNodeIdx = (*cellNodes)(oldNghbCellID+nodes1D-2+(nodes1D-2)*m_elementDivision,2);
		  (*cellNodes)(idxFirstNewElement+iKsi+(iEta)*m_elementDivision,0) = existingNodeIdx;
		  nodeFound = true;
		}
	      }
              if (faceNeighbors[0] != -1 && !nodeFound)
	      {
		if (isUpdated[faceNeighbors[0]])
		{
		  const CFuint oldNghbCellID = nbNewCellsPerOldCell*faceNeighbors[0];
		  const CFuint existingNodeIdx = (*cellNodes)(oldNghbCellID+0+(nodes1D-2)*m_elementDivision,3);
		  (*cellNodes)(idxFirstNewElement+iKsi+(iEta)*m_elementDivision,0) = existingNodeIdx;
		  nodeFound = true;
		}
	      }
	      if (faceNeighbors[3] != -1 && !nodeFound)
	      {
		if (isUpdated[faceNeighbors[3]])
		{
		  const CFuint oldNghbCellID = nbNewCellsPerOldCell*faceNeighbors[3];
		  const CFuint existingNodeIdx = (*cellNodes)(oldNghbCellID+nodes1D-2+(0)*m_elementDivision,1);
		  (*cellNodes)(idxFirstNewElement+iKsi+(iEta)*m_elementDivision,0) = existingNodeIdx;
		  nodeFound = true;
		}
	      }
	      // node has not yet been created
	      if (!nodeFound) 
	      {
		// insert in nodes to cells map
	        mapNode2Cells.insert(nodeIdx, idxFirstNewElement+iKsi+(iEta)*m_elementDivision);
		(*cellNodes)(idxFirstNewElement+iKsi+(iEta)*m_elementDivision,0) = nodeIdx;
		
		m_newToOldNodeID[nodeIdx].push_back(initialCellNodes(globalIdx,0));
		  
		// create node 
	        getCFmeshData().createNode(nodeIdx,newNodeCoords[globalIdx][iKsi+nodes1D*iEta]);
	        nodes[nodeIdx]->setGlobalID(nodeIdx);
                nodes[nodeIdx]->setIsOnMesh(true);
                nodes[nodeIdx]->setIsOwnedByState(false);
	        ++nodeIdx;
	      } 
	    }
	    else if (iKsi == nodes1D-1 && iEta == 0)
	    {
	      bool nodeFound = false;
	      // insert in nodes to cells map
	      if (nodeNeighbors[1] != -1)
	      {
		if (isUpdated[nodeNeighbors[1]])
		{
		  const CFuint oldNghbCellID = nbNewCellsPerOldCell*nodeNeighbors[1];
		  const CFuint existingNodeIdx = (*cellNodes)(oldNghbCellID+0+(nodes1D-2)*m_elementDivision,3);
		  (*cellNodes)(idxFirstNewElement+iKsi-1+(iEta)*m_elementDivision,1) = existingNodeIdx;
		  nodeFound = true;
		}
	      }
	      if (faceNeighbors[0] != -1 && !nodeFound)
	      {
		if (isUpdated[faceNeighbors[0]])
		{
		  const CFuint oldNghbCellID = nbNewCellsPerOldCell*faceNeighbors[0];
		  const CFuint existingNodeIdx = (*cellNodes)(oldNghbCellID+nodes1D-2+(nodes1D-2)*m_elementDivision,2);
		  (*cellNodes)(idxFirstNewElement+iKsi-1+(iEta)*m_elementDivision,1) = existingNodeIdx;
		  nodeFound = true;
		}
	      }
	      if (faceNeighbors[1] != -1 && !nodeFound)
	      {
	        if (isUpdated[faceNeighbors[1]])
		{
		  const CFuint oldNghbCellID = nbNewCellsPerOldCell*faceNeighbors[1];
		  const CFuint existingNodeIdx = (*cellNodes)(oldNghbCellID+0+(0)*m_elementDivision,0);
		  (*cellNodes)(idxFirstNewElement+iKsi-1+(iEta)*m_elementDivision,1) = existingNodeIdx;
		  nodeFound = true;
		}
	      }
	      // node has not yet been created
	      if (!nodeFound) 
	      {
		// insert in nodes to cells map
	        mapNode2Cells.insert(nodeIdx, idxFirstNewElement+iKsi-1+(iEta)*m_elementDivision);
		(*cellNodes)(idxFirstNewElement+iKsi-1+(iEta)*m_elementDivision,1) = nodeIdx;
		
		m_newToOldNodeID[nodeIdx].push_back(initialCellNodes(globalIdx,1));
		  
		// create node 
	        getCFmeshData().createNode(nodeIdx,newNodeCoords[globalIdx][iKsi+nodes1D*iEta]);
	        nodes[nodeIdx]->setGlobalID(nodeIdx);
                nodes[nodeIdx]->setIsOnMesh(true);
                nodes[nodeIdx]->setIsOwnedByState(false);
	        ++nodeIdx;
	      } 
	    }
	    else if (iKsi == nodes1D-1 && iEta == nodes1D-1)
	    {
	      bool nodeFound = false;
	      // insert in nodes to cells map
	      if (nodeNeighbors[2] != -1)
	      {
		if (isUpdated[nodeNeighbors[2]])
		{
		  const CFuint oldNghbCellID = nbNewCellsPerOldCell*nodeNeighbors[2];
		  const CFuint existingNodeIdx = (*cellNodes)(oldNghbCellID+0+(0)*m_elementDivision,0);
		  (*cellNodes)(idxFirstNewElement+iKsi-1+(iEta-1)*m_elementDivision,2) = existingNodeIdx;
		  nodeFound = true;
		}
	      }
	      if (faceNeighbors[1] != -1 && !nodeFound)
	      {
		if (isUpdated[faceNeighbors[1]])
		{
		  const CFuint oldNghbCellID = nbNewCellsPerOldCell*faceNeighbors[1];
		  const CFuint existingNodeIdx = (*cellNodes)(oldNghbCellID+0+(nodes1D-2)*m_elementDivision,3);
		  (*cellNodes)(idxFirstNewElement+iKsi-1+(iEta-1)*m_elementDivision,2) = existingNodeIdx;
		  nodeFound = true;
		}
	      }
	      if (faceNeighbors[2] != -1 && !nodeFound)
	      {
		if (isUpdated[faceNeighbors[2]])
		{
		  const CFuint oldNghbCellID = nbNewCellsPerOldCell*faceNeighbors[2];
		  const CFuint existingNodeIdx = (*cellNodes)(oldNghbCellID+nodes1D-2+(0)*m_elementDivision,1);
		  (*cellNodes)(idxFirstNewElement+iKsi-1+(iEta-1)*m_elementDivision,2) = existingNodeIdx;
		  nodeFound = true;
		}
	      }
	      // node has not yet been created
	      if (!nodeFound) 
	      {
		// insert in nodes to cells map
	        mapNode2Cells.insert(nodeIdx, idxFirstNewElement+iKsi-1+(iEta-1)*m_elementDivision);
		(*cellNodes)(idxFirstNewElement+iKsi-1+(iEta-1)*m_elementDivision,2) = nodeIdx;
		
		m_newToOldNodeID[nodeIdx].push_back(initialCellNodes(globalIdx,2));
		  
		// create node 
	        getCFmeshData().createNode(nodeIdx,newNodeCoords[globalIdx][iKsi+nodes1D*iEta]);
	        nodes[nodeIdx]->setGlobalID(nodeIdx);
                nodes[nodeIdx]->setIsOnMesh(true);
                nodes[nodeIdx]->setIsOwnedByState(false);
	        ++nodeIdx;
	      } 
	    }
	    else if (iKsi == 0 && iEta == nodes1D-1)
	    {
	      bool nodeFound = false;
	      // insert in nodes to cells map
	      if (nodeNeighbors[3] != -1)
	      {
		if (isUpdated[nodeNeighbors[3]])
		{
		  const CFuint oldNghbCellID = nbNewCellsPerOldCell*nodeNeighbors[3];
		  const CFuint existingNodeIdx = (*cellNodes)(oldNghbCellID+nodes1D-2+(0)*m_elementDivision,1);
		  (*cellNodes)(idxFirstNewElement+iKsi+(iEta-1)*m_elementDivision,3) = existingNodeIdx;
		  nodeFound = true;
		}
	      }
	      if (faceNeighbors[2] != -1 && !nodeFound)
	      {
		if (isUpdated[faceNeighbors[2]])
		{
		  const CFuint oldNghbCellID = nbNewCellsPerOldCell*faceNeighbors[2];
		  const CFuint existingNodeIdx = (*cellNodes)(oldNghbCellID+0+(0)*m_elementDivision,0);
		  (*cellNodes)(idxFirstNewElement+iKsi+(iEta-1)*m_elementDivision,3) = existingNodeIdx;
		  nodeFound = true;
		}
	      }
	      if (faceNeighbors[3] != -1 && !nodeFound)
	      {
		if (isUpdated[faceNeighbors[3]])
		{
		  const CFuint oldNghbCellID = nbNewCellsPerOldCell*faceNeighbors[3];
		  const CFuint existingNodeIdx = (*cellNodes)(oldNghbCellID+nodes1D-2+(nodes1D-2)*m_elementDivision,2);
		  (*cellNodes)(idxFirstNewElement+iKsi+(iEta-1)*m_elementDivision,3) = existingNodeIdx;
		  nodeFound = true;
		}
	      }
	      // node has not yet been created
	      if (!nodeFound) 
	      {
		// insert in nodes to cells map
	        mapNode2Cells.insert(nodeIdx, idxFirstNewElement+iKsi+(iEta-1)*m_elementDivision);
		(*cellNodes)(idxFirstNewElement+iKsi+(iEta-1)*m_elementDivision,3) = nodeIdx;
		
		m_newToOldNodeID[nodeIdx].push_back(initialCellNodes(globalIdx,3));
		  
		// create node 
	        getCFmeshData().createNode(nodeIdx,newNodeCoords[globalIdx][iKsi+nodes1D*iEta]);
	        nodes[nodeIdx]->setGlobalID(nodeIdx);
                nodes[nodeIdx]->setIsOnMesh(true);
                nodes[nodeIdx]->setIsOwnedByState(false);
	        ++nodeIdx;
	      } 
	    }
	  }
	}

	// create the states, these will be deleted later and made properly, when the solution order is upgraded.
        for (CFuint iEta = 0; iEta < m_elementDivision; ++iEta)
	{
	  for (CFuint iKsi = 0; iKsi < m_elementDivision; ++iKsi, ++firstFreeIdx)
	  {
	    for (CFuint iNode = 0; iNode < nbNodesPerCell; ++iNode, ++stateIdx)
	    { 
	      (*cellStates)(firstFreeIdx,iNode) = stateIdx;
	      RealVector stateData (nbeq);
              getCFmeshData().createState(stateIdx,stateData);
              //states[nodeIdx]->setParUpdatable(updatables[globalIdx]);
              states[stateIdx]->setGlobalID(stateIdx);
	    }
	  }
	}
	
        break;
	
      case CFGeoShape::PRISM:

        break;

      case CFGeoShape::TETRA:

        break;

      case CFGeoShape::HEXA:

        break;

      default:
        std::string shape =
          CFGeoShape::Convert::to_str(elemGeoShape);
        std::string msg = std::string("Element type not implemented: ") + shape;
        throw Common::NotImplementedException (FromHere(),msg);
      }
      isUpdated[globalIdx] = true;
    }
    // set the new nb of elements
    type_itr->setNbElems(nbNewCellsPerOldCell*nbrCellsPerType);
    type_itr2->setNbElems(nbNewCellsPerOldCell*nbrCellsPerType);
    //CFLog(VERBOSE,"Nb of new elems of this type: " << nbNewCellsPerOldCell*nbrCellsPerType << "\n");
  }
  
  // set the oversized node data handle to the correct size, same for m_newToOldNodeID
  nodes.resize(nodeIdx);
  m_newToOldNodeID.resize(nodeIdx);

  // update the CFmeshData
  getCFmeshData().setNbElements(newNbElements);
  getCFmeshData().setNbUpdatableNodes(nodeIdx);
  getCFmeshData().setNbUpdatableStates(nbNewStates);
  getCFmeshData().setNbNonUpdatableNodes(0);
  getCFmeshData().setNbNonUpdatableStates(0);
  
}

//////////////////////////////////////////////////////////////////////////////

void MeshUpgradeBuilder::upgradeStateConnectivity()
{
  CFAUTOTRACE;
  CFLog(VERBOSE, "upgradeStateConnectivity start\n");
  // get the element type data
  SafePtr< vector< ElementTypeData > > elementType = getCFmeshData().getElementTypeData();
  vector< ElementTypeData >::iterator type_itr;

  // get the connectivity that we will change
  SafePtr<MeshData::ConnTable> cellStates = MeshDataStack::getActive()->getConnectivity("cellStates_InnerCells");
  
  MeshData::ConnTable cellStatesOld = *cellStates;
  
  // get the old states
  DataHandle < Framework::State*, Framework::GLOBAL > oldStates = getCFmeshData().getStatesHandle();



  const CFuint nbElems = getNbElements();
  
  const CFuint oldNbStatesPerElem = std::max((oldStates.size())/nbElems,(CFuint) 1);
  
  cout<<" oldStates.size() "<< oldStates.size() << endl;
 cout<<" oldNbStatesPerElem below "<< oldNbStatesPerElem << endl;

  // temporary vector to store which old states are parallel updatable
  vector< bool > updatablesElems;
  updatablesElems.resize(nbElems);
  vector< vector< CFuint > > globalIDsTemp;
  globalIDsTemp.resize(nbElems);
  m_updatables.resize(0);
  m_globalIDs.resize(0);
  m_elemFirstStateLocalID.resize(0);

  const CFuint oldMaxGlobalID = oldStates.getGlobalSize();

  // loop on all the elements to store which old states are parallel updatable
  for (type_itr = elementType->begin(); type_itr != elementType->end(); ++type_itr)
  {
    // get the number of elements
    const CFuint nbrElems = type_itr->getNbElems();

    // get start index of this element type in global element list
    CFuint globalIdx = type_itr->getStartIdx();

    // loop over elements of this type
    for (CFuint iElem = 0; iElem < nbrElems; ++iElem, ++globalIdx)
    {
      globalIDsTemp[globalIdx].resize(oldNbStatesPerElem);
      CFuint stateID = (*cellStates)(globalIdx,0);
      m_elemFirstStateLocalID.push_back(stateID);
      updatablesElems[globalIdx] = oldStates[stateID]->isParUpdatable();
      
      for (CFuint jState = 0; jState < oldNbStatesPerElem; ++jState)
      {
        globalIDsTemp[globalIdx][jState] = oldStates[(*cellStates)(globalIdx,jState)]->getGlobalID();
      }
    }
  }

  // we will create a mesh with equal order on all the elements
  // and we reset the connectivity accordingly, using a std::valarray
  // first create the std::valarray
  std::valarray<CFuint> columnPattern(nbElems);
  m_elemLocalIDOfState.resize(0);
  m_elemIDOfState.resize(0);

  for (type_itr = elementType->begin(); type_itr != elementType->end(); ++type_itr)
  {
    // get the number of control volumes (states) in this element type
    const CFuint nbStatesPerElem = getNbrOfStatesInSDCellType(type_itr->getGeoShape(),m_solPolyOrder);

    // get the number of elements
    const CFuint nbrElems = type_itr->getNbElems();

    // get start index of this element type in global element list
    CFuint globalIdx = type_itr->getStartIdx();

    // loop over elements of this type
    for (CFuint iElem = 0; iElem < nbrElems; ++iElem, ++globalIdx)
    {

      columnPattern[globalIdx] = nbStatesPerElem;
    }
  }
  // then resize the connectivity
  cellStates->resize(columnPattern);
  SafePtr< vector< RealVector > > newSolPntCoords;
  std::vector< std::vector< CFreal > > solPolyValsOld;

  if (m_upgradeInit) 
  {
    FluxReconstructionElementData* frElemData;
    FluxReconstructionElementData* frElemData2;
    
    m_newStatesVal.resize(cellStates->size());
    const CFGeoShape::Type cellShape = ((*elementType)[0]).getGeoShape();
    CFuint oldOrder;
    CFPolyOrder::Type oldOrder2;
  
    switch (cellShape)
    {
      case CFGeoShape::LINE:
      {
        throw Common::NotImplementedException (FromHere(),"Flux Reconstruction has not been implemented for 1D");
      } break;
      case CFGeoShape::QUAD:
      {
	      oldOrder = static_cast<CFuint> (sqrt(oldNbStatesPerElem)+0.5) - 1;
	      oldOrder2 = CFPolyOrder::Convert::to_enum(oldOrder);
        frElemData = new QuadFluxReconstructionElementData(oldOrder2);
	      frElemData2 = new QuadFluxReconstructionElementData(m_solPolyOrder);
      } break;
      case CFGeoShape::HEXA:
      {
	  oldOrder = static_cast<CFuint>  (pow(oldNbStatesPerElem,1./3.)+0.5) - 1;
	  oldOrder2 = CFPolyOrder::Convert::to_enum(oldOrder);
        frElemData = new HexaFluxReconstructionElementData(oldOrder2);
	frElemData2 = new HexaFluxReconstructionElementData(m_solPolyOrder);
      } break;
      case CFGeoShape::TRIAG:
      {
	      oldOrder = static_cast<CFuint>  ((-3 + sqrt(1+8*oldNbStatesPerElem))/2 +0.5)   ; 
	      oldOrder2 = CFPolyOrder::Convert::to_enum(oldOrder);
        cout<< " MeshUpgradeBuilder::oldNbStatesPerElem   " << oldNbStatesPerElem << endl;
        cout<< " MeshUpgradeBuilder::oldOrder2  " << oldOrder2 << endl;
        cout<< " MeshUpgradeBuilder::m_solPolyOrder  "<< m_solPolyOrder << endl;
        frElemData = new TriagFluxReconstructionElementData(oldOrder2);
	      frElemData2 = new TriagFluxReconstructionElementData(m_solPolyOrder);
      } break;
      case CFGeoShape::TETRA:
      {
  oldOrder = static_cast<CFuint> (0.5 -2. + pow((27.*oldNbStatesPerElem + sqrt(-3. + 729.*pow(oldNbStatesPerElem,2))),(1./3.)) / pow(3.,(2./3.)) + 1./pow((81.*oldNbStatesPerElem + 3.*sqrt(-3. + 729.*pow(oldNbStatesPerElem,2.))),(1./3.)));
	oldOrder2 = CFPolyOrder::Convert::to_enum(oldOrder);
        frElemData = new TetraFluxReconstructionElementData(oldOrder2);
	frElemData2 = new TetraFluxReconstructionElementData(m_solPolyOrder);
      } break;
      default:
      {
        throw Common::ShouldNotBeHereException (FromHere(),"Unsupported cell shape");
      }
    }
    newSolPntCoords = frElemData2->getSolPntsLocalCoords();
    cout<< " ----------------------------- newSolPntCoords   " << newSolPntCoords->size()<<endl;
    solPolyValsOld = frElemData->getSolPolyValsAtNode(*newSolPntCoords);
    
    delete frElemData;
    delete frElemData2;
  }

  // loop on all the elements and reassign the IDs of the states
  CFuint stateID = 0;
  for (type_itr = elementType->begin(); type_itr != elementType->end(); ++type_itr)
  {
    // get the number of control volumes (states) in this element type
    const CFuint nbStatesPerElem = getNbrOfStatesInSDCellType(type_itr->getGeoShape(),m_solPolyOrder);

    // get the number of elements
    const CFuint nbrElems = type_itr->getNbElems();

    // get start index of this element type in global element list
    CFuint globalIdx = type_itr->getStartIdx();

    cout<< " MeshUpgradeBuilder::solPolyValsOld.size()  "<< solPolyValsOld.size()<<endl;
    cout<< "  MeshUpgradeBuilder::nbStatesPerElem   "<< nbStatesPerElem << endl;
    // loop over elements of this type
    for (CFuint iElem = 0; iElem < nbrElems; ++iElem, ++globalIdx)
    {

      // if needed upgrade the low-order initialization
      if (m_upgradeInit)
      {
	cf_assert(solPolyValsOld.size() == nbStatesPerElem);
	
	const CFuint nbeq = PhysicalModelStack::getActive()->getNbEq();
	
	for (CFuint iNewState = 0; iNewState < nbStatesPerElem; ++iNewState)
        {
	  m_newStatesVal[stateID+iNewState].resize(nbeq);
	  m_newStatesVal[stateID+iNewState] = 0.0;
          for (CFuint iOldState = 0; iOldState < oldNbStatesPerElem; ++iOldState)
          {
	    const CFuint oldID = (cellStatesOld)(globalIdx,iOldState);
	    RealVector tempState = *(oldStates[oldID]);  

	    m_newStatesVal[stateID+iNewState] += solPolyValsOld[iNewState][iOldState]*tempState;
	    
	    cf_assert((solPolyValsOld[iNewState]).size() == oldNbStatesPerElem);
          }
        }
      }

      for (CFuint jState = 0; jState < nbStatesPerElem; ++jState, ++stateID)
      {
        (*cellStates)(globalIdx,jState) = stateID;

	m_elemIDOfState.push_back(globalIdx);
	m_elemLocalIDOfState.push_back(jState);
	m_updatables.push_back(updatablesElems[globalIdx]);

	if (jState == 0)
	{
	  m_globalIDs.push_back(globalIDsTemp[globalIdx][0]);
	}
	else
	{
	  m_globalIDs.push_back(oldMaxGlobalID+globalIDsTemp[globalIdx][0]*(nbStatesPerElem-1)+jState-1);//(globalIDsTemp[globalIdx][0]/oldNbStatesPerElem)*nbStatesPerElem+jState
	}
      }
    }
  }
  updatablesElems.resize(0);
  globalIDsTemp.resize(0);
  cf_assert(stateID == columnPattern.sum());
  CFLog(VERBOSE, "upgradeStateConnectivity end\n");
}

//////////////////////////////////////////////////////////////////////////////

void MeshUpgradeBuilder::recreateStates()
{
  CFAUTOTRACE;
CFLog(VERBOSE, "recreateStates start\n");
  SafePtr<MeshData::ConnTable> cellStates = MeshDataStack::getActive()->getConnectivity("cellStates_InnerCells");

  //DataHandle < Framework::State*, Framework::GLOBAL > states = getCFmeshData().getStatesHandle();
  DataHandle < Framework::State*, Framework::GLOBAL > states =
    MeshDataStack::getActive()->getStateDataSocketSink().getDataHandle();

  const CFuint newNbStates = cellStates->size();
  const CFuint oldNbStates = states.size();
  //const CFuint oldNbGlobalStates = states.getGlobalSize();

  // delete the existing states
  for (CFuint i = 0; i < oldNbStates; ++i)
  {
    //CFLog(VERBOSE, "Global ID: " << states[i]->getGlobalID() << "\n");
    deletePtr(states[i]);
  }
  IndexList<State>::getList().reset();

  // for now only serial upgrading is supported
  const bool isPara = false; //PE::GetPE().IsParallel();
  // Resize the datahandle for the states
  states.resize(newNbStates);

  // allocate the new states
  const CFuint nbeq = PhysicalModelStack::getActive()->getNbEq();
  
  for (CFuint iState = 0; iState < newNbStates; ++iState)
  {
    CFuint globalID = m_globalIDs[iState];
    CFuint localID = 0;

    if (isPara)
    {
      if (m_elemLocalIDOfState[iState] != 0)
      {
        if (m_updatables[iState])
        {
          localID = states.addLocalPoint(globalID);
        } 
        else
        {
          localID = states.addGhostPoint(globalID);
        }
      }
      else
      {
        localID = m_elemFirstStateLocalID[m_elemIDOfState[iState]];
      }
      (*cellStates)(m_elemIDOfState[iState],m_elemLocalIDOfState[iState]) = localID;
    }
    else
    {
      localID = iState;
    }

    if (!m_upgradeInit)
    {
      RealVector stateData (nbeq);
      getCFmeshData().createState(localID,stateData);
    }
    else
    {
      getCFmeshData().createState(localID,m_newStatesVal[localID]);
    }
    states[localID]->setParUpdatable(m_updatables[iState]);
    states[localID]->setGlobalID(globalID);
    
  }

  // AL: .resize(0) crashes with some compilers on some systems, better to use clear()
  m_updatables.clear();
  m_globalIDs.clear();
  m_elemLocalIDOfState.clear();
  m_elemIDOfState.clear();
  m_elemFirstStateLocalID.clear();
  m_newStatesVal.clear();
    
  //   CFuint localID = 0;
//     bool isGhost = false;
//     bool isFound = false;
//     if (hasEntry(m_localStateIDs, iState)) {
//       countLocals++;
//       localID = states.addLocalPoint (iState);
//       cf_assert(localID < nbLocalStates);
//       isFound = true;
//     }
//     else if (hasEntry(m_ghostStateIDs, iState)) {
//       countLocals++;
//       localID = states.addGhostPoint (iState);
//       cf_assert(localID < nbLocalStates);
//       isGhost = true;
//       isFound = true;
//     }
// 
//     if (isFound) {
//       State* newState = getReadData().createState
//   (localID, states.getGlobalData(localID), tmpState, !isGhost);
//       newState->setGlobalID(iState);
// 
//       
//     }
  
  CFLog(VERBOSE, "recreateStates end\n");
}

//////////////////////////////////////////////////////////////////////////////

void MeshUpgradeBuilder::upgradeNodeConnectivity()
{
  CFAUTOTRACE;

  // ADD NODE INDEXES AND UPGRADE CELLS TO NODES CONNECTIVITY
  // get the element type data
  SafePtr< vector< ElementTypeData > > elementType = getCFmeshData().getElementTypeData();
  vector< ElementTypeData >::iterator type_itr;
  const CFuint nbrElemTypes = elementType->size();

  // get cell-nodes connectivity
  SafePtr<MeshData::ConnTable> cellNodes = MeshDataStack::getActive()->getConnectivity("cellNodes_InnerCells");

  // temporary storage for original (lower order) nodes:
  const MeshData::ConnTable cellNodesLowOrder = (*cellNodes);

  // number of cells
  const CFuint nbrCells = getNbElements();

  // dimensionality
  const CFuint dim = getCFmeshData().getDimension();

  // create a mesh with equal order on all the elements
  // reset the connectivity accordingly, using a std::valarray
  std::valarray<CFuint> newPattern(nbrCells);
  CFuint nodeToCellsMapSize = 0;
  for (type_itr = elementType->begin(); type_itr != elementType->end(); ++type_itr)
  {
    // get the number of nodes in this cell type
    const CFuint newNbrNodesPerCell = type_itr->getNbNodes();

    // get the number of nodes in this cell type originally
    const CFuint nbNodesPerLowOrderCell = getNbrOfNodesInCellType(type_itr->getGeoShape(),m_prevGeoPolyOrder);

    // get the number of elements
    const CFuint nbrCellsPerType = type_itr->getNbElems();

    // loop over elements of this type
    CFuint globalIdx = type_itr->getStartIdx();
    for (CFuint iCell = 0; iCell < nbrCellsPerType; ++iCell, ++globalIdx)
    {
      newPattern[globalIdx] = newNbrNodesPerCell;
      nodeToCellsMapSize += nbNodesPerLowOrderCell;
    }
  }

  // resize cell-nodes connectivity
  cellNodes->resize(newPattern);

  // variable marking whether a cell has been updated
  // the element type is also stored here, though this is not really necessary since there is only one type of element
  vector<CFuint> isUpdated(nbrCells);

  // nodes to cells multimap
  CFMultiMap<CFuint,CFuint> mapNode2Cells;
  typedef CFMultiMap<CFuint, CFuint>::MapIterator MapIterator;
  mapNode2Cells.reserve(nodeToCellsMapSize);

  // loop over all the cells to reassign the IDs of the nodes, create nodes to cells map en initialize isUpdated variable
  m_totNbrNodes = 0;
  for (type_itr = elementType->begin(); type_itr != elementType->end(); ++type_itr)
  {
    // get the number of nodes in this cell type originally
    const CFuint nbNodesPerLowOrderCell = getNbrOfNodesInCellType(type_itr->getGeoShape(),m_prevGeoPolyOrder);

    // get the number of cells
    const CFuint nbrCellsPerType = type_itr->getNbElems();

    // loop over cells of this type
    /// @warning KVDA: here it is assumed that the original nodes are P1 nodes... (P2 --> P3: not all P2 nodes are kept)
    CFuint globalIdx = type_itr->getStartIdx();
    for (CFuint iCell = 0; iCell < nbrCellsPerType; ++iCell, ++globalIdx)
    {
      // set the original nodes and add elements to mapNode2Cells
      for (CFuint jNode = 0; jNode < nbNodesPerLowOrderCell; ++jNode)
      {
        // set nodes
        (*cellNodes)(globalIdx,jNode) = cellNodesLowOrder(globalIdx,jNode);

        // compute max node index
        m_totNbrNodes = std::max((*cellNodes)(globalIdx,jNode),m_totNbrNodes);

        // insert in nodes to cells map
        mapNode2Cells.insert((*cellNodes)(globalIdx,jNode), globalIdx);
      }

      // initialize isUpdated[globalIdx]
      isUpdated[globalIdx] = std::numeric_limits<CFuint>::max();
    }
  }
  ++m_totNbrNodes;

  // sort the nodes to cells map
  mapNode2Cells.sortKeys();

  // loop over all the cells to add the node indexes
  type_itr = elementType->begin();
  for (CFuint iType = 0; iType < nbrElemTypes; ++iType, ++type_itr)
  {
    // get geometrical shape
    const CFGeoShape::Type shape = type_itr->getGeoShape();

    // get the number of nodes in this cell type now
    const CFuint newNbrNodesPerCell = type_itr->getNbNodes();

    // get face to node and edge to node connectivity tables
    Table<CFuint>* face2NodeOld =
        LocalConnectionData::getInstance().getFaceDofLocal(shape,m_prevGeoPolyOrder,NODE,CFPolyForm::LAGRANGE);
    Table<CFuint>* face2NodeNew =
        LocalConnectionData::getInstance().getFaceDofLocal(shape,m_geoPolyOrder    ,NODE,CFPolyForm::LAGRANGE);

    Table<CFuint>* edge2NodeOld = CFNULL;
    Table<CFuint>* edge2NodeNew = CFNULL;
    CFuint nbrEdgesPerCell = 0;
    if (dim == 3)
    {
      edge2NodeOld =
          LocalConnectionData::getInstance().getEdgeDofLocal(shape,m_prevGeoPolyOrder,NODE,CFPolyForm::LAGRANGE);
      edge2NodeNew =
          LocalConnectionData::getInstance().getEdgeDofLocal(shape,m_geoPolyOrder    ,NODE,CFPolyForm::LAGRANGE);
      nbrEdgesPerCell = edge2NodeOld->nbRows();
    }

    // get the number of faces in this cell type
    const CFuint nbrFacesPerCell = LocalConnectionData::getInstance().getNbFacesInShape(shape);

    // get cell internal nodes local IDs (they do not belong to any face)
    vector<CFuint> cellIntNodeLocalIDs(0);
    for (CFuint iNode = 0; iNode < newNbrNodesPerCell; ++iNode)
    {
      bool isFaceNode = false;
      for (CFuint iFace = 0; iFace < nbrFacesPerCell && !isFaceNode; ++iFace)
      {
        const CFuint nbrNodesInFace = face2NodeNew->nbCols(iFace);
        for (CFuint iFaceNode = 0; iFaceNode < nbrNodesInFace && !isFaceNode; ++iFaceNode)
        {
          isFaceNode = iNode == (*face2NodeNew)(iFace,iFaceNode);
        }
      }
      if (!isFaceNode)
      {
        cellIntNodeLocalIDs.push_back(iNode);
      }
    }

    // number of edge nodes if 3D
    const CFuint nbrEdgeNodesNew = m_geoPolyOrder + 1;

    // number of internal nodes in cell
    const CFuint nbrIntNodesPerCell = getNbrOfInternalNodesInCellType(shape,m_geoPolyOrder);
    cf_assert(cellIntNodeLocalIDs.size() == nbrIntNodesPerCell);

    // number of internal nodes per face
    vector<CFuint> nbrIntNodesPerFace(nbrFacesPerCell);
    for (CFuint iFace = 0; iFace < nbrFacesPerCell; ++iFace)
    {
      const CFGeoShape::Type faceShape = LocalConnectionData::getInstance().getFaceShape(shape,iFace);
      nbrIntNodesPerFace[iFace] = getNbrOfInternalNodesInFaceType(faceShape,m_geoPolyOrder);
    }

    // get the number of cells
    const CFuint nbrCellsPerType = type_itr->getNbElems();

    // loop over cells of this type
    /// @warning KVDA: here it is assumed that the original nodes are P1 nodes... (P2 --> P3: not all P2 nodes are kept)
    CFuint globalIdx = type_itr->getStartIdx();
    for (CFuint iCell = 0; iCell < nbrCellsPerType; ++iCell, ++globalIdx)
    {
      // in 3D, loop over edges to add high-order edge nodes
      if (dim == 3)
      {
        // loop over edges to add edge nodes
        for (CFuint iEdge = 0; iEdge < nbrEdgesPerCell; ++iEdge)
        {
          cf_assert(edge2NodeNew->nbCols(iEdge) == nbrEdgeNodesNew);

          // find (updated) cell that shares this edge
          CFuint updEdgeCellIdx = std::numeric_limits<CFuint>::max();
          const CFuint firstNodeLocalID = (*edge2NodeOld)(iEdge,0);
          const CFuint firstNodeID      = (*cellNodes)(globalIdx,firstNodeLocalID);
          
	  bool fo = false;
	  pair<MapIterator, MapIterator> nghbrCellsOfNode1 = mapNode2Cells.find(firstNodeID,fo);
          cf_assert(fo);
	  
	  const CFuint secondNodeLocalID = (*edge2NodeOld)(iEdge,1);
          const CFuint secondNodeID      = (*cellNodes)(globalIdx,secondNodeLocalID);
	  
	  bool fo1 = false;
	  pair<MapIterator, MapIterator> nghbrCellsOfNode2 = mapNode2Cells.find(secondNodeID, fo1);
	  cf_assert(fo1);
	  
          for (MapIterator nghbrCellItr1 = nghbrCellsOfNode1.first;
               nghbrCellItr1 != nghbrCellsOfNode1.second && updEdgeCellIdx == std::numeric_limits<CFuint>::max();
               ++nghbrCellItr1)
          {
            const CFuint currCellIdx = nghbrCellItr1->second;
            if(currCellIdx != globalIdx)
            {
              for (MapIterator nghbrCellItr2 = nghbrCellsOfNode2.first;
                   nghbrCellItr2 != nghbrCellsOfNode2.second;
                   ++nghbrCellItr2)
              {
                const CFuint otherCellIdx = nghbrCellItr2->second;
                if(otherCellIdx == currCellIdx && isUpdated[currCellIdx] !=  std::numeric_limits<CFuint>::max())
                {
                  updEdgeCellIdx = currCellIdx;
                  break;
                }
              }
            }
          }

          // check if updated edge-neighbour cell was found
          /// @warning KVDA: here it is assumed that the original nodes are P1 nodes (--> 2 nodes per edge)
          if (updEdgeCellIdx == std::numeric_limits<CFuint>::max())
          // did not find updated edge-neighbour cell
          {
            // add new nodes
            for (CFuint iNode = 2; iNode < nbrEdgeNodesNew; ++iNode)
            {
              const CFuint localNodeID = (*edge2NodeNew)(iEdge,iNode);
              (*cellNodes)(globalIdx,localNodeID) = m_totNbrNodes;
//               CFLog(NOTICE,"Edge node (new):\n");
//               CF_DEBUG_OBJ(localNodeID);
//               CF_DEBUG_OBJ((*cellNodes)(globalIdx,localNodeID));
              ++m_totNbrNodes;
            }
          }
          else
          // found updated edge-neighbour cell
          {
            // get neighbour cell type
            const CFuint nghbrCellType = isUpdated[updEdgeCellIdx];

            // get neighbour geometrical shape
            const CFGeoShape::Type nghbrShape = (*elementType)[nghbrCellType].getGeoShape();

            // get neighbour cell old local connectivity
            Table<CFuint>* edge2NodeNghbr =
                LocalConnectionData::getInstance().getEdgeDofLocal(nghbrShape,m_geoPolyOrder,NODE,CFPolyForm::LAGRANGE);

            // find edge-neighbour cell edge local ID and add node indexes to cell-node connectivity
            bool foundEdge = false;
            const CFuint nbrNghbrEdges = edge2NodeNghbr->nbRows();
            for (CFuint iEdge2 = 0; iEdge2 < nbrNghbrEdges && !foundEdge; ++iEdge2)
            {
              const CFuint firstNghbrNodeLocalID = (*edge2NodeNghbr)(iEdge2        ,0                    );
              const CFuint firstNghbrNodeID      = (*cellNodes     )(updEdgeCellIdx,firstNghbrNodeLocalID);
              if (firstNghbrNodeID == firstNodeID)
              {
                const CFuint secondNghbrNodeLocalID = (*edge2NodeNghbr)(iEdge2        ,1                     );
                const CFuint secondNghbrNodeID      = (*cellNodes     )(updEdgeCellIdx,secondNghbrNodeLocalID);
                if (secondNghbrNodeID == secondNodeID)
                {
                  cf_assert(edge2NodeNghbr->nbCols(iEdge2) == nbrEdgeNodesNew);
                  foundEdge = true;
                  CFuint nodeEdgeIdx = 2;
                  for (CFuint iNode = 2; iNode < nbrEdgeNodesNew; ++iNode, ++nodeEdgeIdx)
                  {
                    const CFuint nghbrLocalNodeID = (*edge2NodeNghbr)(iEdge2,nodeEdgeIdx);
                    const CFuint nodeID = (*cellNodes)(updEdgeCellIdx,nghbrLocalNodeID);
                    const CFuint localNodeID = (*edge2NodeNew)(iEdge,iNode);
                    (*cellNodes)(globalIdx,localNodeID) = nodeID;
//                     CFLog(NOTICE,"Edge node:\n");
//                     CF_DEBUG_OBJ(localNodeID);
//                     CF_DEBUG_OBJ((*cellNodes)(globalIdx,localNodeID));
                  }
                }
              }
              else if (firstNghbrNodeID == secondNodeID)
              {
                const CFuint secondNghbrNodeLocalID = (*edge2NodeNghbr)(iEdge2        ,1                     );
                const CFuint secondNghbrNodeID      = (*cellNodes     )(updEdgeCellIdx,secondNghbrNodeLocalID);
                if (secondNghbrNodeID == firstNodeID)
                {
                  cf_assert(edge2NodeNghbr->nbCols(iEdge2) == nbrEdgeNodesNew);
                  foundEdge = true;
                  CFuint nodeEdgeIdx = nbrEdgeNodesNew-1;
                  for (CFuint iNode = 2; iNode < nbrEdgeNodesNew; ++iNode, --nodeEdgeIdx)
                  {
                    const CFuint nghbrLocalNodeID = (*edge2NodeNghbr)(iEdge2,nodeEdgeIdx);
                    const CFuint nodeID = (*cellNodes)(updEdgeCellIdx,nghbrLocalNodeID);
                    const CFuint localNodeID = (*edge2NodeNew)(iEdge,iNode);
                    (*cellNodes)(globalIdx,localNodeID) = nodeID;
//                     CFLog(NOTICE,"Edge node:\n");
//                     CF_DEBUG_OBJ(localNodeID);
//                     CF_DEBUG_OBJ((*cellNodes)(globalIdx,localNodeID));
                  }
                }
              }
            }
            cf_assert(foundEdge);
          }
        }
      }

      // loop over faces to add `internal' face nodes
      for (CFuint iFace = 0; iFace < nbrFacesPerCell; ++iFace)
      {
        if (nbrIntNodesPerFace[iFace] > 0)
        {
          // number of face nodes
          const CFuint nbrNodesPerFaceOld = face2NodeOld->nbCols(iFace);
          const CFuint nbrNodesPerFaceNew = face2NodeNew->nbCols(iFace);

          // find neighbouring cell to this face
          CFuint nghbrCellIdx = std::numeric_limits<CFuint>::max();
          const CFuint firstNodeLocalID = (*face2NodeOld)(iFace, 0);
          const CFuint firstNodeID = (*cellNodes)(globalIdx,firstNodeLocalID);
	  
	  bool fo = false;
          pair<MapIterator, MapIterator> nghbrCellsOfNode0 = mapNode2Cells.find(firstNodeID, fo);
	  cf_assert(fo);
	  
          for (MapIterator nghbrCellItr1 = nghbrCellsOfNode0.first;
               nghbrCellItr1 != nghbrCellsOfNode0.second && nghbrCellIdx == std::numeric_limits<CFuint>::max();
               ++nghbrCellItr1)
          {
            CFuint nbCommonNodes = 1;
            const CFuint currCellIdx = nghbrCellItr1->second;
            if(currCellIdx != globalIdx)
            {
              for (CFuint iNode = 1; iNode < nbrNodesPerFaceOld; ++iNode)
              {
                const CFuint localNodeID = (*face2NodeOld)(iFace, iNode);
                const CFuint nodeID = (*cellNodes)(globalIdx, localNodeID);
                
		bool fo = false;
		pair<MapIterator, MapIterator> nghbrCellsOfNode = mapNode2Cells.find(nodeID, fo);
		cf_assert(fo);
		
                for (MapIterator nghbrCellItr2 = nghbrCellsOfNode.first;
                     nghbrCellItr2 != nghbrCellsOfNode.second;
                     ++nghbrCellItr2)
                {
                  const CFuint otherCellIdx = nghbrCellItr2->second;
                  if(otherCellIdx == currCellIdx)
                  {
                    ++nbCommonNodes;
                    break;
                  }
                }
              }
              if(nbCommonNodes == nbrNodesPerFaceOld)
              {
                nghbrCellIdx = currCellIdx;
              }
            }
          }

          // check if neighbouring cell was found or current face is a boundary face
          const bool isBndFace = nghbrCellIdx == std::numeric_limits<CFuint>::max();
          if (isBndFace)
          // boundary face
          {
            // add new nodes
            CFuint iNode = nbrNodesPerFaceNew - nbrIntNodesPerFace[iFace];
            for (CFuint iIntNode = 0; iIntNode < nbrIntNodesPerFace[iFace]; ++iIntNode, ++iNode)
            {
              const CFuint localNodeID = (*face2NodeNew)(iFace,iNode);
              (*cellNodes)(globalIdx,localNodeID) = m_totNbrNodes;
//               CFLog(NOTICE,"Boundary face node (new):\n");
//               CF_DEBUG_OBJ(localNodeID);
//               CF_DEBUG_OBJ((*cellNodes)(globalIdx,localNodeID));
              ++m_totNbrNodes;
            }

            // save boundary face connectivity
            vector<CFuint*> currBndFaceNodes;
            for (iNode = 0; iNode < nbrNodesPerFaceNew; ++iNode)
            {
              const CFuint localNodeID = (*face2NodeNew)(iFace,iNode);
              currBndFaceNodes.push_back(&(*cellNodes)(globalIdx,localNodeID));
            }
            m_bndFacesNodes.push_back(currBndFaceNodes);
          }
          else
          // internal face
          {
            // check if neighbouring cell has already been updated
            const bool nghbrCellIsUpdated = isUpdated[nghbrCellIdx] !=  std::numeric_limits<CFuint>::max();
            if (nghbrCellIsUpdated)
            {
              // get neighbour cell type
              const CFuint nghbrCellType = isUpdated[nghbrCellIdx];

              // get neighbour geometrical shape
              const CFGeoShape::Type nghbrShape = (*elementType)[nghbrCellType].getGeoShape();

              // get neighbour cell old local connectivity
              Table<CFuint>* face2NodeNghbr =
                  LocalConnectionData::getInstance().getFaceDofLocal(nghbrShape,m_prevGeoPolyOrder,NODE,CFPolyForm::LAGRANGE);

              // find neighbour cell face that matches the current face of the current cell
              CFuint nghbrFaceLocalID = std::numeric_limits<CFuint>::max();
              const CFuint nbrNghbrFaces = face2NodeNghbr->nbRows();
              for(CFuint iNghbrFace = 0;
                  iNghbrFace < nbrNghbrFaces && nghbrFaceLocalID == std::numeric_limits<CFuint>::max();
                  ++iNghbrFace)
              {
                const CFuint nbrNodesInNghbrFaceOld = face2NodeNghbr->nbCols(iNghbrFace);
                if (nbrNodesInNghbrFaceOld == nbrNodesPerFaceOld)
                {
                  CFuint nbrMatchNodes = 0;
                  for(CFuint iNodeNghbr = 0; iNodeNghbr < nbrNodesInNghbrFaceOld; ++iNodeNghbr)
                  {
                    const CFuint nghbrNodeLocalID = (*face2NodeNghbr)(iNghbrFace,iNodeNghbr);
                    const CFuint nghbrNodeID = (*cellNodes)(nghbrCellIdx,nghbrNodeLocalID);
                    for(CFuint iNodeCurr = 0; iNodeCurr < nbrNodesPerFaceOld; ++iNodeCurr)
                    {
                      const CFuint currNodeLocalID = (*face2NodeOld)(iFace,iNodeCurr);
                      const CFuint currNodeID = (*cellNodes)(globalIdx,currNodeLocalID);
                      if(currNodeID == nghbrNodeID)
                      {
                        ++nbrMatchNodes;
                        break;
                      }
                    }
                  }
                  if(nbrMatchNodes == nbrNodesPerFaceOld)
                  {
                    nghbrFaceLocalID = iNghbrFace;
                  }
                }
              }
              cf_assert(nghbrFaceLocalID != std::numeric_limits<CFuint>::max());

              // share `internal' face nodes with neighbouring cell face
              /// @warning KVDA: here, a P2 geometrical element is assumed
              /// --> only one `internal' face node, which is last in the local face to node connectivity
              cf_assert(nbrIntNodesPerFace[iFace] == 1);

              // get neighbour cell new local connectivity
              Table<CFuint>* face2NodeNghbrNew =
                  LocalConnectionData::getInstance().getFaceDofLocal(nghbrShape,m_geoPolyOrder,NODE,CFPolyForm::LAGRANGE);

              // get node global ID
              const CFuint nghbrLocalNodeID = (*face2NodeNghbrNew)(nghbrFaceLocalID,nbrNodesPerFaceNew-1);
              const CFuint nodeID = (*cellNodes)(nghbrCellIdx,nghbrLocalNodeID);

              // set current face `internal' node
              const CFuint localNodeID = (*face2NodeNew)(iFace,nbrNodesPerFaceNew-1);
              (*cellNodes)(globalIdx,localNodeID) = nodeID;
//               CFLog(NOTICE,"Face node:\n");
//               CF_DEBUG_OBJ(localNodeID);
//               CF_DEBUG_OBJ((*cellNodes)(globalIdx,localNodeID));
            }
            else
            {
              // add new nodes
              CFuint iNode = nbrNodesPerFaceNew - nbrIntNodesPerFace[iFace];
              for (CFuint iIntNode = 0; iIntNode < nbrIntNodesPerFace[iFace]; ++iIntNode, ++iNode)
              {
                const CFuint localNodeID = (*face2NodeNew)(iFace,iNode);
                (*cellNodes)(globalIdx,localNodeID) = m_totNbrNodes;
//                 CFLog(NOTICE,"Face node (new):\n");
//                 CF_DEBUG_OBJ(localNodeID);
//                 CF_DEBUG_OBJ((*cellNodes)(globalIdx,localNodeID));
                ++m_totNbrNodes;
              }
            }
          }
        }
      }

      // add internal cell node indexes
      for (CFuint iNode = 0; iNode < nbrIntNodesPerCell; ++iNode)
      {
        (*cellNodes)(globalIdx,cellIntNodeLocalIDs[iNode]) = m_totNbrNodes;
//         CFLog(NOTICE,"Cell node (new):\n");
//         CF_DEBUG_OBJ(cellIntNodeLocalIDs[iNode]);
//         CF_DEBUG_OBJ((*cellNodes)(globalIdx,cellIntNodeLocalIDs[iNode]));
        ++m_totNbrNodes;
      }

      // set isUpdated[globalIdx]
      isUpdated[globalIdx] = iType;
    }
  }

  // UPGRADE BOUNDARY FACES TO NODES CONNECTIVITY
  // total number of boundary faces
  const CFuint totNbrBndFaces = m_bndFacesNodes.size();

  // get number of Topological Region Sets
  const CFuint nbTRSs = getCFmeshData().getNbTRSs();

  // get number of Topological Regions per TRS
  SafePtr< vector<CFuint> > nbTRs = getCFmeshData().getNbTRs();

  // get number of geometrical entities per TR
  SafePtr< vector< vector<CFuint> > > nbGeomEntsPerTR = getCFmeshData().getNbGeomEntsPerTR();

  // loop over TRSs
  for (CFuint iTRS = 0; iTRS < nbTRSs; ++iTRS)
  {
    // get current TRS TRGeoConn
    TRGeoConn& currTRSGeoConnPerTR = getCFmeshData().getTRGeoConn(iTRS);

    // loop over TRs
    const CFuint nbrTRs = (*nbTRs)[iTRS];
    for (CFuint iTR = 0; iTR < nbrTRs; ++iTR)
    {
      // loop over faces
      const CFuint nbrFaces = (*nbGeomEntsPerTR)[iTRS][iTR];
      for (CFuint iFace = 0; iFace < nbrFaces; ++iFace)
      {
        // get face to nodes connectivity
        GeoConnElementPart& face2Nodes = currTRSGeoConnPerTR[iTR][iFace].first;
        const CFuint nbrFaceNodesOld = face2Nodes.size();

        // find face in m_bndFacesNodes
        CFuint bndFaceIdx = std::numeric_limits<CFuint>::max();
        for (CFuint iFace2 = 0;
             iFace2 < totNbrBndFaces && bndFaceIdx == std::numeric_limits<CFuint>::max();
             ++iFace2)
        {
          if (m_bndFacesNodes[iFace2].size() >= nbrFaceNodesOld)
          {
            CFuint nbrMatchNodes = 0;
            for (CFuint iNode = 0; iNode < nbrFaceNodesOld; ++iNode)
            {
              const CFuint nodeID = face2Nodes[iNode];
              for (CFuint iNode2 = 0; iNode2 < nbrFaceNodesOld; ++iNode2)
              {
                const CFuint nodeID2 = *m_bndFacesNodes[iFace2][iNode2];
                if (nodeID == nodeID2)
                {
                  ++nbrMatchNodes;
                  break;
                }
              }
            }
            if (nbrMatchNodes == nbrFaceNodesOld)
            {
              bndFaceIdx = iFace2;
            }
          }
        }
        cf_assert(bndFaceIdx != std::numeric_limits<CFuint>::max());

        // change face to node connectivity
        const CFuint nbrFaceNodesNew = m_bndFacesNodes[bndFaceIdx].size();
        face2Nodes.resize(nbrFaceNodesNew);
        for (CFuint iNode = 0; iNode < nbrFaceNodesNew; ++iNode)
        {
          face2Nodes[iNode] = *m_bndFacesNodes[bndFaceIdx][iNode];
        }
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void MeshUpgradeBuilder::recreateNodes()
{
  CFAUTOTRACE;

  // get (upgraded) cell to node connectivity
  SafePtr<MeshData::ConnTable> cellNodes = MeshDataStack::getActive()->getConnectivity("cellNodes_InnerCells");

  // get the nodes
  DataHandle < Framework::Node*, Framework::GLOBAL > nodes = getCFmeshData().getNodesHandle();

  // old and new number of nodes
  const CFuint newNbrNodes = m_totNbrNodes;
  const CFuint oldNbrNodes = nodes.size();

  // backup the existing nodes
  vector<Node*> oldNodes(oldNbrNodes);
  for (CFuint iNode = 0; iNode < oldNbrNodes; ++iNode)
  {
    oldNodes[iNode] = nodes[iNode];
  }

  // resize the datahandle for the nodes and copy the old nodes back
  nodes.resize(newNbrNodes);
  for (CFuint iNode = 0; iNode < oldNbrNodes; ++iNode)
  {
    nodes[iNode] = oldNodes[iNode];
    nodes[iNode]->setIsOwnedByState(false);
    nodes[iNode]->setIsOnMesh(true);
  }

  // number of nodes that need to be added
  const CFuint nbrNewNodes = newNbrNodes - oldNbrNodes;

  // booleans to keep track of the already created nodes
  vector<bool> isCreated(nbrNewNodes, false);
  vector<bool> secondTimeFace(newNbrNodes,false);

  // get the element type data
  SafePtr< vector< ElementTypeData > > elementType = getCFmeshData().getElementTypeData();
  vector< ElementTypeData >::iterator type_itr;
  const CFuint nbrElemTypes = elementType->size();

  // dimensionality
  const CFuint dim = getCFmeshData().getDimension();

  // variable for coordinates
  RealVector coords(dim);

  // loop over all the cells to add the node
  type_itr = elementType->begin();
  for (CFuint iType = 0; iType < nbrElemTypes; ++iType, ++type_itr)
  {
    // get geometrical shape
    const CFGeoShape::Type shape = type_itr->getGeoShape();

    // get the number of nodes in this cell type originally
    const CFuint nbrNodesPerCellOld = getNbrOfNodesInCellType(shape,m_prevGeoPolyOrder);

    // get the number of nodes in this cell type now
    const CFuint newNbrNodesPerCell = type_itr->getNbNodes();

    // get face to node and edge to node connectivity tables
    Table<CFuint>* face2NodeOld =
        LocalConnectionData::getInstance().getFaceDofLocal(shape,m_prevGeoPolyOrder,NODE,CFPolyForm::LAGRANGE);
    Table<CFuint>* face2NodeNew =
        LocalConnectionData::getInstance().getFaceDofLocal(shape,m_geoPolyOrder    ,NODE,CFPolyForm::LAGRANGE);
    Table<CFuint>* edge2NodeOld = CFNULL;
    Table<CFuint>* edge2NodeNew = CFNULL;
    CFuint nbrEdgesPerCell = 0;
    if (dim == 3)
    {
      edge2NodeOld =
          LocalConnectionData::getInstance().getEdgeDofLocal(shape,m_prevGeoPolyOrder,NODE,CFPolyForm::LAGRANGE);
      edge2NodeNew =
          LocalConnectionData::getInstance().getEdgeDofLocal(shape,m_geoPolyOrder    ,NODE,CFPolyForm::LAGRANGE);
      nbrEdgesPerCell = edge2NodeOld->nbRows();
    }

    // get the number of faces in this cell type
    const CFuint nbrFacesPerCell = LocalConnectionData::getInstance().getNbFacesInShape(shape);

    // get cell internal nodes local IDs (they do not belong to any face)
    vector<CFuint> cellIntNodeLocalIDs(0);
    for (CFuint iNode = 0; iNode < newNbrNodesPerCell; ++iNode)
    {
      bool isFaceNode = false;
      for (CFuint iFace = 0; iFace < nbrFacesPerCell && !isFaceNode; ++iFace)
      {
        const CFuint nbrNodesInFace = face2NodeNew->nbCols(iFace);
        for (CFuint iFaceNode = 0; iFaceNode < nbrNodesInFace && !isFaceNode; ++iFaceNode)
        {
          isFaceNode = iNode == (*face2NodeNew)(iFace,iFaceNode);
        }
      }
      if (!isFaceNode)
      {
        cellIntNodeLocalIDs.push_back(iNode);
      }
    }

    // number of edge nodes if 3D
    const CFuint nbrEdgeNodesNew = m_geoPolyOrder + 1;

    // number of internal nodes in cell
    const CFuint nbrIntNodesPerCell = getNbrOfInternalNodesInCellType(shape,m_geoPolyOrder);
    cf_assert(cellIntNodeLocalIDs.size() == nbrIntNodesPerCell);

    // number of internal nodes per face
    vector<CFuint> nbrIntNodesPerFace(nbrFacesPerCell);
    for (CFuint iFace = 0; iFace < nbrFacesPerCell; ++iFace)
    {
      const CFGeoShape::Type faceShape = LocalConnectionData::getInstance().getFaceShape(shape,iFace);
      nbrIntNodesPerFace[iFace] = getNbrOfInternalNodesInFaceType(faceShape,m_geoPolyOrder);
    }

    // get the number of cells
    const CFuint nbrCellsPerType = type_itr->getNbElems();

    // loop over cells of this type
    /// @warning KVDA: here it is assumed that the original nodes are P1 nodes... (P2 --> P3: not all P2 nodes are kept)
    CFuint globalIdx = type_itr->getStartIdx();
    for (CFuint iCell = 0; iCell < nbrCellsPerType; ++iCell, ++globalIdx)
    {
      // in 3D, loop over edges to add high-order edge nodes
      if (dim == 3)
      {
        // loop over edges to add edge nodes
        for (CFuint iEdge = 0; iEdge < nbrEdgesPerCell; ++iEdge)
        {
          cf_assert(edge2NodeNew->nbCols(iEdge) == nbrEdgeNodesNew);

          /// @warning KVDA: here, a P2 geometrical element is assumed
          /// --> only three nodes on each edge, last in the local edge to node connectivity
          /// and the third node in the middle of the edge
          /// extension to higher order: local mapped coordinate of the high-order nodes on the edge are needed
          cf_assert(nbrEdgeNodesNew == 3);

          // get current node ID
          const CFuint nodeLocalID = (*edge2NodeNew)(iEdge    ,2          );
          const CFuint nodeID      = (*cellNodes   )(globalIdx,nodeLocalID);

          // add current edge middle node if not created yet
          if (!isCreated[nodeID-oldNbrNodes])
          {
            // compute new node coordinates (for P2 element!!!)
            const CFuint vNode1LocalID = (*edge2NodeNew)(iEdge    ,0            );
            const CFuint vNode1ID      = (*cellNodes   )(globalIdx,vNode1LocalID);
            const CFuint vNode2LocalID = (*edge2NodeNew)(iEdge    ,1            );
            const CFuint vNode2ID      = (*cellNodes   )(globalIdx,vNode2LocalID);
            coords = 0.5*(*nodes[vNode1ID] + *nodes[vNode2ID]);

            // create the node in the mesh data
            getCFmeshData().createNode(nodeID,coords);
            nodes[nodeID]->setGlobalID(nodeID);
            nodes[nodeID]->setIsOnMesh(true);
            nodes[nodeID]->setIsOwnedByState(false);

            // set isCreated[nodeID-oldNbrNodes]
            isCreated[nodeID-oldNbrNodes] = true;
          }
        }
      }

      // loop over faces to add `internal' face nodes
      for (CFuint iFace = 0; iFace < nbrFacesPerCell; ++iFace)
      {
        if (nbrIntNodesPerFace[iFace] > 0)
        {
          // number of face nodes
          const CFuint nbrNodesPerFaceOld = face2NodeOld->nbCols(iFace);
          const CFuint nbrNodesPerFaceNew = face2NodeNew->nbCols(iFace);

          /// @warning KVDA: here, a P2 geometrical element is assumed
          /// --> only one `internal' face node, which is last in the local face to node connectivity
          /// and in the center of the face
          /// extension to higher order: local mapped coordinates of the high-order nodes in the face are needed
          cf_assert(nbrIntNodesPerFace[iFace] == 1);

          // get current node ID
          const CFuint nodeLocalID = (*face2NodeNew)(iFace    ,nbrNodesPerFaceNew-1);
          const CFuint nodeID      = (*cellNodes   )(globalIdx,nodeLocalID         );

          // add current face `internal' node if not created yet
          if (!isCreated[nodeID-oldNbrNodes])
          {
            // compute new node coordinates (for P2 element!!!)
            coords = 0.0;
            for (CFuint iNode = 0; iNode < nbrNodesPerFaceOld; ++iNode)
            {
              const CFuint vNodeLocalID = (*face2NodeOld)(iFace,iNode);
              const CFuint vNodeID      = (*cellNodes)(globalIdx,vNodeLocalID);
              coords += *nodes[vNodeID];
            }
            coords /= nbrNodesPerFaceOld;

            // create the node in the mesh data
            getCFmeshData().createNode(nodeID,coords);
            nodes[nodeID]->setGlobalID(nodeID);
            nodes[nodeID]->setIsOnMesh(true);
            nodes[nodeID]->setIsOwnedByState(false);

            // set isCreated[nodeID-oldNbrNodes]
            isCreated[nodeID-oldNbrNodes] = true;
          }
        }
      }

      // add internal cell nodes
      /// @warning KVDA: here, a P2 geometrical element is assumed
      /// --> only one internal cell node, which is in the center of the cell
      /// extension to higher order: local mapped coordinates of the high-order nodes in the face are needed
      cf_assert(1 <= nbrIntNodesPerCell);
      for (CFuint iNode = 0; iNode < nbrIntNodesPerCell; ++iNode)
      {
        const CFuint nodeID = (*cellNodes)(globalIdx,cellIntNodeLocalIDs[iNode]);

        // no check for already created necessary, only the current cell has this (internal) node
        // compute new node coordinates (for P2 element!!!)
        coords = 0.0;
        for (CFuint iNode = 0; iNode < nbrNodesPerCellOld; ++iNode)
        {
          const CFuint nodeID = (*cellNodes)(globalIdx,iNode);
          coords += *nodes[nodeID];
        }
        coords /= nbrNodesPerCellOld;

        // create the node in the mesh data
        getCFmeshData().createNode(nodeID,coords);
        nodes[nodeID]->setGlobalID(nodeID);
        nodes[nodeID]->setIsOnMesh(true);
        nodes[nodeID]->setIsOwnedByState(false);

        // set isCreated[nodeID-oldNbrNodes]
        isCreated[nodeID-oldNbrNodes] = true;
      }
    }
  }

//   CFuint nbrNotCreated = 0;
//   CFuint nbrCreated = 0;
//   for (CFuint iNode = 0; iNode < isCreated.size(); ++iNode)
//   {
//     if (!isCreated[iNode])
//     {
//       ++nbrNotCreated;
//     }
//     else
//     {
//       ++nbrCreated;
//     }
// //     cf_assert(isCreated[iNode]);
//   }
//   CF_DEBUG_OBJ(isCreated.size());
//   CF_DEBUG_OBJ(nbrNotCreated);
//   CF_DEBUG_OBJ(nbrCreated);
//   CF_DEBUG_EXIT;
}

//////////////////////////////////////////////////////////////////////////////

vector< RealVector > MeshUpgradeBuilder::getNewNodesCoords(CFGeoShape::Type shape, CFuint solOrder, CFuint cellIdx, MeshData::ConnTable& nodesConn)
{
  // output variable
  vector< RealVector > nodeMappedCoords;
  
  // get the nodes
  DataHandle < Framework::Node*, Framework::GLOBAL > nodes = getCFmeshData().getNodesHandle();

  // for zeroth order solution polynomial, use same points as for first order
  if (solOrder == 0)
  {
    solOrder = 1;
  }

  // number of points needed for representing a polynomial of order degree solOrder
  const CFuint nbrNodes1D = solOrder + 1;
  
  const CFreal nbFaceParts = solOrder;

  switch (shape)
  {
    case CFGeoShape::LINE:
    {
      throw Common::NotImplementedException (FromHere(),"MeshUpgradeBuilder::getOutputPntsMappedCoords() for LINE");
//       for (CFuint iKsi = 0; iKsi < nbrNodes1D; ++iKsi)
//       {
//         RealVector coords(1);
//         coords[KSI] = -1.0 + iKsi*2.0/solOrder;
//         nodeMappedCoords.push_back(coords);
//       }
    } break;
    case CFGeoShape::TRIAG:
    {
      throw Common::NotImplementedException (FromHere(),"MeshUpgradeBuilder::getOutputPntsMappedCoords() for PYRAM");
//       for (CFuint iKsi = 0; iKsi < nbrNodes1D; ++iKsi)
//       {
//         const CFreal ksi = iKsi*1.0/solOrder;
//         for (CFuint iEta = 0; iEta < nbrNodes1D-iKsi; ++iEta)
//         {
//           RealVector coords(2);
//           coords[KSI] = ksi;
//           coords[ETA] = iEta*1.0/solOrder;
//           nodeMappedCoords.push_back(coords);
//         }
//       }
    } break;
    case CFGeoShape::QUAD:
    {
      RealVector nodeIDs(4);
      vector< Node* > nodeCoordsTemp;
      vector< Node* > nodeCoords;
      nodeIDs[0] = nodesConn(cellIdx,0);
      nodeCoordsTemp.push_back(nodes[nodeIDs[0]]);
      CFreal maxSumXY = (*nodeCoordsTemp[0])[0]+(*nodeCoordsTemp[0])[1];
      CFreal minSumXY = (*nodeCoordsTemp[0])[0]+(*nodeCoordsTemp[0])[1];
      CFreal minDiffYX = (*nodeCoordsTemp[0])[1]-(*nodeCoordsTemp[0])[0];
      CFuint minXY = 0;
      CFuint maxXY = 0;
      CFuint minDiff = 0;
      RealVector newNodeOrder(4);
      for (CFuint iNode = 1; iNode < 4; ++iNode)
      {
        nodeIDs[iNode] = nodesConn(cellIdx,iNode);
	nodeCoordsTemp.push_back(nodes[nodeIDs[iNode]]);
	if ((*nodeCoordsTemp[iNode])[0]+(*nodeCoordsTemp[iNode])[1] > maxSumXY)
	{
	  maxXY = iNode;
	  maxSumXY = (*nodeCoordsTemp[iNode])[0]+(*nodeCoordsTemp[iNode])[1];
	}
	if ((*nodeCoordsTemp[iNode])[0]+(*nodeCoordsTemp[iNode])[1] < minSumXY)
	{
	  minXY = iNode;
	  minSumXY = (*nodeCoordsTemp[iNode])[0]+(*nodeCoordsTemp[iNode])[1];
	}
	if ((*nodeCoordsTemp[iNode])[1]-(*nodeCoordsTemp[iNode])[0] < minDiffYX)
	{
	  minDiff = iNode;
	  minDiffYX = (*nodeCoordsTemp[iNode])[1]-(*nodeCoordsTemp[iNode])[0];
	}
      }

      for (CFuint iNode = 0; iNode < 4; ++iNode)
      {
	if (iNode == maxXY)
	{
	  newNodeOrder[2] = iNode;
	}
	else if (iNode == minXY)
	{
	  newNodeOrder[0] = iNode;
	}
	else if (iNode == minDiff)
	{
	  newNodeOrder[1] = iNode;
	}
	else 
	{
	  newNodeOrder[3] = iNode;
	}
      }

      for (CFuint iNode = 0; iNode < 4; ++iNode)
      {
        nodesConn(cellIdx,iNode) = nodeIDs[newNodeOrder[iNode]];
	nodeCoords.push_back(nodes[nodeIDs[newNodeOrder[iNode]]]);
	CFLog(VERBOSE,"nodes: " << (*nodeCoords[iNode]) << "\n");
      }

      for (CFuint iEta = 0; iEta < nbrNodes1D; ++iEta)
      {
        const CFreal xLeft = (*nodeCoords[0])[0]+iEta/nbFaceParts*((*nodeCoords[3])[0]-(*nodeCoords[0])[0]);
	const CFreal xRight = (*nodeCoords[1])[0]+iEta/nbFaceParts*((*nodeCoords[2])[0]-(*nodeCoords[1])[0]);
        const CFreal yLeft = (*nodeCoords[0])[1]+iEta/nbFaceParts*((*nodeCoords[3])[1]-(*nodeCoords[0])[1]);
	const CFreal yRight = (*nodeCoords[1])[1]+iEta/nbFaceParts*((*nodeCoords[2])[1]-(*nodeCoords[1])[1]);
	
        for (CFuint iKsi = 0; iKsi < nbrNodes1D; ++iKsi)
        { 
	  const CFreal xDown = (*nodeCoords[0])[0]+iKsi/nbFaceParts*((*nodeCoords[1])[0]-(*nodeCoords[0])[0]);
	  const CFreal xUp = (*nodeCoords[3])[0]+iKsi/nbFaceParts*((*nodeCoords[2])[0]-(*nodeCoords[3])[0]);
	  const CFreal yDown = (*nodeCoords[0])[1]+iKsi/nbFaceParts*((*nodeCoords[1])[1]-(*nodeCoords[0])[1]);
	  const CFreal yUp = (*nodeCoords[3])[1]+iKsi/nbFaceParts*((*nodeCoords[2])[1]-(*nodeCoords[3])[1]);
	  const CFreal mLR = (yRight-yLeft)/(xRight-xLeft);
	  CFreal x;
	  if ((xUp-xDown) == 0)
	  {
	    x = xDown;
	  }
	  else 
	  {
	    const CFreal mUD = (yUp-yDown)/(xUp-xDown);
            x = (mLR*xLeft-mUD*xDown+yDown-yLeft)/(mLR-mUD);
	  }
	  const CFreal y = mLR*x+yLeft-xLeft*mLR;
	  //CFLog(VERBOSE,"x = " << x << ", y = " << y << ", for xL = " << xLeft << ", xR = " << xRight << ", yL = " << yLeft << ", yR = " << yRight << ", xD = " << xDown << ", xU = " << xUp << ", yD = " << yDown << ", yU = " << yUp << "\n");
	  RealVector currCoords(2);
	  currCoords[0] = x;
	  currCoords[1] = y;
	  nodeMappedCoords.push_back(currCoords);
        }
      }
    } break;
    case CFGeoShape::TETRA:
    {
      throw Common::NotImplementedException (FromHere(),"MeshUpgradeBuilder::getOutputPntsMappedCoords() for TETRA");
      /// @warn: for tetra this is only implemented for P1
//       RealVector coords(3);
//       coords[KSI] = 0.0;
//       coords[ETA] = 0.0;
//       coords[ZTA] = 0.0;
//       nodeMappedCoords.push_back(coords);
//       coords[KSI] = 1.0;
//       coords[ETA] = 0.0;
//       coords[ZTA] = 0.0;
//       nodeMappedCoords.push_back(coords);
//       coords[KSI] = 0.0;
//       coords[ETA] = 1.0;
//       coords[ZTA] = 0.0;
//       nodeMappedCoords.push_back(coords);
//       coords[KSI] = 0.0;
//       coords[ETA] = 0.0;
//       coords[ZTA] = 1.0;
//       nodeMappedCoords.push_back(coords);
    } break;
    case CFGeoShape::PYRAM:
    {
      throw Common::NotImplementedException (FromHere(),"MeshUpgradeBuilder::getOutputPntsMappedCoords() for PYRAM");
    } break;
    case CFGeoShape::PRISM:
    {
      throw Common::NotImplementedException (FromHere(),"MeshUpgradeBuilder::getOutputPntsMappedCoords() for PRISM");
    } break;
    case CFGeoShape::HEXA:
    {
      throw Common::NotImplementedException (FromHere(),"MeshUpgradeBuilder::getOutputPntsMappedCoords() for HEXA");
//       for (CFuint iKsi = 0; iKsi < nbrNodes1D; ++iKsi)
//       {
//         const CFreal ksi = -1.0 + iKsi*2.0/solOrder;
//         for (CFuint iEta = 0; iEta < nbrNodes1D; ++iEta)
//         {
//           const CFreal eta = -1.0 + iEta*2.0/solOrder;
//           for (CFuint iZta = 0; iZta < nbrNodes1D; ++iZta)
//           {
//             RealVector coords(3);
//             coords[KSI] = ksi;
//             coords[ETA] = eta;
//             coords[ZTA] = -1.0 + iZta*2.0/solOrder;
//             nodeMappedCoords.push_back(coords);
//           }
//         }
//       }
    } break;
    default:
    {
      throw Common::ShouldNotBeHereException (FromHere(),"Unknown GeoShape");
    }
  }

  return nodeMappedCoords;
}

//////////////////////////////////////////////////////////////////////////////

void MeshUpgradeBuilder::createBoundaryFacesTRS()
{
  CFAUTOTRACE;

  // get number of Topological Region Sets
  const CFuint nbTRSs = getCFmeshData().getNbTRSs();

  // get number of Topological Regions per TRS
  SafePtr< vector<CFuint> > nbTRs = getCFmeshData().getNbTRs();

  // get number of geometrical entities per TR
  SafePtr< vector< vector<CFuint> > > nbGeomEntsPerTR = getCFmeshData().getNbGeomEntsPerTR();
  
  // get TRS names
  SafePtr< vector<std::string> > nameTRS = getCFmeshData().getNameTRS();

  // get number of boundary faces + partition faces
  const CFuint nbBPlusPartitionFaces = m_bLocalGeoIDs.size();
  CFLog(VERBOSE,"Bnd Face: nb bnd+part faces: " << nbBPlusPartitionFaces << "\n");

  // flag telling if the face is a partition face
  m_isPartitionFace.resize(nbBPlusPartitionFaces);
  m_isPartitionFace = true;
  

  for(CFuint iTRS = 0; iTRS < nbTRSs; ++iTRS)
  {
    cf_assert((*nameTRS)[iTRS] != "InnerCells");
    cf_assert((*nameTRS)[iTRS] != "InnerFaces");

    CFLogDebugMed("Nb TRs in TRS: " << (*nbTRs)[iTRS] << "\n");
    
    for(CFuint iTR = 0; iTR < ((*nbGeomEntsPerTR)[iTRS]).size(); ++iTR)
    {
      (*nbGeomEntsPerTR)[iTRS][iTR] *= m_elementDivision;
    }

    Common::SafePtr<Framework::TopologicalRegionSet> ptrs =
    createTopologicalRegionSet((*nbGeomEntsPerTR)[iTRS],
                               (*nameTRS)[iTRS],
                               getCFmeshData().getTRGeoConn(iTRS),
                               iTRS);
    
    ptrs->attachTag("writable");

    CFLog(NOTICE, "Built TRS named " << (*nameTRS)[iTRS] << "\n");
    
    
  }
  // build partition boundary faces
  createPartitionBoundaryFacesTRS();

  CFLog(NOTICE, "Built PartitionFaces TRS \n");

}

//////////////////////////////////////////////////////////////////////////////

Common::SafePtr<Framework::TopologicalRegionSet>
MeshUpgradeBuilder::createTopologicalRegionSet
(const vector<CFuint>& nbFacesPerTR,
 const std::string& name,
 TRGeoConn& trGeoConn,
 const CFuint iTRS)
{
  CFAUTOTRACE;
  
  CFLog(VERBOSE,"MeshUpgradeBuilder::createTopologicalRegionSet\n");

  // get P1 number of face nodes
  SafePtr<vector<ElementTypeData> > elementType = getCFmeshData().getElementTypeData();
  cf_assert(elementType->size() > 0);
  const CFGeoShape::Type cellShape = (*elementType)[0].getGeoShape();
  vector< vector< CFuint > > bFaceOrientations = getBFaceOrientations(cellShape);
  const CFuint nbFaceNodesP1 = bFaceOrientations[0].size();

  // create the TopologicalRegion storage
  const CFuint nbTRs = nbFacesPerTR.size();
  vector<TopologicalRegion*>* storeTR = new vector<TopologicalRegion*>(nbTRs);

  // get total number of boundary faces + partition faces
  const CFuint nbBPlusPartitionFaces = m_bLocalGeoIDs.size();

  // compute the total number of faces in this TRS
  const CFuint totalNbFaces = accumulate(nbFacesPerTR.begin(), nbFacesPerTR.end(), static_cast<CFuint>(0));
        

  // fill in two arrays specifying the number of nodes and neighbour cells
  // per geometric entity in this TR
  std::valarray<CFuint> nbFaceNodes(totalNbFaces);
  std::valarray<CFuint> nbFaceCells(totalNbFaces);
  if (m_elementDivision == 1)
  {
    setNbNodesNbStatesInGeo(nbFacesPerTR, trGeoConn, nbFaceNodes, nbFaceCells);
  }
  else 
  {
    switch (cellShape)
  {
    case CFGeoShape::LINE:
    {
      throw Common::NotImplementedException (FromHere(),"MeshUpgradeBuilder::divideElements for LINE");
    } break;
    case CFGeoShape::TRIAG:
    {
      throw Common::NotImplementedException (FromHere(),"MeshUpgradeBuilder::divideElements for PYRAM");
    } break;
    case CFGeoShape::QUAD:
    {
      // for now, geoPolyOrder is assumed to be P1
      for (CFuint iFace = 0; iFace < totalNbFaces; ++iFace)
      {
        nbFaceNodes[iFace] = 2;
      }
    } break;
    case CFGeoShape::TETRA:
    {
      throw Common::NotImplementedException (FromHere(),"MeshUpgradeBuilder::divideElements for TETRA");
    } break;
    case CFGeoShape::PYRAM:
    {
      throw Common::NotImplementedException (FromHere(),"MeshUpgradeBuilder::divideElements for PYRAM");
    } break;
    case CFGeoShape::PRISM:
    {
      throw Common::NotImplementedException (FromHere(),"MeshUpgradeBuilder::divideElements for PRISM");
    } break;
    case CFGeoShape::HEXA:
    {
      throw Common::NotImplementedException (FromHere(),"MeshUpgradeBuilder::divideElements for HEXA");
    } break;
    default:
    {
      throw Common::ShouldNotBeHereException (FromHere(),"Unknown GeoShape");
    }
  }
  }

  // override the value of nbFaceCells set in setNbNodesNbStatesInGeo
  for (CFuint iFace = 0; iFace < totalNbFaces; ++iFace)
  {
    nbFaceCells[iFace] = 1;
  }

  // create face-node and face-cell connectivity tables
  ConnTable* faceNodes = new ConnTable(nbFaceNodes);
  ConnTable* faceCells = new ConnTable(nbFaceCells);

  // store the connectivities in the MeshData
  /// @todo think about having SharedPtr<ConnTable > and not put
  /// the TRS connectivities in MeshData
  MeshDataStack::getActive()->storeConnectivity(name + "Nodes", faceNodes);
  MeshDataStack::getActive()->storeConnectivity(name + "-Faces2Cells", faceCells);


  // array with all the IDs of all the geometric entities in this TRS
  // ownership of this array will be given to the TR
  vector<CFuint>* localFaceIDs  = new vector<CFuint>(totalNbFaces);
  vector<CFuint>* globalFaceIDs = new vector<CFuint>(totalNbFaces);
  vector<CFuint>* faceTypes = new vector<CFuint>(totalNbFaces);

  // name of the face provider
  const std::string faceProviderName = "Face";

  // get the global face IDs in this TRS
  SafePtr<vector<vector<vector<CFuint> > > > trsGlobalIDs = MeshDataStack::getActive()->getGlobalTRSGeoIDs();

  const bool hasGlobalIDs = (trsGlobalIDs->size() > 0) && (m_elementDivision == 1);

  CFuint nbProcessedFaces = 0;
  // let's create the required number of TR's
  for (CFuint iTR = 0; iTR < nbTRs; ++iTR)
  {
    CFLog(VERBOSE, "iTR= " << iTR << ", nbProcessedFaces= " << nbProcessedFaces << "\n");

    // number of faces in this TR
    const CFuint nbFaces = nbFacesPerTR[iTR]/m_elementDivision;

    // get the geoconnectivity (faces --> nodes and --> states (latter is not used))
    GeoConn& geoConn = trGeoConn[iTR];

    // if the cells are divided, resize the geoConn
    if (m_elementDivision != 1)
    {
      geoConn.resize(nbFacesPerTR[iTR]);
      ///@todo is this necessary?
//       for (CFuint iNewFace = nbFaces; iNewFace < nbFacesPerTR[iTR]; ++iNewFace)
//       {
// 	GeoConnElementPart temp1;
// 	GeoConnElementPart temp2;
// 	GeoConnElement temp(temp1,temp2);
// 	geoConn[iNewFace] = temp;
//       }
    }
    
    // allocate the TopologicalRegion
    (*storeTR)[iTR] = new TopologicalRegion();
    (*storeTR)[iTR]->setLocalNbGeoEnts(nbFacesPerTR[iTR]);
    cf_assert(nbFaces <= totalNbFaces);

    // check that the number of faces in this TR is > 0
    // if not just skip all this
    if (nbFaces > 0)
    {
      cf_assert((nbProcessedFaces < localFaceIDs->size()));

      (*storeTR)[iTR]->setGeoEntsLocalIdx(&(*localFaceIDs)[nbProcessedFaces], nbProcessedFaces);

      // set the connectivity tables of the TRS
      for (CFuint iFace = 0; iFace < nbFaces; ++iFace, ++nbProcessedFaces)
      {
        // face-node connectivity (not really needed during actual run)
        // number of nodes in this face
        const CFuint nbFaceNodes = faceNodes->nbCols(nbProcessedFaces*m_elementDivision);

        // boolean to check if face has been found
        bool faceFound = false;

        // loop over boundary faces
        CFuint faceIdx;
	
	// check if the cells need to be divided
	if (m_elementDivision == 1)
	{
          for (faceIdx = 0; faceIdx < nbBPlusPartitionFaces; ++faceIdx)
          {
            // if all the face nodes corresponding to faceIdx match (even if not in order),
            // then choose this faceIdx as the right one
            CFuint nbMatchingNodes = 0;
	  
            for (CFuint iNode = 0; iNode < nbFaceNodes; ++iNode)
            {
              const CFuint nodeID = (*m_bFaceNodes)(faceIdx, iNode);
              for (CFuint jNode = 0; jNode < nbFaceNodes; ++jNode)
              {
                if (geoConn[iFace].first[jNode] == nodeID)
                {
                  ++nbMatchingNodes;
                  break;
                } // end if
              } // end loop over boundary face nodes
            } // end loop over face nodes

            if (nbMatchingNodes == nbFaceNodes)
            {
              // the current faceIdx is the right one, this face is surely not on the partition boundary
              cf_assert(faceIdx < m_isPartitionFace.size());
              m_isPartitionFace[faceIdx] = false;
              faceFound = true;
              break;
            }
	  } // end loop over boundary faces
	  
	  // check if the face has been found
          if (!faceFound)
          {
            throw BadFormatException (FromHere(),"boundary face not found");
          }

          // store face-node connectivity
          for (CFuint iNode = 0; iNode < nbFaceNodes; ++iNode)
          {
            const CFuint nodeID =  (*m_bFaceNodes)(faceIdx, iNode);
	    
            (*faceNodes)(nbProcessedFaces, iNode) = nodeID;

            if (nodeID >= MeshDataStack::getActive()->getNbNodes() )
            {
              CFLogDebugMax( "NodeID: " << nodeID << "\n");
              throw BadFormatException (FromHere(),"CFmesh had bad node index in GeometricEntity");
            }
          }

          // face-cell connectivity
          cf_assert(faceCells->nbCols(nbProcessedFaces) == 1);
          (*faceCells)(nbProcessedFaces,0) = m_bFaceCell[faceIdx];

          // assign the local ID of the current face
          (*localFaceIDs)[nbProcessedFaces] = m_bLocalGeoIDs[faceIdx];

          // assign the global ID of the current geometric entity
          (*globalFaceIDs)[nbProcessedFaces] = (hasGlobalIDs) ? (*trsGlobalIDs)[iTRS][iTR][iFace] : iFace;

          // face geotype name
          /// @note getGeoShape() assumes a P1 entity. Here the number of nodes is passed for a P1 entity,
          /// even though the entity may be higher order. Another way to fix this is to extend getGeoShape()
          /// by passing also the dimensionality and the geometric order of the entity.
          const CFuint dim = PhysicalModelStack::getActive()->getDim();
          const std::string faceGeoTypeName = makeGeomEntName (getGeoShape(CFGeoEnt::FACE, dim, nbFaceNodesP1),
                                                               getGeometricPolyType(),
                                                               getGeometricPolyOrder(),
                                                               getSolutionPolyType(),
                                                               getSolutionPolyOrder());

          // face provider name
          const std::string geoProviderName = faceProviderName + faceGeoTypeName;
          (*faceTypes)[nbProcessedFaces] = m_mapGeoProviderNameToType.find(geoProviderName);
	}
	else 
	{
	  // new faces for a given old face
	  vector< CFuint > newFaces;
	  
	  // loop over faces to get the correct face Idx in this TR for the new faces 
	  for (faceIdx = 0; faceIdx < nbBPlusPartitionFaces; ++faceIdx)
          {
            // if all the face nodes corresponding to faceIdx match (even if not in order),
            // then choose this faceIdx as the right one
            CFuint nbMatchingNodes = 0;
	    // required amount of old nodes to be in the TRS for a new face to be in this TRS  
	    CFuint reqMatchingNodes = 0;
	    
	    // loop over the nodes of the face
            for (CFuint iNode = 0; iNode < nbFaceNodes; ++iNode)
            {
	      // localID of the new node
              const CFuint newNodeID = (*m_bFaceNodes)(faceIdx, iNode);
	      
	      // add to the required amount of nodes found to be a correct face
	      reqMatchingNodes += (m_newToOldNodeID[newNodeID]).size();
	      
	      cf_assert((m_newToOldNodeID[newNodeID]).size() > 0);
	      
              for (CFuint jNode = 0; jNode < nbFaceNodes; ++jNode)
              {
		for (CFuint iOldNode = 0; iOldNode < (m_newToOldNodeID[newNodeID]).size(); ++iOldNode)
		{
                  if (geoConn[iFace].first[jNode] == m_newToOldNodeID[newNodeID][iOldNode])
                  {
		    // matching old node found
                    ++nbMatchingNodes;
                  } // end if
		}
              } // end loop over boundary face nodes
            } // end loop over face nodes

            // if the nb of matching nodes is the required one, the faceIdx is a part of this TR
            if (nbMatchingNodes == reqMatchingNodes)
            {
              // the current faceIdx is the right one, this face is surely not on the partition boundary
              cf_assert(faceIdx < m_isPartitionFace.size());
              m_isPartitionFace[faceIdx] = false;
              newFaces.push_back(faceIdx);
            }
	  } // end loop over boundary faces
	  
	  // check if the faces have been found
          cf_assert(newFaces.size() == m_elementDivision);
	  
	  // create the entries of the new faces in the face node connectivity
	  for (CFuint iNewFace = 0; iNewFace < m_elementDivision; ++iNewFace)
	  {
            // store face-node connectivity
            for (CFuint iNode = 0; iNode < nbFaceNodes; ++iNode)
            {
	      // new node ID
              const CFuint nodeID =  (*m_bFaceNodes)(newFaces[iNewFace], iNode);
	      
	      // update the geoConn
	      geoConn[iFace+iNewFace*nbFaces].first[iNode] = nodeID;
	      
	      // update the face node conn
              (*faceNodes)(nbProcessedFaces*m_elementDivision+iNewFace, iNode) = nodeID;

              if (nodeID >= MeshDataStack::getActive()->getNbNodes() )
              {
                CFLogDebugMax( "NodeID: " << nodeID << "\n");
                throw BadFormatException (FromHere(),"CFmesh had bad node index in GeometricEntity");
              }
            }
            
            // face-cell connectivity
            cf_assert(faceCells->nbCols(nbProcessedFaces*m_elementDivision+iNewFace) == 1);
            (*faceCells)(nbProcessedFaces*m_elementDivision+iNewFace,0) = m_bFaceCell[newFaces[iNewFace]];

            // assign the local ID of the current face
            (*localFaceIDs)[nbProcessedFaces*m_elementDivision+iNewFace] = m_bLocalGeoIDs[newFaces[iNewFace]];

            // assign the global ID of the current geometric entity
            (*globalFaceIDs)[nbProcessedFaces*m_elementDivision+iNewFace] = m_bLocalGeoIDs[newFaces[iNewFace]];

            // face geotype name
            /// @note getGeoShape() assumes a P1 entity. Here the number of nodes is passed for a P1 entity,
            /// even though the entity may be higher order. Another way to fix this is to extend getGeoShape()
            /// by passing also the dimensionality and the geometric order of the entity.
            const CFuint dim = PhysicalModelStack::getActive()->getDim();
            const std::string faceGeoTypeName = makeGeomEntName (getGeoShape(CFGeoEnt::FACE, dim, nbFaceNodesP1),
                                                                 getGeometricPolyType(),
                                                                 getGeometricPolyOrder(),
                                                                 getSolutionPolyType(),
                                                                 getSolutionPolyOrder());

            // face provider name
            const std::string geoProviderName = faceProviderName + faceGeoTypeName;
            (*faceTypes)[nbProcessedFaces*m_elementDivision+iNewFace] = m_mapGeoProviderNameToType.find(geoProviderName);
	  }
        } 

      } // end for loop over faces

    } // end if statement

  } // end for loop over TRs

  // Create TopologicalRegionSet
  TopologicalRegionSet* ptrs = new TopologicalRegionSet (name, storeTR);

  // this is a boundary TRS
  ptrs->attachTag("boundary");
  // this is a TRS of faces
  ptrs->attachTag("face");

  // set local GeometricEntity IDs
  ptrs->setGeoEntsLocalIdx(localFaceIDs);
  // set global GeometricEntity IDs
  ptrs->setGeoEntsGlobalIdx(globalFaceIDs);
  // set the connectivity GeometricEntity to nodeIDs
  ptrs->setGeo2NodesConn(faceNodes);
  // set the GeometricEntity types
  ptrs->setGeoTypes(faceTypes);
  
  // don't set the connectivity GeometricEntity to stateIDs !!! there are no states in the boundary trs
  // set empty connectivity?
  // create empty face-state connectivity table
  std::valarray<CFuint> nbFaceStates(static_cast<CFuint>(0),totalNbFaces);
  ConnTable* faceStates = new ConnTable(nbFaceStates);
  ptrs->setGeo2StatesConn(faceStates);
  // create the cached list of state indexes in the all TRS
  ptrs->createStatesList();

  // put it into MeshData TRS list
  MeshDataStack::getActive()->addTrs(ptrs);

  cf_assert(nbProcessedFaces*m_elementDivision == totalNbFaces);

  CFLogDebugMin( "MeshUpgradeBuilder::setTopologicalRegions() : nameTRS "
                 << name
                 << ", nb boundary faces detected : "
                 << nbProcessedFaces << "\n");
  return ptrs;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluxReconstructionMethod

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

