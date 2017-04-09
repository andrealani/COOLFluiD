#include "Common/CFLog.hh"
#include "Environment/ObjectProvider.hh"

#include "Common/BadValueException.hh"
#include "Framework/MeshData.hh"
#include "Framework/PhysicalModel.hh"

#include "FluxReconstructionMethod/FluxReconstruction.hh"
#include "FluxReconstructionMethod/MeshUpgradeBuilder.hh"

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
  options.addConfigOption< CFuint >("DivideElements","Divide elements on equal parts to form new cells. This number is equal to te number of element adjescend to an old element face");
}

//////////////////////////////////////////////////////////////////////////////

MeshUpgradeBuilder::MeshUpgradeBuilder(const std::string& name) :
  FluxReconstructionBuilder(name),
  m_solPolyOrder(),
  m_geoPolyOrder(),
  m_elementDivision(),
  m_prevGeoPolyOrder(),
  m_bndFacesNodes(),
  m_globalIDs(),
  m_elemLocalIDOfState(),
  m_elemIDOfState(),
  m_elemFirstStateLocalID(),
  m_updatables()
{
  addConfigOptionsTo(this);

  m_solPolyOrderStr = "P1";
  setParameter( "PolynomialOrder", &m_solPolyOrderStr);

  m_geoPolyOrderStr = "P1";
  setParameter( "GeoPolynomialOrder", &m_geoPolyOrderStr);
  
  m_elementDivision = 1;
  setParameter( "DivideElements", &m_elementDivision);
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
        case CFPolyOrder::ORDER2:
        {
          return 0;
        }
        default:
        {
          throw Common::NotImplementedException (FromHere(),"Face type not implemented...");
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
  SafePtr< vector<ElementTypeData> > elementType =
    getCFmeshData().getElementTypeData();
    
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
  const CFuint nbNewCellsPerOldCell = pow(m_elementDivision,getCFmeshData().getDimension());
  
  const CFuint newNbElements = nbElements*nbNewCellsPerOldCell;
  
  // temporary storage for original (lower order) nodes:
  const MeshData::ConnTable initialCellNodes = (*cellNodes);
  const MeshData::ConnTable initialCellStates = (*cellStates);
  
  vector< vector< RealVector > > newNodeCoords;
  vector< bool > updatables;
  
  m_pattern.resize(newNbElements);
  
  vector< ElementTypeData >::iterator type_itr;
  
  // set the correct number of nodes per element in m_pattern
  CFuint firstFreeIdx = 0;
  CFuint nodeToCellsMapSize = 0;
  CFuint nbNewNodes = 0;
  // nodes to cells multimap
  CFMultiMap<CFuint,CFuint> mapNode2CellsOld;
  mapNode2CellsOld.reserve(nbNodes);
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
      
      nodeToCellsMapSize += nbNodesPerCell*nbNewCellsPerOldCell;
      
      newNodeCoords.push_back(getNewNodesMappedCoords(elemGeoShape,m_elementDivision,globalIdx,initialCellNodes));
      
      updatables.push_back(states[(*cellStates)(globalIdx,0)]->isParUpdatable());
      
      for (CFuint iNode = 0; iNode < nbNodesPerCell; ++iNode)
      {
        mapNode2CellsOld.insert((*cellNodes)(globalIdx,iNode),globalIdx);
      }

      switch(elemGeoShape) 
      {

      case CFGeoShape::TRIAG:
	for (CFuint iNewCell = 0; iNewCell < nbNewCellsPerOldCell; ++iNewCell,++firstFreeIdx)
	{
	  m_pattern[firstFreeIdx] = 3;
	  nbNewNodes += 3;
	}
        break;

      case CFGeoShape::QUAD:
        for (CFuint iNewCell = 0; iNewCell < nbNewCellsPerOldCell; ++iNewCell,++firstFreeIdx)
	{
	  m_pattern[firstFreeIdx] = 4;
	  nbNewNodes += 4;
	}
        break;
	
      case CFGeoShape::HEXA:
        for (CFuint iNewCell = 0; iNewCell < nbNewCellsPerOldCell; ++iNewCell,++firstFreeIdx)
	{
	  m_pattern[firstFreeIdx] = 8;
	  nbNewNodes += 8;
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
  
  cellNodes->resize(m_pattern);
  cellStates->resize(m_pattern);
  
  // resize the states data handle to put dummy states in it
  states.resize(nbNewNodes);
  
  nodes.resize(nbNewNodes);
  
  firstFreeIdx = 0;
  CFuint nodeIdx = 0;
  const CFuint nodes1D = m_elementDivision + 1;
  
  // resize the states data handle to put dummy states in it
  states.resize(nbNewNodes);
  
  // variable marking whether a cell has been updated
  // the element type is also stored here, though this is not really necessary since there is only one type of element
  vector<CFuint> isUpdated(nbElements);

  // nodes to cells multimap
  CFMultiMap<CFuint,CFuint> mapNode2Cells;
  typedef CFMultiMap<CFuint, CFuint>::MapIterator MapIterator;
  mapNode2Cells.reserve(nodeToCellsMapSize);
  
  const CFuint nbeq = PhysicalModelStack::getActive()->getNbEq();
  
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
      
      const CFuint idxFirstNewElement = firstFreeIdx;
      
      
      vector< pair<MapIterator, MapIterator> > neighbors;

      switch(elemGeoShape) 
      {
      case CFGeoShape::TRIAG:

        break;

      case CFGeoShape::QUAD:
	
	
	
	for (CFuint iNode = 0; iNode < 4; ++iNode)
	{
	  const CFuint id = initialCellNodes(globalIdx,iNode);

	  bool fo = false;
	  neighbors.push_back(mapNode2Cells.find(id,fo));
	  cf_assert(fo);
	}
	
	
	for (CFuint iKsi = 0; iKsi < nodes1D; ++iKsi)
	{
	  for (CFuint iEta = 0; iEta < nodes1D; ++iEta)
	  {
	    if (iKsi != 0 && iKsi != nodes1D-1 && iEta != 0 && iEta != nodes1D-1)
	    {
	      // insert in nodes to cells map
              mapNode2Cells.insert(nodeIdx, idxFirstNewElement+iEta-1+(iKsi-1)*m_elementDivision);
	      mapNode2Cells.insert(nodeIdx, idxFirstNewElement+iEta+(iKsi-1)*m_elementDivision);
	      mapNode2Cells.insert(nodeIdx, idxFirstNewElement+iEta-1+(iKsi)*m_elementDivision);
	      mapNode2Cells.insert(nodeIdx, idxFirstNewElement+iEta+(iKsi)*m_elementDivision);
	      ++nodeIdx;
	    }
	  }
	}

        for (CFuint iKsi = 0; iKsi < m_elementDivision; ++iKsi)
	{
	  for (CFuint iEta = 0; iEta < m_elementDivision; ++iEta, ++firstFreeIdx)
	  {
	    for (CFuint iNode = 0; iNode < nbNodesPerCell; ++iNode)
	    {
	      (*cellNodes)(firstFreeIdx,iNode) = nodeIdx;
	      // create the node in the mesh data
	      if (iNode == 0)
	      {
                getCFmeshData().createNode(nodeIdx,newNodeCoords[globalIdx][iEta+nodes1D*iKsi]);
	      }
	      else if (iNode == 1)
	      {
                getCFmeshData().createNode(nodeIdx,newNodeCoords[globalIdx][iEta+nodes1D*iKsi+1]);
	      }
	      else if (iNode == 2)
	      {
                getCFmeshData().createNode(nodeIdx,newNodeCoords[globalIdx][iEta+nodes1D*(iKsi+1)+1]);
	      }
	      else
	      {
                getCFmeshData().createNode(nodeIdx,newNodeCoords[globalIdx][iEta+nodes1D*(iKsi+1)]);
	      }
              nodes[nodeIdx]->setGlobalID(nodeIdx);
              nodes[nodeIdx]->setIsOnMesh(true);
              nodes[nodeIdx]->setIsOwnedByState(false);
	      
	      (*cellStates)(firstFreeIdx,iNode) = nodeIdx;
	      RealVector stateData (nbeq);
              getCFmeshData().createState(nodeIdx,stateData);
              states[nodeIdx]->setParUpdatable(updatables[globalIdx]);
              states[nodeIdx]->setGlobalID(nodeIdx);
	    }
	  }
	}
        break;
	
      case CFGeoShape::HEXA:

        break;

      default:
        std::string shape =
          CFGeoShape::Convert::to_str(elemGeoShape);
        std::string msg = std::string("Element type not implemented: ") + shape;
        throw Common::NotImplementedException (FromHere(),msg);
      }
    }
    type_itr->setNbElems(nbNewCellsPerOldCell*nbrCellsPerType);
  }
  
  getCFmeshData().setNbElements(newNbElements);
  
}

//////////////////////////////////////////////////////////////////////////////

void MeshUpgradeBuilder::upgradeStateConnectivity()
{
  CFAUTOTRACE;

  // get the element type data
  SafePtr< vector< ElementTypeData > > elementType = getCFmeshData().getElementTypeData();
  vector< ElementTypeData >::iterator type_itr;

  // get the connectivity that we will change
  SafePtr<MeshData::ConnTable> cellStates = MeshDataStack::getActive()->getConnectivity("cellStates_InnerCells");
  
  // get the old states
  DataHandle < Framework::State*, Framework::GLOBAL > oldStates = getCFmeshData().getStatesHandle();

  const CFuint nbElems = getNbElements();
  
  const CFuint oldNbStatesPerElem = (oldStates.size())/nbElems;
  
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

    // loop over elements of this type
    for (CFuint iElem = 0; iElem < nbrElems; ++iElem, ++globalIdx)
    {
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
}

//////////////////////////////////////////////////////////////////////////////

void MeshUpgradeBuilder::recreateStates()
{
  CFAUTOTRACE;

  SafePtr<MeshData::ConnTable> cellStates = MeshDataStack::getActive()->getConnectivity("cellStates_InnerCells");

  //DataHandle < Framework::State*, Framework::GLOBAL > states = getCFmeshData().getStatesHandle();
  DataHandle < Framework::State*, Framework::GLOBAL > states =
    MeshDataStack::getActive()->getStateDataSocketSink().getDataHandle();

  const CFuint newNbStates = cellStates->size();
  const CFuint oldNbStates = states.size();
  const CFuint oldNbGlobalStates = states.getGlobalSize();

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
  RealVector stateData (nbeq);
  
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

    getCFmeshData().createState(localID,stateData);
    states[localID]->setParUpdatable(m_updatables[iState]);
    states[localID]->setGlobalID(globalID);
    
  }

  // AL: .resize(0) crashes with some compilers on some systems, better to use clear()
  m_updatables.clear();
  m_globalIDs.clear();
  m_elemLocalIDOfState.clear();
  m_elemIDOfState.clear();
  m_elemFirstStateLocalID.clear();
    
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

vector< RealVector > MeshUpgradeBuilder::getNewNodesMappedCoords(CFGeoShape::Type shape, CFuint solOrder, CFuint cellIdx, const MeshData::ConnTable nodesConn)
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

  switch (shape)
  {
    case CFGeoShape::LINE:
    {
      throw Common::NotImplementedException (FromHere(),"ParaWriterData::getOutputPntsMappedCoords() for LINE");
//       for (CFuint iKsi = 0; iKsi < nbrNodes1D; ++iKsi)
//       {
//         RealVector coords(1);
//         coords[KSI] = -1.0 + iKsi*2.0/solOrder;
//         nodeMappedCoords.push_back(coords);
//       }
    } break;
    case CFGeoShape::TRIAG:
    {
      throw Common::NotImplementedException (FromHere(),"ParaWriterData::getOutputPntsMappedCoords() for PYRAM");
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
      vector< Node* > nodeCoords;
      for (CFuint iNode = 0; iNode < 4; ++iNode)
      {
        nodeIDs[iNode] = nodesConn(cellIdx,iNode);
	nodeCoords.push_back(nodes[nodeIDs[iNode]]);
      }

      for (CFuint iKsi = 0; iKsi < nbrNodes1D; ++iKsi)
      {
        const CFreal xDown = (*nodeCoords[0])[0]+iKsi/solOrder*((*nodeCoords[1])[0]-(*nodeCoords[0])[0]);
	const CFreal xUp = (*nodeCoords[3])[0]+iKsi/solOrder*((*nodeCoords[2])[0]-(*nodeCoords[3])[0]);
	const CFreal yDown = (*nodeCoords[0])[1]+iKsi/solOrder*((*nodeCoords[1])[1]-(*nodeCoords[0])[1]);
	const CFreal yUp = (*nodeCoords[3])[1]+iKsi/solOrder*((*nodeCoords[2])[1]-(*nodeCoords[3])[1]);
        for (CFuint iEta = 0; iEta < nbrNodes1D; ++iEta)
        {
	  const CFreal xLeft = (*nodeCoords[0])[0]+iEta/solOrder*((*nodeCoords[3])[0]-(*nodeCoords[0])[0]);
	  const CFreal xRight = (*nodeCoords[1])[0]+iEta/solOrder*((*nodeCoords[2])[0]-(*nodeCoords[1])[0]);
          const CFreal yLeft = (*nodeCoords[0])[1]+iEta/solOrder*((*nodeCoords[3])[1]-(*nodeCoords[0])[1]);
	  const CFreal yRight = (*nodeCoords[1])[1]+iEta/solOrder*((*nodeCoords[2])[1]-(*nodeCoords[1])[1]);
	  const CFreal mLR = (yRight-yLeft)/(xRight-xLeft);
	  const CFreal mUD = (yUp-yDown)/(xUp-xDown);
          const CFreal x = (mLR*xLeft-mUD*xDown+yDown-yLeft)/(mLR-mUD);
	  const CFreal y = mLR*x+yLeft-xLeft*mLR;
	  RealVector currCoords(2);
	  currCoords[0] = x;
	  currCoords[1] = y;
	  nodeMappedCoords.push_back(currCoords);
        }
      }
    } break;
    case CFGeoShape::TETRA:
    {
      throw Common::NotImplementedException (FromHere(),"ParaWriterData::getOutputPntsMappedCoords() for TETRA");
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
      throw Common::NotImplementedException (FromHere(),"ParaWriterData::getOutputPntsMappedCoords() for PYRAM");
    } break;
    case CFGeoShape::PRISM:
    {
      throw Common::NotImplementedException (FromHere(),"ParaWriterData::getOutputPntsMappedCoords() for PRISM");
    } break;
    case CFGeoShape::HEXA:
    {
      throw Common::NotImplementedException (FromHere(),"ParaWriterData::getOutputPntsMappedCoords() for HEXA");
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

    } // namespace FluxReconstructionMethod

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

