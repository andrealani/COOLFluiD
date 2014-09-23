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
#include "FluctSplit/FluctSplit.hh"
#include "FluctSplit/RDS_HighOrderMeshUpdater.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



    namespace FluctSplit{

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<RDS_HighOrderMeshUpdater,
                            MeshDataBuilder,
                            FluctSplitModule,
                            1>
RDS_HighOrderMeshUpdaterProvider("FluctSplitHO");

//////////////////////////////////////////////////////////////////////////////

void RDS_HighOrderMeshUpdater::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOption< CFuint >("GeoOrder","Order of geometry of the new mesh");
   options.addConfigOption< CFuint >("SolOrder","Order of representation of numerical solution");
   options.addConfigOption< bool >("UpdateGeometry","Flag to update nodes of the mesh");
   options.addConfigOption< bool >("UpdateSolution","Flag to update states of the mesh");
}
//////////////////////////////////////////////////////////////////////////////

RDS_HighOrderMeshUpdater::RDS_HighOrderMeshUpdater(const std::string& name) :
  RDS_MeshDataBuilder(name),
  socket_nodes("nodes")
{
  addConfigOptionsTo(this);

  m_newGeoPolyOrder_int = 1;
  setParameter("GeoOrder",&m_newGeoPolyOrder_int);

//   m_newGeoPolyOrder = static_cast<CFPolyOrder::Type>(m_newGeoPolyOrder_int);

  m_newSolPolyOrder_int = 2;
  setParameter("SolOrder",&m_newSolPolyOrder_int);

//   m_newSolPolyOrder = static_cast<CFPolyOrder::Type>(m_newSolPolyOrder_int);

  m_updateGeometry = true;
  setParameter("UpdateGeometry",&m_updateGeometry);

  m_updateSolution = true;
  setParameter("UpdateSolution",&m_updateSolution);


}

//////////////////////////////////////////////////////////////////////////////

RDS_HighOrderMeshUpdater::~RDS_HighOrderMeshUpdater()
{
}

//////////////////////////////////////////////////////////////////////////////

void RDS_HighOrderMeshUpdater::configure ( Config::ConfigArgs& args )
{

  RDS_MeshDataBuilder::configure(args);

  m_newGeoPolyOrder = static_cast<CFPolyOrder::Type>(m_newGeoPolyOrder_int);
  m_newSolPolyOrder = static_cast<CFPolyOrder::Type>(m_newSolPolyOrder_int);

//   CF_DEBUG_OBJ(m_newGeoPolyOrder);
//   CF_DEBUG_OBJ(m_newGeoPolyOrder_int);
//   CF_DEBUG_OBJ(m_newSolPolyOrder);
//   CF_DEBUG_OBJ(m_newSolPolyOrder_int);


}
//////////////////////////////////////////////////////////////////////////////


void RDS_HighOrderMeshUpdater::releaseMemory()
{
  RDS_MeshDataBuilder::releaseMemory();
}

//////////////////////////////////////////////////////////////////////////////

void RDS_HighOrderMeshUpdater::computeGeoTypeInfo()
{
  CFAUTOTRACE;

//   CFout << "\n-----------------------------------------------------\n\n";
//   CFout << "\tExecuting computeGeoTypeInfo()\n";
//   CFout << "\n-----------------------------------------------------\n";


  //cf_assert(m_newGeoPolyOrder.isNotNull());

  // get the element type data
  SafePtr< std::vector< ElementTypeData > > elementType = getCFmeshData().getElementTypeData();

  std::vector< ElementTypeData >::iterator type_itr = elementType->begin();
  for (; type_itr != elementType->end(); ++type_itr)
  {

//     CFout << "m_newGeoPolyOrder = " << m_newGeoPolyOrder << "\n";
//     CFout << "m_newSolPolyOrder = " << m_newSolPolyOrder << "\n";

    // get the number of nodes in this element type: first number of nodes that lie on the faces of the elements,
    // then those nodes that are inside the element
    const CFuint nbNodesPerElem = getNbDofInNewType(type_itr->getGeoShape(),m_newGeoPolyOrder);
    const CFuint nbInteriorNodesPerElem = getNbInteriorDofInNewType(type_itr->getGeoShape(),m_newGeoPolyOrder);


    // get the number of control volumes (states) in this element type
    const CFuint nbStatesPerElem = getNbDofInNewType(type_itr->getGeoShape(),m_newSolPolyOrder);
    const CFuint nbInteriorStatesPerElem = getNbInteriorDofInNewType(type_itr->getGeoShape(),m_newSolPolyOrder);

    type_itr->setNbNodes(nbNodesPerElem+nbInteriorNodesPerElem);
    type_itr->setGeoOrder(m_newGeoPolyOrder_int);

    type_itr->setNbStates(nbStatesPerElem+nbInteriorStatesPerElem);
    type_itr->setSolOrder(m_newSolPolyOrder_int);



    CFout << "GEOMETRY:\n";
//     CFout << "\tNumber of nodes per element = " << type_itr->getNbNodes() << "\n";
    CFout << "\tNumber of nodes per element = " << nbNodesPerElem << " + " << nbInteriorNodesPerElem << "\n";
    CFout << "\tOrder of polynomial: " << type_itr->getGeoOrder() << "\n\n";
    CFout << "SOLUTION:\n";
    CFout << "\tNumber of states per element = " << nbStatesPerElem << " + " << nbInteriorStatesPerElem << "\n";
    CFout << "\tOrder of polynomial: " << type_itr->getSolOrder() << "\n\n";

  }


//   CFout << "Geometric polyOrder = " << getGeometricPolyOrder() << "\n";
  getCFmeshData().setGeometricPolyOrder(m_newGeoPolyOrder);
//   CFout << "Geometric polyOrder = " << getGeometricPolyOrder() << "\n";


  /// @todo since only meshes with the same order are currently supported
  ///       we have to change the order in the CFmeshData and
  ///       not only on the ElementTypeData
  getCFmeshData().setSolutionPolyOrder(m_newSolPolyOrder);

//  CFout << "------------------------------------------------------------\n";
//  CFout << "Printing LocalConnectionData before adding new type:\n";
//  CFout << "------------------------------------------------------------\n";
//  LocalConnectionData::getInstance().print();


//  LocalConnectionData::getInstance().addType(TRIAG,m_newGeoPolyOrder,CFPolyForm::LAGRANGE,m_newSolPolyOrder,CFPolyForm::LAGRANGE,NODE);
//  LocalConnectionData::getInstance().addType(TRIAG,CFPolyOrder::ORDER1,CFPolyForm::LAGRANGE,CFPolyOrder::ORDER1,CFPolyForm::LAGRANGE,NODE);
//  LocalConnectionData::getInstance().addType(TRIAG,getGeometricPolyOrder(),CFPolyForm::LAGRANGE,getSolutionPolyOrder(),CFPolyForm::LAGRANGE,NODE);




  // continue with the standard algorithm
  RDS_MeshDataBuilder::computeGeoTypeInfo();



//   CFout << "------------------------------------------------------------\n";
//  CFout << "Printing LocalConnectionData after adding new type:\n";
//  CFout << "------------------------------------------------------------\n";
//
//  LocalConnectionData::getInstance().print();


}


//////////////////////////////////////////////////////////////////////////////


CFuint RDS_HighOrderMeshUpdater::getNbDofInNewType(const CFGeoShape::Type geoShape, const CFPolyOrder::Type polyOrder)
{

//   CFout << "\n-----------------------------------------------------\n\n";
//   CFout << "\tExecuting getNbDofInNewType()\n";
//   CFout << "\n-----------------------------------------------------\n";



  //cf_assert(polyOrder == 2);

  switch(polyOrder)
  {
  case CFPolyOrder::ORDER1:
      switch (geoShape)
     {
      case CFGeoShape::LINE:
        return 2;
      case CFGeoShape::TRIAG:
        return 3;
      case CFGeoShape::QUAD:
        return 4;
      case CFGeoShape::TETRA:
        return 4;
      case CFGeoShape::PRISM:
        throw Common::NotImplementedException(FromHere(),"Element type PRISM not implemented...");
      case CFGeoShape::PYRAM:
        throw Common::NotImplementedException(FromHere(),"Element type PYRAM not implemented...");
      case CFGeoShape::HEXA:
        return 8;
      default:
        throw Common::NotImplementedException(FromHere(),"Element type not found...");
    }
  break;
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
        throw Common::NotImplementedException(FromHere(),"Element type PRISM not implemented...");
      case CFGeoShape::PYRAM:
        throw Common::NotImplementedException(FromHere(),"Element type PYRAM not implemented...");
      case CFGeoShape::HEXA:
        return 20;
      default:
        throw Common::NotImplementedException(FromHere(),"Element type not found...");
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
        throw Common::NotImplementedException(FromHere(),"Element type QUAD not implemented...");
      case CFGeoShape::TETRA:
        throw Common::NotImplementedException(FromHere(),"Element type TETRA not implemented...");
      case CFGeoShape::PRISM:
        throw Common::NotImplementedException(FromHere(),"Element type PRISM not implemented...");
      case CFGeoShape::PYRAM:
        throw Common::NotImplementedException(FromHere(),"Element type PYRAM not implemented...");
      case CFGeoShape::HEXA:
        throw Common::NotImplementedException(FromHere(),"Element type HEXA not implemented...");
      default:
        throw Common::NotImplementedException(FromHere(),"Element type not found...");
    }
  break;

//Default: order = 2
  default:
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
        throw Common::NotImplementedException(FromHere(),"Element type PRISM not implemented...");
      case CFGeoShape::PYRAM:
        throw Common::NotImplementedException(FromHere(),"Element type PYRAM not implemented...");
      case CFGeoShape::HEXA:
        return 20;
      default:
        throw Common::NotImplementedException(FromHere(),"Element type not found...");
    }


  }
}



//////////////////////////////////////////////////////////////////////////////

CFuint RDS_HighOrderMeshUpdater::getNbInteriorDofInNewType(const CFGeoShape::Type geoShape, const CFPolyOrder::Type polyOrder)
{

//   CFout << "\n-----------------------------------------------------\n\n";
//   CFout << "\tExecuting getNbInteriorDofInNewType()\n";
//   CFout << "\n-----------------------------------------------------\n";



switch(polyOrder)
  {
  case CFPolyOrder::ORDER1:
  return 0;
  break;

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
      throw Common::NotImplementedException(FromHere(),"Element type PRISM not implemented...");
    case CFGeoShape::PYRAM:
      throw Common::NotImplementedException(FromHere(),"Element type PYRAM not implemented...");
    case CFGeoShape::HEXA:
      return 0;
    default:
      throw Common::NotImplementedException(FromHere(),"Element type not found...");
  }
  break;
  case CFPolyOrder::ORDER3:
  switch (geoShape)
  {
    case CFGeoShape::LINE:
      return 2;
    case CFGeoShape::TRIAG:
      return 1;
     case CFGeoShape::QUAD:
        throw Common::NotImplementedException(FromHere(),"Element type QUAD not implemented...");
      case CFGeoShape::TETRA:
        throw Common::NotImplementedException(FromHere(),"Element type TETRA not implemented...");
      case CFGeoShape::PRISM:
        throw Common::NotImplementedException(FromHere(),"Element type PRISM not implemented...");
      case CFGeoShape::PYRAM:
        throw Common::NotImplementedException(FromHere(),"Element type PYRAM not implemented...");
      case CFGeoShape::HEXA:
        throw Common::NotImplementedException(FromHere(),"Element type HEXA not implemented...");
      default:
        throw Common::NotImplementedException(FromHere(),"Element type not found...");
  }
  break;

  //default: order = 2

  default:
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
      throw Common::NotImplementedException(FromHere(),"Element type PRISM not implemented...");
    case CFGeoShape::PYRAM:
      throw Common::NotImplementedException(FromHere(),"Element type PYRAM not implemented...");
    case CFGeoShape::HEXA:
      return 0;
    default:
      throw Common::NotImplementedException(FromHere(),"Element type not found...");
  }


  }
}

//////////////////////////////////////////////////////////////////////////////

CFuint RDS_HighOrderMeshUpdater::getNewNbDofInBoundaryGeo(const CFuint oldNbDof, const CFPolyOrder::Type polyOrder)
{

//   CFout << "\n-----------------------------------------------------\n\n";
//   CFout << "\tExecuting getNbStatesInBoundaryGeo()\n";
//   CFout << "\n-----------------------------------------------------\n";



  const CFuint dim = getCFmeshData().getDimension();

  switch(polyOrder)
{
  case CFPolyOrder::ORDER2:
  if(dim == DIM_2D)
    {
      cf_assert(oldNbDof == 2);
      return 3;
    }
    else
    {
      cf_assert(dim == DIM_3D);
      switch (oldNbDof)
     {
        case 3: // we have a triangular face
          return 6;
        case 4: // we have a quad face
          return 8;
        default:
          throw Common::NotImplementedException(FromHere(),"Face type not found...");
      }
    }
  break;
  case CFPolyOrder::ORDER3:
  if(dim == DIM_2D)
    {
      cf_assert(oldNbDof == 2);
      return 4;
    }
    else
    {
      cf_assert(dim == DIM_3D);
      switch (oldNbDof)
     {
        case 3: // we have a triangular face
          return 9;
        case 4: // we have a quad face
         throw Common::NotImplementedException(FromHere(),"Element type PRISM not implemented...");
        default:
          throw Common::NotImplementedException(FromHere(),"Face type not found...");
      }
    }
  break;


  //default: CFPolyOrder::ORDER2
  default:
  if(dim == DIM_2D)
    {
      cf_assert(oldNbDof == 2);
      return 3;
    }
    else
    {
      cf_assert(dim == DIM_3D);
      switch (oldNbDof)
     {
        case 3: // we have a triangular face
          return 6;
        case 4: // we have a quad face
          return 8;
        default:
          throw Common::NotImplementedException(FromHere(),"Face type not found...");
      }
    }



}
}

//////////////////////////////////////////////////////////////////////////////

void RDS_HighOrderMeshUpdater::createTopologicalRegionSets()
{
  CFAUTOTRACE;

//   CFout << "\n-----------------------------------------------------\n\n";
//   CFout << "\tExecuting createTopologicalRegionSets()\n";
//   CFout << "\n-----------------------------------------------------\n";


  CFout << "Converting the P1P1 mesh to a PmPn mesh\n\n";

  SafePtr< std::vector< ElementTypeData > > elementType = getCFmeshData().getElementTypeData();
  const CFuint nbElemTypes = elementType->size();

  ///Store the element local connectivity in an easy way

  /// Local connectivity face-node for each element type, P1 geometry and solution
  /// This is used as reference when adding new nodes and states to faces
  m_faceNodeP1Element.resize(nbElemTypes);
  m_faceStateP1Element.resize(nbElemTypes);

  for (CFuint iType = 0; iType < nbElemTypes; ++iType) {
    m_faceNodeP1Element[iType] = LocalConnectionData::getInstance().getFaceDofLocal
      ((*elementType)[iType].getGeoShape(),
       CFPolyOrder::ORDER1,
       NODE,
       CFPolyForm::LAGRANGE);
    m_faceStateP1Element[iType] = LocalConnectionData::getInstance().getFaceDofLocal
      ((*elementType)[iType].getGeoShape(),
       CFPolyOrder::ORDER1,
       STATE,
       CFPolyForm::LAGRANGE);
  }


m_faceNodePnElement.resize(nbElemTypes);

switch(m_newGeoPolyOrder) {

  case CFPolyOrder::ORDER2: {
    CFout << "Creating local P2 connectivity table for nodes\n\n";
    /// Local connectivity face-node for each element type, P2 geometry
      for (CFuint iType = 0; iType < nbElemTypes; ++iType) {
        m_faceNodePnElement[iType] = LocalConnectionData::getInstance().getFaceDofLocal
          ((*elementType)[iType].getGeoShape(),
          CFPolyOrder::ORDER2,
          NODE,
          CFPolyForm::LAGRANGE);
      }

    }
    break;

  case CFPolyOrder::ORDER3: {

    CFout << "Creating local P3 connectivity table for nodes\n\n";
    /// Local connectivity face-node for each element type, P2 geometry
      for (CFuint iType = 0; iType < nbElemTypes; ++iType) {
        m_faceNodePnElement[iType] = LocalConnectionData::getInstance().getFaceDofLocal
          ((*elementType)[iType].getGeoShape(),
          CFPolyOrder::ORDER3,
          NODE,
          CFPolyForm::LAGRANGE);
      }

//     for(CFuint iface = 0; iface < 6; ++iface) {
//     CFout << "[";
//       for(CFuint inode = 0; inode < 3; ++inode) {
//         CFout << (*m_faceNodePnElement[0])(iface,inode) << ", ";
//       }
//         CFout << (*m_faceNodePnElement[0])(iface,3) << "]\n";
//     }

    } //order 3

    break;


  default:
  {
    CFout << "Creating default local higher order connectivity table for nodes: P1\n\n";
    /// local connectivity face-node for each element type, default: P1 geometry
      for (CFuint iType = 0; iType < nbElemTypes; ++iType) {
        m_faceNodePnElement[iType] = LocalConnectionData::getInstance().getFaceDofLocal
          ((*elementType)[iType].getGeoShape(),
          CFPolyOrder::ORDER1,
          NODE,
          CFPolyForm::LAGRANGE);
      }
   }

} //switch


m_faceStatePnElement.resize(nbElemTypes);

switch(m_newSolPolyOrder) {

  case CFPolyOrder::ORDER2: {
    CFout << "Creating local P2 connectivity table for states\n\n";
    /// Local connectivity face-state for each element type, P2 (solution)
      for (CFuint iType = 0; iType < nbElemTypes; ++iType) {
        m_faceStatePnElement[iType] = LocalConnectionData::getInstance().getFaceDofLocal
          ((*elementType)[iType].getGeoShape(),
          CFPolyOrder::ORDER2,
          STATE,
          CFPolyForm::LAGRANGE);
      }
    }
    break;

  case CFPolyOrder::ORDER3: {

    CFout << "Creating local P3 connectivity table for states\n\n";
    /// Local connectivity face-node for each element type, P2 geometry
      for (CFuint iType = 0; iType < nbElemTypes; ++iType) {
        m_faceStatePnElement[iType] = LocalConnectionData::getInstance().getFaceDofLocal
          ((*elementType)[iType].getGeoShape(),
          CFPolyOrder::ORDER3,
          STATE,
          CFPolyForm::LAGRANGE);
      }
    }
    break;

  default:
  {
    CFout << "Creating default local higher order connectivity table for states: P1\n\n";
    /// local connectivity face-node for each element type, default: P1 solution
      for (CFuint iType = 0; iType < nbElemTypes; ++iType) {
        m_faceStatePnElement[iType] = LocalConnectionData::getInstance().getFaceDofLocal
          ((*elementType)[iType].getGeoShape(),
          CFPolyOrder::ORDER1,
          STATE,
          CFPolyForm::LAGRANGE);
      }
   }

} //switch



  //transform the cell-nodes connectivity
  // from P1 geometry to P2 or P3 geometry
  if((m_newGeoPolyOrder > 1) && m_updateGeometry) upgradeNodalConnectivity();

  // first transform the cell-states connectivity
  // from P1 to P2 or P3
  if((m_newSolPolyOrder > 1) && m_updateSolution) upgradeStateConnectivity();

  //recreate the nodes
  if ((m_newGeoPolyOrder > 1) && m_updateGeometry) recreateNodes();

  // recreate the states
  if ((m_newSolPolyOrder > 1) && m_updateSolution) recreateStates();

  // modify the TRS states
  updateTRSData();

  // continue with the standard algorithm
  RDS_MeshDataBuilder::createTopologicalRegionSets();


}


//////////////////////////////////////////////////////////////////////////////


void RDS_HighOrderMeshUpdater::upgradeNodalConnectivity()
{
 CFAUTOTRACE;

//  CFout << "\n-----------------------------------------------------\n\n";
//  CFout << "\tExecuting upgradeNodalConnectivity()\n";
//  CFout << "\n-----------------------------------------------------\n";

  CFout << "Updating existing nodes connectivity\n";
  // get the element type data
  SafePtr< std::vector< ElementTypeData > > elementType = getCFmeshData().getElementTypeData();
  std::vector< ElementTypeData >::iterator type_itr;

  // get the nodal connectivity
  Common::SafePtr<MeshData::ConnTable> cellNodes = MeshDataStack::getActive()->getConnectivity("cellNodes_InnerCells");

  // temporary storage for original P1 nodes:
  const MeshData::ConnTable cellNodesP1 = (*cellNodes);
//   const CFuint nbP1Nodes = cellNodesP1.size();

  //CFout << "There are total of " << nbP1Nodes << " P1 cell nodes\n";

  const CFuint nbElems = getNbElements();

//   CFout << "There are " << nbElems << " cells in the mesh \n\n";

  nbP1CellNodes.resize(nbElems);
//   nbP1FaceNodes.resize(nbElems);


// we will create a mesh with equal order on all the elements
  // and we reset the connectivity accordingly, using a std::valarray
  // first create the std::valarray
  std::valarray<CFuint> newPattern(nbElems);
  for (type_itr = elementType->begin(); type_itr != elementType->end(); ++type_itr)
  {

    // get the number of nodes in this element type
    // const CFuint newNbNodesPerElem = type_itr->getNbNodes();
    const CFuint newNbNodesPerElem = getNbDofInNewType(type_itr->getGeoShape(),m_newGeoPolyOrder) + \
                                     getNbInteriorDofInNewType(type_itr->getGeoShape(),m_newGeoPolyOrder);

    //get the number of nodes in this element type if it were just P1 (geometry)
    const CFuint nbNodesPerP1Elem = getNbDofInNewType(type_itr->getGeoShape(),CFPolyOrder::ORDER1);

    // get the number of elements
    const CFuint nbElemPerType = type_itr->getNbElems();

    // get start index of this element type in global element list
    CFuint globalIdx = type_itr->getStartIdx();

    // loop over elements of this type
    for (CFuint iElem = 0; iElem < nbElemPerType; ++iElem, ++globalIdx)
    {
      newPattern[globalIdx] = newNbNodesPerElem;
      nbP1CellNodes[globalIdx] = nbNodesPerP1Elem;
    }
  }
  // then resize the connectivity
  cellNodes->resize(newPattern);

  // loop on all the elements and reassign the IDs of the nodal states
  CFuint maxNodeID = 0;
  for (type_itr = elementType->begin(); type_itr != elementType->end(); ++type_itr)
  {

    // get the number of nodes in this element type IF IT WERE JUST P1 ELEMENT IN
    // TERMS OF GEOMETRY
    // This info is needed to copy P1 nodes form the temporary storage cellNodesP1 to
    // resized cellNodes
    const CFuint nbNodesPerP1Elem = getNbDofInNewType(type_itr->getGeoShape(),CFPolyOrder::ORDER1);
    //CFout << "The original number of nodes per type is " << nbNodesP1Type << "\n";


    // get the number of elements
    const CFuint nbElemPerType = type_itr->getNbElems();

    // get start index of this element type in global element list
    CFuint globalIdx = type_itr->getStartIdx();

    // loop over elements of this type and set the P1P1 nodes
    for (CFuint iElem = 0; iElem < nbElemPerType; ++iElem, ++globalIdx)
    {
//       CFout << "Cell: " << globalIdx << "\n";
      for (CFuint jNode = 0; jNode < nbNodesPerP1Elem; ++jNode)
      {
        (*cellNodes)(globalIdx,jNode) = cellNodesP1(globalIdx,jNode);
//  CFout << "Recopying the node: " << cellNodesP1(globalIdx, jNode) <<"\n";
  maxNodeID = std::max((*cellNodes)(globalIdx,jNode),maxNodeID);
  // CFout << "Checking nodal state: " << (*cellStates)(globalIdx,jState) << "\n";
      }
//       CFout << "\n";
    }
  }

  maxNodeID++;

  setExtraNodes(maxNodeID);

  _totalNewNbNodes = maxNodeID;
//  cf_assert(maxStateID == newPattern.sum());

}



//////////////////////////////////////////////////////////////////////////////

void RDS_HighOrderMeshUpdater::upgradeStateConnectivity()
{
  CFAUTOTRACE;

  CFout << "\n-----------------------------------------------------\n\n";
  CFout << "\tExecuting upgradeStateConnectivity()\n";
  CFout << "\n-----------------------------------------------------\n";


//CF_DEBUG_POINT;

CFout << "Updating existing states connectivity\n";
  // get the element type data
  SafePtr< std::vector< ElementTypeData > > elementType = getCFmeshData().getElementTypeData();
  std::vector< ElementTypeData >::iterator type_itr;

  // get the state connectivity
  Common::SafePtr<MeshData::ConnTable> cellStates = MeshDataStack::getActive()->getConnectivity("cellStates_InnerCells");

  //temporary storage for original P1 states
  const MeshData::ConnTable cellStatesP1 = (*cellStates);

  const CFuint nbElems = getNbElements();

  nbP1CellStates.resize(nbElems);

  // we will create a mesh with equal order on all the elements
  // and we reset the connectivity accordingly, using a std::valarray
  // first create the std::valarray
  std::valarray<CFuint> newPattern(nbElems);
  for (type_itr = elementType->begin(); type_itr != elementType->end(); ++type_itr)
  {
    // get the number of control volumes (states) in this element type
    const CFuint newNbStatesPerElem = getNbDofInNewType(type_itr->getGeoShape(),m_newSolPolyOrder) + \
                                        getNbInteriorDofInNewType(type_itr->getGeoShape(),m_newSolPolyOrder);
//     CFout << "\nNew number of of states in this element type is " << newNbStatesPerElem << "\n\n";

    //get the number of states in this element type if it were just P1 (in terms of solution)
    const CFuint nbStatesPerP1Elem = getNbDofInNewType(type_itr->getGeoShape(),CFPolyOrder::ORDER1);

    // get the number of elements
    const CFuint nbElemPerType = type_itr->getNbElems();

    // get start index of this element type in global element list
    CFuint globalIdx = type_itr->getStartIdx();

    // loop over elements of this type
    for (CFuint iElem = 0; iElem < nbElemPerType; ++iElem, ++globalIdx)
    {
      newPattern[globalIdx] = newNbStatesPerElem;
      nbP1CellStates[globalIdx] = nbStatesPerP1Elem;
    }
  }
  // then resize the connectivity
  cellStates->resize(newPattern);

  // loop on all the elements and reassign the IDs of the nodal states
  CFuint maxStateID = 0;
  for (type_itr = elementType->begin(); type_itr != elementType->end(); ++type_itr)
  {

    // get the number of states in this element type IF IT WERE JUST P1 ELEMENT IN
    // TERMS OF SOLUTION
    // This info is needed to copy P1 states from the temorary storage cellStatesP1 to
    // resized cellStates
    // get the number of nodes in this element type
    const CFuint nbStatesPerP1Elem = getNbDofInNewType(type_itr->getGeoShape(),CFPolyOrder::ORDER1);
//     CFout << "The original number of states per type is " << nbStatesPerP1Elem << "\n";

    // get the number of elements
    const CFuint nbElemPerType = type_itr->getNbElems();

    // get start index of this element type in global element list
    CFuint globalIdx = type_itr->getStartIdx();

    // loop over elements of this type and set the P1P1 states
    for (CFuint iElem = 0; iElem < nbElemPerType; ++iElem, ++globalIdx)
    {
//       CFout << "Cell: " << globalIdx << "\n";
      for (CFuint jState = 0; jState < nbStatesPerP1Elem; ++jState)
      {
        (*cellStates)(globalIdx,jState) = cellStatesP1(globalIdx,jState);
//       CFout << "Recopying the state: " << cellStatesP1(globalIdx, jState) <<"\n";
        maxStateID = std::max((*cellStates)(globalIdx,jState),maxStateID);
// CFout << "Checking nodal state: " << (*cellStates)(globalIdx,jState) << "\n";
      }
//           CFout << "\n";
    }
  }

  maxStateID++;
  setExtraStates(maxStateID);

  _totalNewNbStates = maxStateID;
//  cf_assert(maxStateID == newPattern.sum());
}


//////////////////////////////////////////////////////////////////////////////

void RDS_HighOrderMeshUpdater::setExtraNodes(CFuint& nextNodeID)
{

  CFout << "\n-----------------------------------------------------\n\n";
  CFout << "\tExecuting setExtraNodes()\n";
  CFout << "\n-----------------------------------------------------\n";



// get the element type data
  SafePtr< std::vector< ElementTypeData > > elementType = getCFmeshData().getElementTypeData();
  std::vector< ElementTypeData >::iterator type_itr;
  const CFuint nbElemTypes = elementType->size();
//   const CFuint dim = getCFmeshData().getDimension();

  // get the nodal connectivity
  Common::SafePtr<MeshData::ConnTable> cellNodes = MeshDataStack::getActive()->getConnectivity("cellNodes_InnerCells");

  // get nodal coordinates
//   CF_DEBUG_POINT;
  DataHandle<Node*, GLOBAL> nodeCoord = MeshDataStack::getActive()->getNodeDataSocketSink().getDataHandle();
//   CF_DEBUG_POINT;

  const CFuint nbElements = getNbElements();


  const CFuint nbNodesP1 = getNbNodes();
  CFout << "Number of P1 nodes is " << nbNodesP1 << "\n";
//   CFout << "Number of P1 nodes = " << nbNodesP1 << "\n\n";
//   CFout << "These are the nodal coordinates:\n\n";
//         for(CFuint iNode = 0; iNode < nbNodesP1; ++iNode) {
//    const Node& currNode = *nodeCoord[iNode];
//    CFout << "[" << currNode[XX] << "," << currNode[YY] << "]\n";
//         }
//


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

//    CFout << "Number of elements = " << nbElements << "\n";
//    CFout << "Number of cells = " << nbCells << "\n";

//     const CFuint nbNodesInCell = cellNodes->nbCols(iCell);
    const CFuint nbP1NodesInCell = nbP1CellNodes[iCell];
    totalMapSize += nbP1NodesInCell;
  }

  //Create a CFMultiMap that will store for each node, the related cells
  Common::CFMultiMap<CFuint,CFuint> mapNode2Cells;

  //Reserve the memory of totalSize
  mapNode2Cells.reserve(totalMapSize);


  //Loop over the cells
 for(CFuint iCell = 0; iCell < nbCells; ++iCell)
  {
//     const CFuint nbNodesInCell = cellNodes->nbCols(iCell);
    const CFuint nbNodesInP1Cell = nbP1CellNodes[iCell];

//     CFout << "Number of nodes in cell is " << nbNodesInCell << "\n";
//     CFout << "Number of nodes in P1 cell is " << nbP1CellNodes[iCell] << "\n";

    //Loop over the nodes of the cell and insert
    for(CFuint iNode = 0; iNode < nbNodesInP1Cell; ++iNode)
    {
      mapNode2Cells.insert((*cellNodes)(iCell,iNode), iCell);

//       CFout << "Inserting pair [" << (*cellNodes)(iCell,iNode) << "," << iCell << "]\n";
    }
//       CFout << "\n\n";
  }
  mapNode2Cells.sortKeys();


//   CF_DEBUG_POINT;


// CF_DEBUG_POINT;
// const CFuint nbFacesPerP1Element = m_faceNodeP1Element[0]->nbRows();
// CFout << "Number of faces per P1 element = " << nbFacesPerP1Element << "\n\n";
// const CFuint nbFacesPerPnElement = m_faceNodePnElement[0]->nbRows();
// CFout << "Number of faces per Pn element = " << nbFacesPerPnElement << "\n\n";

// CFout << "\n P1 Element local connectivity face->node:\n";
// CFout << (*m_faceNodeP1Element[0]);
// CFout << "\n\n";
//
// CFout << "\n Pn Element local connectivity face->node:\n";
// CFout << (*m_faceNodePnElement[0]);
// CFout << "\n\n";
// CF_DEBUG_POINT;


CFout << "Adding extra nodes on the faces\n";
  //Reset the isUpdated counter
  for(CFuint iElem=0; iElem < nbElements; ++iElem)
  {
    isUpdated[iElem] = std::numeric_limits<CFuint>::max();
  }

  //Loop over cells
  std::deque<CFuint> nextNodeIDReserve(0);
  std::deque<CFuint> lastCreatedNodeCell(0);
  std::deque<CFuint> lastCreatedNodeLocalID(0);

  type_itr = elementType->begin();


  for (CFuint iType = 0; iType < nbElemTypes; ++iType, ++type_itr)
  {
    const CFuint nbNodesInPnElement = (*elementType)[iType].getNbNodes();
    const CFuint nbElemPerType = (*elementType)[iType].getNbElems();
    const CFuint nbNodesInP1Element = getNbDofInNewType(type_itr->getGeoShape(),CFPolyOrder::ORDER1);

//  CF_DEBUG_OBJ(nbNodesPerPnElem);
//  CF_DEBUG_OBJ(nbElemPerType);
//  CF_DEBUG_OBJ(nbNodesPerP1Elem);

    // get start index of this element type in global element list
    CFuint globalIdx = (*elementType)[iType].getStartIdx();


/*
        for(CFuint globalNodeID = 0; globalNodeID < 9;++globalNodeID) {

        CFout << "Neighbors of node " << globalNodeID << ":\n";

  typedef CFMultiMap<CFuint, CFuint>::MapIterator MapIteratorMine;
        pair<MapIteratorMine, MapIteratorMine> neighCellsOfNode0Mine = mapNode2Cells.find(globalNodeID);
  for (MapIteratorMine neighCellItr1 = neighCellsOfNode0Mine.first;
               neighCellItr1 != neighCellsOfNode0Mine.second;
               ++neighCellItr1)
        {
          const CFuint neighCellIDMine = neighCellItr1->second;
    CFout << "\t \t" << neighCellIDMine << "\n";
  }

  }
*/


    // loop over elements of this type and set the P1P1 nodes
    for (CFuint iElem = 0; iElem < nbElemPerType; ++iElem, ++globalIdx)
    {
      const CFuint elemID = globalIdx;
      const CFuint nbFaces = LocalConnectionData::getInstance().getNbFacesInShape
        ((*elementType)[iType].getGeoShape());
      vector <bool> isAlreadyCreatedNode(nbNodesInPnElement, false);

//       CFout << "Number of nodes per element = " << nbNodesPerPnElem << "\n\n";


//  CFout << "\nElement ID = " << elemID << "\n";

      //-------------------------------------------------------------
      //Loop over the element faces to insert the face extra nodes
      //-------------------------------------------------------------
      for(CFuint iFace = 0; iFace < nbFaces; ++iFace)
      {

        const CFuint nbP1NodesPerFace = m_faceNodeP1Element[iType]->nbCols(iFace);
        const CFuint nbPnNodesPerFace = m_faceNodePnElement[iType]->nbCols(iFace);
//         CFout << "Number of nodes per P1 face is " << nbP1NodesPerFace << "\n";
//         CFout << "Number of nodes per Pn face is " << nbPnNodesPerFace << "\n\n";


        CFuint matchingCellID = std::numeric_limits<CFuint>::max();

        // Get the neighborCell of the face
        const CFuint firstNodeID = (*cellNodes)(elemID, (*m_faceNodeP1Element[iType])(iFace, 0));
        typedef CFMultiMap<CFuint, CFuint>::MapIterator MapIterator;

        bool fo4 = false;
	pair<MapIterator, MapIterator> neighCellsOfNode0 = mapNode2Cells.find(firstNodeID, fo4);
	if (!fo4) cout << "FEM_HighOrderMeshUpdater::setExtraNodes() => node " << firstNodeID << " not found!\n"; 
	
//  CFout << "\t Face " << iFace << "\n";
//  CFout << "\t\tfirstNodeID = " << firstNodeID << "\n";

        for (MapIterator neighCellItr1 = neighCellsOfNode0.first;
               neighCellItr1 != neighCellsOfNode0.second;
               ++neighCellItr1)
        {
          CFuint nbCommonNodes = 1;
          const CFuint neighCellID = neighCellItr1->second;
          if(neighCellID != elemID)
          {
            for (CFuint iNode = 1; iNode < nbP1NodesPerFace; ++iNode)
            {
              const CFuint localNodeID = (*m_faceNodeP1Element[iType])(iFace, iNode);
              const CFuint nodeID = (*cellNodes)(elemID, localNodeID);

	      bool fo5 = false;
              pair<MapIterator, MapIterator> neighCellsOfNode = mapNode2Cells.find(nodeID, fo5);
	      if (!fo5) cout << "FEM_HighOrderMeshUpdater::setExtraNodes() => node " << nodeID << " not found!\n"; 
	      
              for (MapIterator neighCellItr2 = neighCellsOfNode.first;
                     neighCellItr2 != neighCellsOfNode.second;
                     ++neighCellItr2)
              {
                const CFuint otherNeighCellID = neighCellItr2->second;
                if(otherNeighCellID == neighCellID) ++nbCommonNodes;
              }
            }
            if(nbCommonNodes == nbP1NodesPerFace) matchingCellID = neighCellID;
          }
        }



        // Check if the neighborCell has already been updated
        bool highOrder = false;
        if((matchingCellID != std::numeric_limits<CFuint>::max()) &&
           (isUpdated[matchingCellID] !=  std::numeric_limits<CFuint>::max())) highOrder = true;

        // If the neighbour is still P1, we should create nodes on the common face
        if(!highOrder)
        {

//    CFout << "\t\tMatching cell " << matchingCellID << " is still P1\n";

          for (CFuint iNode = 0; iNode < nbPnNodesPerFace; ++iNode)
          {
            const CFuint localNodeID = (*m_faceNodePnElement[iType])(iFace, iNode);
      //const CFuint nbNodesInP1Element = nbP1CellNodes[elemID];

            bool vertexNode(false);
            /// @todo find a safer way to check that it is a vertexNode or not...
            if(localNodeID < nbNodesInP1Element) vertexNode = true;

            if((isAlreadyCreatedNode[localNodeID] == false) && (!vertexNode))
            {
              if(nextNodeIDReserve.size() > 0){
                (*cellNodes)(elemID, localNodeID) = nextNodeIDReserve[0];
                nextNodeIDReserve.pop_front();
              }
              else{
                (*cellNodes)(elemID, localNodeID) = nextNodeID;

//    CFout << "\t\t\t Adding node: \n";
//    CFout << "\t\t\t Local node ID = " << localNodeID << "\n";
//    CFout << "\t\t\t Global node ID = " << nextNodeID << "\n\n\n";


                nextNodeID++;
                lastCreatedNodeCell.push_back(elemID);
                lastCreatedNodeLocalID.push_back(localNodeID);

                if(lastCreatedNodeCell.size() > 10) lastCreatedNodeCell.pop_front();
                if(lastCreatedNodeLocalID.size() > 10) lastCreatedNodeLocalID.pop_front();
              }
              isAlreadyCreatedNode[localNodeID] = true;
            }
          }
        }
        // if the neighbour is already P2, we should share the nodes on the common face
        else
        {

//  CFout << "\t\tMatching cell " << matchingCellID << " is already updated to P" << m_newGeoPolyOrder << "\n";

          // check which face of the neighbor is matching with the current face
          CFuint matchingNeighborFaceID = std::numeric_limits<CFuint>::max();

          const CFuint neighborType = isUpdated[matchingCellID];
          const CFuint nbNeighborFaces = m_faceNodeP1Element[neighborType]->nbRows();
          for(CFuint neighFace = 0; neighFace < nbNeighborFaces; ++neighFace)
          {
            const CFuint nbP1NodesInNeighborFace = m_faceNodeP1Element[neighborType]->nbCols(neighFace);
            CFuint nbMatchingNodes = 0;
            for(CFuint iNodeNeighFace = 0; iNodeNeighFace < nbP1NodesInNeighborFace; ++iNodeNeighFace)
            {
              const CFuint neighNodeID = (*m_faceNodeP1Element[neighborType])(neighFace, iNodeNeighFace);
              const CFuint neighNode = (*cellNodes)(matchingCellID, neighNodeID);
              for(CFuint iNodeFace = 0; iNodeFace < nbP1NodesPerFace; ++iNodeFace)
              {
                const CFuint currentNodeID = (*m_faceNodeP1Element[iType])(iFace, iNodeFace);
                const CFuint currentNode = (*cellNodes)(elemID, currentNodeID);
                if(currentNode == neighNode) nbMatchingNodes++;
              }
            }
            if(nbMatchingNodes == nbP1NodesPerFace) matchingNeighborFaceID = neighFace;
          }

          cf_assert(matchingNeighborFaceID != std::numeric_limits<CFuint>::max());

          //we have the matching face of the matching neighbor
          // for each node, we have to check if it is a vertex node
          //if not we can add it
          const CFuint nbPnNodesInNeighborFace = m_faceNodePnElement[neighborType]->nbCols(matchingNeighborFaceID);
          CFuint idx = 1;

          for (CFuint iNeighNode = 0; iNeighNode < nbPnNodesInNeighborFace; ++iNeighNode)
          {
            bool alreadyThere(false);

            const CFuint neighLocalNodeID = (*m_faceNodePnElement[iType])(matchingNeighborFaceID, iNeighNode);
            const CFuint neighNodeID = (*cellNodes)(matchingCellID, neighLocalNodeID);
//             CFout << "Should we add the state: " << (*cellNodes)(matchingCellID, neighLocalNodeID) <<"\n";

            for (CFuint iNode = 0; iNode < nbPnNodesPerFace; ++iNode)
            {
              const CFuint currentLocalNodeID = (*m_faceNodePnElement[iType])(iFace, iNode);
              const CFuint currentNodeID = (*cellNodes)(elemID, currentLocalNodeID);
              if( neighNodeID == currentNodeID) alreadyThere = true;
            }

            //the state localNeighNodeID should be added to the face iFace
            // but to which node in the face??
            if(!alreadyThere)
            {
//                CFout << "YES it is not yet there...\n";
//  CFout << "idx: " << nbPnNodesPerFace - idx <<" \n";
              CFuint newNodeLocalID = (*m_faceNodePnElement[iType])(iFace, nbPnNodesPerFace - idx);

              ///@todo here we are not sure about the order of the states to be added
              ///->to be modified or make sure that it is ok!!!!
//   CFout << "Trying to insert it at location: " << newNodeLocalID << " in Cell: " <<elemID << " \n";
//   CFout << "CellNodes connectivity of cell: " << elemID << " has size: " << cellNodes->nbCols(elemID) <<"\n";
              if(isAlreadyCreatedNode[newNodeLocalID] == true)
              {
                nextNodeIDReserve.push_back((*cellNodes)(elemID, newNodeLocalID));
              }

              (*cellNodes)(elemID, newNodeLocalID) = neighNodeID;
//        CFout << "\t\tAdding node: " << elemID << " Local ID = " << newNodeLocalID << " Global ID = " << neighNodeID << "\n";
              ++idx;
            }
          }
        } //if neighbor is already P2, we share the nodes on common face

        //if there is no neighbour sharing the face iFace,
        //then it is a boundary face, so save its connectivity for later
        if(matchingCellID == std::numeric_limits<CFuint>::max())
        {
// CFout << "This is a boundary face!!!\n";
          BoundaryFaceConnect faceNodeConn;
          faceNodeConn.resize(nbPnNodesPerFace);

          for (CFuint iNode = 0; iNode < nbPnNodesPerFace; ++iNode)
          {
            const CFuint localNodeID = (*m_faceNodePnElement[iType])(iFace, iNode);
            faceNodeConn[iNode] = &(*cellNodes)(elemID, localNodeID);
          }

          _boundaryFacesNodes.push_back(faceNodeConn);
        }
      } //Loop over faces



        //
        if((elemID == nbCells-1) && (nextNodeIDReserve.size() > 0))
        {
            CFout << "There are " << nextNodeIDReserve.size() <<"remaining unused nodeID\n";

cf_assert(lastCreatedNodeCell.size() > nextNodeIDReserve.size());
cf_assert(lastCreatedNodeLocalID.size() > nextNodeIDReserve.size());
for(CFuint iNode = 0 ; iNode < nextNodeIDReserve.size(); ++iNode)
{
const CFuint cellID2Replace = lastCreatedNodeCell.size()-1-iNode;
const CFuint localNodeID2Replace = lastCreatedNodeLocalID.size()-1-iNode;

// CFout << "Last created node is in cell: " << lastCreatedNodeCell[cellID2Replace]<< " at localID: " << lastCreatedNodeLocalID[localNodeID2Replace]<<"\n";
// CFout << "Replacing: "<<(*cellNodes)(lastCreatedNodeCell[cellID2Replace],lastCreatedNodeLocalID[localNodeID2Replace])<<"\n";

(*cellNodes)(lastCreatedNodeCell[cellID2Replace],lastCreatedNodeLocalID[localNodeID2Replace]) = nextNodeIDReserve[iNode];

//             CFout << "by: "<<(*cellNodes)(lastCreatedNodeCell[cellID2Replace],lastCreatedNodeLocalID[localNodeID2Replace])<<"\n";
  nextNodeID -= 1;
}
// CFout << "Total nb nodes: "<< nextNodeID<<"\n";
nextNodeIDReserve.resize(0);
        }



      const CFuint nbInteriorNodes = getNbInteriorDofInNewType((*elementType)[iType].getGeoShape(),m_newGeoPolyOrder);
//       CFout << "Number of Pn nodes = " << nbNodesInPnElement << "\n";
//       CFout << "Number of interior nodes = " << nbInteriorNodes << "\n";

      for(CFuint iInteriorNode=0; iInteriorNode < nbInteriorNodes ; ++iInteriorNode)
      {
        (*cellNodes)(elemID,nbNodesInPnElement-nbInteriorNodes) = nextNodeID;
//         CFout << "Creating the central node: " << nextNodeID <<"\n";
        nextNodeID++;
      }



      isUpdated[elemID] = iType;


    } //loop over elements

    } //loop over element types



//     CFout << "\n\n\Listing boundary faces for nodes:\n\n";
//       for(CFuint iFace = 0; iFace < _boundaryFacesNodes.size(); ++iFace)
//       {
//         for(CFuint iNode = 0; iNode < _boundaryFacesNodes[iNode].size(); ++iNode)
//           {
//             CFout << *((_boundaryFacesNodes[iFace])[iNode]) << "  ";
//
//
//           }
//           CFout << "\n";
//       }




///---------------------------------------------------------------------------------
///
///   Print complete info about nodal connectivity
///
///---------------------------------------------------------------------------------


// CFout << "\n\n--------------------------------------------------------------------------------\n";
// CFout << "\t\tPRINTING INFORMATION ABOUT NODAL CONNECTIVITY:\n";
// CFout << "--------------------------------------------------------------------------------\n\n";
//
//
// type_itr = elementType->begin();
//
// for (CFuint iType = 0; iType < nbElemTypes; ++iType, ++type_itr) {
//
// //   const CFuint nbNodesPerPnElem = (*elementType)[iType].getNbNodes();
//
//  // get start index of this element type in global element list
//      CFuint globalIdx = (*elementType)[iType].getStartIdx();
//  const CFuint nbElemPerType = (*elementType)[iType].getNbElems();
//
//
//  // loop over elements of this type and print connectivity info
//      for (CFuint iElem = 0; iElem < nbElemPerType; ++iElem, ++globalIdx) {
//        const CFuint elemID = globalIdx;
//        const CFuint nbFaces = LocalConnectionData::getInstance().getNbFacesInShape
//           ((*elementType)[iType].getGeoShape());
//
//    CFout << "Element " << globalIdx << "\n";
//
//    for(CFuint iFace = 0; iFace < nbFaces; ++iFace) {
//            const CFuint nbPnNodesPerFace = m_faceNodePnElement[iType]->nbCols(iFace);
//
//      CFout << "\t";
//      for (CFuint iNode = 0; iNode < nbPnNodesPerFace; ++iNode) {
//        const CFuint localNodeID = (*m_faceNodePnElement[iType])(iFace, iNode);
//        const CFuint nodeID = (*cellNodes)(elemID, localNodeID);
//        CFout << nodeID << " ";
//      }
//      CFout << "\n";
//
//    }
//    CFout << "\n";
//  }
//
// }
//

// CFout << "\n############ NODAL CONNECTIVITY ###########\n";
// CFout << *cellNodes << "\n";




//   CFout << "\n-----------------------------------------------------\n\n";
//   CFout << "\tFinished executing setExtraNodes()\n";
//   CFout << "\n-----------------------------------------------------\n";
//
//   CFout << "\n\n\n";


}  //End of setExtraNodes(CFuint& nextNodeID)


//////////////////////////////////////////////////////////////////////////////

void RDS_HighOrderMeshUpdater::setExtraStates(CFuint& nextStateID)
{

  CFout << "\n-----------------------------------------------------\n\n";
  CFout << "\tExecuting setExtraStates()\n";
  CFout << "\n-----------------------------------------------------\n";

  // get the element type data
  SafePtr< std::vector< ElementTypeData > > elementType = getCFmeshData().getElementTypeData();
  std::vector< ElementTypeData >::iterator type_itr;
  const CFuint nbElemTypes = elementType->size();
//   const CFuint dim = getCFmeshData().getDimension();

  // get the state connectivity
  Common::SafePtr<MeshData::ConnTable> cellStates = MeshDataStack::getActive()->getConnectivity("cellStates_InnerCells");

  const CFuint nbElements = getNbElements();

  //marking if it is updated and the element type
  std::vector< CFuint > isUpdated(nbElements);
  for(CFuint iElem=0; iElem < nbElements; ++iElem)
  {
    isUpdated[iElem] = std::numeric_limits<CFuint>::max();
  }

CFout << "Creating neighbour connectivity\n";
  //Create neighbour connectivity
  const CFuint nbCells = cellStates->nbRows();
  CFuint totalMapSize = 0;
  for(CFuint iCell = 0; iCell < nbCells; ++iCell)
  {
    const CFuint nbP1StatesInCell = nbP1CellStates[iCell];
    totalMapSize += nbP1StatesInCell;
  }

  //Create a CFMultiMap that will store for each state, the related cells
  Common::CFMultiMap<CFuint,CFuint> mapState2Cells;

  //Reserve the memory of totalSize
  mapState2Cells.reserve(totalMapSize);

  //Loop over the cells
  for(CFuint iCell = 0; iCell < nbCells; ++iCell)
  {
    const CFuint nbStatesInP1Cell = nbP1CellStates[iCell];
    //Loop over the nodes of the cell and insert
    for(CFuint iState = 0; iState < nbStatesInP1Cell; ++iState)
    {
      mapState2Cells.insert((*cellStates)(iCell,iState), iCell);
//       CFout << "Inserting pair [" << (*cellStates)(iCell,iState) << "," << iCell << "]\n";
    }
//       CFout << "\n\n";
  }
  mapState2Cells.sortKeys();



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

  type_itr = elementType->begin();

  for (CFuint iType = 0; iType < nbElemTypes; ++iType, ++type_itr)
  {
    const CFuint nbStatesInPnElement = (*elementType)[iType].getNbStates();
    const CFuint nbElemPerType = (*elementType)[iType].getNbElems();
    const CFuint nbStatesInP1Element = getNbDofInNewType(type_itr->getGeoShape(),CFPolyOrder::ORDER1);

//     CF_DEBUG_OBJ(nbStatesInPnElement);
//     CF_DEBUG_OBJ(nbElemPerType);
//     CF_DEBUG_OBJ(nbStatesInP1Element);


    // get start index of this element type in global element list
    CFuint globalIdx = (*elementType)[iType].getStartIdx();

    // loop over elements of this type and set the P1P1 states
    for (CFuint iElem = 0; iElem < nbElemPerType; ++iElem, ++globalIdx)
    {
      const CFuint elemID = globalIdx;
      const CFuint nbFaces = LocalConnectionData::getInstance().getNbFacesInShape
        ((*elementType)[iType].getGeoShape());
      vector <bool> isAlreadyCreatedState(nbStatesInPnElement, false);


//       CFout << "\nElement ID = " << elemID << "\n";

      //--------------------------------------------------------------
      //Loop over the element faces to insert the face extra states
      //--------------------------------------------------------------
      for(CFuint iFace = 0; iFace < nbFaces; ++iFace)
      {
//  CFout << "Face: " << iFace <<"\n";

        const CFuint nbP1StatesPerFace = m_faceStateP1Element[iType]->nbCols(iFace);
        const CFuint nbPnStatesPerFace = m_faceStatePnElement[iType]->nbCols(iFace);

        CFuint matchingCellID = std::numeric_limits<CFuint>::max();

        // Get the neighborCell of the face
        const CFuint firstStateID = (*cellStates)(elemID, (*m_faceStateP1Element[iType])(iFace, 0));
        typedef CFMultiMap<CFuint, CFuint>::MapIterator MapIterator;

	bool fo = false;
        pair<MapIterator, MapIterator> neighCellsOfState0 = mapState2Cells.find(firstStateID, fo);
	if (!fo) cout << "RDS_HighOrderMeshUpdater::setExtraStates() => state " << firstStateID << " not found!\n";
	
	//         CFout << "\t Face " << iFace << "\n";
	
        for (MapIterator neighCellItr1 = neighCellsOfState0.first;
               neighCellItr1 != neighCellsOfState0.second;
               ++neighCellItr1)
        {
          CFuint nbCommonStates = 1;
          const CFuint neighCellID = neighCellItr1->second;
          if(neighCellID != elemID)
          {
            for (CFuint iState = 1; iState < nbP1StatesPerFace; ++iState)
            {
              const CFuint localStateID = (*m_faceStateP1Element[iType])(iFace, iState);
              const CFuint stateID = (*cellStates)(elemID, localStateID);
	      
	      bool fo1 = false;
              pair<MapIterator, MapIterator> neighCellsOfState = mapState2Cells.find(stateID, fo1);
	      if (!fo1) cout << "RDS_HighOrderMeshUpdater::setExtraStates() => state " << stateID << " not found!\n";
	      
              for (MapIterator neighCellItr2 = neighCellsOfState.first;
                     neighCellItr2 != neighCellsOfState.second;
                     ++neighCellItr2)
              {
                const CFuint otherNeighCellID = neighCellItr2->second;
                if(otherNeighCellID == neighCellID) ++nbCommonStates;
              }
            }
            if(nbCommonStates == nbP1StatesPerFace) matchingCellID = neighCellID;
          }
        }

        // Check if the neighborCell has already been updated
        bool highOrder = false;
        if((matchingCellID != std::numeric_limits<CFuint>::max()) &&
           (isUpdated[matchingCellID] !=  std::numeric_limits<CFuint>::max())) highOrder = true;

        // If the neighbour is still P1, we should create states on the common face
        if(!highOrder)
        {

//           CFout << "\t\tMatching cell " << matchingCellID << " is still P1\n";

          for (CFuint iState = 0; iState < nbPnStatesPerFace; ++iState)
          {
            const CFuint localStateID = (*m_faceStatePnElement[iType])(iFace, iState);

            bool nodalState(false);
            /// @todo find a safer way to check that it is a nodalState or not...
            if(localStateID < nbStatesInP1Element) nodalState = true;

            if((isAlreadyCreatedState[localStateID] == false) && (!nodalState))
            {
              if(nextStateIDReserve.size() > 0){
                (*cellStates)(elemID, localStateID) = nextStateIDReserve[0];
                nextStateIDReserve.pop_front();
              }
              else{
                (*cellStates)(elemID, localStateID) = nextStateID;

//                 CFout << "\t\t\t Adding state: \n";
//                 CFout << "\t\t\t Local state ID = " << localStateID << "\n";
//                 CFout << "\t\t\t Global state ID = " << nextStateID << "\n\n\n";
//

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

//         CFout << "\t\tMatching cell " << matchingCellID << " is already updated to P" << m_newSolPolyOrder << "\n";

          // check which face of the neighbor is matching with the current face
          CFuint matchingNeighborFaceID = std::numeric_limits<CFuint>::max();

          const CFuint neighborType = isUpdated[matchingCellID];
          const CFuint nbNeighborFaces = m_faceStateP1Element[neighborType]->nbRows();
          for(CFuint neighFace = 0; neighFace < nbNeighborFaces; ++neighFace)
          {
            const CFuint nbP1StatesInNeighborFace = m_faceStateP1Element[neighborType]->nbCols(neighFace);
            CFuint nbMatchingStates = 0;
            for(CFuint iStateNeighFace = 0; iStateNeighFace < nbP1StatesInNeighborFace; ++iStateNeighFace)
            {
              const CFuint neighStateID = (*m_faceStateP1Element[neighborType])(neighFace, iStateNeighFace);
              const CFuint neighState = (*cellStates)(matchingCellID, neighStateID);
              for(CFuint iStateFace = 0; iStateFace < nbP1StatesPerFace; ++iStateFace)
              {
                const CFuint currentStateID = (*m_faceStateP1Element[iType])(iFace, iStateFace);
                const CFuint currentState = (*cellStates)(elemID, currentStateID);
                if(currentState == neighState) nbMatchingStates++;
              }
            }
            if(nbMatchingStates == nbP1StatesPerFace) matchingNeighborFaceID = neighFace;
          }

          cf_assert(matchingNeighborFaceID != std::numeric_limits<CFuint>::max());

          //we have the matching face of the matching neighbor
          // for each state, we have to check if it is a nodal state
          //if not we can add it
          const CFuint nbPnStatesInNeighborFace = m_faceStatePnElement[neighborType]->nbCols(matchingNeighborFaceID);
          CFuint idx = 1;

          for (CFuint iNeighState = 0; iNeighState < nbPnStatesInNeighborFace; ++iNeighState)
          {
            bool alreadyThere(false);

            const CFuint neighLocalStateID = (*m_faceStatePnElement[iType])(matchingNeighborFaceID, iNeighState);
            const CFuint neighStateID = (*cellStates)(matchingCellID, neighLocalStateID);
//             CFout << "Should we add the state: " << (*cellStates)(matchingCellID, neighLocalStateID) <<"\n";

            for (CFuint iState = 0; iState < nbPnStatesPerFace; ++iState)
            {
              const CFuint currentLocalStateID = (*m_faceStatePnElement[iType])(iFace, iState);
              const CFuint currentStateID = (*cellStates)(elemID, currentLocalStateID);
              if( neighStateID == currentStateID) alreadyThere = true;
            }

            //the state localNeighStateID should be added to the face iFace
            // but to which state in the face??
            if(!alreadyThere)
            {
//                CFout << "YES it is not yet there...\n";
//  CFout << "idx: " << nbStatesPerFace - idx <<" \n";
              CFuint newStateLocalID = (*m_faceStatePnElement[iType])(iFace, nbPnStatesPerFace - idx);

              ///@todo here we are not sure about the order of the states to be added
              ///->to be modified or make sure that it is ok!!!!
//   CFout << "Trying to insert it at location: " << newStateLocalID << " in Cell: " <<elemID << " \n";
//   CFout << "CellStates connectivity of cell: " << elemID << " has size: " << cellStates->nbCols(elemID) <<"\n";
              if(isAlreadyCreatedState[newStateLocalID] == true)
              {
                nextStateIDReserve.push_back((*cellStates)(elemID, newStateLocalID));
              }

              (*cellStates)(elemID, newStateLocalID) = neighStateID;
//               CFout << "\t\tAdding state: " << elemID << " Local ID = " << newStateLocalID << " GlobalID = " << neighStateID << "\n";
              ++idx;
            }
          }
        }

        //if there is no neighbour sharing the face iFace,
        //then it is a boundary face, so save its connectivity for later
        if(matchingCellID == std::numeric_limits<CFuint>::max())
        {
// CFout << "This is a boundary face!!!\n";
          BoundaryFaceConnect faceStateConn;
          faceStateConn.resize(nbPnStatesPerFace);

          for (CFuint iState = 0; iState < nbPnStatesPerFace; ++iState)
          {
            const CFuint localStateID = (*m_faceStatePnElement[iType])(iFace, iState);
            faceStateConn[iState] = &(*cellStates)(elemID, localStateID);
          }

          _boundaryFacesStates.push_back(faceStateConn);
        }
      }




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

      const CFuint nbInteriorStates = getNbInteriorDofInNewType((*elementType)[iType].getGeoShape(),m_newSolPolyOrder);
      for(CFuint iInteriorState=0; iInteriorState < nbInteriorStates ; ++iInteriorState)
      {
        (*cellStates)(elemID,nbStatesInPnElement-nbInteriorStates) = nextStateID;
// CFout << "Creating the central state: " << nextStateID <<"\n";
        nextStateID++;
      }

      isUpdated[elemID] = iType;
    }
  }


//     CFout << "\n\n\Listing boundary faces for states:\n\n";
//       for(CFuint iFace = 0; iFace < _boundaryFacesStates.size(); ++iFace)
//       {
//         for(CFuint iState = 0; iState < _boundaryFacesStates[iState].size(); ++iState)
//           {
//             CFout << *((_boundaryFacesStates[iFace])[iState]) << "  ";
//
//
//           }
//           CFout << "\n";
//       }



///---------------------------------------------------------------------------------
///
///   Print complete info about state connectivity
///
///---------------------------------------------------------------------------------


// CFout << "\n\n--------------------------------------------------------------------------------\n";
// CFout << "\t\tPRINTING INFORMATION ABOUT STATE CONNECTIVITY:\n";
// CFout << "--------------------------------------------------------------------------------\n\n";
//
//
// type_itr = elementType->begin();
//
// for (CFuint iType = 0; iType < nbElemTypes; ++iType, ++type_itr) {
//
// //   const CFuint nbNodesPerPnElem = (*elementType)[iType].getNbNodes();
//
//  // get start index of this element type in global element list
//      CFuint globalIdx = (*elementType)[iType].getStartIdx();
//  const CFuint nbElemPerType = (*elementType)[iType].getNbElems();
//
//
//  // loop over elements of this type and print connectivity info
//      for (CFuint iElem = 0; iElem < nbElemPerType; ++iElem, ++globalIdx) {
//        const CFuint elemID = globalIdx;
//        const CFuint nbFaces = LocalConnectionData::getInstance().getNbFacesInShape
//           ((*elementType)[iType].getGeoShape());
//
//    CFout << "Element " << globalIdx << "\n";
//
//    for(CFuint iFace = 0; iFace < nbFaces; ++iFace) {
//            const CFuint nbPnStatesPerFace = m_faceStatePnElement[iType]->nbCols(iFace);
//
//      CFout << "\t";
//      for (CFuint iState = 0; iState < nbPnStatesPerFace; ++iState) {
//        const CFuint localStateID = (*m_faceStatePnElement[iType])(iFace, iState);
//        const CFuint stateID = (*cellStates)(elemID, localStateID);
//        CFout << stateID << " ";
//      }
//      CFout << "\n";
//
//    }
//    CFout << "\n";
//  }
//
// }
//


// Common::SafePtr<MeshData::ConnTable> cellNodes = MeshDataStack::getActive()->getConnectivity("cellNodes_InnerCells");
//
// CFout << "\n############ NODAL CONNECTIVITY ###########\n";
// CFout << *cellNodes << "\n";
// CFout << "\n############ STATE CONNECTIVITY ###########\n";
// CFout << *cellStates << "\n";
//
// CF_DEBUG_EXIT;
/*
  CFout << "\n-----------------------------------------------------\n\n";
  CFout << "\tFinished executing setExtraStates()\n";
  CFout << "\n-----------------------------------------------------\n";

  CFout << "\n\n\n";
*/




}

//////////////////////////////////////////////////////////////////////////////

void RDS_HighOrderMeshUpdater::updateTRSData()
{
  CFAUTOTRACE;

  CFout << "\n-----------------------------------------------------\n\n";
  CFout << "\tExecuting updateTRSData()\n";
  CFout << "\n-----------------------------------------------------\n";


  //The array _boundaryFacesNodes stores all boundary faces with updated
  //connectivity these faces are not stored per TRS, but 'randomly'
  //we loop over TRS, find corresponding faces in these TRS and update them
  //using information from _boundaryFacesNodes
  //The same is then done for states

  CFout << "Converting the TRS data to higher order ...\n";


  SafePtr< vector<CFuint> > nbTRs = getCFmeshData().getNbTRs();
  SafePtr< vector<vector<CFuint> > > nbGeomEntsPerTR = getCFmeshData().getNbGeomEntsPerTR();

  CFuint matchingFaceID;

  // loop on TRS's
  for(CFuint iTRS = 0; iTRS < getCFmeshData().getNbTRSs(); ++iTRS )
  {

//    CFout << "############  TRS : [" << iTRS << "]\n";

    // loop on each TR on the TRS
    TRGeoConn& geoConn = getCFmeshData().getTRGeoConn(iTRS);
    for(CFuint iTR = 0; iTR < (*nbTRs)[iTRS]; ++iTR )
    {
//       CFout << "############  TR : [" << iTR << "]\n";
      GeoConn& geoC = geoConn[iTR];
      const CFuint nbGeos = geoC.size();

      // loop on each GeometricEntity on the TR
      for (CFuint iGeo = 0; iGeo < nbGeos; ++iGeo)
      {
//         CFout << "############  GEO : [" << iGeo << "]\n";
        std::valarray<CFuint>& nodeCon  = geoC[iGeo].first;


// for (int i = 0; i < nodeCon.size(); i++) CFout << " " << nodeCon[i]; CFout << "\n";

          std::valarray<CFuint>& stateCon = geoC[iGeo].second;

          const CFuint oldNbNodesInGeo  =  nodeCon.size();
          const CFuint oldNbStatesInGeo = stateCon.size();


///////// FROM HERE (NODES)


//         CFout << "\n\n---------------------------------------------------------------\n\n";
//         CFout << "\t\t Updating nodes ...";
//         CFout << "\n\n---------------------------------------------------------------\n\n";

if((m_newGeoPolyOrder > 1) && m_updateGeometry) {

        const CFuint newNbNodesInGeo  = getNewNbDofInBoundaryGeo(oldNbNodesInGeo, m_newGeoPolyOrder);

// for (int i = 0; i < nodeCon.size(); i++) CFout << " " << nodeCon[i]; CFout << "\n";

        matchingFaceID = std::numeric_limits<CFuint>::max();
        const CFuint nbBoundaryFacesNodes = _boundaryFacesNodes.size();

//         CFout << "Trying to get the corresponding boundaryFace among "<< nbBoundaryFacesNodes << " faces \n";

        // Look for the matching neighbor face by matching the nodes
        for(CFuint iFace = 0; iFace < nbBoundaryFacesNodes; ++iFace)
        {
//           CFout << "\tTrying to match face: "<< iFace<< "\n";
          vector<CFuint*>& faceNodes = _boundaryFacesNodes[iFace];
          CFuint matchingNodes = 0;
          for(CFuint iNode = 0; iNode < nodeCon.size(); iNode++)
          {
//           CFout << "\tNode: "<< nodeCon[iNode] << "\n";

            for(CFuint jNode = 0; jNode < faceNodes.size(); jNode++)
            {
//             CFout << "\t\tBFace node: "<< *(faceNodes[jNode]) << "\n";
              if(*(faceNodes[jNode]) == nodeCon[iNode]) matchingNodes++;
            }
          }

          if (matchingNodes == nodeCon.size())
          {
            matchingFaceID = iFace;
//             CFout << "\t\t\t Matching face ID = " << matchingFaceID << "\n";
            break; // dont need to keep searching
          }
        }
// CF_DEBUG_POINT;
        // found matching face and now we update is connectivity
        cf_assert(matchingFaceID != std::numeric_limits<CFuint>::max());
        vector<CFuint*>& faceNodes = _boundaryFacesNodes[matchingFaceID];
        nodeCon.resize(newNbNodesInGeo);
        for(CFuint iNode = 0; iNode < faceNodes.size(); ++iNode)
        {
          nodeCon[iNode] = *(faceNodes[iNode]);
//           CFout << nodeCon[iNode] << " ";
        }
//             CFout << "\n";
// CF_DEBUG_POINT;
///////// TO HERE (NODES)

}



///////// FROM HERE (STATES)

if((m_newSolPolyOrder > 1) && m_updateSolution) {

//         CFout << "\n\n---------------------------------------------------------------\n\n";
//         CFout << "\t\t Updating states ...";
//         CFout << "\n\n---------------------------------------------------------------\n\n";


       const CFuint newNbStatesInGeo = getNewNbDofInBoundaryGeo(oldNbStatesInGeo, m_newSolPolyOrder);

        matchingFaceID = std::numeric_limits<CFuint>::max();
        const CFuint nbBoundaryFacesStates = _boundaryFacesStates.size();

//         CFout << "Trying to get the corresponding boundaryFace among "<< nbBoundaryFaces << " faces \n";

        // Look for the matching neighbor face by matching the nodes
        for(CFuint iFace = 0; iFace < nbBoundaryFacesStates; ++iFace)
        {
//           CFout << "\tTrying to match face: "<< iFace<< "\n";
          vector<CFuint*>& faceStates = _boundaryFacesStates[iFace];
          CFuint matchingStates = 0;
          for(CFuint iState = 0; iState < stateCon.size(); iState++)
          {
//           CFout << "\tNode: "<< nodeCon[iNode] << "\n";

            for(CFuint jState = 0; jState < faceStates.size(); jState++)
            {
//             CFout << "\t\tBFace node: "<< *(faceNodes[jNode]) << "\n";
              if(*(faceStates[jState]) == stateCon[iState]) matchingStates++;
            }
          }

          if (matchingStates == stateCon.size())
          {
            matchingFaceID = iFace;
//             CFout << "\t\t\t Matching face ID = " << matchingFaceID << "\n";
            break; // dont need to keep searching
          }
        }
// CF_DEBUG_POINT;
        // found matching face and now we update is connectivity
        cf_assert(matchingFaceID != std::numeric_limits<CFuint>::max());
        vector<CFuint*>& faceStates = _boundaryFacesStates[matchingFaceID];
        stateCon.resize(newNbStatesInGeo);
        for(CFuint iState = 0; iState < faceStates.size(); ++iState)
        {
          stateCon[iState] = *(faceStates[iState]);
//           CFout << "nodes: " << nodeCon[ijNode] <<"\n";
        }
// CF_DEBUG_POINT;
///////// TO HERE (STATES)

}

      } // loop on geoents of the TR

    } // loop on the TR's of the TRS

  } // loop on the TRS's

}

//////////////////////////////////////////////////////////////////////////////

void RDS_HighOrderMeshUpdater::recreateNodes()
{
  CFAUTOTRACE;

  CFout << "\n-----------------------------------------------------\n\n";
  CFout << "\tExecuting recreateNodes()\n";
  CFout << "\n-----------------------------------------------------\n";



  Common::SafePtr<MeshData::ConnTable> cellNodes = MeshDataStack::getActive()->getConnectivity("cellNodes_InnerCells");

  DataHandle < Framework::Node*, Framework::GLOBAL > nodes = getCFmeshData().getNodesHandle();
  //   DataHandle < Framework::Node*, Framework::GLOBAL > nodes =
  // MeshDataStack::getActive()->getNodeDataSocketSink().getDataHandle();

  const CFuint newNbNodes = _totalNewNbNodes;
  const CFuint oldNbNodes = nodes.size();
  CFout << "Updating the nb of nodes from: " << oldNbNodes << "  to  " << newNbNodes <<"\n";

  std::vector<Node*> oldNodes(oldNbNodes);

  // backup the existing nodes
  for (CFuint i = 0; i < oldNbNodes; ++i)
  {
    oldNodes[i] = nodes[i];
  }


  // Resize the datahandle for the states
  nodes.resize(newNbNodes);

  // copy back the original P1 nodes
  for (CFuint i = 0; i < oldNbNodes; ++i)
  {
    nodes[i] = oldNodes[i];
    nodes[i]->setIsOwnedByState(false);
    nodes[i]->setIsOnMesh(true);
  }


//  CFout << "Node info ... before recreating the nodal data\n";
//           for(CFuint iNode = 0; iNode < oldNbNodes; ++iNode) {
//                 CFout << "Node " << nodes[iNode]->getGlobalID() << "\t[" << (*nodes[iNode])[XX] << "," << (*nodes[iNode])[YY] << "]\t is on mesh: "
//                       << nodes[iNode]->isOnMesh() << "\t and is owned by state: " << nodes[iNode]->isOwnedByState() << "\n";
//           }



//   // allocate memory for the new nodes
//   for (CFuint i = oldNbNodes; i < newNbNodes; ++i)
//   {
//     nodes[i] = new Node();
//   }


  //Set coordinates of new nodes by interpolation between P1 nodes

  //Number of nodes that need to be added: all non-vertex nodes:
  const CFuint nbPnNodes = newNbNodes - oldNbNodes;


  //Indexes of Pn nodes whose coordinates have been updated:
  vector<bool> isUpdated(nbPnNodes);

  for(CFuint iNode = 0; iNode < nbPnNodes; ++iNode)  isUpdated[iNode] = false;


  SafePtr< std::vector< ElementTypeData > > elementType = getCFmeshData().getElementTypeData();
  std::vector< ElementTypeData >::iterator type_itr;
  const CFuint nbElemTypes = elementType->size();

  ///@todo: check dimension. If dim = 3, loop over element edges first and set nodal coordinates
  ///       then set coordinates of nodes that are inside faces or inside elements

/*
try
  {
    throw Common::NotImplementedException(FromHere(),"RDS_HighOrderMeshUpdater::dim = DIM_3D");
  }
  catch (NotImplementedException& e)
  {
  }
*/


  RealVector nodeCoord(2);

  type_itr = elementType->begin();

  for (CFuint iType = 0; iType < nbElemTypes; ++iType, ++type_itr) {

  // get element shape
  const CFGeoShape::Type elmtShape = type_itr->getGeoShape();

  // get number of nodes in this type

//   const CFuint nbNodesPerElem = getNbDofInNewType(type_itr->getGeoShape(),m_newGeoPolyOrder) + getNbInteriorDofInNewType(type_itr->getGeoShape(),m_newGeoPolyOrder);

  // get start index of this element type in global element list
      CFuint globalIdx = (*elementType)[iType].getStartIdx();
  const CFuint nbElemPerType = (*elementType)[iType].getNbElems();


  // loop over elements of this type and print connectivity info
      for (CFuint iElem = 0; iElem < nbElemPerType; ++iElem, ++globalIdx) {

        const CFuint elemID = globalIdx;
        const CFuint nbFaces = LocalConnectionData::getInstance().getNbFacesInShape
         ((*elementType)[iType].getGeoShape());


    for(CFuint iFace = 0; iFace < nbFaces; ++iFace) {
            const CFuint nbPnNodesPerFace = m_faceNodePnElement[iType]->nbCols(iFace);

                        //2D: the first and second node indexes in the list of nodes per each face
                        //are the indexes of the vertex nodes, all following nodes indexes refer to
                        //nodes between these two
                        const CFuint localVertexNodeID1 = (*m_faceNodePnElement[iType])(iFace,0);
                        const CFuint localVertexNodeID2 = (*m_faceNodePnElement[iType])(iFace,1);

                        const CFuint vertexNodeID1 = (*cellNodes)(elemID,localVertexNodeID1);
                        const CFuint vertexNodeID2 = (*cellNodes)(elemID,localVertexNodeID2);

                        const CFreal dx = ((*nodes[vertexNodeID2])[XX] - (*nodes[vertexNodeID1])[XX])/(nbPnNodesPerFace-1);
                        const CFreal dy = ((*nodes[vertexNodeID2])[YY] - (*nodes[vertexNodeID1])[YY])/(nbPnNodesPerFace-1);

      for (CFuint iNode = 2; iNode < nbPnNodesPerFace; ++iNode) {
        const CFuint localNodeID = (*m_faceNodePnElement[iType])(iFace, iNode);
        const CFuint nodeID = (*cellNodes)(elemID, localNodeID);

                                if (!isUpdated[nodeID-oldNbNodes]) {

                                     nodeCoord[XX] = (*nodes[vertexNodeID1])[XX] + (iNode-1)*dx;
                                     nodeCoord[YY] = (*nodes[vertexNodeID1])[YY] + (iNode-1)*dy;

                                    getCFmeshData().createNode(nodeID,nodeCoord);

                                    isUpdated[nodeID-oldNbNodes] = true;
                                    nodes[nodeID]->setGlobalID(nodeID);
                                    nodes[nodeID]->setIsOnMesh(true);
                                    nodes[nodeID]->setIsOwnedByState(false);
                                }
      }
    } //Loop over faces - new nodes on faces inserted

      //Insert inner nodes, at least in case of P3 triangles:

      if(m_newGeoPolyOrder_int == 3) {
        switch(elmtShape) {

          case CFGeoShape::TRIAG: {

            const CFuint locID0 = (*m_faceNodePnElement[iType])(0, 0);
            const CFuint locID1 = (*m_faceNodePnElement[iType])(1, 0);
            const CFuint locID2 = (*m_faceNodePnElement[iType])(2, 0);

            const CFuint ID0 = (*cellNodes)(elemID,locID0);
            const CFuint ID1 = (*cellNodes)(elemID,locID1);
            const CFuint ID2 = (*cellNodes)(elemID,locID2);

            ///Only node with local index 9 is an inner node:
            const CFuint updateNodeID = (*cellNodes)(elemID, 9);

            if (!isUpdated[updateNodeID-oldNbNodes]) {

              const CFreal onethird = 1.0/3.0;
              nodeCoord[XX] = onethird * ((*nodes[ID0])[XX] + (*nodes[ID1])[XX] + (*nodes[ID2])[XX]);
              nodeCoord[YY] = onethird * ((*nodes[ID0])[YY] + (*nodes[ID1])[YY] + (*nodes[ID2])[YY]);

              getCFmeshData().createNode(updateNodeID,nodeCoord);

              isUpdated[updateNodeID-oldNbNodes] = true;
              nodes[updateNodeID]->setGlobalID(updateNodeID);
              nodes[updateNodeID]->setIsOnMesh(true);
              nodes[updateNodeID]->setIsOwnedByState(false);
           }
           break;

          }

         default: { }
        };//switch

      } //if geometric polyorder == 3

  } //Loop over elements

  } //Loop over element types


//     CFout << "Listing nodal coordinates ... after recreating the nodal data\n";
//           for(CFuint iNode = 0; iNode < newNbNodes; ++iNode) {
//      const Node& currNode = *nodes[iNode];
//      CFout << "Node " << nodes[iNode]->getGlobalID() << "\t[" << currNode[XX] << "," << currNode[YY] << "] is on mesh: "
//                       << nodes[iNode]->isOnMesh() << "\t and is owned by state: " << nodes[iNode]->isOwnedByState() << "\n";
//           }


  }

//////////////////////////////////////////////////////////////////////////////


void RDS_HighOrderMeshUpdater::recreateStates()
{
  CFAUTOTRACE;

  CFout << "\n-----------------------------------------------------\n\n";
  CFout << "\tExecuting recreateStates()\n";
  CFout << "\n-----------------------------------------------------\n";



  Common::SafePtr<MeshData::ConnTable> cellStates = MeshDataStack::getActive()->getConnectivity("cellStates_InnerCells");

  //   DataHandle < Framework::State*, Framework::GLOBAL > states = getCFmeshData().getStatesHandle();
  DataHandle < Framework::State*, Framework::GLOBAL > states =
    MeshDataStack::getActive()->getStateDataSocketSink().getDataHandle();
  const CFuint newNbStates = _totalNewNbStates;
  const CFuint oldNbStates = states.size();
  CFout << "Updating the nb of states from: " << oldNbStates << "  to  " << newNbStates <<"\n";

  Common::SafePtr<MeshData::ConnTable> cellNodes =
    MeshDataStack::getActive()->getConnectivity("cellNodes_InnerCells");

  DataHandle < Framework::Node*, Framework::GLOBAL > nodes =
    MeshDataStack::getActive()->getNodeDataSocketSink().getDataHandle();
  // DataHandle < Framework::Node*, Framework::GLOBAL > nodes = getCFmeshData().getNodesHandle();

  std::vector<State*> oldStates(oldNbStates);
  // backup the existing states
  for (CFuint i = 0; i < oldNbStates; ++i)
  {
    oldStates[i] = states[i];
  }

  // Resize the datahandle for the states
  states.resize(newNbStates);

  // restore P1 states from backup
  /// We suppose that P1 nodes and P1 states have the same global IDs and the same coordinates!
  for (CFuint i = 0; i < oldNbStates; ++i)
  {
    states[i] = oldStates[i];
    states[i]->setSpaceCoordinates(nodes[i]);
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


   ///Finally, the new states may sit 'on top' of some nodes - set the state coordinates
   ///Make sure nodes know they are owned by some states

  //Number of states that need to be modified: all non-vertex states:
  const CFuint nbPnStates = newNbStates - oldNbStates;


  //Indexes of Pn states whose coordinates have been updated:
  vector<bool> isUpdated(nbPnStates);
  for(CFuint iState = 0; iState < nbPnStates; ++iState)  isUpdated[iState] = false;


  SafePtr< std::vector< ElementTypeData > > elementType = getCFmeshData().getElementTypeData();
  std::vector< ElementTypeData >::iterator type_itr;
  const CFuint nbElemTypes = elementType->size();

  ///@todo: check dimension. If dim = 3, loop over element edges first and set nodal coordinates
  ///       then set coordinates of nodes that are inside faces or inside elements

  RealVector nodeCoord(2);
  Node* newNode;
  CFuint nodeSkip, stateSkip;

  type_itr = elementType->begin();

  for (CFuint iType = 0; iType < nbElemTypes; ++iType, ++type_itr) {

  // get element shape
  const CFGeoShape::Type elmtShape = type_itr->getGeoShape();

  // get start index of this element type in global element list
      CFuint globalIdx = (*elementType)[iType].getStartIdx();
  const CFuint nbElemPerType = (*elementType)[iType].getNbElems();

    // loop over elements of this type and print connectivity info
    for (CFuint iElem = 0; iElem < nbElemPerType; ++iElem, ++globalIdx) {

        const CFuint elemID = globalIdx;
        const CFuint nbFaces = LocalConnectionData::getInstance().getNbFacesInShape
         ((*elementType)[iType].getGeoShape());

//               CFout << "\n\n Element " << elemID << "\n";

  for(CFuint iFace = 0; iFace < nbFaces; ++iFace) {
          const CFuint nbPnNodesPerFace = m_faceNodePnElement[iType]->nbCols(iFace);
                const CFuint nbPnStatesPerFace = m_faceStatePnElement[iType]->nbCols(iFace);

                //The coordinates of the first and last nodes of the face
                const CFuint localVertexNodeID1 = (*m_faceNodePnElement[iType])(iFace,0);
                const CFuint localVertexNodeID2 = (*m_faceNodePnElement[iType])(iFace,1);

                const CFuint vertexNodeID1 = (*cellNodes)(elemID,localVertexNodeID1);
                const CFuint vertexNodeID2 = (*cellNodes)(elemID,localVertexNodeID2);

                const CFreal dx = ((*nodes[vertexNodeID2])[XX] - (*nodes[vertexNodeID1])[XX])/(nbPnStatesPerFace-1);
                const CFreal dy = ((*nodes[vertexNodeID2])[YY] - (*nodes[vertexNodeID1])[YY])/(nbPnStatesPerFace-1);

                     //Do some nodes and states have the same location?
                     //(i.e for example every second node has the same location as some state on the face)
                     const bool isMultiple = ( (nbPnStatesPerFace-1) % (nbPnNodesPerFace-1) == 0) || \
                                            ( (nbPnNodesPerFace-1) % (nbPnStatesPerFace-1) == 0);


                     if(isMultiple) {

                      if(nbPnNodesPerFace < nbPnStatesPerFace) {
                              nodeSkip = ((nbPnStatesPerFace-1) / (nbPnNodesPerFace-1));
                              stateSkip = 1; }
                      else if (nbPnNodesPerFace > nbPnStatesPerFace) {
                              nodeSkip = ((nbPnNodesPerFace-1) / (nbPnStatesPerFace-1));
                              stateSkip = 1; }
                      else {
                              nodeSkip = 1;
                              stateSkip = 1;
                            }

                          //Loop over states of the face, skip the first and last states (they are P1 states - updated already)
                         for(CFuint iState = 1; iState < (nbPnStatesPerFace-1); ++iState) {


                          //Note that the local connectivity table for each edge in 2D stores the first and
                          //the last local stateID first, all intermediate nodes follow. That's why we say iState+1 to retrieve
                          //the localStateID of the first state that does not coincide with the first/last node of the edge

                              const CFuint localStateID = (*m_faceStatePnElement[iType])(iFace,iState+1);
                              const CFuint stateID = (*cellStates)(elemID,localStateID);

                           if (!isUpdated[stateID-oldNbStates]) {


                              const CFuint iNode = iState / nodeSkip * stateSkip;
                              const CFuint localNodeID = (*m_faceNodePnElement[iType])(iFace,iNode+1);
                              const CFuint nodeID = (*cellNodes)(elemID,localNodeID);

                            if ((iState*stateSkip) % nodeSkip == 0) {
                                  //This state should have the same coordinates as an existing node in the mesh
                                  //Let this state use the coordinates of corresponding node
//                                CFout << "Nodes:  [" << vertexNodeID1 << "]  | " << nodeID << " |  [" << vertexNodeID2 << "]\n";
//                                CFout << "States: [" << vertexNodeID1 << "]  | " << stateID << " |  [" << vertexNodeID2 << "]\n\n";

                              states[stateID]->setSpaceCoordinates(nodes[nodeID]);
                              isUpdated[stateID-oldNbStates] = true;
                            }
                            //Else create a new node owned by the state
                            else {
//                               CFout << "StateID = " << stateID << " , NodeID = " << nodeID << "\n";
                              nodeCoord[XX] = (*nodes[vertexNodeID1])[XX] + iState*dx;
                              nodeCoord[YY] = (*nodes[vertexNodeID1])[YY] + iState*dy;

                              newNode = new Node(nodeCoord,false);

                              const CFuint localStateID = (*m_faceStatePnElement[iType])(iFace,iState+1);
                              const CFuint stateID = (*cellStates)(elemID,localStateID);
                              states[stateID]->setSpaceCoordinates(newNode);
                              newNode->setIsOwnedByState(true); //This node is not a meshpoint
                              newNode->setIsOnMesh(false);
                              isUpdated[stateID-oldNbStates] = true;

                            }

                           }//Check if the state is updated

                           } //Loop over states

                     }//if(isMultiple)

                     ///If no nodal and state coordinates coincide, all states on this face have
                     ///to have their own nodes that are not in the mesh
                     else {

//                              CFout << "Create own coordinates\n";
                             for(CFuint iState = 1; iState < (nbPnStatesPerFace-1); ++iState) {

                              nodeCoord[XX] = (*nodes[vertexNodeID1])[XX] + iState*dx;
                              nodeCoord[YY] = (*nodes[vertexNodeID1])[YY] + iState*dy;

                              newNode = new Node(nodeCoord,false);

                              const CFuint localStateID = (*m_faceStatePnElement[iType])(iFace,iState+1);
                              const CFuint stateID = (*cellStates)(elemID,localStateID);
                              states[stateID]->setSpaceCoordinates(newNode);
                              newNode->setIsOwnedByState(true); //This node is not a meshpoint
                              newNode->setIsOnMesh(false);
                              isUpdated[stateID-oldNbStates] = true;

                            }
                     }


        } //Loop over faces

      //Insert inner nodes, at least in case of P3 triangles:

      if(m_newSolPolyOrder_int == 3) {
        switch(elmtShape) {

          case CFGeoShape::TRIAG: {

            ///Only node with local index 9 is an inner state:
            const CFuint updateStateID = (*cellStates)(elemID, 9);

            if (!isUpdated[updateStateID-oldNbStates]) {

            //This state coincides with certain node
            //The state will have the coordinates given by its node
            if(m_newGeoPolyOrder_int ==3) {

              const CFuint updateNodeID = (*cellNodes)(elemID,9);
              nodeCoord[XX] = (*nodes[updateNodeID])[XX];
              nodeCoord[YY] = (*nodes[updateNodeID])[YY];

              states[updateStateID]->setSpaceCoordinates(nodes[updateNodeID]);

            }
            else { //let's create a node for this state

              const CFuint locID0 = (*m_faceNodePnElement[iType])(0, 0);
              const CFuint locID1 = (*m_faceNodePnElement[iType])(1, 0);
              const CFuint locID2 = (*m_faceNodePnElement[iType])(2, 0);

              const CFuint ID0 = (*cellNodes)(elemID,locID0);
              const CFuint ID1 = (*cellNodes)(elemID,locID1);
              const CFuint ID2 = (*cellNodes)(elemID,locID2);

              const CFreal onethird = 1.0/3.0;
              nodeCoord[XX] = onethird * ((*nodes[ID0])[XX] + (*nodes[ID1])[XX] + (*nodes[ID2])[XX]);
              nodeCoord[YY] = onethird * ((*nodes[ID0])[YY] + (*nodes[ID1])[YY] + (*nodes[ID2])[YY]);

              newNode = new Node(nodeCoord,false);

              states[updateStateID]->setSpaceCoordinates(newNode);
              newNode->setIsOwnedByState(true); //This node is not a meshpoint
              newNode->setIsOnMesh(false);

            } //creating a new node for this state
          } //if the state is already updated

          ///Finally, mark this state as updated and set its global id:
          isUpdated[updateStateID-oldNbStates] = true;
          states[updateStateID]->setGlobalID(updateStateID);


           }
           break;


         default: { }
        };//switch

      } //if geometric polyorder == 3



    } //Loop over elements


  } //Loop over element types


//Now loop over all states and print info

//     CFout << "\n\nPRINTING INFORMATION ABOUT STATES ...\n\n";
//     CFout << "Number of states = " << newNbStates << "\n";
//
//           for(CFuint iState = 0; iState < newNbStates; ++iState) {
//                 CFout << "State " << states[iState]->getGlobalID() << " has coordinates " << states[iState]->getCoordinates() << "\t";
//                 CFout << "Node is on mesh: " << nodes[iState]->isOnMesh() << "\n\n";
//
//           }
//
//     CF_DEBUG_EXIT;

}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
