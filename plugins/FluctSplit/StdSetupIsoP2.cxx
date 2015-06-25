#include "Common/NoSuchValueException.hh"
#include "Common/CFMultiMap.hh"

// #include "MathTools/RealMatrix.hh"

#include "Framework/MeshData.hh"
#include "Framework/SubSystemStatus.hh"
#include "Framework/ComputeNormals.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/MapGeoToTrsAndIdx.hh"
#include "Framework/NamespaceSwitcher.hh"

#include "FluctSplit/InwardNormalsData.hh"
#include "FluctSplit/ComputeInwardNormals.hh"
#include "FluctSplit/FluctSplit.hh"
#include "FluctSplit/StdSetupIsoP2.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<StdSetupIsoP2, FluctuationSplitData, FluctSplitModule> stdSetupIsoP2Provider("StdSetupIsoP2");

//////////////////////////////////////////////////////////////////////////////

StdSetupIsoP2::StdSetupIsoP2(const std::string& name) :
  FluctuationSplitCom(name),
  m_nodeIdToStateId(),
  m_dynamicSockets(),
  socket_states("states"),
  socket_nodes("nodes"),
  socket_normals("normals"),
  socket_nstatesProxy("nstatesProxy"),
  socket_isUpdated("isUpdated"),
  socket_isBState("isBState"),
  socket_discardTimeJacob("discardTimeJacob"),
  socket_normalsData("normalsData"),
  socket_faceNeighCell("faceNeighCell"),
  socket_tempSize("tempSize"),
  socket_volumes("volumes")
{ }

//////////////////////////////////////////////////////////////////////////////

StdSetupIsoP2::~StdSetupIsoP2()
{
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSource> >
StdSetupIsoP2::providesSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSource> > result;

  result.push_back(&socket_normals);
  result.push_back(&socket_nstatesProxy);
  result.push_back(&socket_isUpdated);
  result.push_back(&socket_isBState);
  result.push_back(&socket_discardTimeJacob);
  result.push_back(&socket_normalsData);
  result.push_back(&socket_faceNeighCell);
  result.push_back(&socket_tempSize);
  result.push_back(&socket_volumes);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
StdSetupIsoP2::needsSockets()
{

  std::vector<Common::SafePtr<BaseDataSocketSink> > result = m_dynamicSockets.getAllSinkSockets();

  result.push_back(&socket_states);
  result.push_back(&socket_nodes);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

void StdSetupIsoP2::configure ( Config::ConfigArgs& args )
{
  CFAUTOTRACE;

  FluctuationSplitCom::configure(args);

  //Loop over the TRS's and add the "TRSName" + "-boundaryNormals" datasocketsink to the m_dynamicSockets
  const std::string name = getMethodData().getNamespace();
  Common::SafePtr<Namespace> nsp = NamespaceSwitcher::getInstance
    (SubSystemStatusStack::getCurrentName()).getNamespace(name);
  Common::SafePtr<MeshData> meshData = MeshDataStack::getInstance().getEntryByNamespace(nsp);

  vector<std::string> trsList = meshData->getTRSNameList();
  const CFint nbTRSs = trsList.size();

// CF_DEBUG_OBJ(nbTRSs);
// CF_DEBUG_OBJ(trsList[0]);

  for (CFint iTRS = 0; iTRS < nbTRSs; ++iTRS) {
    const std::string trsName = trsList[iTRS];

    if (trsName != "PartitionFaces" &&
        trsName != "InnerCells" &&
        trsName != "InnerFaces")
    {
      const std::string socketName = trsName + "-boundaryNormals";
// CF_DEBUG_OBJ(socketName);
      const bool isEssential = false;
      m_dynamicSockets.createSocketSink<const CFreal*>(socketName,isEssential);
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void StdSetupIsoP2::execute()
{
  CFAUTOTRACE;

  Common::SafePtr<TopologicalRegionSet> cells =
    MeshDataStack::getActive()->getTrs("InnerCells");
  const CFuint nbCells = cells->getLocalNbGeoEnts();

  DataHandle<InwardNormalsData*> normals =
    socket_normals.getDataHandle();
  normals.resize(nbCells);

  computeNormalsData();

  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();
  DataHandle < Framework::Node*, Framework::GLOBAL > nodes  = socket_nodes.getDataHandle();

  const CFuint nbStates = states.size();
  const CFuint nbNodes = nodes.size();

  DataHandle<ProxyDofIterator<RealVector>*> nstatesProxy =
    socket_nstatesProxy.getDataHandle();
  nstatesProxy.resize(1);

  // set the mapping from nodeID to stateID
  m_nodeIdToStateId.resize(nbNodes);
  for (CFuint i = 0; i < nbStates; ++i)
  {
    const bool indexed = states[i]->getCoordinates().isIndexed();
    if(indexed)
    {
      const CFuint nodeID = states[i]->getCoordinates().getLocalID();
      cf_assert (nodeID < nbStates);
      m_nodeIdToStateId[nodeID] = i;
    }
  }
  nstatesProxy[0] =
    new DofDataHandleIterator<RealVector, State, GLOBAL>(states, &m_nodeIdToStateId);
  
  DataHandle<bool> isUpdated = socket_isUpdated.getDataHandle();
  isUpdated.resize(nbStates);
  isUpdated = false;

  DataHandle<bool> isBState = socket_isBState.getDataHandle();
  isBState.resize(nbStates);
  isBState = false;

  // the following DataHandle has to be filled in by the single implicit BCs
  DataHandle<vector<bool> > discardTimeJacob = socket_discardTimeJacob.getDataHandle();
  discardTimeJacob.resize(nbStates);

  // set the neighbor cell for each boundary face
  setFaceNeighCell();

  // flag the states that are on the boundary
  flagBoundaryStates();

  // store the boundaryNormals
  storeBoundaryNormals();

  // compute volumes
  computeCellVolumes();
}

//////////////////////////////////////////////////////////////////////////////

void StdSetupIsoP2::storeBoundaryNormals()
{

  //Loop over the TRS's and add the "TRSName" + "-boundaryNormals" datasocketsink to the m_dynamicSockets
  vector<Common::SafePtr<TopologicalRegionSet> > trs = MeshDataStack::getActive()->getTrsList();
  const CFint nbTRSs = trs.size();

  for (CFint iTRS = 0; iTRS < nbTRSs; ++iTRS)
  {
    SafePtr<TopologicalRegionSet> currTrs = trs[iTRS];

    if (!currTrs->hasTag("inner"))
    {
// CF_DEBUG_POINT;
      const std::string socketName = currTrs->getName() + "-boundaryNormals";
      Common::SafePtr<DataSocketSink< const CFreal*> > boundaryNormals =
	m_dynamicSockets.getSocketSink<const CFreal*>(socketName);
// CF_DEBUG_POINT;

      if(boundaryNormals->isConnected())
      {
        Common::SafePtr<GeometricEntityPool<StdTrsGeoBuilder> >
          geoBuilder = getMethodData().getStdTrsGeoBuilder();

        StdTrsGeoBuilder::GeoData& geoData = geoBuilder->getDataGE();

        geoData.trs = currTrs;

        const CFuint nbFaces = currTrs->getLocalNbGeoEnts();

        // get the face neighbor cell data
        DataHandle<Common::Trio<CFuint, CFuint, Common::SafePtr<Framework::TopologicalRegionSet> > >
          faceNeighCell = socket_faceNeighCell.getDataHandle();

        // get the normals
        DataHandle< InwardNormalsData*> normals = socket_normals.getDataHandle();

        // get the boundaryNormals
        DataHandle< const CFreal*> boundaryNormalsHandle = boundaryNormals->getDataHandle();

        //resize the boundaryNormals datahandle
        boundaryNormalsHandle.resize(nbFaces);

        for (CFuint iFace = 0; iFace < nbFaces; ++iFace)
        {
          geoData.idx = iFace;

          GeometricEntity *const currFace = geoBuilder->buildGE();

          CFuint faceID = currFace->getID();
          const CFuint cellTrsID = faceNeighCell[faceID].first;
          const CFuint iFaceLocal = faceNeighCell[faceID].second;
          Common::SafePtr<TopologicalRegionSet> cellTrs = faceNeighCell[faceID].third;
          const CFuint cellLocalID = cellTrs->getLocalGeoID(cellTrsID);

          boundaryNormalsHandle[iFace] = normals[cellLocalID]->getFaceNormalPtr(iFaceLocal);

          // release the face
          geoBuilder->releaseGE();
        }
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void StdSetupIsoP2::computeNormalsData()
{
  CFAUTOTRACE;
// CF_DEBUG_POINT;
  DataHandle<CFreal> normalsData = socket_normalsData.getDataHandle();

  SafePtr<vector<ElementTypeData> > elemTypes =
    MeshDataStack::getActive()->getElementTypeData();

  const CFuint nbElemTypes = elemTypes->size();

// CF_DEBUG_OBJ(nbElemTypes);

  // resize the map between state ID and opposite face ID in the
  // local connectivity (this is useful if you have triangles or
  // tetras)
  InwardNormalsData::resizeStateToFaceIDMap(nbElemTypes);

  const CFuint dim = PhysicalModelStack::getActive()->getDim();
  CFuint sizeFaceData = 0;
  CFuint sizeNodeData = 0;
  CFuint sizeFaceArea = 0;
  CFuint sizeNodeArea = 0;

  /// @note this is hardcoded to P2P2 curved elements
  const CFuint nb_subfaces_elem = 9;

  // compute the size of the normals data storages
   for (CFuint iType = 0; iType < nbElemTypes; ++iType)
   {
    const ElementTypeData& eData = (*elemTypes)[iType];

	const CFuint nbelems = eData.getNbElems();
// unused //	const CFuint nbfaces = eData.getNbFaces();
// 	CF_DEBUG_OBJ(nbelems);
// 	CF_DEBUG_OBJ(nbfaces);
	//Check if the element is a P2-P2 triangle:
	if (eData.getGeoShape() == CFGeoShape::TRIAG && eData.getGeoOrder() == CFPolyOrder::ORDER2)
        {
		sizeFaceData += nbelems*nb_subfaces_elem*dim;
		sizeFaceArea += nbelems*nb_subfaces_elem;
	}
	else
	{
    	sizeFaceData += eData.getNbElems()*eData.getNbFaces()*dim; //Normals to faces
    	sizeFaceArea += eData.getNbElems()*eData.getNbFaces();     //One CFreal per face - area size
	}

	if (eData.getGeoShape() != CFGeoShape::TRIAG && eData.getGeoShape() != CFGeoShape::TETRA) {
	sizeNodeData += eData.getNbElems()*eData.getNbNodes()*dim;
	sizeNodeArea += eData.getNbElems()*eData.getNbNodes();
	}
  }

  // resize normals data storages
  normalsData.resize(sizeFaceData + sizeFaceArea +
                     sizeNodeData + sizeNodeArea);

  // this storage keeps track of the temporary sizes of the
  // normals related data while they are being filled in
  // when push_back will be available in GrowArray, this will
  // not be needed anymore
  DataHandle<CFuint> tempSize = socket_tempSize.getDataHandle();

  tempSize.resize(1);
  tempSize= 0;

  for (CFuint iType = 0; iType < nbElemTypes; ++iType)
  {

    const CFuint geoOrder = (*elemTypes)[iType].getGeoOrder();
    const std::string elemName = (*elemTypes)[iType].getShape() + CFPolyOrder::Convert::to_str( geoOrder );

    SelfRegistPtr<ComputeNormals> computeFaceNormals =
      Environment::Factory<ComputeNormals>::getInstance().getProvider
      ("Inward" + elemName)->create();

    SelfRegistPtr<ComputeInwardNormals> dPtr = computeFaceNormals.d_castTo<ComputeInwardNormals>();

    // set the data sockets
    DataSocketSink< InwardNormalsData*> socket_normalsSink(socket_normals);
    DataSocketSink< CFreal> socket_normalsDataSink(socket_normalsData);
    DataSocketSink< CFuint> socket_tempSizeSink(socket_tempSize);

    dPtr->setNormalsSockets(&socket_normalsSink);
    dPtr->setNormalsDataSockets(&socket_normalsDataSink);
    dPtr->setTempSizeSockets(&socket_tempSizeSink);

    const CFuint firstElem = (*elemTypes)[iType].getStartIdx();
    const CFuint lastElem  = (*elemTypes)[iType].getEndIdx();

    (*dPtr)(firstElem, lastElem, iType);

  }
}

//////////////////////////////////////////////////////////////////////////////

void StdSetupIsoP2::flagBoundaryStates()
{
  CFAUTOTRACE;

  DataHandle< bool> isBoundaryState = socket_isBState.getDataHandle();

  // set the list of faces
  vector<Common::SafePtr<TopologicalRegionSet> > trs = MeshDataStack::getActive()->getTrsList();
  vector< Common::SafePtr<TopologicalRegionSet> >::iterator itrs;
  for (itrs = trs.begin(); itrs != trs.end(); ++itrs) {
    if ((*itrs)->getName() != "InnerCells") {
      Common::SafePtr< vector<CFuint> > const statesIdxInTRS =
        (*itrs)->getStatesInTrs();

      for (CFuint iState = 0; iState < statesIdxInTRS->size(); ++iState) {
        isBoundaryState[(*statesIdxInTRS)[iState]] = true;
//         CFout << (*statesIdxInTRS)[iState] << " is boundary state\n";
      }
    }
  }
 }

//////////////////////////////////////////////////////////////////////////////

void StdSetupIsoP2::setFaceNeighCell()
{
  CFAUTOTRACE;

  vector<Common::SafePtr<TopologicalRegionSet> > trs =
    MeshDataStack::getActive()->getTrsList();

  DataHandle < Framework::Node*, Framework::GLOBAL > nodes = socket_nodes.getDataHandle();

  const CFuint totNbNodes = nodes.size();
  vector<bool> isBNode(totNbNodes, false);

  //count the number of boundary faces
  CFuint nbBFaces = 0;
  // create mapping between boundary nodes and faces
  CFMultiMap<CFuint, CFuint> mapBNode2BFace;

  for (CFuint i = 0; i < trs.size(); ++i) {
    SafePtr<TopologicalRegionSet> currTrs = trs[i];
    if (currTrs->getName() != "InnerCells") {
      const CFuint nbTrsFaces = currTrs->getLocalNbGeoEnts();
      for (CFuint iFace = 0; iFace < nbTrsFaces; ++iFace) {
        const CFuint faceID = currTrs->getLocalGeoID(iFace);
        const CFuint nbNodesInFace = currTrs->getNbNodesInGeo(iFace);
        for (CFuint iNode = 0; iNode < nbNodesInFace; ++iNode) {
          const CFuint nodeID = currTrs->getNodeID(iFace,iNode);
//           CFout << "Inserting node " << nodeID << ", face " << faceID << "\n";
          mapBNode2BFace.insert(nodeID, faceID);
          // flag the boundary nodes
          cf_assert (nodeID < isBNode.size());
          isBNode[nodeID] = true;
        }
//           CFout << "Node " << nodeID << " is on boundary\n";
      }
      nbBFaces += nbTrsFaces;
    }
  }
  mapBNode2BFace.sortKeys();

  Common::SafePtr<MapGeoToTrsAndIdx> mapGeoToTrs =
    MeshDataStack::getActive()->getMapGeoToTrs("MapFacesToTrs");

  // face neighbor cell
  DataHandle<Common::Trio<CFuint, CFuint, Common::SafePtr<Framework::TopologicalRegionSet> > >
    faceNeighCell = socket_faceNeighCell.getDataHandle();
  faceNeighCell.resize(nbBFaces);

  std::vector< SafePtr<TopologicalRegionSet> > cellList =
    MeshDataStack::getActive()->getFilteredTrsList("inner");

  CFuint countF = 0;

  for(CFuint iTrs=0; iTrs < cellList.size(); ++iTrs){

    SafePtr<TopologicalRegionSet> cells = cellList[iTrs];

    SafePtr<vector<ElementTypeData> > elementType =
      MeshDataStack::getActive()->getElementTypeData(cells->getName());

    const CFuint nbElemTypes = elementType->size();

    // local connectivity face-node for each element type
    vector<Table<CFuint>*> faceNodeElement(nbElemTypes);
    for (CFuint iType = 0; iType < nbElemTypes; ++iType) {
      faceNodeElement[iType] = LocalConnectionData::getInstance().getFaceDofLocal
        ((*elementType)[iType].getGeoShape(),
        static_cast<CFPolyOrder::Type>((*elementType)[iType].getGeoOrder()),
        NODE,
        CFPolyForm::LAGRANGE);
    }

    // set the face shapes per element type
    vector< vector<CFGeoShape::Type> > faceShapesPerElemType(nbElemTypes);
    LocalConnectionData::getInstance().setFaceShapesPerElemType(faceShapesPerElemType);

    // state list
    DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();

    // atomic number to indicate the maximum possible number
    // of nodes in a face (this allow to avoid frequent reallocations
    // of the vector nodesInFace)
    const CFuint maxNbNodesInFace = 100;
    vector<CFuint> nodesInFace(maxNbNodesInFace);

    // loop over the elements and construct faceIDs
    CFuint elemID = 0;
    typedef CFMultiMap<CFuint, CFuint>::MapIterator MapIterator;

    // loop over the types
    for (CFuint iType = 0; iType < nbElemTypes; ++iType)
    {
      /// @todo for now all geoents have same geometric and solution polyorder
      const CFuint nbElemFaces = faceShapesPerElemType[iType].size();

      // loop over the elements of this type
      const CFuint nbElemPerType = (*elementType)[iType].getNbElems();

      for (CFuint iElem = 0; iElem < nbElemPerType; ++iElem, ++elemID)
      {

        // loop over the faces in the current element
        for (CFuint iFace = 0; iFace < nbElemFaces; ++iFace)
         {
            const CFuint nbNodesPerFace = faceNodeElement[iType]->nbCols(iFace);
//   CF_DEBUG_OBJ( nbNodesPerFace );
            // construct sets of nodes that make the corresponding face
            // in this element
            CFuint countFaceNodesOnBoundary = 0;
            for (CFuint iNode = 0; iNode < nbNodesPerFace; ++iNode) {
              const CFuint localNodeID = (*faceNodeElement[iType])(iFace, iNode);
              // cout << "localNodeID = "<< localNodeID << endl;
              const CFuint nodeID = cells->getNodeID(elemID, localNodeID);
              nodesInFace[iNode] = nodeID;
              if (isBNode[nodeID]) {
                countFaceNodesOnBoundary++;
              }
              else break;
            }
//   CF_DEBUG_OBJ( countFaceNodesOnBoundary );

            // face is on the boundary if ALL its nodes are on the boundary
            if (countFaceNodesOnBoundary == nbNodesPerFace)
            {
// CF_DEBUG_POINT;
              // consider the first node belonging to the current face
              // check if you find a face ID shared between all the other
              // nodes building a face
              const CFuint nodeID = nodesInFace[0];
              
	      bool fo = false;
	      pair<MapIterator, MapIterator> faces = mapBNode2BFace.find(nodeID, fo);
	      cf_assert(fo);
	      
              // loop over all the faceIDs referenced by the first node to see if
              // all the other nodes reference the same face
              bool faceFound = false;
              CFuint currFaceID = 0;
              for (MapIterator faceInMapItr = faces.first;
                  faceInMapItr != faces.second && (faceFound == false);
                  ++faceInMapItr) {

                currFaceID = faceInMapItr->second;
                SafePtr<TopologicalRegionSet> faceTrs =
                  mapGeoToTrs->getTrs(currFaceID);

                const CFuint faceIdx = mapGeoToTrs->getIdxInTrs(currFaceID);
                const CFuint nbNodesInCurrFace = faceTrs->getNbNodesInGeo(faceIdx);

                CFuint countNodes = 0;
                // first check if the number of nodes of the face is equal
                if (nbNodesInCurrFace == nbNodesPerFace) {
                  for (CFuint iNode = 0; iNode < nbNodesPerFace; ++iNode) {
                    const CFuint newNodeID = nodesInFace[iNode];
                    for (CFuint jNode = 0; jNode < nbNodesPerFace; ++jNode) {
                      if (faceTrs->getNodeID(faceIdx,jNode) == newNodeID) {
                        countNodes++;
                        break;
                      }
                    }
                  }
                }

                if (countNodes == nbNodesPerFace) {
                  // set the neighbor cell ID
                  faceNeighCell[currFaceID].first = elemID;
                  // set the local face ID in the corresponding cell
                  faceNeighCell[currFaceID].second = iFace;
                  //set the inner trs to which the cell belongs
                  faceNeighCell[currFaceID].third = cells;

                  faceFound = true;
                  countF++;
                }
              }

              // For "corner" cells, it could happen a face whose all
              // nodes are on the boundary but belong to different TRs
              // (possibly even within the same TRS).
              // In this case the face is NOT a boundary face and so
              // at this point faceFound would remain "false"
              if (!faceFound) {
                /// @todo a more rigorous check could be made to verify
                /// that this exceptional situation is occurring
                CFLogDebugMin("corner inner face with all boundary nodes found\n");
                CFLogDebugMin("face nodes are: ");
                for (CFuint in = 0; in < nbNodesPerFace; ++in) {
                  CFLogDebugMin(nodesInFace[in] << " ");
                }
                CFLogDebugMin("\n");
              }
            }

         } // end loop faces
       } // end loop elements
    }  // loop element types
  }  //end of loop over inner trs

  if (countF != nbBFaces)
  {
    throw Common::NoSuchValueException
      (FromHere(), "Not all the face neighbor cell have been detected correctly!!");
  }
}

//////////////////////////////////////////////////////////////////////////////

void StdSetupIsoP2::computeCellVolumes()
{
  SafePtr<TopologicalRegionSet> cells =
    MeshDataStack::getActive()->getTrs("InnerCells");

  const CFuint nbCells = cells->getLocalNbGeoEnts();
  DataHandle<CFreal> volumes = socket_volumes.getDataHandle();
  volumes.resize(nbCells);

  Common::SafePtr<GeometricEntityPool<StdTrsGeoBuilder> >
    geoBuilder = getMethodData().getStdTrsGeoBuilder();

  StdTrsGeoBuilder::GeoData& geoData = geoBuilder->getDataGE();
  geoData.trs = cells;

  for(CFuint iGeoEnt = 0; iGeoEnt < nbCells; ++iGeoEnt) {
    // build the GeometricEntity
    geoData.idx = iGeoEnt;
    GeometricEntity& cell = *geoBuilder->buildGE();

    volumes[iGeoEnt] = cell.computeVolume();

    //release the GeometricEntity
    geoBuilder->releaseGE();
  }
}

//////////////////////////////////////////////////////////////////////////////

 } // namespace FluctSplit



} // namespace COOLFluiD
