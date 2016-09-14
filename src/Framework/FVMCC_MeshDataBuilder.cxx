#include <set>
#include <numeric>

#include "Common/SwapEmpty.hh"

#include "Common/CFLog.hh"
#include "Environment/ObjectProvider.hh"

#include "Framework/MeshData.hh"
#include "Framework/LocalConnectionData.hh"
#include "Framework/BaseGeometricEntityProvider.hh"
#include "Framework/Face.hh"
#include "Framework/Cell.hh"
#include "Framework/PhysicalModel.hh"
#include "Framework/SetElementStateCoord.hh"
#include "Framework/BadFormatException.hh"
#include "Framework/MapGeoToTrsAndIdx.hh"
#include "Framework/Framework.hh"
#include "Framework/FVMCC_MeshDataBuilder.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::MathTools;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<FVMCC_MeshDataBuilder,
			    MeshDataBuilder,
			    FrameworkLib,
			    1>
connFVMCCProvider("FVMCC");

//////////////////////////////////////////////////////////////////////////////

FVMCC_MeshDataBuilder::FVMCC_MeshDataBuilder(const std::string& name) :
   MeshDataBuilder(name),
   m_nbFaces(0),
   m_inGeoTypes(CFNULL),
   m_inLocalGeoIDs(CFNULL),
   m_bGeoType(),
   m_bLocalGeoIDs(),
   m_faceNodeElement(),
   m_nbFacesPerElem(),
   m_isBFace(),
   m_geoTypeIDs(),
   m_nbInFacesNodes(),
   m_nbBFacesNodes(),
   m_bFaceStateID(),
   m_bFaceNodes(CFNULL),
   m_mapStateIdToFaceIdx(),
   m_isPartitionFace()
 {
 }

//////////////////////////////////////////////////////////////////////////////

FVMCC_MeshDataBuilder::~FVMCC_MeshDataBuilder()
{
}

//////////////////////////////////////////////////////////////////////////////

void FVMCC_MeshDataBuilder::releaseMemory()
{
  MeshDataBuilder::releaseMemory();
  
  SwapEmpty(m_bGeoType);
  SwapEmpty(m_bLocalGeoIDs);
  SwapEmpty(m_faceNodeElement);
  m_nbFacesPerElem.resize(0);
  SwapEmpty(m_isBFace);
  SwapEmpty(m_geoTypeIDs);
  m_nbInFacesNodes.resize(0);
  m_nbBFacesNodes.resize(0);
  SwapEmpty(m_bFaceStateID);
  deletePtr(m_bFaceNodes);
  m_mapStateIdToFaceIdx.clear();
  m_isPartitionFace.resize(0);
}

//////////////////////////////////////////////////////////////////////////////

void FVMCC_MeshDataBuilder::createTopologicalRegionSets()
{
  CFAUTOTRACE;

  // first create the cells and the renumber them
  // as if it would be a cell-vertex mesh
  CFLog(VERBOSE, "FVMCC_MeshDataBuilder::createTopologicalRegionSets() => 1\n");
  createInnerCells();
 
  CFLog(VERBOSE, "FVMCC_MeshDataBuilder::createTopologicalRegionSets() => 2\n");
  renumberCells();
  
  CFLog(VERBOSE, "FVMCC_MeshDataBuilder::createTopologicalRegionSets() => 3\n");
  // put the coordinates in the states
  setCoordInCellStates();
  
  CFLog(VERBOSE, "FVMCC_MeshDataBuilder::createTopologicalRegionSets() => 4\n");
  // create the cell-face connectivity
  createCellFaces();
  
  CFLog(VERBOSE, "FVMCC_MeshDataBuilder::createTopologicalRegionSets() => 5\n");
  // create the inner faces TRS
  createInnerFacesTRS();
  
  CFLog(VERBOSE, "FVMCC_MeshDataBuilder::createTopologicalRegionSets() => 6\n");
  // create all the boundary faces TRSs
  createBoundaryFacesTRS();
  
  CFLog(VERBOSE, "FVMCC_MeshDataBuilder::createTopologicalRegionSets() => 7\n");
  // set the mappping allowing to get Face info from TRS
  // by geometric entity local (in the processor) ID
  setMapGeoToTrs();
  
  CFLog(VERBOSE, "FVMCC_MeshDataBuilder::createTopologicalRegionSets() => 8\n");
}
    
//////////////////////////////////////////////////////////////////////////////

void FVMCC_MeshDataBuilder::createCellFaces()
{
  CFAUTOTRACE;

  const CFuint nbElem      = getNbElements();
  const CFuint nbElemTypes = getNbElementTypes();

  // local connectivity face-node for each element type
  m_faceNodeElement.resize(nbElemTypes);

  SafePtr<vector<ElementTypeData> > elementType =
    getCFmeshData().getElementTypeData();

  cf_assert(nbElemTypes == elementType->size());

  for (CFuint iType = 0; iType < nbElemTypes; ++iType) {
    m_faceNodeElement[iType] = LocalConnectionData::getInstance().getFaceDofLocal
      ((*elementType)[iType].getGeoShape(),
       getGeometricPolyOrder(),
       NODE, CFPolyForm::LAGRANGE);
  }

  // array storing the number of faces per element
  m_nbFacesPerElem.resize(nbElem);
  m_nbFacesPerElem = 0;
  LocalConnectionData::getInstance().setNbFacesPerElement(m_nbFacesPerElem);

  // set the face shapes per element type
  vector< vector<CFGeoShape::Type> > faceShapesPerElemType(nbElemTypes);
  LocalConnectionData::getInstance().setFaceShapesPerElemType(faceShapesPerElemType);

  // create a table to store the connectivity element-face ID locally to this processor
  // and store the connectivity in MeshData
  ConnTable* cellFaces = new ConnTable(m_nbFacesPerElem);
  MeshDataStack::getActive()->storeConnectivity("cellFaces", cellFaces);

  const CFuint totNbNodes = MeshDataStack::getActive()->getNbNodes();

  // allocate a table mapping node-face ID
  vector < vector<CFuint> > mapNodeFace(totNbNodes);

  // atomic number to indicate the maximum possible number
  // of nodes in a face
  // allows to avoid frequent reallocations of the vector nodesInFace
  const CFuint maxNbNodesInFace = 100;
  vector<CFuint> nodesInFace(maxNbNodesInFace);

  const std::string faceProviderName = "Face";

  // number of boundary faces in TRS data read from mesh file
  // this DOESN'T include partition faces !!
  const CFuint nbBoundaryFaces = getNbBoundaryFaces();
  const CFuint sumNbFacesPerElem = m_nbFacesPerElem.sum();
  // max possible number of inner faces (>= actual value in SERIAL simulation)
  // const CFuint maxNbInnerFaces = (m_nbFacesPerElem.sum() - nbBoundaryFaces)/2;
  // the max number of total TRSs is always slightly overestimated
  // ONLY in serial run totalNbFaces = nbBoundaryFaces + maxNbInnerFaces
  const CFuint maxTotalNbFaces = sumNbFacesPerElem;

  CFLog(INFO, "FVMCC BoundaryFaces [" << nbBoundaryFaces << "]\n");
  CFLog(INFO, "FVMCC Max Total NbFaces [" << maxTotalNbFaces << "]\n");

  // the following arrays are oversized
  m_geoTypeIDs.resize(maxTotalNbFaces);
  m_isBFace.resize(maxTotalNbFaces);
  for (CFuint i = 0; i < maxTotalNbFaces; ++i) { m_isBFace[i] = true; }

  vector<CFuint> nbFaceNodes;
  nbFaceNodes.reserve(maxTotalNbFaces);

  // loop over the elements and construct faceIDs
  CFuint elemID = 0;
  CFuint countInFaces = 0;
  std::string providerName = "";

  // during the first big loop the following is done:
  // 1. set the element-faces connectivity
  // 2. set the geometric entity type IDs of each face
  // 3. select which are the boundary faces and which are internal ones

  // loop over the types
  for (CFuint iType = 0; iType < nbElemTypes; ++iType)
  {
    /// @todo for now all geoents have same geometric and solution polyorder
    const CFuint nbElemFaces = faceShapesPerElemType[iType].size();
    vector<CFuint> faceGeoTypeID(nbElemFaces);

    // create the GeometricEntityProviders corresponding to each
    // single face of this element type
    for (CFuint iFace = 0; iFace < nbElemFaces; ++iFace)
    {
      providerName = faceProviderName +
        makeGeomEntName(faceShapesPerElemType[iType][iFace],
         getGeometricPolyType(),
         getGeometricPolyOrder(),
         getSolutionPolyType(),
         getSolutionPolyOrder());

      CFLogDebugMin("FVMCC Face provider [" << providerName << "]\n");

      faceGeoTypeID[iFace] = m_mapGeoProviderNameToType.find(providerName);
    }

    // loop over the elements of this type
    const CFuint nbElemPerType = (*elementType)[iType].getNbElems();
    for (CFuint iElem = 0; iElem < nbElemPerType; ++iElem, ++elemID)
    {
      // loop over the faces in the current element
      const CFuint nbElemFaces = m_nbFacesPerElem[elemID];
      for (CFuint iFace = 0; iFace < nbElemFaces; ++iFace)
      {
	// construct sets of nodes that make the corresponding face in this element
	const CFuint nbNodesPerFace = m_faceNodeElement[iType]->nbCols(iFace);
	for (CFuint iNode = 0; iNode < nbNodesPerFace; ++iNode)
	{
	  const CFuint localNodeID = (*m_faceNodeElement[iType])(iFace, iNode); 
          nodesInFace[iNode] = getCFmeshData().getElementNode(elemID, localNodeID);
	}

        // consider the first node belonging to the current face
        // check if you find a face ID shared between all the other
        // nodes building a face
        const CFuint nodeID = nodesInFace[0];
	
        // number of faceIDs referencing node
        const CFuint nbFaceIDsRefNode = mapNodeFace[nodeID].size();
	// loop over all the faceIDs referenced by the first node to see if
        // all the other nodes reference the same face
        bool faceFound = false;
        CFuint countNodes = 0;
        CFuint currFaceID = 0;
        for (CFuint iFaceID = 0;
             iFaceID < nbFaceIDsRefNode && (faceFound == false);
             ++iFaceID)
        {
          currFaceID = mapNodeFace[nodeID][iFaceID];
          countNodes = 1;
	  
          if (nbNodesPerFace > 1) {
	    for (CFuint iNode = 1; iNode < nbNodesPerFace; ++iNode) {
	      const CFuint newNodeID = nodesInFace[iNode];
	      // number of faceIDs referencing node
	      const CFuint nbFaceIDsRefThisNode = mapNodeFace[newNodeID].size();
	      for (CFuint jFaceID = 0; jFaceID < nbFaceIDsRefThisNode; ++jFaceID) {
		if (mapNodeFace[newNodeID][jFaceID] == currFaceID) {
		  ++countNodes;
		  break;
		}
	      }
	      
	      if (countNodes == nbNodesPerFace)
	      {
		// the corresponding faceID already exists, meaning
		// that the face is an internal one, shared by two elements
		// here you set the second element (==state) neighbor of the face
		faceFound = true;
		(*cellFaces)(elemID, iFace) = currFaceID;
		
		// since it has two neighbor cells,
		// this face is surely NOT a boundary face
		m_isBFace[currFaceID] = false;
		
		// increment number of inner faces (they always have 2 states)
		countInFaces++;
		break;
	      }
	    }
	  }
	  else {
	    // this only applies to the 1D case 
	    cf_assert(nbNodesPerFace == 1);
	    
	    // the corresponding faceID already exists, meaning
	    // that the face is an internal one, shared by two elements
	    // here you set the second element (==state) neighbor of the face
	    faceFound = true;
	    (*cellFaces)(elemID, iFace) = currFaceID;
	    
	    // since it has two neighbor cells,
	    // this face is surely NOT a boundary face
	    m_isBFace[currFaceID] = false;
	    
	    // increment number of inner faces (they always have 2 states)
	    countInFaces++;
	  }
	}
	
        if (!faceFound || nbFaceIDsRefNode == 0) {
          // a new face has been found
          // add the ID of the new face in the corresponding nodes
          // referencing it
          for (CFuint i = 0; i < nbNodesPerFace; ++i) {
            const CFuint currNodeID = nodesInFace[i];
            mapNodeFace[currNodeID].push_back(m_nbFaces);
          }

          // store the geometric entity type for the current face
          m_geoTypeIDs[m_nbFaces] = faceGeoTypeID[iFace];

          (*cellFaces)(elemID, iFace) = m_nbFaces;
          nbFaceNodes.push_back(nbNodesPerFace);

          // increment the number of faces
          m_nbFaces++;
        }
      }
    }
  }

  cf_assert(m_nbFaces <= maxTotalNbFaces);
  cf_assert(countInFaces <= maxTotalNbFaces);

  const CFuint totalNbFaces = m_nbFaces;
  const CFuint nbInnerFaces = countInFaces;

  CFLog(INFO, "FVMCC Total nb faces [" << totalNbFaces << "]\n");
  CFLog(INFO, "FVMCC Inner nb faces [" << nbInnerFaces << "]\n");

  // total number of boundary + partition boundary faces
  const CFuint nbBPlusPartitionFaces = totalNbFaces - nbInnerFaces;
  CFLog(INFO, "FVMCC Boundary and Partition faces [" << nbBPlusPartitionFaces << "]\n");

  m_nbInFacesNodes.resize(nbInnerFaces);
  m_nbBFacesNodes.resize(nbBPlusPartitionFaces);

  // set the number of nodes in faces
  CFuint iBFace = 0;
  CFuint iInFace = 0;
  for (CFuint i = 0; i < nbFaceNodes.size(); ++i) {
    if (!m_isBFace[i]) {
      m_nbInFacesNodes[iInFace++] = nbFaceNodes[i];
    }
    else {
      m_nbBFacesNodes[iBFace++] = nbFaceNodes[i];
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void FVMCC_MeshDataBuilder::createInnerFacesTRS()
{
  CFAUTOTRACE;

  const CFuint nbInnerFaces = m_nbInFacesNodes.size();
  const CFuint nbBPlusPartitionFaces = m_nbBFacesNodes.size();
  m_bFaceStateID.resize(nbBPlusPartitionFaces);

  ConnTable* innerFaceNodes = new ConnTable(m_nbInFacesNodes);
  m_bFaceNodes = new ConnTable(m_nbBFacesNodes);

  // reset the patterns
  m_nbInFacesNodes = 2;
  ConnTable* innerFaceStates = new ConnTable(m_nbInFacesNodes);

  std::valarray<CFint> idxState(-1, m_nbFaces);
  
  if (m_inGeoTypes != CFNULL) {deletePtr(m_inGeoTypes);}
  if (m_inLocalGeoIDs != CFNULL) {deletePtr(m_inLocalGeoIDs);}
  
  m_inGeoTypes = new vector<CFuint>(nbInnerFaces);
  m_inLocalGeoIDs = new vector<CFuint>(nbInnerFaces);
  m_bGeoType.resize(nbBPlusPartitionFaces);
  m_bLocalGeoIDs.resize(nbBPlusPartitionFaces);

  SafePtr<ConnTable> cellFaces = MeshDataStack::getActive()->
    getConnectivity("cellFaces");

  SafePtr<vector<ElementTypeData> > elementType =
    getCFmeshData().getElementTypeData();

  const CFuint nbElemTypes = elementType->size();
  // reset the number of faces
  CFint countBFaces = -1;
  CFint countInFaces = -1;
  CFuint elemID = 0;

  // loop over the types
  for (CFuint iType = 0; iType < nbElemTypes; ++iType) {
    // loop over the elements of this type
    const CFuint nbElemPerType = (*elementType)[iType].getNbElems();
    for (CFuint iElem = 0; iElem < nbElemPerType; ++iElem, ++elemID) {
      const CFuint nbElemFaces = m_nbFacesPerElem[elemID];

      // loop over the faces in the current element
      for (CFuint iFace = 0; iFace < nbElemFaces; ++iFace)
      {
        const CFuint nbNodesPerFace = m_faceNodeElement[iType]->nbCols(iFace);
        const CFuint faceID = (*cellFaces)(elemID, iFace);

        // construct sets of nodes that make the corresponding face in this element
        if (!m_isBFace[faceID]) // is not boundary face
        {
          if (idxState[faceID] == -1) {
            ++countInFaces;
            // first time that this internal face is detected
            for (CFuint iNode = 0; iNode < nbNodesPerFace; ++iNode) {
              const CFuint localNodeID =
                (*m_faceNodeElement[iType])(iFace, iNode);

              const CFuint nodeID =
                getCFmeshData().getElementNode(elemID, localNodeID);
              (*innerFaceNodes)(countInFaces, iNode) = nodeID;
            }

            // set the local geo ID and the geometric entity type
            (*m_inLocalGeoIDs)[countInFaces] = faceID;
            (*m_inGeoTypes)[countInFaces] = m_geoTypeIDs[faceID];
            // set the first state
            (*innerFaceStates)(countInFaces, 0) = elemID;
            // set the local idx in the inner faces TRS for this face
            idxState[faceID] = static_cast<CFint>(countInFaces);
          }
          else {
            // set the second state in the inner face
            cf_assert(idxState[faceID] != -1);
            (*innerFaceStates)(idxState[faceID],1) = elemID;
          }
        }
        else // is boundary face
        {
          countBFaces++;
          for (CFuint iNode = 0; iNode < nbNodesPerFace; ++iNode)
          {
            const CFuint localNodeID = (*m_faceNodeElement[iType])(iFace, iNode);
            const CFuint nodeID = getCFmeshData().getElementNode(elemID, localNodeID);
            (*m_bFaceNodes)(countBFaces, iNode) = nodeID;
          }

          // set the local geo ID and the geometric entity type
          m_bLocalGeoIDs[countBFaces] = faceID;
          m_bGeoType[countBFaces] = m_geoTypeIDs[faceID];
          // set the first face state
          m_bFaceStateID[countBFaces] = elemID;
        }
      }
    }
  }

  // creation of the inner faces storage
  // list of TopologicalRegions
  vector<TopologicalRegion*>* innerFaceTrList = new vector<TopologicalRegion*>(1);
  (*innerFaceTrList)[0] = new TopologicalRegion();

  TopologicalRegionSet* innerFaceTRS = new TopologicalRegionSet("InnerFaces", innerFaceTrList);
  innerFaceTRS->attachTag("inner");
  // this is a trs with faces
  innerFaceTRS->attachTag("face");

  // set the TRS
  innerFaceTRS->setGeoTypes(m_inGeoTypes);
  innerFaceTRS->setGeo2NodesConn(innerFaceNodes);
  innerFaceTRS->setGeo2StatesConn(innerFaceStates);
  innerFaceTRS->setGeoEntsLocalIdx(m_inLocalGeoIDs);

  // check this in parallel
  // innerFaceTRS->setGlobalNbGeoEnts(cellGlobalIDs->size());
  // innerFaceTRS->setGeoEntsGlobalIdx(cellGlobalIDs);
  innerFaceTRS->createStatesList();

  // set the TR
  (*innerFaceTrList)[0]->setLocalNbGeoEnts(nbInnerFaces);
  (*innerFaceTrList)[0]->setGeoEntsLocalIdx(&(*m_inLocalGeoIDs)[0], 0);

  MeshDataStack::getActive()->addTrs(innerFaceTRS);
  MeshDataStack::getActive()->storeConnectivity("InnerFacesNodes", innerFaceNodes);
  MeshDataStack::getActive()->storeConnectivity("InnerFacesStates", innerFaceStates);
  MeshDataStack::getActive()->Statistics().setNbFaces(m_nbFaces);

  CFLog(NOTICE,"FVM Faces created : " <<  m_nbFaces << "\n");
}

//////////////////////////////////////////////////////////////////////////////

void FVMCC_MeshDataBuilder::createBoundaryFacesTRS()
{
  CFAUTOTRACE;

  const CFuint nbTRSs = getCFmeshData().getNbTRSs();

  SafePtr< vector<CFuint> > nbTRs =
    getCFmeshData().getNbTRs();

  SafePtr< vector< vector<CFuint> > > nbGeomEntsPerTR =
    getCFmeshData().getNbGeomEntsPerTR();

  SafePtr< vector<std::string> > nameTRS =
    getCFmeshData().getNameTRS();


  const CFuint nbBPlusPartitionFaces = m_bLocalGeoIDs.size();

  // mapping between the neighbor state ID of the boundary face and
  // its idx in the list of boundary faces
  m_mapStateIdToFaceIdx.reserve(nbBPlusPartitionFaces);
  for (CFuint i = 0; i < nbBPlusPartitionFaces; ++i) {
    m_mapStateIdToFaceIdx.insert(m_bFaceStateID[i], i);
  }
  m_mapStateIdToFaceIdx.sortKeys();

  // flag telling if the face is a partition face
  m_isPartitionFace.resize(nbBPlusPartitionFaces);
  m_isPartitionFace = true;

  for(CFuint iTRS = 0; iTRS < nbTRSs; ++iTRS) {

    cf_assert((*nameTRS)[iTRS] != "InnerCells");
    cf_assert((*nameTRS)[iTRS] != "InnerFaces");

    CFLogDebugMed("Nb TRs in TRS: " << (*nbGeomEntsPerTR)[iTRS].size() << "\n");

    Common::SafePtr<TopologicalRegionSet> ptrs =
    createTopologicalRegionSet((*nbGeomEntsPerTR)[iTRS],
                               (*nameTRS)[iTRS],
                               getCFmeshData().getTRGeoConn(iTRS),
                               iTRS);

    ptrs->attachTag("writable");
    ptrs->attachTag("boundary");

    CFLog(NOTICE, "Built TRS named " << (*nameTRS)[iTRS] << "\n");
  }

  // build partition boundary faces
  createPartitionBoundaryFacesTRS();

  CFLog(NOTICE, "Built PartitionFaces TRS \n");
}

//////////////////////////////////////////////////////////////////////////////

void FVMCC_MeshDataBuilder::setMaxNbStatesInCell()
{
  const CFuint maxNbStates = 1;
  MeshDataStack::getActive()->Statistics().setMaxNbStatesInCell(maxNbStates);
}

//////////////////////////////////////////////////////////////////////////////

void FVMCC_MeshDataBuilder::setMaxNbNodesInCell()
{
  SafePtr<vector<ElementTypeData> > elementType =
    getCFmeshData().getElementTypeData();

  CFuint maxNbNodes = 0;
  for (CFuint iType = 0; iType < getNbElementTypes(); ++iType) {
    const CFuint nbNodes =  (*elementType)[iType].getNbNodes();
    if (nbNodes > maxNbNodes) {
      maxNbNodes = nbNodes;
    }
  }

  MeshDataStack::getActive()->Statistics().setMaxNbNodesInCell(maxNbNodes);
}

//////////////////////////////////////////////////////////////////////////////

void FVMCC_MeshDataBuilder::setMaxNbFacesInCell()
{
  SafePtr<vector<ElementTypeData> > elementType =
    getCFmeshData().getElementTypeData();

  CFuint maxNbFaces = 0;
  for (CFuint iType = 0; iType < getNbElementTypes(); ++iType) {
    const CFuint nbFaces = LocalConnectionData::getInstance().getNbFacesInShape
      ((*elementType)[iType].getGeoShape());
    if (nbFaces > maxNbFaces) {
      maxNbFaces = nbFaces;
    }
  }

  MeshDataStack::getActive()->Statistics().setMaxNbFacesInCell(maxNbFaces);
}

//////////////////////////////////////////////////////////////////////////////

Common::SafePtr<Framework::TopologicalRegionSet>
FVMCC_MeshDataBuilder::createTopologicalRegionSet
(const vector<CFuint>& nbGeomEntsPerTR,
 const std::string& name,
 const TRGeoConn& trGeoConn,
 const CFuint iTRS)
{
  CFAUTOTRACE;

  // Create the TopologicalRegion storage
  const CFuint nbTRs = nbGeomEntsPerTR.size();
  vector<TopologicalRegion*>* storeTR = new vector<TopologicalRegion*>(nbTRs);

  DataHandle < Framework::State*, Framework::GLOBAL > states =
    MeshDataStack::getActive()->getStateDataSocketSink().getDataHandle();
  DataHandle < Framework::Node*, Framework::GLOBAL > nodes =
    MeshDataStack::getActive()->getNodeDataSocketSink().getDataHandle();

  // compute the total number of geometric entities
  const CFuint totalNbGeos = std::accumulate(nbGeomEntsPerTR.begin(),
                                        nbGeomEntsPerTR.end(),
                                        static_cast<CFuint>(0));

  // fill in two arrays specifying the number of states and nodes
  // per geometric entity in this TR
  std::valarray<CFuint> nbNodesInGeo(totalNbGeos);
  std::valarray<CFuint> nbStatesInGeo(totalNbGeos);
  setNbNodesNbStatesInGeo(nbGeomEntsPerTR, trGeoConn,
                          nbNodesInGeo, nbStatesInGeo);

  // overridden value (left and right-ghost states for the boundary faces)
  nbStatesInGeo = 2;

  // allocate connectivity tables (you could allocate just 1 if you were in a case
  // when the node-connectivity is the same as the state connectivity)
  ConnTable* geoNodes = new ConnTable(nbNodesInGeo);
  ConnTable* geoStates = new ConnTable(nbStatesInGeo);

  // store the connectivities in the MeshData
  /// @todo think about having SharedPtr<ConnTable > and not put
  /// the TRS connectivities in MeshData
  MeshDataStack::getActive()->storeConnectivity(name + "Nodes", geoNodes);
  MeshDataStack::getActive()->storeConnectivity(name + "States", geoStates);

  // array with all the IDs of all the geometric entities in this TRS
  // ownership of this array will be given to the TR
  vector<CFuint>* localGeoIDs = new vector<CFuint>(totalNbGeos);
  vector<CFuint>* globalGeoIDs = new vector<CFuint>(totalNbGeos);
  vector<CFuint>* geoTypes = new vector<CFuint>(totalNbGeos);
  CFuint nbProcessedGeoEnts = 0;

  const std::string faceProviderName = "Face";

  typedef CFMultiMap<CFuint, CFuint>::MapIterator MapIterator;

  SafePtr<vector<vector<vector<CFuint> > > > trsGlobalIDs =
    MeshDataStack::getActive()->getGlobalTRSGeoIDs();

  const bool hasGlobalIDs = (trsGlobalIDs->size() > 0);

  // let's create the required number of TR's
  for (CFuint iTR = 0; iTR < nbTRs; ++iTR)
  {
    const GeoConn& geoConn = trGeoConn[iTR];
    const CFuint nbGeomEnts = nbGeomEntsPerTR[iTR];

    CFLog(VERBOSE, "TopologicalRegion [" << iTR << "] : GeoEnts [" << nbProcessedGeoEnts << "] \n");

    // allocate the TopologicalRegion
    (*storeTR)[iTR] = new TopologicalRegion();
    (*storeTR)[iTR]->setLocalNbGeoEnts(nbGeomEnts);
    cf_assert(nbGeomEnts <= totalNbGeos);


// debugging block
// when there is a state that is not found with an adjacent face this will indicate which.
//
//     MapIterator mitr = m_mapStateIdToFaceIdx.begin();
//     for ( ; mitr != m_mapStateIdToFaceIdx.end(); ++mitr )
//       CFLog(INFO, "[" << mitr->first << "] [" << mitr->second << "]\n" );

    // check that the number of GeometricEntities in this TR is > 0
    // if not just skip all this
    if (nbGeomEnts > 0)
    {
      cf_assert((nbProcessedGeoEnts < localGeoIDs->size()) || 
		((nbProcessedGeoEnts == localGeoIDs->size()) && (localGeoIDs->size())== 0));
      
      (*storeTR)[iTR]->setGeoEntsLocalIdx(&(*localGeoIDs)[nbProcessedGeoEnts], nbProcessedGeoEnts);

      // set the connectivity tables of the TRS
      for (CFuint iGeo = 0; iGeo < nbGeomEnts; ++iGeo, ++nbProcessedGeoEnts)
      {

        CFLogDebugMed(" +++ iGeo [" << iGeo << "]\n");

        cf_assert(nbProcessedGeoEnts < geoStates->nbRows());
	// const CFuint nbGeoStates = geoStates->nbCols(nbProcessedGeoEnts);
        cf_assert(geoStates->nbCols(nbProcessedGeoEnts) == 2);
        cf_assert(geoConn[iGeo].second.size() == 1);
        const CFuint stateID = geoConn[iGeo].second[0];
        (*geoStates)(nbProcessedGeoEnts, 0) = stateID;

        if (stateID >= MeshDataStack::getActive()->getNbStates() )
        {
          CFLogDebugMax( "StateID: " << stateID << "\n");
          throw BadFormatException (FromHere(),"CFmesh had bad state index in GeometricEntity");
        }

        cf_assert(nbProcessedGeoEnts < geoNodes->nbRows());
        const CFuint nbGeoNodes = geoNodes->nbCols(nbProcessedGeoEnts);

// debugging block
//
// the typical error is a thrown exception from the CFMultiMap with a key not found
// This usually occurs when the mesh has TRS's which do not bound the domain
// and therefore are floating around as result of unclean CAD
//
//       CFout << name << " ";
//       CFout << "stateID =" << stateID << " ";
//       CFout << "global stateID =" << states[stateID]->getGlobalID() << "\n";
//       CFout << "nodeID =" << geoConn[iGeo].first[0] << " " << geoConn[iGeo].first[1] << " ";
//       CFout << "global nodeID =" << nodes[geoConn[iGeo].first[0]]->getGlobalID()
//             << " " << nodes[geoConn[iGeo].first[1]]->getGlobalID() << "\n";
	
	bool fo = false;
        pair<MapIterator, MapIterator> faces = m_mapStateIdToFaceIdx.find(stateID, fo);
	cf_assert(fo);
	
        CFuint faceIdx = 0;
        bool faceFound = false;
        for (MapIterator faceInMapItr = faces.first;
            faceInMapItr != faces.second; ++faceInMapItr)
        {

          // consider the current face local idx candidate
          faceIdx = faceInMapItr->second;

          // if all the face nodes corresponding to faceIdx
          // match (even if not in order), then choose this faceIdx as
          // the right one
          CFuint nbMatchingNodes = 0;
          if (nbGeoNodes == m_bFaceNodes->nbCols(faceIdx))
          {
            for (CFuint iNode = 0; iNode < nbGeoNodes; ++iNode)
            {
              const CFuint nodeID =  (*m_bFaceNodes)(faceIdx, iNode);
              for (CFuint jNode = 0; jNode < nbGeoNodes; ++jNode)
              {
                if (geoConn[iGeo].first[jNode] == nodeID)
                {
                  ++nbMatchingNodes;
                  break;
                }
              }
            }
          }

          if (nbMatchingNodes == nbGeoNodes)
          {
            // the current faceIdx is the right one
            // this face is surely not on the partition boundary
            cf_assert(faceIdx < m_isPartitionFace.size());
            m_isPartitionFace[faceIdx] = false;
            faceFound = true;
            break;
          }
        }

        if (!faceFound) throw BadFormatException (FromHere(),"boundary face not found");

        for (CFuint iNode = 0; iNode < nbGeoNodes; ++iNode) {
          const CFuint nodeID =  (*m_bFaceNodes)(faceIdx, iNode);
          (*geoNodes)(nbProcessedGeoEnts, iNode) = nodeID;

          if (nodeID >= MeshDataStack::getActive()->getNbNodes() ) {
            CFLogDebugMax( "NodeID: " << nodeID << "\n");
            throw BadFormatException (FromHere(),"CFmesh had bad node index in GeometricEntity");
          }
        }

        // assign the local ID of the current geometric entity
        (*localGeoIDs)[nbProcessedGeoEnts] = m_bLocalGeoIDs[faceIdx];

        // assign the global ID of the current geometric entity
        (*globalGeoIDs)[nbProcessedGeoEnts] = (hasGlobalIDs) ?
          (*trsGlobalIDs)[iTRS][iTR][iGeo] : iGeo;

        const CFuint dim = PhysicalModelStack::getActive()->getDim();
        const std::string faceGeoTypeName = makeGeomEntName (getGeoShape(CFGeoEnt::FACE, dim, nbGeoNodes),
                      getGeometricPolyType(),
                      getGeometricPolyOrder(),
                      getSolutionPolyType(),
                      getSolutionPolyOrder());

        const std::string geoProviderName = faceProviderName + faceGeoTypeName;
        (*geoTypes)[nbProcessedGeoEnts] = m_mapGeoProviderNameToType.find(geoProviderName);

      }  // loop geos
    }  // nb geos > 0
  }  // end loop over TR's

  // Create TopologicalRegionSet
  TopologicalRegionSet* ptrs = new TopologicalRegionSet (name, storeTR);
  // this is a trs with faces
  ptrs->attachTag("face");

  // set local GeometricEntity IDs
  ptrs->setGeoEntsLocalIdx(localGeoIDs);
  // set global GeometricEntity IDs
  ptrs->setGeoEntsGlobalIdx(globalGeoIDs);
  // set the connectivity GeometricEntity to nodeIDs
  ptrs->setGeo2NodesConn(geoNodes);
  // set the connectivity GeometricEntity to stateIDs
  ptrs->setGeo2StatesConn(geoStates);
  // set the GeometricEntity types
  ptrs->setGeoTypes(geoTypes);
  // create the cached list of state indexes in the all TRS
  ptrs->createStatesList();

  // put it into MeshData TRS list
  MeshDataStack::getActive()->addTrs(ptrs);

  cf_assert(nbProcessedGeoEnts == totalNbGeos);

  CFLogDebugMin( "FVMCC_MeshDataBuilder::setTopologicalRegions() : nameTRS "
                 << name
                 << ", nb boundary faces detected : "
                 << nbProcessedGeoEnts << "\n");

  return ptrs;
}

//////////////////////////////////////////////////////////////////////////////

void FVMCC_MeshDataBuilder::createPartitionBoundaryFacesTRS()
{
  CFAUTOTRACE;

  DataHandle < Framework::State*, Framework::GLOBAL > states =
    MeshDataStack::getActive()->getStateDataSocketSink().getDataHandle();
  DataHandle < Framework::Node*, Framework::GLOBAL > nodes =
    MeshDataStack::getActive()->getNodeDataSocketSink().getDataHandle();

  // compute the total number of geometric entities in this TRS
  const CFuint totalNbGeos = m_isPartitionFace.size() - getNbBoundaryFaces();

  // fill in two arrays specifying the number of states and nodes
  // per geometric entity in this TR
  std::valarray<CFuint> nbNodesInGeo(static_cast<CFuint>(0), totalNbGeos);
  std::valarray<CFuint> nbStatesInGeo(static_cast<CFuint>(0), totalNbGeos);

  CFuint idx = 0;
  for (CFuint i = 0; i < m_isPartitionFace.size(); ++i) {
    if (m_isPartitionFace[i]) {
      nbNodesInGeo[idx] = m_bFaceNodes->nbCols(i);
      nbStatesInGeo[idx] = 2;
      idx++;
    }
  }

  // allocate connectivity tables (you could allocate just 1 if you were in a case
  // when the node-connectivity is the same as the state connectivity)
  ConnTable* geoNodes = new ConnTable(nbNodesInGeo);
  ConnTable* geoStates = new ConnTable(nbStatesInGeo);

  // store the connectivities in the MeshData
  /// @todo think about having SharedPtr<ConnTable> and not put
  /// the TRS connectivities in MeshData
  const std::string name = "PartitionFaces";
  MeshDataStack::getActive()->storeConnectivity(name + "Nodes", geoNodes);
  MeshDataStack::getActive()->storeConnectivity(name + "States", geoStates);

  // array with all the IDs of all the geometric entities in this TRS
  // ownership of this array will be given to the TR
  vector<CFuint>* localGeoIDs = new vector<CFuint>(totalNbGeos);
  vector<CFuint>* geoTypes = new vector<CFuint>(totalNbGeos);
  const std::string faceProviderName = "Face";

  // let's create the required number of TR's
  // allocate the TopologicalRegion
  // Create the TopologicalRegion storage
  vector<TopologicalRegion*>* storeTR = new vector<TopologicalRegion*>(1);
  (*storeTR)[0] = new TopologicalRegion();
  if (totalNbGeos > 0) {
    (*storeTR)[0]->setGeoEntsLocalIdx(&(*localGeoIDs)[0], 0);
  }
  else {
    (*storeTR)[0]->setGeoEntsLocalIdx(CFNULL, 0);
  }
  (*storeTR)[0]->setLocalNbGeoEnts(totalNbGeos);

  // set the connectivity tables of the TRS
  CFuint nbProcessedFaces = 0;
  for (CFuint ibface = 0; ibface < m_isPartitionFace.size(); ++ibface) {

    if (m_isPartitionFace[ibface]) {
      const CFuint stateID = m_bFaceStateID[ibface];

      (*geoStates)(nbProcessedFaces, 0) = stateID;

      if (stateID >= states.size() ) {
        CFLogDebugMax( "StateID: " << stateID << "\n");
        throw BadFormatException (FromHere(),"CFmesh had bad state index in GeometricEntity");
      }

      const CFuint nbGeoNodes = m_bFaceNodes->nbCols(ibface);
      for (CFuint iNode = 0; iNode < nbGeoNodes; ++iNode) {
	assert(iNode < m_bFaceNodes->nbCols(ibface));
	const CFuint nodeID =  (*m_bFaceNodes)(ibface, iNode);
	assert(iNode < geoNodes->nbCols(nbProcessedFaces));
	(*geoNodes)(nbProcessedFaces, iNode) = nodeID;

        if (nodeID >= nodes.size() ) {
          CFLogDebugMax( "NodeID: " << nodeID << "\n");
          throw BadFormatException (FromHere(),"CFmesh had bad node index in GeometricEntity");
        }
      }

      // assign the local ID of the current geometric entity
      (*localGeoIDs)[nbProcessedFaces] = m_bLocalGeoIDs[ibface];

      const CFuint dim = PhysicalModelStack::getActive()->getDim();
      const std::string faceGeoTypeName = makeGeomEntName (getGeoShape(CFGeoEnt::FACE, dim, nbGeoNodes),
							getGeometricPolyType(),
							getGeometricPolyOrder(),
							getSolutionPolyType(),
							getSolutionPolyOrder());

      const std::string geoProviderName = faceProviderName + faceGeoTypeName;
      (*geoTypes)[nbProcessedFaces] = m_mapGeoProviderNameToType.find(geoProviderName);

      // increment the counter of partition faces
      nbProcessedFaces++;
    }
  }

  if (nbProcessedFaces != totalNbGeos) {
    throw BadFormatException (FromHere(),"Wrong number of partition faces");
  }

  CFLog(NOTICE,"Number of partition faces detected = " << totalNbGeos << "\n");

  // Create TopologicalRegionSet
  TopologicalRegionSet* ptrs = new TopologicalRegionSet(name, storeTR);
  ptrs->attachTag("partition");
  // this is a trs with faces
  ptrs->attachTag("face");

  // set local GeometricEntity IDs
  ptrs->setGeoEntsLocalIdx(localGeoIDs);
  // set the connectivity GeometricEntity to nodeIDs
  ptrs->setGeo2NodesConn(geoNodes);
  // set the connectivity GeometricEntity to stateIDs
  ptrs->setGeo2StatesConn(geoStates);
  // set the GeometricEntity types
  ptrs->setGeoTypes(geoTypes);
  // create the cached list of state indexes in the all TRS
  /// @todo is this needed here ???
  ptrs->createStatesList();

  // put it into MeshData TRS list
  MeshDataStack::getActive()->addTrs(ptrs);

#ifdef CF_WRITE_FV_PARTITION
  // ----------------------------------------------//
  // write the partition faces to tecplot!!!       //
  // ----------------------------------------------//

   const CFuint dim = PhysicalModelStack::getActive()->getDim();
   ofstream fout("partitionFaces.plt");

   //  Tecplot Header
   fout << "TITLE      = Boundary data" << "\n";
   fout << "VARIABLES  = ";
   for (CFuint i = 0; i < dim; ++i) {
     fout << " \"x" << i << '\"';
   }
   fout << endl;

   Common::SafePtr<TopologicalRegionSet> currTrs =
     MeshDataStack::getActive()->getTrs("PartitionFaces");
   const CFuint nbTRs = currTrs->getNbTRs();
   for (CFuint iTR = 0; iTR < nbTRs; ++iTR) {
     SafePtr<TopologicalRegion> tr = currTrs->getTopologicalRegion(iTR);
     const CFuint nbBFaces = tr->getLocalNbGeoEnts();
     if (nbBFaces > 0) {
       // get all the unique nodes in the TR
       vector<CFuint> trNodes;
       tr->putNodesInTR(trNodes);

       const CFuint nbAllNodes = trNodes.size();
       cout <<"nbAllNodes = " << nbAllNodes << endl;
       CFMap<CFuint,CFuint> mapNodesID(nbAllNodes);
       for (CFuint i = 0; i < nbAllNodes; ++i) {
         mapNodesID.insert(trNodes[i],i+1);
       }
       mapNodesID.sortKeys();

       std::string elemShape;
       if (dim == 2) {
         elemShape = "LINESEG";
       }
       else if (dim == 3) {
         /// @TODO AL: this assumes all the faces have the same shape
         ///           within the same TR (triangle or quad) !!!
         /// think about getting the maximum size and then adding an extra virtual
         /// node if TETRA and QUADRILATERAL are both present
         const CFuint nbNodesInFace = tr->getNbNodesInGeo(0);
         elemShape = (nbNodesInFace == 3) ? "TRIANGLE" : "QUADRILATERAL";
         cout << "nbNodesInFace = " << nbNodesInFace << endl;
       }

       // Tecplot doesn't read zones with 0 elements !!! (this was a problem in parallel)
       if (tr->getLocalNbGeoEnts() > 0) {
         // print zone header
         // one zone per TR
         fout << "ZONE N=" << nbAllNodes
  	    << ", T=\"" << currTrs->getName() << ", TR " << iTR << "\""
  	    << ", E=" << tr->getLocalNbGeoEnts()
  	    << ", F=FEPOINT"
  	    << ", ET=" << elemShape
  	    << "\n";

         vector<CFuint>::const_iterator itr;
         // print  nodal coordinates and stored nodal variables
         for (itr = trNodes.begin(); itr != trNodes.end(); ++itr) {
  	 // node has to be printed with the right length
  	 const CFuint nodeID = *itr;
  	 fout << *nodes[nodeID] << endl;
         }

         for (CFuint iGeo = 0; iGeo < nbBFaces; ++iGeo) {
  	 const CFuint nbNodesInGeo = tr->getNbNodesInGeo(iGeo);
  	 for (CFuint in = 0; in < nbNodesInGeo; ++in) {
  	   fout << mapNodesID.find(tr->getNodeID(iGeo,in)) << " ";
  	 }
  	 fout << "\n";
  	 fout.flush();
         }
       }
     }
   }

   fout.close();
#endif
}

//////////////////////////////////////////////////////////////////////////////

void FVMCC_MeshDataBuilder::setMapGeoToTrs()
{
  MapGeoToTrsAndIdx* mapGeoToTrs = new MapGeoToTrsAndIdx();
  const CFuint nbFaces = MeshDataStack::getActive()->Statistics().getNbFaces();
  mapGeoToTrs->resize(nbFaces);

  vector<Common::SafePtr<TopologicalRegionSet> > trs = MeshDataStack::getActive()->getTrsList();

  for (CFuint i = 0; i < trs.size(); ++i) {
    SafePtr<TopologicalRegionSet> currTrs = trs[i];
    if (currTrs->hasTag("face"))
    {
      const bool isOnBoundary = (currTrs->hasTag("partition") || currTrs->hasTag("boundary"));
      const CFuint nbTrsFaces = currTrs->getLocalNbGeoEnts();
      for (CFuint iFace = 0; iFace < nbTrsFaces; ++iFace) {
        const CFuint faceID = currTrs->getLocalGeoID(iFace);
        mapGeoToTrs->setMappingData
          (faceID, currTrs, iFace, isOnBoundary);
      }
    }
  }

  MeshDataStack::getActive()->storeMapGeoToTrs("MapFacesToTrs", mapGeoToTrs);
}

//////////////////////////////////////////////////////////////////////////////

void FVMCC_MeshDataBuilder::renumberCells()
{
  CFLog(VERBOSE, "FVMCC_MeshDataBuilder::renumberCells() => START\n");
  
  SafePtr<TopologicalRegionSet> cells =
    MeshDataStack::getActive()->getTrs("InnerCells");

  DataHandle < Framework::State*, Framework::GLOBAL > states =
    MeshDataStack::getActive()->getStateDataSocketSink().getDataHandle();

  const CFuint nbCells = cells->getLocalNbGeoEnts();
  bool doRenumber = true;

  for (CFuint i = 0; i < nbCells; ++i) {
    const CFuint globalStateID = states[cells->getStateID(i,0)]->getGlobalID();
    if (globalStateID != cells->getGlobalGeoID(i)) {
      doRenumber = false;
      break;
    }
  } 
  
  CFLog(INFO, "nbCells = " << nbCells << ", nbStates = " << states.size() << "\n");
  
  if (doRenumber) {
    // build a mapping global cellID to local cellID
    CFMap<CFuint, CFuint> mapGlobalToLocalCellID;
    mapGlobalToLocalCellID.reserve(nbCells);
    cf_assert(nbCells == states.size());
    
    for (CFuint i = 0; i < nbCells; ++i) {
      mapGlobalToLocalCellID.insert(cells->getGlobalGeoID(i),i);
    }
    mapGlobalToLocalCellID.sortKeys();

    // back up the cell-node connectivity
    ConnectivityTable<CFuint> oldCellNodes
      (*MeshDataStack::getActive()->getConnectivity("cellNodes_InnerCells"));

    // perform a local renumbering of the cells so that their local IDs match
    // the local state IDs
    for (CFuint iState = 0; iState < states.size(); ++iState) {
      const CFuint globalStateID = states[iState]->getGlobalID();
      const CFuint oldLocalCellID = mapGlobalToLocalCellID.find(globalStateID);
      const CFuint nbNodesInCell = cells->getNbNodesInGeo(oldLocalCellID);
      // renumber the node connectivity of the current cell
      for (CFuint iNode = 0; iNode < nbNodesInCell; ++iNode) {
	const CFuint nodeID = oldCellNodes(oldLocalCellID, iNode);
	cells->setNodeID(iState, iNode, nodeID);
      }
      // renumber the state connectivity of the current cell
      cells->setStateID(iState, 0, iState);
      // set the global ID of the current cell to be equal to the
      // corresponding global state ID
      cells->setGlobalGeoID(iState, globalStateID);
    }

    // redo the consistency check
    for (CFuint i = 0; i < nbCells; ++i) {
      const CFuint globalStateID = states[cells->getStateID(i,0)]->getGlobalID();
      if (globalStateID != cells->getGlobalGeoID(i)) {
	CFLog(ERROR, "ERROR: globalStateID " << globalStateID
	      << " != globalCellID " << cells->getGlobalGeoID(i) << "\n");
	abort();
      }
    }
  }

  CFLog(VERBOSE, "FVMCC_MeshDataBuilder::renumberCells() => END\n");
}

//////////////////////////////////////////////////////////////////////////////

} // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

