#include <set>
#include <algorithm>
#include <numeric>

#include "Common/SwapEmpty.hh"

#include "Common/CFLog.hh"
#include "Environment/ObjectProvider.hh"

#include "Framework/CFSide.hh"
#include "Framework/MeshData.hh"
#include "Framework/LocalConnectionData.hh"
#include "Framework/BaseGeometricEntityProvider.hh"
#include "Framework/Face.hh"
#include "Framework/Cell.hh"
#include "Framework/PhysicalModel.hh"
#include "Framework/SetElementStateCoord.hh"
#include "Framework/BadFormatException.hh"
#include "Framework/MapGeoToTrsAndIdx.hh"

#include "FluxReconstructionMethod/HexaFluxReconstructionElementData.hh"
#include "FluxReconstructionMethod/QuadFluxReconstructionElementData.hh"
#include "FluxReconstructionMethod/TriagFluxReconstructionElementData.hh"
#include "FluxReconstructionMethod/TetraFluxReconstructionElementData.hh"
#include "FluxReconstructionMethod/FluxReconstruction.hh"
#include "FluxReconstructionMethod/FluxReconstructionBuilder.hh"
#include "FluxReconstructionMethod/FluxReconstructionElementData.hh"


//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<FluxReconstructionBuilder,
			    MeshDataBuilder,
			    FluxReconstructionModule,
			    1>
StdFRBuilderProvider("StdBuilder");

//////////////////////////////////////////////////////////////////////////////

FluxReconstructionBuilder::FluxReconstructionBuilder(const std::string& name) :
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
  m_bFaceCell(),
  m_bFaceNodes(CFNULL),
  m_isPartitionFace(),
  m_maxNbrStatesInCell()
 {
 }

//////////////////////////////////////////////////////////////////////////////

FluxReconstructionBuilder::~FluxReconstructionBuilder()
{
}

//////////////////////////////////////////////////////////////////////////////

void FluxReconstructionBuilder::configure ( Config::ConfigArgs& args )
{
  MeshDataBuilder::configure(args);
}

//////////////////////////////////////////////////////////////////////////////

void FluxReconstructionBuilder::releaseMemory()
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
  deletePtr(m_bFaceNodes);
  m_isPartitionFace.resize(0);
}

//////////////////////////////////////////////////////////////////////////////

void FluxReconstructionBuilder::setMapGeoToTrs()
{
}

//////////////////////////////////////////////////////////////////////////////

void FluxReconstructionBuilder::createTopologicalRegionSets()
{
  CFAUTOTRACE;

  // first create the cells
  createInnerCells();

  // put the coordinates in the states
  /// @todo is the following line still necessary?
  setCoordInCellStates();

  // create the cell-face connectivity
  createCellFaces();

  // create the inner faces TRS
  createInnerFacesTRS();

  // reorder the faces list
  reorderInnerFacesTRS();

  // create all the boundary faces TRSs
  createBoundaryFacesTRS();

  // reorder the boundary faces list
  reorderBoundaryFacesTRS();
  
  // reorder the partition faces list
  reorderPartitionFacesTRS();

  // needed for the sparsity of the matrix
  setMapGeoToTrs();
}

//////////////////////////////////////////////////////////////////////////////

void FluxReconstructionBuilder::computeGeoTypeInfo()
{
  CFAUTOTRACE;

  // set polynomial type
  getCFmeshData().setSolutionPolyType(CFPolyForm::FLUXRECONSTRUCTION);
  CFLog(INFO, "solPolyType = " << getCFmeshData().getSolutionPolyType() << "\n");

  // continue with the standard algorithm
  Framework::MeshDataBuilder::computeGeoTypeInfo();
}

//////////////////////////////////////////////////////////////////////////////

void FluxReconstructionBuilder::createCellFaces()
{
  CFAUTOTRACE;

  const CFuint nbElem      = getNbElements();
  const CFuint nbElemTypes = getNbElementTypes();
  CFLog(INFO, "nbElem = " << nbElem << "\n");
  CFLog(VERBOSE, "nbElemTypes = " << nbElemTypes << "\n");


  // local connectivity face-node for each element type
  // m_faceNodeElem is a vector of pointers to tables containing CFuints
  m_faceNodeElement.resize(nbElemTypes);

  // get the elementTypeData
  SafePtr<vector<ElementTypeData> > elementType = getCFmeshData().getElementTypeData();

  cf_assert(nbElemTypes == elementType->size());

  
  // loop over the types of elements
  for (CFuint iType = 0; iType < nbElemTypes; ++iType)
  {
    m_faceNodeElement[iType] = LocalConnectionData::getInstance().getFaceDofLocal((*elementType)[iType].getGeoShape(),
                                                                    getGeometricPolyOrder(),
                                                                    NODE,CFPolyForm::LAGRANGE);
  }

  // array storing the number of faces per element
  m_nbFacesPerElem.resize(nbElem);
  m_nbFacesPerElem = 0;
  LocalConnectionData::getInstance().setNbFacesPerElement(m_nbFacesPerElem); // pass m_nbFacesPerElem by reference

  // set the face shapes per element type
  vector< vector<CFGeoShape::Type> > faceShapesPerElemType(nbElemTypes);
  LocalConnectionData::getInstance().setFaceShapesPerElemType(faceShapesPerElemType); // pass faceShapesPerElemType by reference

  // table storing the connectivity element-local faceID
  // (in this processor)
  ConnTable* cellFaces = new ConnTable(m_nbFacesPerElem);
  // store the connectivity in MeshData
  MeshDataStack::getActive()->storeConnectivity("InnerCells-Cells2Faces", cellFaces);

  const CFuint totNbNodes = MeshDataStack::getActive()->getNbNodes();
  CFLog(INFO, "nbNodes = " << totNbNodes << "\n");

  // allocate a table mapping node - faceID
  vector < vector<CFuint> > mapNodeFace(totNbNodes);

  // atomic number to indicate the maximum possible number
  // of nodes in a face (this allow to avoid frequent reallocations
  // of the vector nodesInFace)
  const CFuint maxNbNodesInFace = 100;
  vector<CFuint> nodesInFace(maxNbNodesInFace);
  const std::string faceProviderName = "Face";

  // number of boundary faces in TRS data read from mesh file
  // this DOESN'T include partition faces
  const CFuint nbBoundaryFaces = getNbBoundaryFaces();
  const CFuint sumNbFacesPerElem = m_nbFacesPerElem.sum();

  // max possible number of inner faces (>= actual value in SERIAL simulation)
  // const CFuint maxNbInnerFaces = (m_nbFacesPerElem.sum() - nbBoundaryFaces)/2;
  // the max number of total TRSs is always slightly overestimated
  // ONLY in serial run totalNbFaces = nbBoundaryFaces + maxNbInnerFaces
  const CFuint maxTotalNbFaces = sumNbFacesPerElem;

  CFLog(INFO, "nbBoundaryFaces = " << nbBoundaryFaces << "\n");
  CFLog(VERBOSE, "maxTotalNbFaces = " << maxTotalNbFaces << "\n");

  // the following arrays are oversized
  m_geoTypeIDs.resize(maxTotalNbFaces);
  m_isBFace.resize(maxTotalNbFaces);
  for (CFuint i = 0; i < maxTotalNbFaces; ++i)
  {
    m_isBFace[i] = true;
  }

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

      CFLog(INFO, "Face provider name = " << providerName << "\n");

      faceGeoTypeID[iFace] = m_mapGeoProviderNameToType.find(providerName);
    }

    // loop over the elements of this type
    const CFuint nbElemPerType = (*elementType)[iType].getNbElems();
    
    CFLog(INFO, "nbElemTypes: " << nbElemPerType << "\n");

    for (CFuint iElem = 0; iElem < nbElemPerType; ++iElem, ++elemID)
    {
      const CFuint nbElemFaces = m_nbFacesPerElem[elemID];
      
      // loop over the faces in the current element
      for (CFuint iFace = 0; iFace < nbElemFaces; ++iFace)
      {
        const CFuint nbNodesPerFace = m_faceNodeElement[iType]->nbCols(iFace);

        // construct sets of nodes that make the corresponding face
        // in this element
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
        for (CFuint iFaceID = 0; iFaceID < nbFaceIDsRefNode; ++iFaceID)
        {
          currFaceID = mapNodeFace[nodeID][iFaceID];

          countNodes = 1;
          for (CFuint iNode = 1; iNode < nbNodesPerFace; ++iNode)
          {
            const CFuint newNodeID = nodesInFace[iNode];
            // number of faceIDs referencing node
            const CFuint nbFaceIDsRefThisNode = mapNodeFace[newNodeID].size();
            for (CFuint jFaceID = 0; jFaceID < nbFaceIDsRefThisNode; ++jFaceID)
            {
              if (mapNodeFace[newNodeID][jFaceID] == currFaceID)
              {
                ++countNodes;

                // break from the loop over faces referencing the new node
                break;
              }
            }
          }
          if (countNodes == nbNodesPerFace)
          {
            // the corresponding faceID already exists, meaning
            // that the face is an internal one, shared by two elements
            // here you set the second element neighbor of the face
            faceFound = true;
            (*cellFaces)(elemID, iFace) = currFaceID;

            // since it has two neighbor cells,
            // this face is surely NOT a boundary face
            m_isBFace[currFaceID] = false;

            // increment number of inner faces (they always have 2 states)
            ++countInFaces;

            // break from loop over faces referencing the first node
            break;
          }
        }

        if (!faceFound)
        {
          // a new face has been found
          // add the ID of the new face in the corresponding nodes referencing it
          for (CFuint iNode = 0; iNode < nbNodesPerFace; ++iNode)
          {
            const CFuint currNodeID = nodesInFace[iNode];
            mapNodeFace[currNodeID].push_back(m_nbFaces);
          }

          // store the geometric entity type for the current face
          m_geoTypeIDs[m_nbFaces] = faceGeoTypeID[iFace];

          (*cellFaces)(elemID, iFace) = m_nbFaces;
          nbFaceNodes.push_back(nbNodesPerFace);

          // increment the number of faces
          ++m_nbFaces;
        }
      }
    }
  }

  cf_assert(m_nbFaces <= maxTotalNbFaces);
  cf_assert(countInFaces <= maxTotalNbFaces);
  
    //CFLog(VERBOSE, "cellFaces = \n" << *cellFaces << "\n");

  const CFuint totalNbFaces = m_nbFaces;
  const CFuint nbInnerFaces = countInFaces;

  CFLog(INFO, "totalNbFaces = " << totalNbFaces << "\n");
  CFLog(INFO, "nbInnerFaces = " << nbInnerFaces << "\n");

  // total number of boundary + partition boundary faces
  const CFuint nbBPlusPartitionFaces = totalNbFaces - nbInnerFaces;
  // CFout << "nbBPlusPartitionFaces = " << nbBPlusPartitionFaces << "\n";
  CFLog(INFO, "nbBPlusPartitionFaces = " << nbBPlusPartitionFaces << "\n");

  m_nbInFacesNodes.resize(nbInnerFaces);
  m_nbBFacesNodes.resize(nbBPlusPartitionFaces);

  // set the number of nodes in faces
  CFuint iBFace = 0;
  CFuint iInFace = 0;
  for (CFuint i = 0; i < nbFaceNodes.size(); ++i)
  {
    if (!m_isBFace[i])
    {
      m_nbInFacesNodes[iInFace++] = nbFaceNodes[i];
    }
    else
    {
      m_nbBFacesNodes[iBFace++] = nbFaceNodes[i];
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void FluxReconstructionBuilder::createInnerFacesTRS()
{
  CFAUTOTRACE;

  const CFuint nbInnerFaces = m_nbInFacesNodes.size();
  const CFuint nbBPlusPartitionFaces = m_nbBFacesNodes.size();

  // Create connectivity table inner face - nodes
  ConnTable* innerFaceNodes = new ConnTable(m_nbInFacesNodes);

  // Create connectivity table boundary face - nodes
  m_bFaceNodes = new ConnTable(m_nbBFacesNodes);

  // Create connectivity table inner face - cell
  std::valarray<CFuint> nbInFacesCells(2,nbInnerFaces);
  ConnTable* innerFaceCells = new ConnTable(nbInFacesCells);

  // Create connectivity boundary face - cell (is a vector!!)
  m_bFaceCell.resize(nbBPlusPartitionFaces);

  // Get the cell - face connectivity from the MeshData
  SafePtr<ConnTable> cellFaces = MeshDataStack::getActive()->getConnectivity("InnerCells-Cells2Faces");

  // Get the ElemenTypeData from CFMeshData
  SafePtr<vector<ElementTypeData> > elementType = getCFmeshData().getElementTypeData();

  // std::valarray containing the face local index in InnerFaces TRS
  std::valarray<CFint> idxFace(-1, m_nbFaces);

  // Connectivity inner face idx - geotype
  m_inGeoTypes = new vector<CFuint>(nbInnerFaces);
  // Connectivity inner face idx - ID in local processor
  m_inLocalGeoIDs = new vector<CFuint>(nbInnerFaces);
  // Connectivity boundary face local idx - geotype
  m_bGeoType.resize(nbBPlusPartitionFaces);
  // Connectivity boundary face local idx - ID in local processor
  m_bLocalGeoIDs.resize(nbBPlusPartitionFaces);

  // reset the number of faces
  CFint countBFaces = -1;
  CFint countInFaces = -1;

  // loop over the types (there should be only one type!)
  CFuint elemID = 0;
  const CFuint nbElemTypes = elementType->size();
  cf_assert(nbElemTypes == 1);
  for (CFuint iType = 0; iType < nbElemTypes; ++iType)
  {
    // loop over the elements of this type
    const CFuint nbElemPerType = (*elementType)[iType].getNbElems();
    for (CFuint iElem = 0; iElem < nbElemPerType; ++iElem, ++elemID)
    {
      // loop over the faces in the current element
      const CFuint nbElemFaces = m_nbFacesPerElem[elemID];
      for (CFuint iFace = 0; iFace < nbElemFaces; ++iFace)
      {
        const CFuint nbNodesPerFace = m_faceNodeElement[iType]->nbCols(iFace);
        const CFuint faceID = (*cellFaces)(elemID, iFace);

        // construct sets of nodes that make the corresponding face in this element
        if (!m_isBFace[faceID])
        {
          if (idxFace[faceID] == -1)
          {
            ++countInFaces;

            // first time that this internal face is detected
            for (CFuint iNode = 0; iNode < nbNodesPerFace; ++iNode)
            {
              const CFuint localNodeID = (*m_faceNodeElement[iType])(iFace, iNode);

              const CFuint nodeID = getCFmeshData().getElementNode(elemID, localNodeID);

              (*innerFaceNodes)(countInFaces, iNode) = nodeID;
            }

            // set the local geo ID and the geometric entity type
            (*m_inLocalGeoIDs)[countInFaces] = faceID;
            (*m_inGeoTypes)[countInFaces] = m_geoTypeIDs[faceID];

            // set the first neighbouring element
            (*innerFaceCells)(countInFaces, LEFT) = elemID;

            // set the local idx in the inner faces TRS for this face
            idxFace[faceID] = static_cast<CFint>(countInFaces);
          }
          else
          {
            // set the second neighbouring element in the inner face
            cf_assert(idxFace[faceID] != -1);
            (*innerFaceCells)(idxFace[faceID], RIGHT) = elemID;
          }
        }
        else
        {
          ++countBFaces;
          for (CFuint iNode = 0; iNode < nbNodesPerFace; ++iNode)
          {
            const CFuint localNodeID = (*m_faceNodeElement[iType])(iFace, iNode);
            const CFuint nodeID = getCFmeshData().getElementNode(elemID, localNodeID);
            (*m_bFaceNodes)(countBFaces, iNode) = nodeID;
          }

          // set the local geo ID and the geometric entity type
          m_bLocalGeoIDs[countBFaces] = faceID;
          m_bGeoType[countBFaces] = m_geoTypeIDs[faceID];
          // set the neighbouring element
          m_bFaceCell[countBFaces] = elemID;
        }
      }
    }
  }

  // creation of the TR's for "InnerFaces"
  vector<TopologicalRegion*>* innerFaceTrList = new vector<TopologicalRegion*>(1);
  (*innerFaceTrList)[0] = new TopologicalRegion();

  // creation of the "InnerFaces" TRS
  TopologicalRegionSet* innerFaceTRS = new TopologicalRegionSet ("InnerFaces", innerFaceTrList);
  // this is a inner trs
  innerFaceTRS->attachTag("inner");
  // this is a trs with faces
  innerFaceTRS->attachTag("face");

  // set the TRS
  innerFaceTRS->setGeoTypes(m_inGeoTypes);
  innerFaceTRS->setGeo2NodesConn(innerFaceNodes);
  innerFaceTRS->setGeoEntsLocalIdx(m_inLocalGeoIDs);

  // check this in parallel
  // innerFaceTRS->setGlobalNbGeoEnts(cellGlobalIDs->size());
  // innerFaceTRS->setGeoEntsGlobalIdx(cellGlobalIDs);

  // set the TR
  (*innerFaceTrList)[0]->setLocalNbGeoEnts(nbInnerFaces);
  (*innerFaceTrList)[0]->setGeoEntsLocalIdx(&(*m_inLocalGeoIDs)[0], 0);

  // Add to MeshData
  MeshDataStack::getActive()->addTrs(innerFaceTRS);
  MeshDataStack::getActive()->storeConnectivity("InnerFacesNodes", innerFaceNodes);
  MeshDataStack::getActive()->storeConnectivity("InnerFaces-Faces2Cells", innerFaceCells);
  MeshDataStack::getActive()->Statistics().setNbFaces(m_nbFaces);

  CFLog(NOTICE,"Created FR TRS InnerFaces with "
        << innerFaceTRS->getLocalNbGeoEnts() << " local faces\n");
}

//////////////////////////////////////////////////////////////////////////////

void FluxReconstructionBuilder::reorderInnerFacesTRS()
{
  CFAUTOTRACE;
  
  CFLog(VERBOSE,"Reordering inner faces TRS\n");
  // Pointer to list containing the GeoType of the face ==> class variable of FluxReconstructionBuilder, we can access this variable directly through m_inGeoTypes
  // Same for list containing the local (this processor) index of the face, variable m_inLocalGeoIDs

  // Get connectivities
  SafePtr< ConnTable > cellNodes      = MeshDataStack::getActive()->getConnectivity("cellNodes_InnerCells");
  SafePtr< ConnTable > innerFaceNodes = MeshDataStack::getActive()->getConnectivity("InnerFacesNodes");
  SafePtr< ConnTable > innerFaceCells = MeshDataStack::getActive()->getConnectivity("InnerFaces-Faces2Cells");

  // Number of inner faces
  const CFuint nbInnerFaces = m_inGeoTypes->size();

  // Vector to hold the already detected face GeoTypes
  vector< CFuint > faceGeoTypes;
  // Create connectivity table face - inner face idx and orientation or bc type
  std::valarray<CFuint> dataSize(m_nbFaces);
  for (CFuint iFace = 0; iFace < m_nbFaces; ++iFace)
  {
    dataSize[iFace] = m_isBFace[iFace] ? 1 : 2;
  }
  ConnTable* faceToInFaceIdxOrientOrBCIdx = new ConnTable(dataSize);

  // Get the number of different face types in mesh (there should only be one type for FluxReconstruction for now!!)
  for (CFuint iFace = 0; iFace < nbInnerFaces; ++iFace)
  {
    bool isNewType = true;
    const CFuint faceGeoType = (*m_inGeoTypes)[iFace];
    for (CFuint iType = 0; iType < faceGeoTypes.size(); ++iType)
      if (faceGeoTypes[iType] == faceGeoType) isNewType = false;

    // Add new type to list of detected geo types
    if (isNewType)
      faceGeoTypes.push_back(faceGeoType);
  }
  CFLog(INFO, "nbFaceTypes: " << faceGeoTypes.size() << "\n");
  
  // Get vector containing the nodes connectivities for each orientation
  SafePtr<vector<ElementTypeData> > elementType = getCFmeshData().getElementTypeData();
  cf_assert(elementType->size() == 1); // there should be only one element type
  const CFuint nbElemTypes = elementType->size();
  const CFGeoShape::Type cellShape = (*elementType)[0].getGeoShape();
  const vector< vector < vector < CFuint > > > nodeConnPerOrientation =
      getNodeConnPerOrientation(cellShape);
  const CFuint nbrNodesPerFace = nodeConnPerOrientation[0][0].size();
  const vector < vector < CFuint > > faceConnPerOrientation =
      getFaceConnPerOrientation(cellShape);
      
  // Rearrange the faces
  const CFuint nbrOrient = nodeConnPerOrientation.size();

  /// @warning I'm cheating a little here, by using a connectivity table to store the
  /// start indexes of the faces with a certain orientation. Is there a better way?
  std::valarray<CFuint> nbrOrientValAr(1,nbrOrient+1);
  ConnTable* innerFacesStartIdxs = new ConnTable(nbrOrientValAr);
  CFuint faceIdx = 0;
  for (CFuint iOrient = 0; iOrient < nbrOrient; ++iOrient)
  {
    // assign start index of first range of faces
    (*innerFacesStartIdxs)(iOrient,0) = faceIdx;
    for (CFuint iFace = faceIdx; iFace < nbInnerFaces; ++iFace)
    {
      // Get face cell neighbour local IDs
      CFuint leftCellID  = (*innerFaceCells)(iFace,LEFT );
      CFuint rightCellID = (*innerFaceCells)(iFace,RIGHT);

      // Vectors for cell nodes IDs
      vector< CFuint > leftCellFaceNodes (nbrNodesPerFace);
      vector< CFuint > rightCellFaceNodes(nbrNodesPerFace);

      // Check if face has this orientation
      // First try for current face orientation
      // Get neighbour cell face nodes local (this processor) IDs
      for (CFuint iNode = 0; iNode < nbrNodesPerFace; ++iNode)
      {
        const CFuint leftCellNodeLocalID  = nodeConnPerOrientation[iOrient][LEFT ][iNode];
        const CFuint rightCellNodeLocalID = nodeConnPerOrientation[iOrient][RIGHT][iNode];

        leftCellFaceNodes [iNode] = (*cellNodes)(leftCellID ,leftCellNodeLocalID );
        rightCellFaceNodes[iNode] = (*cellNodes)(rightCellID,rightCellNodeLocalID);
      }
      // Check if nodes match
      bool hasThisOrient = true;
      for (CFuint iNode = 0; iNode < nbrNodesPerFace && hasThisOrient; ++iNode)
      {
        if (leftCellFaceNodes[iNode] != rightCellFaceNodes[iNode])
        {
          hasThisOrient = false;
          break;
        }
      }
      // If no match was found, try for opposite face orientation
      if (!hasThisOrient)
      {
        // Change face orientation
        const CFuint swap = leftCellID;
        leftCellID = rightCellID;
        rightCellID = swap;

        // Get neighbour cell face nodes local (this processor) IDs
        for (CFuint iNode = 0; iNode < nbrNodesPerFace; ++iNode)
        {
          const CFuint leftCellNodeLocalID  = nodeConnPerOrientation[iOrient][LEFT ][iNode];
          const CFuint rightCellNodeLocalID = nodeConnPerOrientation[iOrient][RIGHT][iNode];

          leftCellFaceNodes[iNode]  = (*cellNodes)(leftCellID ,leftCellNodeLocalID );
          rightCellFaceNodes[iNode] = (*cellNodes)(rightCellID,rightCellNodeLocalID);
        }

        // Check if nodes match
        hasThisOrient = true;
        for (CFuint iNode = 0; iNode < nbrNodesPerFace && hasThisOrient; ++iNode)
        {
          if (leftCellFaceNodes[iNode] != rightCellFaceNodes[iNode])
          {
            hasThisOrient = false;
            break;
          }
        }
      }
      // If a match was found, move this face to the next place in the InnerFaces list
      if (hasThisOrient)
      {
        // swap list entries if needed
        if (iFace != faceIdx)
        {
          CFuint swap;

          // swap inner face-local (this processor) ID connectivity
          swap = (*m_inLocalGeoIDs)[iFace];
          (*m_inLocalGeoIDs)[iFace] = (*m_inLocalGeoIDs)[faceIdx];
          (*m_inLocalGeoIDs)[faceIdx] = swap;

          /// @todo swap inner face-global ID connectivity necessary?

          // swap inner face-geotype connectivity
          swap = (*m_inGeoTypes)[iFace];
          (*m_inGeoTypes)[iFace] = (*m_inGeoTypes)[faceIdx];
          (*m_inGeoTypes)[faceIdx] = swap;
        }
        // swap inner face-cell connectivity (outside loop to ensure proper normal orientation)
        (*innerFaceCells)(iFace,LEFT ) = (*innerFaceCells)(faceIdx,LEFT );
        (*innerFaceCells)(iFace,RIGHT) = (*innerFaceCells)(faceIdx,RIGHT);
        (*innerFaceCells)(faceIdx,LEFT ) = leftCellID ;
        (*innerFaceCells)(faceIdx,RIGHT) = rightCellID;

        // find left neighbouring cell type
        CFuint leftCellType;
        for (leftCellType = 0; leftCellType < nbElemTypes; ++leftCellType)
        {
          if ((*elementType)[leftCellType].getStartIdx() <= leftCellID &&
              (*elementType)[leftCellType].getEndIdx()   >  leftCellID   )
          {
            break;
          }
        }
        // set inner face-node connectivity (to ensure correct face orientation for normal computation)
        // no need for swapping here, innerFaceNodes is not used in the identification of the face orientations
        // this will give problems if the faces at faceIdx and at iFace have a different number of nodes,
        // namely when e.g. faceIdx is P1 and iFace is P2.
        const CFuint faceLocalIDLeftCell = faceConnPerOrientation[iOrient][LEFT];
        const CFuint nbrCurrFaceNodes = innerFaceNodes->nbCols(iFace);
        for (CFuint iNode = 0; iNode < nbrCurrFaceNodes; ++iNode)
        {
          const CFuint nodeLocalID = (*m_faceNodeElement[leftCellType])(faceLocalIDLeftCell,iNode);
          (*innerFaceNodes)(faceIdx,iNode) = (*cellNodes)(leftCellID,nodeLocalID);
        }
        // set inner face index
        (*faceToInFaceIdxOrientOrBCIdx)((*m_inLocalGeoIDs)[faceIdx],0) = faceIdx;

        // set inner face orientation
        (*faceToInFaceIdxOrientOrBCIdx)((*m_inLocalGeoIDs)[faceIdx],1) = iOrient;

        // Increase face counter
        ++faceIdx;
      }
    }
  }

  // store connectivity global face ID - inner face ID or BC type
  MeshDataStack::getActive()->storeConnectivity("faceToInFaceIdxOrientOrBCIdx", faceToInFaceIdxOrientOrBCIdx);

  // assign end index of last range of faces and store the start indexes
  (*innerFacesStartIdxs)(nbrOrient,0) = faceIdx;
  MeshDataStack::getActive()->storeConnectivity("innerFacesStartIdxs", innerFacesStartIdxs);
  CFLog(INFO, "nbInnerFaces: " << nbInnerFaces << ", last face idx: " << faceIdx << "\n");
  // Check if all faces have been assigned an orientation
  //cf_assert(faceIdx == nbInnerFaces);
  if (faceIdx != nbInnerFaces) CFLog(INFO, "WARNING: MESH PROBLEM DETECTED!\n");
  
  CFLog(VERBOSE,"Reordering inner faces TRS complete\n");
}

//////////////////////////////////////////////////////////////////////////////

vector< vector < vector < CFuint > > > FluxReconstructionBuilder::getNodeConnPerOrientation(const CFGeoShape::Type shape)
{
  CFAUTOTRACE;

  FluxReconstructionElementData* frElemData;
  switch (shape)
  {
    case CFGeoShape::LINE:
    {
      throw Common::NotImplementedException (FromHere(),"Flux Reconstruction has not been implemented for 1D");
    } break;
    case CFGeoShape::QUAD:
    {
      frElemData = new QuadFluxReconstructionElementData(CFPolyOrder::ORDER0);
    } break;
    case CFGeoShape::HEXA:
    {
      frElemData = new HexaFluxReconstructionElementData(CFPolyOrder::ORDER0);
    } break;
    case CFGeoShape::TRIAG:
    {
      frElemData = new TriagFluxReconstructionElementData(CFPolyOrder::ORDER0);
    } break;
    case CFGeoShape::TETRA:
    {
      frElemData = new TetraFluxReconstructionElementData(CFPolyOrder::ORDER0);
    } break;
    default:
    {
      throw Common::ShouldNotBeHereException (FromHere(),"Unsupported cell shape");
    }
  }
  
  // get vector containing the nodes connectivities for each orientation
  vector< vector < vector < CFuint > > > nodeConnPerOrientation
  = *frElemData->getFaceNodeConnPerOrient();

  delete frElemData;

  return nodeConnPerOrientation;
}

//////////////////////////////////////////////////////////////////////////////

vector < vector < CFuint > > FluxReconstructionBuilder::getFaceConnPerOrientation(const CFGeoShape::Type shape)
{
  CFAUTOTRACE;

  FluxReconstructionElementData* frElemData;
  switch (shape)
  {
    case CFGeoShape::LINE:
    {
      throw Common::NotImplementedException (FromHere(),"Flux Reconstruction has not been implemented for 1D");
    } break;
    case CFGeoShape::QUAD:
    {
      frElemData = new QuadFluxReconstructionElementData(CFPolyOrder::ORDER0);
    } break;
    case CFGeoShape::HEXA:
    {
      frElemData = new HexaFluxReconstructionElementData(CFPolyOrder::ORDER0);
    } break;
    case CFGeoShape::TRIAG:
    {
      frElemData = new TriagFluxReconstructionElementData(CFPolyOrder::ORDER0);
    } break;
    case CFGeoShape::TETRA:
    {
      frElemData = new TetraFluxReconstructionElementData(CFPolyOrder::ORDER0);
    } break;
    default:
    {
      throw Common::ShouldNotBeHereException (FromHere(),"Unsupported cell shape");
    }
  }

  // get vector containing the nodes connectivities for each orientation
  vector < vector < CFuint > > faceConnPerOrientation
  = *frElemData->getFaceConnPerOrient();
  delete frElemData;

  return faceConnPerOrientation;
}

//////////////////////////////////////////////////////////////////////////////

void FluxReconstructionBuilder::setMaxNbStatesInCell()
{
  if (m_maxNbrStatesInCell == 0)
  {
    // get cell-state connectivity
    Common::SafePtr<MeshData::ConnTable> cellStates = MeshDataStack::getActive()->getConnectivity("cellStates_InnerCells");

    const CFuint nbElems = getNbElements();
    for (CFuint iCell = 0; iCell < nbElems; ++iCell)
    {
      const CFuint nbrStatesInCell = cellStates->nbCols(iCell);
      m_maxNbrStatesInCell = m_maxNbrStatesInCell > nbrStatesInCell ? m_maxNbrStatesInCell : nbrStatesInCell;
    }
  }

  // set the maximum number of states in a cell
  MeshDataStack::getActive()->Statistics().setMaxNbStatesInCell(m_maxNbrStatesInCell);
}

//////////////////////////////////////////////////////////////////////////////

void FluxReconstructionBuilder::setMaxNbNodesInCell()
{
  SafePtr<vector<ElementTypeData> > elementType = getCFmeshData().getElementTypeData();

  CFuint maxNbNodes = 0;
  for (CFuint iType = 0; iType < getNbElementTypes(); ++iType)
  {
    const CFuint nbNodes =  (*elementType)[iType].getNbNodes();
    if (nbNodes > maxNbNodes)
    {
      maxNbNodes = nbNodes;
    }
  }

  MeshDataStack::getActive()->Statistics().setMaxNbNodesInCell(maxNbNodes);
}

//////////////////////////////////////////////////////////////////////////////

void FluxReconstructionBuilder::setMaxNbFacesInCell()
{
  SafePtr<vector<ElementTypeData> > elementType = getCFmeshData().getElementTypeData();

  CFuint maxNbFaces = 0;
  for (CFuint iType = 0; iType < getNbElementTypes(); ++iType)
  {
    const CFuint nbFaces = LocalConnectionData::getInstance().getNbFacesInShape((*elementType)[iType].getGeoShape());
    if (nbFaces > maxNbFaces)
    {
      maxNbFaces = nbFaces;
    }
  }

  MeshDataStack::getActive()->Statistics().setMaxNbFacesInCell(maxNbFaces);
}

//////////////////////////////////////////////////////////////////////////////

void FluxReconstructionBuilder::createBoundaryFacesTRS()
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

  // flag telling if the face is a partition face
  m_isPartitionFace.resize(nbBPlusPartitionFaces);
  m_isPartitionFace = true;
  

  for(CFuint iTRS = 0; iTRS < nbTRSs; ++iTRS)
  {
    cf_assert((*nameTRS)[iTRS] != "InnerCells");
    cf_assert((*nameTRS)[iTRS] != "InnerFaces");

    CFLogDebugMed("Nb TRs in TRS: " << (*nbTRs)[iTRS] << "\n");

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
FluxReconstructionBuilder::createTopologicalRegionSet
(const vector<CFuint>& nbFacesPerTR,
 const std::string& name,
 const TRGeoConn& trGeoConn,
 const CFuint iTRS)
{
  CFAUTOTRACE;
  
  CFLog(VERBOSE,"FRBuilder::createTopologicalRegionSet\n");

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
  setNbNodesNbStatesInGeo(nbFacesPerTR, trGeoConn, nbFaceNodes, nbFaceCells);

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

  const bool hasGlobalIDs = (trsGlobalIDs->size() > 0);

  CFuint nbProcessedFaces = 0;
  // let's create the required number of TR's
  for (CFuint iTR = 0; iTR < nbTRs; ++iTR)
  {
    CFLog(VERBOSE, "iTR= " << iTR << ", nbProcessedFaces= " << nbProcessedFaces << "\n");

    // number of faces in this TR
    const CFuint nbFaces = nbFacesPerTR[iTR];

    // get the geoconnectivity (faces --> nodes and --> states (latter is not used))
    const GeoConn& geoConn = trGeoConn[iTR];

    // allocate the TopologicalRegion
    (*storeTR)[iTR] = new TopologicalRegion();
    (*storeTR)[iTR]->setLocalNbGeoEnts(nbFaces);
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
        //CFLog(VERBOSE, "iFace= " << iFace << "\n");

        // face-node connectivity (not really needed during actual run)
        // number of nodes in this face
        const CFuint nbFaceNodes = faceNodes->nbCols(nbProcessedFaces);

        // boolean to check if face has been found
        bool faceFound = false;

        // loop over boundary faces
        CFuint faceIdx;
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

  cf_assert(nbProcessedFaces == totalNbFaces);

  CFLogDebugMin( "FluxReconstructionBuilder::setTopologicalRegions() : nameTRS "
                 << name
                 << ", nb boundary faces detected : "
                 << nbProcessedFaces << "\n");
  return ptrs;
}

////////////////////////////////////////////////////////////////////////////////////////

void FluxReconstructionBuilder::createPartitionBoundaryFacesTRS()
{
  CFAUTOTRACE;

  // get nodes
  DataHandle < Framework::Node*, Framework::GLOBAL > nodes =
	MeshDataStack::getActive()->getNodeDataSocketSink().getDataHandle();

  // get P1 number of face nodes
  SafePtr<vector<ElementTypeData> > elementType = getCFmeshData().getElementTypeData();
  cf_assert(elementType->size() > 0);
  const CFGeoShape::Type cellShape = (*elementType)[0].getGeoShape();
  vector< vector< CFuint > > bFaceOrientations = getBFaceOrientations(cellShape);
  const CFuint nbFaceNodesP1 = bFaceOrientations[0].size();

  // compute the total number of faces in this TRS
  const CFuint totalNbFaces = m_isPartitionFace.size() - getNbBoundaryFaces();

  // fill in two arrays specifying the number of cells and nodes per face in this TR
  std::valarray<CFuint> nbFaceNodes(static_cast<CFuint>(0), totalNbFaces);
  std::valarray<CFuint> nbFaceCells(static_cast<CFuint>(0), totalNbFaces);

  CFuint idx = 0;
  for (CFuint iFace = 0; iFace < m_isPartitionFace.size(); ++iFace)
  {
    if (m_isPartitionFace[iFace])
    {
      nbFaceNodes[idx] = m_bFaceNodes->nbCols(iFace);
      nbFaceCells[idx] = 1;
      ++idx;
    }
  }
  cf_assert(idx == totalNbFaces);

  // allocate connectivity tables
  ConnTable* faceNodes = new ConnTable(nbFaceNodes);
  ConnTable* faceCells = new ConnTable(nbFaceCells);

  // store the connectivities in the MeshData
  /// @todo think about having SharedPtr<ConnTable> and not put
  /// the TRS connectivities in MeshData
  const std::string name = "PartitionFaces";
  MeshDataStack::getActive()->storeConnectivity(name + "Nodes", faceNodes);
  MeshDataStack::getActive()->storeConnectivity("PartitionFaces-Faces2Cells",faceCells);

  // array with the local IDs and types of all the faces in this TRS
  // ownership of this array will be given to the TR
  vector<CFuint>* localFaceIDs = new vector<CFuint>(totalNbFaces);
  vector<CFuint>* faceTypes    = new vector<CFuint>(totalNbFaces);

  // face provider name
  const std::string faceProviderName = "Face";

  // let's create the required number of TR's
  // allocate the TopologicalRegion
  // Create the TopologicalRegion storage
  vector<TopologicalRegion*>* storeTR = new vector<TopologicalRegion*>(1);
  (*storeTR)[0] = new TopologicalRegion();
  (*storeTR)[0]->setLocalNbGeoEnts(totalNbFaces);
  if (totalNbFaces > 0)
  {
    (*storeTR)[0]->setGeoEntsLocalIdx(&(*localFaceIDs)[0], 0);
  }
  else
  {
    (*storeTR)[0]->setGeoEntsLocalIdx(CFNULL, 0);
  }

  // set the connectivity tables of the TRS
  CFuint nbProcessedFaces = 0;
  for (CFuint iBFace = 0; iBFace < m_isPartitionFace.size(); ++iBFace)
  {
    if (m_isPartitionFace[iBFace])
    {
      // face-cell connectivity
      const CFuint cellID = m_bFaceCell[iBFace];
      (*faceCells)(nbProcessedFaces, 0) = cellID;

//      if (cellID >= /*number of local cells*/ )
//      {
//        CFLogDebugMax( "cellID: " << cellID << "\n");
//        throw BadFormatException (FromHere(),"CFmesh had bad state index in GeometricEntity");
//      }

      const CFuint nbFaceNodes = m_bFaceNodes->nbCols(iBFace);
      for (CFuint iNode = 0; iNode < nbFaceNodes; ++iNode)
      {
        const CFuint nodeID =  (*m_bFaceNodes)(iBFace, iNode);
        (*faceNodes)(nbProcessedFaces, iNode) = nodeID;

        if (nodeID >= nodes.size() )
        {
          CFLogDebugMax( "NodeID: " << nodeID << "\n");
          throw BadFormatException (FromHere(),"CFmesh had bad node index in GeometricEntity");
        }
      }

      // assign the local ID of the current geometric entity
      (*localFaceIDs)[nbProcessedFaces] = m_bLocalGeoIDs[iBFace];

      /// @note getGeoShape() assumes a P1 entity. Here the number of nodes is passed for a P1 entity,
      /// even though the entity may be higher order. Another way to fix this is to extend getGeoShape()
      /// by passing also the dimensionality and the geometric order of the entity.
      const CFuint dim = PhysicalModelStack::getActive()->getDim();
      const std::string faceGeoTypeName = makeGeomEntName(getGeoShape(CFGeoEnt::FACE, dim, nbFaceNodesP1),
                                                       getGeometricPolyType(),
                                                       getGeometricPolyOrder(),
                                                       getSolutionPolyType(),
                                                       getSolutionPolyOrder());

      const std::string geoProviderName = faceProviderName + faceGeoTypeName;
      (*faceTypes)[nbProcessedFaces] = m_mapGeoProviderNameToType.find(geoProviderName);

      // increment the counter of partition faces
      ++nbProcessedFaces;
    }
  }

  // check if all faces have been found
  if (nbProcessedFaces != totalNbFaces)
  {
    throw BadFormatException (FromHere(),"Wrong number of partition faces");
  }

  CFLog(NOTICE,"Number of partition faces detected = " << totalNbFaces << "\n");

  // Create TopologicalRegionSet
  TopologicalRegionSet* ptrs = new TopologicalRegionSet (name, storeTR);
  // this is a partition TRS
  ptrs->attachTag("partition");
  // this is TRS of faces
  ptrs->attachTag("face");

  // set local GeometricEntity IDs
  ptrs->setGeoEntsLocalIdx(localFaceIDs);

  // set the connectivity GeometricEntity to nodeIDs
  ptrs->setGeo2NodesConn(faceNodes);

  // there is no face-state connectivity
  // set empty connectivity?
  // create empty face-state connectivity table
  std::valarray<CFuint> nbFaceStates(static_cast<CFuint>(0),totalNbFaces);
  ConnTable* faceStates = new ConnTable(nbFaceStates);
  ptrs->setGeo2StatesConn(faceStates);

  // set the GeometricEntity types
  ptrs->setGeoTypes(faceTypes);

  // create the cached list of state indexes in the all TRS
  /// @todo is this needed here ???
  ptrs->createStatesList();

  // put it into MeshData TRS list
  MeshDataStack::getActive()->addTrs(ptrs);
}

//////////////////////////////////////////////////////////////////////////////

vector< vector < CFuint > > FluxReconstructionBuilder::getBFaceOrientations(const CFGeoShape::Type shape)
{
  CFAUTOTRACE;

  FluxReconstructionElementData* frElemData;
  switch (shape)
  {
    case CFGeoShape::LINE:
    {
      throw Common::NotImplementedException (FromHere(),"Flux Reconstruction has not been implemented for 1D");
    } break;
    case CFGeoShape::QUAD:
    {
      frElemData = new QuadFluxReconstructionElementData(CFPolyOrder::ORDER0);
    } break;
    case CFGeoShape::HEXA:
    {
      frElemData = new HexaFluxReconstructionElementData(CFPolyOrder::ORDER0);
    } break;
    case CFGeoShape::TRIAG:
    {
      frElemData = new TriagFluxReconstructionElementData(CFPolyOrder::ORDER0);
    } break;
    case CFGeoShape::TETRA:
    {
      frElemData = new TetraFluxReconstructionElementData(CFPolyOrder::ORDER0);
    } break;
    default:
    {
      throw Common::ShouldNotBeHereException (FromHere(),"Unsupported cell shape...");
    }
  }

  // get vector containing the nodes connectivities for each orientation
  vector < vector < CFuint > > bFaceOrientations
  = *frElemData->getFaceNodeConn();
   delete frElemData;

  return bFaceOrientations;
}

////////////////////////////////////////////////////////////////////////////////////////

void FluxReconstructionBuilder::reorderBoundaryFacesTRS()
{
  CFAUTOTRACE;

  // get TRS names
  SafePtr< vector<std::string> > nameTRS = getCFmeshData().getNameTRS();

  // get cell-node connectivity
  SafePtr< ConnTable > cellNodeConn = MeshDataStack::getActive()->getConnectivity("cellNodes_InnerCells");

  // get number of TRSs
  const CFuint nbTRSs = nameTRS->size();

  // get possible face orientations and number of face nodes
  SafePtr<vector<ElementTypeData> > elementType = getCFmeshData().getElementTypeData();
  const CFuint nbElemTypes = elementType->size();
  const CFuint nbCellFaces = (*elementType)[0].getNbFaces();
  const CFGeoShape::Type cellShape = (*elementType)[0].getGeoShape();
  vector< vector< CFuint > > bFaceOrientations = getBFaceOrientations(cellShape);
  const CFuint nbFaceNodes = bFaceOrientations[0].size();

  // number of possible orientations
  const CFuint nbOrients = bFaceOrientations.size();
  cf_assert(nbOrients == nbCellFaces);

  // loop over TRSs
  for (CFuint iTRS = 0; iTRS < nbTRSs; ++iTRS)
  {
    // get current boundary TRS
    SafePtr<TopologicalRegionSet> bndTRS = MeshDataStack::getActive()->getTrs((*nameTRS)[iTRS]);

    // get boundary face - cell connectivity
    SafePtr< ConnTable > bFaceCellConn =
      MeshDataStack::getActive()->getConnectivity((*nameTRS)[iTRS]+"-Faces2Cells");

    // get boundary face - node connectivity
    SafePtr< ConnectivityTable< CFuint > > bFaceNodeConn = bndTRS->getGeo2NodesConn();

    // get local (this processor) face IDs
    SafePtr<vector<CFuint> > localFaceIDs = bndTRS->getGeoEntsLocalIdx();

    // get global face IDs
    SafePtr<vector<CFuint> > globalFaceIDs = bndTRS->getGeoEntsGlobalIdx();

    // number of TRs
    const CFuint nbTRs = bndTRS->getNbTRs();

    /// @warning I'm cheating a little here, by using a connectivity table to store
    /// the start indexes of the faces with a certain orientation. Is there a better way?
    std::valarray<CFuint> nbOrientsValAr(1,nbTRs*(nbOrients+1));
    ConnTable* bndFacesStartIdxs = new ConnTable(nbOrientsValAr);

    // get TR list
    SafePtr<vector<TopologicalRegion*> > bndTRList = bndTRS->getTopologicalRegionList();

    // index over groups of faces with a certain orientation
    CFuint iFaceOrient = 0;

    // loop over TRs
    for (CFuint iTR = 0; iTR < nbTRs; ++iTR)
    {
      // number of boundary faces in this TR
      const CFuint nbBFaces = (*bndTRList)[iTR]->getLocalNbGeoEnts();

      // start face index in this TR
      CFuint faceIdx = (*bndTRList)[iTR]->getGeoIDInTrs(0);

      // stop face index in this TR
      const CFuint stopFaceIdx = (*bndTRList)[iTR]->getGeoIDInTrs(nbBFaces);

      // loop over possible orientations
      for (CFuint iOrient = 0; iOrient < nbOrients; ++iOrient)
      {
        // store start index of faces with this orientation
        (*bndFacesStartIdxs)(iFaceOrient,0) = faceIdx; ++iFaceOrient;

        // loop over boundary faces in this TR
        for (CFuint iBFace = faceIdx; iBFace < stopFaceIdx; ++iBFace)
        {
          // neighbouring cell
          const CFuint cellID = (*bFaceCellConn)(iBFace,0);

          // get cell face nodes
          vector< CFuint > cellFaceNodes(nbFaceNodes);
          for (CFuint iNode = 0; iNode < nbFaceNodes; ++iNode)
          {
            // node local ID in cell
            const CFuint nodeLocalID = bFaceOrientations[iOrient][iNode];

            // node local ID in this processor
            cellFaceNodes[iNode] = (*cellNodeConn)(cellID,nodeLocalID);
          }

          // get face nodes
          vector< CFuint > faceNodes(nbFaceNodes);
          for (CFuint iNode = 0; iNode < nbFaceNodes; ++iNode)
          {
            // node local ID in this processor
            faceNodes[iNode] = (*bFaceNodeConn)(iBFace,iNode);
          }

          // check if face has this orientation
          bool hasThisOrient = false;

          // the following depends on whether the mesh is 3D or not
          if (nbFaceNodes > 2)
          {
            // loop over possible rotations of face
            for (CFuint iRot = 0; iRot < nbFaceNodes && !hasThisOrient; ++iRot)
            {
              hasThisOrient = true;

              // check if face has this orientation
              for (CFuint iNode = 0; iNode < nbFaceNodes && hasThisOrient; ++iNode)
              {
                if (cellFaceNodes[iNode] != faceNodes[iNode])
                {
                  hasThisOrient = false;
                  break;
                }
              }

              // rotate face nodes
              if (!hasThisOrient)
              {
                const CFuint swap = faceNodes[0];
                for (CFuint iNode = 0; iNode < nbFaceNodes-1; ++iNode)
                {
                  faceNodes[iNode] = faceNodes[iNode+1];
                }
                faceNodes[nbFaceNodes-1] = swap;
              }
            }

            // if no match was found, try opposite face orientation
            if (!hasThisOrient)
            {
              // change orientation (only two nodes need to be switched)
              CFuint swap = faceNodes[nbFaceNodes-1];
              faceNodes[nbFaceNodes-1] = faceNodes[nbFaceNodes-2];
              faceNodes[nbFaceNodes-2] = swap;

              // loop over possible rotations of face
              for (CFuint iRot = 0; iRot < nbFaceNodes && !hasThisOrient; ++iRot)
              {
                hasThisOrient = true;

                // check if face has this orientation
                for (CFuint iNode = 0; iNode < nbFaceNodes && hasThisOrient; ++iNode)
                {
                  if (cellFaceNodes[iNode] != faceNodes[iNode])
                  {
                    hasThisOrient = false;
                    break;
                  }
                }

                // rotate face nodes
                if (!hasThisOrient)
                {
                  const CFuint swap = faceNodes[0];
                  for (CFuint iNode = 0; iNode < nbFaceNodes-1; ++iNode)
                  {
                    faceNodes[iNode] = faceNodes[iNode+1];
                  }
                  faceNodes[nbFaceNodes-1] = swap;
                }
              }
            }
          }
          else
          {
            hasThisOrient = true;

            // check if face has this orientation
            for (CFuint iNode = 0; iNode < nbFaceNodes && hasThisOrient; ++iNode)
            {
              if (cellFaceNodes[iNode] != faceNodes[iNode])
              {
                hasThisOrient = false;
                break;
              }
            }

            // if no match was found, try opposite face orientation
            if (!hasThisOrient)
            {
              // change orientation (only two nodes need to be switched)
              CFuint swap = faceNodes[1];
              faceNodes[1] = faceNodes[0];
              faceNodes[0] = swap;

              // check if face has this orientation
              hasThisOrient = true;
              for (CFuint iNode = 0; iNode < nbFaceNodes && hasThisOrient; ++iNode)
              {
                if (cellFaceNodes[iNode] != faceNodes[iNode])
                {
                  hasThisOrient = false;
                  break;
                }
              }
            }
          }

          if (hasThisOrient)
          {
            // swap list entries if needed
            if (iBFace != faceIdx)
            {
              // swap boundary face - cell connectivity
              CFuint swap = (*bFaceCellConn)(faceIdx,0);
              (*bFaceCellConn)(faceIdx,0) = (*bFaceCellConn)(iBFace,0);
              (*bFaceCellConn)(iBFace,0) = swap;

              // swap boundary face local (this processor) ID list
              swap = (*localFaceIDs)[faceIdx];
              (*localFaceIDs)[faceIdx] = (*localFaceIDs)[iBFace];
              (*localFaceIDs)[iBFace] = swap;

              // swap boundary face global ID list
              swap = (*globalFaceIDs)[faceIdx];
              (*globalFaceIDs)[faceIdx] = (*globalFaceIDs)[iBFace];
              (*globalFaceIDs)[iBFace] = swap;

              // swap boundary face-geotype list
              /// @warning The face-geotype list cannot be retrieved from the TRS. There should be
              /// only one face type in the mesh though, so it is not necessary to swap this list.
            }

            // find neighbouring cell type
            CFuint cellType;
            for (cellType = 0; cellType < nbElemTypes; ++cellType)
            {
              if ((*elementType)[cellType].getStartIdx() <= cellID &&
                  (*elementType)[cellType].getEndIdx()   >  cellID   )
              {
                break;
              }
            }

            // swap boundary face - node connectivity
            // (this has to be done outside the if, in case the order of the nodes has to be altered)
            const CFuint nbNodes = m_faceNodeElement[cellType]->nbCols(iOrient); // iOrient == local face ID (in a cell)
            for (CFuint iNode = 0; iNode < nbNodes; ++iNode)
            {              
              (*bFaceNodeConn)(iBFace,iNode) = (*bFaceNodeConn)(faceIdx,iNode);          
              const CFuint nodeLocalID = (*m_faceNodeElement[cellType])(iOrient,iNode);
              (*bFaceNodeConn)(faceIdx,iNode) = (*cellNodeConn)(cellID,nodeLocalID);              
            }

            // increase faceIdx
            ++faceIdx;
          }
        }
      }

      // store stop index of faces with this orientation
      cf_assert(faceIdx == stopFaceIdx);
      (*bndFacesStartIdxs)(iFaceOrient,0) = faceIdx; ++iFaceOrient;
    }

    // store start indexes of faces with a certain orientation in this TRS
    MeshDataStack::getActive()->storeConnectivity((*nameTRS)[iTRS]+"boundaryFacesStartIdxs", bndFacesStartIdxs);
  }
}

////////////////////////////////////////////////////////////////////////////////////////

void FluxReconstructionBuilder::reorderPartitionFacesTRS()
{
  CFAUTOTRACE;

  // get cell-node connectivity
  SafePtr< ConnTable > cellNodeConn = MeshDataStack::getActive()->getConnectivity("cellNodes_InnerCells");

  // get possible face orientations and number of face nodes
  SafePtr<vector<ElementTypeData> > elementType = getCFmeshData().getElementTypeData();
  const CFuint nbElemTypes = elementType->size();
  CFLog(VERBOSE, "size elemType: " << nbElemTypes << "\n");
  const CFuint nbCellFaces = (*elementType)[0].getNbFaces();
  const CFGeoShape::Type cellShape = (*elementType)[0].getGeoShape();
  vector< vector< CFuint > > bFaceOrientations = getBFaceOrientations(cellShape);
  const CFuint nbFaceNodes = bFaceOrientations[0].size();

  // number of possible orientations
  const CFuint nbOrients = bFaceOrientations.size();
  cf_assert(nbOrients == nbCellFaces);

    // get current boundary TRS
    SafePtr<TopologicalRegionSet> bndTRS = MeshDataStack::getActive()->getTrs("PartitionFaces");

    // get boundary face - cell connectivity
    SafePtr< ConnTable > bFaceCellConn =
      MeshDataStack::getActive()->getConnectivity("PartitionFaces-Faces2Cells");

    // get boundary face - node connectivity
    SafePtr< ConnectivityTable< CFuint > > bFaceNodeConn = bndTRS->getGeo2NodesConn();

    // get local (this processor) face IDs
    SafePtr<vector<CFuint> > localFaceIDs = bndTRS->getGeoEntsLocalIdx();

    // get global face IDs
    SafePtr<vector<CFuint> > globalFaceIDs = bndTRS->getGeoEntsGlobalIdx();

    // number of TRs
    const CFuint nbTRs = bndTRS->getNbTRs();
    cf_assert(nbTRs == 1);

    /// @warning I'm cheating a little here, by using a connectivity table to store
    /// the start indexes of the faces with a certain orientation. Is there a better way?
    std::valarray<CFuint> nbOrientsValAr(1,nbTRs*(nbOrients+1));
    ConnTable* bndFacesStartIdxs = new ConnTable(nbOrientsValAr);

    // get TR list
    SafePtr<vector<TopologicalRegion*> > bndTRList = bndTRS->getTopologicalRegionList();

    // index over groups of faces with a certain orientation
    CFuint iFaceOrient = 0;

    // loop over TRs
    for (CFuint iTR = 0; iTR < nbTRs; ++iTR)
    {
      // number of boundary faces in this TR
      const CFuint nbBFaces = (*bndTRList)[iTR]->getLocalNbGeoEnts();

      // start face index in this TR
      CFuint faceIdx = (*bndTRList)[iTR]->getGeoIDInTrs(0);

      // stop face index in this TR
      const CFuint stopFaceIdx = (*bndTRList)[iTR]->getGeoIDInTrs(nbBFaces);

      // loop over possible orientations
      for (CFuint iOrient = 0; iOrient < nbOrients; ++iOrient)
      {

        // store start index of faces with this orientation
        (*bndFacesStartIdxs)(iFaceOrient,0) = faceIdx; ++iFaceOrient;

        // loop over boundary faces in this TR
        for (CFuint iBFace = faceIdx; iBFace < stopFaceIdx; ++iBFace)
        {

          // neighbouring cell
          const CFuint cellID = (*bFaceCellConn)(iBFace,0);

          // get cell face nodes
          vector< CFuint > cellFaceNodes(nbFaceNodes);
          for (CFuint iNode = 0; iNode < nbFaceNodes; ++iNode)
          {
            // node local ID in cell
            const CFuint nodeLocalID = bFaceOrientations[iOrient][iNode];

            // node local ID in this processor
            cellFaceNodes[iNode] = (*cellNodeConn)(cellID,nodeLocalID);
          }

          // get face nodes
          vector< CFuint > faceNodes(nbFaceNodes);
          for (CFuint iNode = 0; iNode < nbFaceNodes; ++iNode)
          {
            // node local ID in this processor
            faceNodes[iNode] = (*bFaceNodeConn)(iBFace,iNode);
          }

          // check if face has this orientation
          bool hasThisOrient = false;

          // the following depends on whether the mesh is 3D or not
          if (nbFaceNodes > 2)
          {
            // loop over possible rotations of face
            for (CFuint iRot = 0; iRot < nbFaceNodes && !hasThisOrient; ++iRot)
            {
              hasThisOrient = true;

              // check if face has this orientation
              for (CFuint iNode = 0; iNode < nbFaceNodes && hasThisOrient; ++iNode)
              {
                if (cellFaceNodes[iNode] != faceNodes[iNode])
                {
                  hasThisOrient = false;
                  break;
                }
              }

              // rotate face nodes
              if (!hasThisOrient)
              {
                const CFuint swap = faceNodes[0];
                for (CFuint iNode = 0; iNode < nbFaceNodes-1; ++iNode)
                {
                  faceNodes[iNode] = faceNodes[iNode+1];
                }
                faceNodes[nbFaceNodes-1] = swap;
              }
            }

            // if no match was found, try opposite face orientation
            if (!hasThisOrient)
            {
              // change orientation (only two nodes need to be switched)
              CFuint swap = faceNodes[nbFaceNodes-1];
              faceNodes[nbFaceNodes-1] = faceNodes[nbFaceNodes-2];
              faceNodes[nbFaceNodes-2] = swap;

              // loop over possible rotations of face
              for (CFuint iRot = 0; iRot < nbFaceNodes && !hasThisOrient; ++iRot)
              {
                hasThisOrient = true;

                // check if face has this orientation
                for (CFuint iNode = 0; iNode < nbFaceNodes && hasThisOrient; ++iNode)
                {
                  if (cellFaceNodes[iNode] != faceNodes[iNode])
                  {
                    hasThisOrient = false;
                    break;
                  }
                }

                // rotate face nodes
                if (!hasThisOrient)
                {
                  const CFuint swap = faceNodes[0];
                  for (CFuint iNode = 0; iNode < nbFaceNodes-1; ++iNode)
                  {
                    faceNodes[iNode] = faceNodes[iNode+1];
                  }
                  faceNodes[nbFaceNodes-1] = swap;
                }
              }
            }
          }
          else
          {
            hasThisOrient = true;
            // check if face has this orientation
            for (CFuint iNode = 0; iNode < nbFaceNodes && hasThisOrient; ++iNode)
            {
              if (cellFaceNodes[iNode] != faceNodes[iNode])
              {
                hasThisOrient = false;
                break;
              }
            }
            // if no match was found, try opposite face orientation
            if (!hasThisOrient)
            {
              // change orientation (only two nodes need to be switched)
              CFuint swap = faceNodes[1];
              faceNodes[1] = faceNodes[0];
              faceNodes[0] = swap;

              // check if face has this orientation
              hasThisOrient = true;
              for (CFuint iNode = 0; iNode < nbFaceNodes && hasThisOrient; ++iNode)
              {
                if (cellFaceNodes[iNode] != faceNodes[iNode])
                {
                  hasThisOrient = false;
                  break;
                }
              }
            }
          }

          if (hasThisOrient)
          {
            // swap list entries if needed
            if (iBFace != faceIdx)
            {
              // swap boundary face - cell connectivity
              CFuint swap = (*bFaceCellConn)(faceIdx,0);
              (*bFaceCellConn)(faceIdx,0) = (*bFaceCellConn)(iBFace,0);
              (*bFaceCellConn)(iBFace,0) = swap;

              // swap boundary face local (this processor) ID list
              swap = (*localFaceIDs)[faceIdx];
              (*localFaceIDs)[faceIdx] = (*localFaceIDs)[iBFace];
              (*localFaceIDs)[iBFace] = swap;

              // swap boundary face global ID list
              //swap = (*globalFaceIDs)[faceIdx];
              //(*globalFaceIDs)[faceIdx] = (*globalFaceIDs)[iBFace];
              //(*globalFaceIDs)[iBFace] = swap;

              // swap boundary face-geotype list
              /// @warning The face-geotype list cannot be retrieved from the TRS. There should be
              /// only one face type in the mesh though, so it is not necessary to swap this list.
            }

            // find neighbouring cell type
            CFuint cellType;
            for (cellType = 0; cellType < nbElemTypes; ++cellType)
            {
              if ((*elementType)[cellType].getStartIdx() <= cellID &&
                  (*elementType)[cellType].getEndIdx()   >  cellID   )
              {
                break;
              }
            }

            // swap boundary face - node connectivity
            // (this has to be done outside the if, in case the order of the nodes has to be altered)
            const CFuint nbNodes = m_faceNodeElement[cellType]->nbCols(iOrient); // iOrient == local face ID (in a cell)
            for (CFuint iNode = 0; iNode < nbNodes; ++iNode)
            {
              (*bFaceNodeConn)(iBFace,iNode) = (*bFaceNodeConn)(faceIdx,iNode);
              const CFuint nodeLocalID = (*m_faceNodeElement[cellType])(iOrient,iNode);
              (*bFaceNodeConn)(faceIdx,iNode) = (*cellNodeConn)(cellID,nodeLocalID);
            }

            // increase faceIdx
            ++faceIdx;
          }
        }
      }

      // store stop index of faces with this orientation
      cf_assert(faceIdx == stopFaceIdx);
      (*bndFacesStartIdxs)(iFaceOrient,0) = faceIdx; ++iFaceOrient;
    }

    // store start indexes of faces with a certain orientation in this TRS
    MeshDataStack::getActive()->storeConnectivity("partitionFacesStartIdxs", bndFacesStartIdxs);
}


////////////////////////////////////////////////////////////////////////////////////////

} // namespace FluxReconstructionMethod

} // namespace COOLFluiD

////////////////////////////////////////////////////////////////////////////////////////

