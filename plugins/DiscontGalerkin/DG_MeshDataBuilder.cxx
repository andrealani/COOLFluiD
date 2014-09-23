#include <ctime>
#include <set>


#include "Common/SwapEmpty.hh"

#include "Common/CFLog.hh"
#include "Environment/ObjectProvider.hh"

#include "Common/BadValueException.hh"
#include "Framework/MeshData.hh"
#include "Framework/LocalConnectionData.hh"
#include "Framework/GeometricEntityProvider.hh"
#include "Framework/Face.hh"
#include "Framework/Cell.hh"
#include "Framework/PhysicalModel.hh"
#include "Framework/MapGeoToTrsAndIdx.hh"
#include "Framework/CFSide.hh"

#include "DiscontGalerkin/DiscontGalerkin.hh"
#include "DiscontGalerkin/DG_MeshDataBuilder.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {


    namespace DiscontGalerkin {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<DG_MeshDataBuilder,
               MeshDataBuilder,
               DiscontGalerkinModule,
               1>
DG_MeshDataBuilderProvider("DG");

//////////////////////////////////////////////////////////////////////////////

DG_MeshDataBuilder::DG_MeshDataBuilder(const std::string& name) :
  MeshDataBuilder(name),
  m_nbFaces(0),
  m_inGeoTypes(CFNULL),
  m_inLocalGeoIDs(CFNULL),
  m_bGeoType(),
  m_bLocalGeoIDs(),
  m_localFace2Node(),
  m_nbFacesPerElem(),
  m_isBFace(),
  m_geoTypeIDs(),
  m_nbInFacesNodes(),
  m_nbBFacesNodes(),
  m_bFaceCell(),
  m_bFaceNodes(CFNULL)
{
}

//////////////////////////////////////////////////////////////////////////////

DG_MeshDataBuilder::~DG_MeshDataBuilder()
{
}

//////////////////////////////////////////////////////////////////////////////

void DG_MeshDataBuilder::releaseMemory()
{
  MeshDataBuilder::releaseMemory();

  SwapEmpty(m_bGeoType);
  SwapEmpty(m_bLocalGeoIDs);
  SwapEmpty(m_localFace2Node);
  m_nbFacesPerElem.resize(0);
  SwapEmpty(m_isBFace);
  SwapEmpty(m_geoTypeIDs);
  m_nbInFacesNodes.resize(0);
  m_nbBFacesNodes.resize(0);
  deletePtr(m_bFaceNodes);
}

//////////////////////////////////////////////////////////////////////////////

void DG_MeshDataBuilder::createTopologicalRegionSets()
{
  CFAUTOTRACE;

  // first create the cells
  createInnerCells();

  // put the coordinates in the states
  setCoordInCellStates();

  createBoundaryTRS();

  // create the cell-face connectivity
  createCellFaces();

  // create the inner faces TRS
  createInnerFacesTRS();

  setMapGeoToTrs();
}

//////////////////////////////////////////////////////////////////////////////

void DG_MeshDataBuilder::setMaxNbStatesInCell()
{
  Common::SafePtr<TopologicalRegionSet> cells =
    MeshDataStack::getActive()->getTrs("InnerCells");

  CFuint maxNbStates = 0;
  const CFuint nbGeos = cells->getLocalNbGeoEnts();
  for (CFuint iGeo = 0; iGeo < nbGeos; ++iGeo) {
    if (cells->getNbStatesInGeo(iGeo) > maxNbStates) {
      maxNbStates = cells->getNbStatesInGeo(iGeo);
    }
  }

  MeshDataStack::getActive()->Statistics().setMaxNbStatesInCell(maxNbStates);
}

//////////////////////////////////////////////////////////////////////////////

void DG_MeshDataBuilder::setMaxNbNodesInCell()
{
  Common::SafePtr<TopologicalRegionSet> cells =
    MeshDataStack::getActive()->getTrs("InnerCells");

  CFuint maxNbNodes = 0;
  const CFuint nbGeos = cells->getLocalNbGeoEnts();
  for (CFuint iGeo = 0; iGeo < nbGeos; ++iGeo) {
    if (cells->getNbNodesInGeo(iGeo) > maxNbNodes) {
      maxNbNodes = cells->getNbNodesInGeo(iGeo);
    }
  }

  MeshDataStack::getActive()->Statistics().setMaxNbNodesInCell(maxNbNodes);
}

//////////////////////////////////////////////////////////////////////////////

void DG_MeshDataBuilder::setMaxNbFacesInCell()
{
  SafePtr<vector<ElementTypeData> > elementType =
    getCFmeshData().getElementTypeData();

  CFuint maxNbFaces = 0;
  for (CFuint iType = 0; iType < getNbElementTypes(); ++iType) {
    CFuint nbFaces = LocalConnectionData::getInstance().getNbFacesInShape
      ((*elementType)[iType].getGeoShape());
    if (nbFaces > maxNbFaces) {
      maxNbFaces = nbFaces;
    }
  }

  MeshDataStack::getActive()->Statistics().setMaxNbFacesInCell(maxNbFaces);
}

//////////////////////////////////////////////////////////////////////////////

void DG_MeshDataBuilder::setMapGeoToTrs()
{
  CFAUTOTRACE;
  // set the number of faces which is the number of boundary faces in MeshData
  MapGeoToTrsAndIdx* mapGeoToTrs = new MapGeoToTrsAndIdx();
  mapGeoToTrs->resize(MeshDataStack::getActive()->Statistics().getNbFaces());

  vector<Common::SafePtr<TopologicalRegionSet> > trs = MeshDataStack::getActive()->getTrsList();
  for (CFuint i = 0; i < trs.size(); ++i)
  {
    SafePtr<TopologicalRegionSet> currTrs = trs[i];
    if (currTrs->hasTag("face"))
    {
      const bool isOnBoundary = (currTrs->hasTag("partition") || currTrs->hasTag("boundary"));
      const CFuint nbTrsFaces = currTrs->getLocalNbGeoEnts();
      for (CFuint iFace = 0; iFace < nbTrsFaces; ++iFace)
      {
        const CFuint faceID = currTrs->getLocalGeoID(iFace);
        mapGeoToTrs->setMappingData(faceID, currTrs, iFace, isOnBoundary);
      }
    }
  }

  MeshDataStack::getActive()->storeMapGeoToTrs("MapFacesToTrs", mapGeoToTrs);
}

//////////////////////////////////////////////////////////////////////////////

void DG_MeshDataBuilder::createCellFaces()
{
  CFAUTOTRACE;

  const CFuint nbElem      = getNbElements();
  const CFuint nbElemTypes = getNbElementTypes();

  // local connectivity face-node for each element type
  // it is a vector of pointers to tables containing CFuints
  m_localFace2Node.resize(nbElemTypes);

  // get the elementTypeData
  SafePtr<vector<ElementTypeData> > elementType = getCFmeshData().getElementTypeData();

  cf_assert(nbElemTypes == elementType->size());

  for (CFuint iType = 0; iType < nbElemTypes; ++iType)
  {
    m_localFace2Node[iType] = LocalConnectionData::getInstance().getFaceDofLocal
      ((*elementType)[iType].getGeoShape(),
       getGeometricPolyOrder(),
       NODE,
       CFPolyForm::LAGRANGE);
  }

  // array storing the number of faces per element
  m_nbFacesPerElem.resize(nbElem);
  m_nbFacesPerElem = 0;
  LocalConnectionData::getInstance().setNbFacesPerElement(m_nbFacesPerElem); // pass m_nbFacesPerElem by reference

  // estimate the maximum number of faces
  const CFuint sumNbFacesPerElem = m_nbFacesPerElem.sum();
  const CFuint maxTotalNbFaces = sumNbFacesPerElem;

  // allocate a table mapping nodes to face ID
  const CFuint totNbNodes = MeshDataStack::getActive()->getNbNodes();
  vector < vector<CFuint> > mapNode2Face(totNbNodes);

  // place the already created boundary faces inside
  // by looping on the boundary TRSs
  vector<Common::SafePtr<TopologicalRegionSet> > trs = MeshDataStack::getActive()->getTrsList();
  for (CFuint i = 0; i < trs.size(); ++i)
  {
    SafePtr<TopologicalRegionSet> currTrs = trs[i];
    if (currTrs->hasTag("face") && currTrs->hasTag("boundary"))
    {
      const CFuint nbTrsFaces = currTrs->getLocalNbGeoEnts();
      for (CFuint iFace = 0; iFace < nbTrsFaces; ++iFace)
      {
        const CFuint faceID = currTrs->getLocalGeoID(iFace);
        const CFuint nbNodes = currTrs->getNbNodesInGeo(iFace);
        for (CFuint iNode = 0; iNode < nbNodes; ++iNode)
        {
          mapNode2Face[currTrs->getNodeID(iFace,iNode)].push_back(faceID);
        }
      }
    }
  }

  // continue from the last index attributed to a face
  const CFuint nbBCFaces = MeshDataStack::getActive()->Statistics().getNbFaces();
  m_nbFaces = nbBCFaces;
  // nb constructed faces until now should be equal to the nb of faces
  // in the boundary read from the file
  cf_assert( nbBCFaces == getNbBoundaryFaces() );


  // set the face shapes per element type
  vector< vector<CFGeoShape::Type> > faceShapesPerElemType(nbElemTypes);
  LocalConnectionData::getInstance().setFaceShapesPerElemType(faceShapesPerElemType); // pass faceShapesPerElemType by reference

  // table storing the connectivity element-local faceID
  // (in this processor)
  ConnTable* cellFaces = new ConnTable(m_nbFacesPerElem);
  // store the connectivity in MeshData
  MeshDataStack::getActive()->storeConnectivity("cellFaces", cellFaces);

  // atomic number to indicate the maximum possible number
  // of nodes in a face which allows to avoid frequent reallocations
  // of the vector nodesInFace
  const CFuint maxNbNodesInFace = 100;
  vector<CFuint> nodesInFace(maxNbNodesInFace);
  const std::string faceProviderName = "Face";

  // the following arrays are oversized
  m_geoTypeIDs.resize(maxTotalNbFaces);

  // loop over the elements and construct faceIDs
  CFuint elemID = 0;
  CFuint countInFaces = 0;
  std::string providerName = "";

  std::vector< std::pair<CFuint,CFuint> > NbNodesInFace;

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
    for (CFuint iElem = 0; iElem < nbElemPerType; ++iElem, ++elemID)
    {
// CFout << "\n" << iElem;

      const CFuint nbElemFaces = m_nbFacesPerElem[elemID];
      // loop over the faces in the current element
      for (CFuint iFace = 0; iFace < nbElemFaces; ++iFace)
      {
// CFout << "  " << iFace;
        const CFuint nbNodesPerFace = m_localFace2Node[iType]->nbCols(iFace);

        // construct sets of nodes that make the corresponding face
        // in this element
        for (CFuint iNode = 0; iNode < nbNodesPerFace; ++iNode)
        {
          const CFuint localNodeID = (*m_localFace2Node[iType])(iFace, iNode);
          nodesInFace[iNode] = getCFmeshData().getElementNode(elemID, localNodeID);
        }
// CFout << "  " << nodesInFace[0] << "  " << nodesInFace[1] << "  " << nodesInFace[2];
        // consider the first node belonging to the current face
        // check if you find a face ID shared between all the other
        // nodes building a face
        const CFuint nodeID = nodesInFace[0];

        // number of faceIDs referencing node
        // if a face with this node has never been added,
        // then this will be zero and the face is immedietly marked as not found
        const CFuint nbFaceIDsRefNode = mapNode2Face[nodeID].size();

        bool faceFound = false;

        // loop over all the faceIDs referenced by the first node to see if
        // all the other nodes reference the same face
        CFuint countNodes = 0;
        CFuint currFaceID = 0;
        for (CFuint iFaceID = 0;
             iFaceID < nbFaceIDsRefNode && (faceFound == false);
             ++iFaceID)
        {
          currFaceID = mapNode2Face[nodeID][iFaceID];
          countNodes = 1;
          for (CFuint iNode = 1; iNode < nbNodesPerFace; ++iNode)
          {
            const CFuint newNodeID = nodesInFace[iNode];
            // number of faceIDs referencing node
            const CFuint nbFaceIDsRefThisNode = mapNode2Face[newNodeID].size();
            for (CFuint jFaceID = 0; jFaceID < nbFaceIDsRefThisNode; ++jFaceID)
            {
              if (mapNode2Face[newNodeID][jFaceID] == currFaceID)
              {
                ++countNodes;
                // break from the loop over faces referencing the new node
                break;
              }
            }
          }

          if (countNodes == nbNodesPerFace)
          {
            // the corresponding faceID already exists
            faceFound = true;
            (*cellFaces)(elemID, iFace) = currFaceID;

// if (currFaceID == 14286) CFout << "**************rightEL: " << elemID << "\n";
// CFout << "(" << currFaceID << "," << nbNodesPerFace << ")";
            NbNodesInFace.push_back(make_pair(currFaceID,nbNodesPerFace));
            // break from loop over faces referencing the first node
            break;
          }
        }
// CFout << " x" << faceFound;
        if (!faceFound)
        {
          // a new face has been found
          // add the ID of the new face in the corresponding nodes referencing it

          countInFaces++;
// if (m_nbFaces == 14286) CFout << "**************leftEL: " << elemID << "\n";
          for (CFuint iNode = 0; iNode < nbNodesPerFace; ++iNode)
          {
            const CFuint currNodeID = nodesInFace[iNode];
            mapNode2Face[currNodeID].push_back(m_nbFaces);
          }

          // store the geometric entity type for the current face
          m_geoTypeIDs[m_nbFaces] = faceGeoTypeID[iFace];
          NbNodesInFace.push_back(make_pair(m_nbFaces,nbNodesPerFace));
// CFout << "(" << m_nbFaces << "," << nbNodesPerFace << ")";
          (*cellFaces)(elemID, iFace) = m_nbFaces;
          // increment the number of faces
          m_nbFaces++;

        }
      }
    }
  }
// for(std::vector< std::pair<CFuint,CFuint> >::iterator el=NbNodesInFace.begin();el!=NbNodesInFace.end();el++)
//   {
//     CFout << "\n" << el->first << " " << el->second << CFendl;
//   }
  // remove duplicated entries
  // because some faces will be counted twice
  std::sort(NbNodesInFace.begin(), NbNodesInFace.end(), LessThan());
// CFout << "\n";
// for(std::vector< std::pair<CFuint,CFuint> >::iterator el=NbNodesInFace.begin();el!=NbNodesInFace.end();el++)
//   {
//     CFout << "\n" << el->first << " " << el->second << CFendl;
//   }

  NbNodesInFace.erase(std::unique(NbNodesInFace.begin(),NbNodesInFace.end(), Equal()),NbNodesInFace.end());
  cf_assert(m_nbFaces <= maxTotalNbFaces);
  cf_assert(countInFaces <= maxTotalNbFaces);

  const CFuint totalNbFaces = m_nbFaces;
  const CFuint nbInnerFaces = countInFaces;

  MeshDataStack::getActive()->Statistics().setNbFaces(totalNbFaces);

  CFLog(INFO, "Total Nb Faces : " << totalNbFaces << "\n");
  CFLog(INFO, "Nb Inner Faces : " << nbInnerFaces << "\n");
  CFLog(INFO, "Nb Bound Faces : " << nbBCFaces << "\n");
  cf_assert(NbNodesInFace.size() == totalNbFaces);

  // by default faces are not in the boundary
  m_isBFace.resize(totalNbFaces);
  for (CFuint i = 0; i < totalNbFaces; ++i)
  {
    m_isBFace[i] = false;
  }
  // mark which faces are in the boundary
  for (CFuint i = 0; i < trs.size(); ++i)
  {
    SafePtr<TopologicalRegionSet> currTrs = trs[i];
    if (currTrs->hasTag("face") && currTrs->hasTag("boundary"))
    {
      const CFuint nbTrsFaces = currTrs->getLocalNbGeoEnts();
      for (CFuint iFace = 0; iFace < nbTrsFaces; ++iFace)
      {
        const CFuint faceID = currTrs->getLocalGeoID(iFace);
        m_isBFace[faceID] = true;
      }
    }
  }

  m_nbInFacesNodes.resize(nbInnerFaces);
  m_nbBFacesNodes.resize(nbBCFaces);

  // set the number of nodes in faces
  CFuint iBFace = 0;
  CFuint iInFace = 0;

  for (CFuint i = 0; i < NbNodesInFace.size(); ++i)
  {
    if (!m_isBFace[i])
    {
      m_nbInFacesNodes[iInFace++] = NbNodesInFace[i].second;
    }
    else
    {
      m_nbBFacesNodes[iBFace++] = NbNodesInFace[i].second;
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void DG_MeshDataBuilder::createInnerFacesTRS()
{
  CFAUTOTRACE;

  const CFuint nbInnerFaces = m_nbInFacesNodes.size();
  const CFuint nbBCFaces    = m_nbBFacesNodes.size();
  const CFuint totalNbFaces = nbInnerFaces + nbBCFaces;

  // Create connectivity table inner face - nodes
  ConnTable* innerFaceNodes = new ConnTable(m_nbInFacesNodes);

  // Create connectivity table boundary face - nodes
  m_bFaceNodes = new ConnTable(m_nbBFacesNodes);

  // Create connectivity table inner face - cell
  std::valarray<CFuint> nbInFacesCells(2,nbInnerFaces);
  ConnTable* innerFaceCells = new ConnTable(nbInFacesCells);

  // Create connectivity boundary face - cell
  m_bFaceCell.resize(totalNbFaces);
  for (CFuint i = 0 ; i< totalNbFaces ; i++)
    m_bFaceCell[i]=-1;

  // Get the cell - face connectivity from the MeshData
  SafePtr<ConnTable> cellFaces = MeshDataStack::getActive()->getConnectivity("cellFaces");

  // Get the ElemenTypeData from CFMeshData
  SafePtr<vector<ElementTypeData> > elementType = getCFmeshData().getElementTypeData();

  // std::valarray containing the face local index in InnerFaces TRS
  std::valarray<CFint> idxFace(-1, totalNbFaces);

  // Connectivity inner face idx - geotype
  m_inGeoTypes = new vector<CFuint>(nbInnerFaces);
  // Connectivity inner face idx - ID in local processor
  m_inLocalGeoIDs = new vector<CFuint>(nbInnerFaces);
  // Connectivity boundary face local idx - geotype
  m_bGeoType.resize(nbBCFaces);
  // Connectivity boundary face local idx - ID in local processor
  m_bLocalGeoIDs.resize(nbBCFaces);

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

        const CFuint nbNodesPerFace = m_localFace2Node[iType]->nbCols(iFace);

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
              const CFuint localNodeID = (*m_localFace2Node[iType])(iFace, iNode);

              const CFuint nodeID = getCFmeshData().getElementNode(elemID, localNodeID);

              (*innerFaceNodes)(countInFaces, iNode) = nodeID;
            }

            // set the local geo ID and the geometric entity type
            (*m_inLocalGeoIDs)[countInFaces] = faceID;
            (*m_inGeoTypes)[countInFaces] = m_geoTypeIDs[faceID];

            // set the first neighbouring element
            (*innerFaceCells)(countInFaces, LEFT) = elemID;
// if (faceID == 14286) CFout << "\n**************ElLeft: "<< elemID << "\n" << CFendl;
            // set the local idx in the inner faces TRS for this face
            idxFace[faceID] = static_cast<CFint>(countInFaces);

          }
          else
          {
            // set the second neighbouring element in the inner face
            cf_assert(idxFace[faceID] != -1);
            (*innerFaceCells)(idxFace[faceID], RIGHT) = elemID;
// if (faceID == 14286) CFout << "\n**************ElRight: "<< elemID << "\n" << CFendl;
          }
        }
        else
        {
          ++countBFaces;
          for (CFuint iNode = 0; iNode < nbNodesPerFace; ++iNode)
          {
            const CFuint localNodeID = (*m_localFace2Node[iType])(iFace, iNode);
            const CFuint nodeID = getCFmeshData().getElementNode(elemID, localNodeID);
            (*m_bFaceNodes)(countBFaces, iNode) = nodeID;
          }

          // set the local geo ID and the geometric entity type
          m_bLocalGeoIDs[countBFaces] = faceID;
          m_bGeoType[countBFaces] = m_geoTypeIDs[faceID];
          // set the neighbouring element
          m_bFaceCell[faceID] = elemID;
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

  CFLog(NOTICE,"Created DG TRS InnerFaces with "
        << innerFaceTRS->getLocalNbGeoEnts() << " local faces\n");


// Add conectivity boundary faces to cell
//  const CFuint nbTRSs = getCFmeshData().getNbTRSs();
  std::vector< Common::SafePtr<TopologicalRegionSet> > m_TrsList = MeshDataStack::getActive()->getTrsList();

  const CFuint number_TRS = m_TrsList.size();

  for(CFuint iTRS = 0; iTRS < number_TRS; ++iTRS)
  {
    if ((*m_TrsList[iTRS]).hasTag("boundary"))
    {
      Common::SafePtr<std::vector<CFuint> > bFaces = m_TrsList[iTRS]->getGeoEntsLocalIdx();

      CFuint numberBoundaryFaces = bFaces->size();

      // Create connectivity table boundary face - cell
      std::valarray<CFuint> nbCellFace(1,numberBoundaryFaces);
      ConnTable* boundaryCell2Face = new ConnTable(nbCellFace);

      for(CFuint iFace=0; iFace < numberBoundaryFaces; iFace++)
      {
        cf_assert( m_bFaceCell[(*bFaces)[iFace]] != -1);
        (*boundaryCell2Face)(iFace, LEFT) = m_bFaceCell[(*bFaces)[iFace]];
      }

      MeshDataStack::getActive()->storeConnectivity((*m_TrsList[iTRS]).getName()+"-Faces2Cells", boundaryCell2Face);
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace DiscontGalerkin


} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
