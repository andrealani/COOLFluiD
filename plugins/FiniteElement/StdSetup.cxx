#include "Common/CFMultiMap.hh"
#include "MathTools/RealMatrix.hh"
#include "Framework/MeshData.hh"

#include "Framework/MapGeoToTrsAndIdx.hh"
#include "Framework/MethodCommandProvider.hh"

#include "FiniteElement/FiniteElement.hh"
#include "FiniteElement/StdSetup.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteElement {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<StdSetup, FiniteElementMethodData, FiniteElementModule> stdSetupProvider("StdSetup");

//////////////////////////////////////////////////////////////////////////////

StdSetup::StdSetup(const std::string& name) :
  FiniteElementMethodCom(name),
  _nodeIdToStateId(),
  socket_nstatesProxy("nstatesProxy"),
  socket_isUpdated("isUpdated"),
  socket_states("states"),
  socket_nodes("nodes"),
  socket_faceNeighCell("faceNeighCell"),
  socket_appliedStrongBC("appliedStrongBC")
{
}

//////////////////////////////////////////////////////////////////////////////

StdSetup::~StdSetup()
{
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSource> >
StdSetup::providesSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSource> > result;

  result.push_back(&socket_nstatesProxy);
  result.push_back(&socket_isUpdated);
  result.push_back(&socket_faceNeighCell);
  result.push_back(&socket_appliedStrongBC);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

void StdSetup::execute()
{
  CFAUTOTRACE;
  
  DataHandle <Framework::State*, Framework::GLOBAL> states  = socket_states.getDataHandle();
  DataHandle <Framework::Node*, Framework::GLOBAL> nodes  = socket_nodes.getDataHandle();
  const CFuint nbStates = states.size();
  const CFuint nbNodes = nodes.size();

  DataHandle<ProxyDofIterator<RealVector>*> nstatesProxy = socket_nstatesProxy.getDataHandle();
  nstatesProxy.resize(1);

  // set the mapping from nodeID to stateID
  _nodeIdToStateId.resize(nbNodes);
  for (CFuint i = 0; i < nbStates; ++i) {
    const bool indexed = states[i]->getCoordinates().isIndexed();
    if(indexed)
    {
      const CFuint nodeID = states[i]->getCoordinates().getLocalID();
      cf_assert(nodeID < nbNodes);
      _nodeIdToStateId[nodeID] = i;
    }
  }
  nstatesProxy[0] =
    new DofDataHandleIterator<RealVector, State, GLOBAL>(states, &_nodeIdToStateId);
  
  DataHandle<bool> isUpdated = socket_isUpdated.getDataHandle();
  isUpdated.resize(nbStates);
  isUpdated = false;

  DataHandle<vector<bool> > appliedStrongBC = socket_appliedStrongBC.getDataHandle();
  appliedStrongBC.resize(nbStates);

  setFaceNeighCell();
}

//////////////////////////////////////////////////////////////////////////////

void StdSetup::setFaceNeighCell()
{
  CFAUTOTRACE;

  vector< SafePtr<TopologicalRegionSet> > trs = MeshDataStack::getActive()->getTrsList();

  DataHandle < Framework::Node*, Framework::GLOBAL > nodes = socket_nodes.getDataHandle();

  const CFuint totNbNodes = nodes.size();
  vector<bool> isBNode(totNbNodes, false);

  //count the number of boundary faces
  CFuint nbBFaces = 0;
  // create mapping between boundary nodes and faces
  CFMultiMap<CFuint, CFuint> mapBNode2BFace;

  for (CFuint i = 0; i < trs.size(); ++i) {
    SafePtr<TopologicalRegionSet> currTrs = trs[i];
    if (!currTrs->hasTag("inner")) {
      const CFuint nbTrsFaces = currTrs->getLocalNbGeoEnts();
      for (CFuint iFace = 0; iFace < nbTrsFaces; ++iFace) {
        const CFuint faceID = currTrs->getLocalGeoID(iFace);
        const CFuint nbNodesInFace = currTrs->getNbNodesInGeo(iFace);
        for (CFuint iNode = 0; iNode < nbNodesInFace; ++iNode) {
          const CFuint nodeID = currTrs->getNodeID(iFace,iNode);
          mapBNode2BFace.insert(nodeID, faceID);
          // flag the boundary nodes
          cf_assert(nodeID < isBNode.size());
          isBNode[nodeID] = true;
        }
      }
      nbBFaces += nbTrsFaces;
    }
  }
  mapBNode2BFace.sortKeys();

  Common::SafePtr<MapGeoToTrsAndIdx> mapGeoToTrs =
    MeshDataStack::getActive()->getMapGeoToTrs("MapFacesToTrs");

  // face neighbor cell
  DataHandle<Common::Trio<CFuint, CFuint, SafePtr<TopologicalRegionSet> > > faceNeighCell =
    socket_faceNeighCell.getDataHandle();
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

      for (CFuint iElem = 0; iElem < nbElemPerType; ++iElem, ++elemID) {
        // loop over the faces in the current element
        for (CFuint iFace = 0; iFace < nbElemFaces; ++iFace) {
          const CFuint nbNodesPerFace = faceNodeElement[iType]->nbCols(iFace);

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

          // face is on the boundary if ALL its nodes are on the boundary
          if (countFaceNodesOnBoundary == nbNodesPerFace) {
            // consider the first node belonging to the current face
            // check if you find a face ID shared between all the other
            // nodes building a face
            const CFuint nodeID = nodesInFace[0];
	    bool fo = false;
            pair<MapIterator, MapIterator> faces = mapBNode2BFace.find(nodeID, fo);
	    
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
          } //  end of if(countFaceNodesOnBoundary == nbNodesPerFace)
        }  //end of loop over face of element
      }  //end of loop over elements
    } // en dof loop over types
  } //end of loop over inner trs

  if (countF != nbBFaces) {
    throw Common::NoSuchValueException
      (FromHere(), "Not all the face neighbor cell have been detected correctly!");
  }
}


//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
StdSetup::needsSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSink> > result;

  result.push_back(&socket_states);
  result.push_back(&socket_nodes);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

 } // namespace FiniteElement

  } // namespace Numerics

} // namespace COOLFluiD
