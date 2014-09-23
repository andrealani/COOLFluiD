
#include "Common/BadValueException.hh"
#include "Framework/MethodCommandProvider.hh"
#include "UFEM/UFEM.hh"
#include "UFEM/ElemAssembler.hh"
#include "UFEM/SetupPLaS.hh"

using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

namespace COOLFluiD {
  namespace UFEM {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider< SetupPLaS,UFEMSolverData,UFEMPlugin > setupPLaSProvider("SetupPLaS");

//////////////////////////////////////////////////////////////////////////////

SetupPLaS::SetupPLaS(const std::string& name) :
  UFEMSolverCom(name),
  s_nodes("nodes"),                   // socket sink
  s_faceNeighCell("faceNeighCell"),   // ...
  s_nvfraction("NodalVoidFraction"),  // socket source
  s_nvolume("NodalVolume"),           // ...
  s_evolume("ElementVolumes"),        // ...
  s_ienormals("ElementNormals"),      // ...
  s_benormals("BoundariesNormals")    // ...
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

SetupPLaS::~SetupPLaS()
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

void SetupPLaS::execute()
{
  CFAUTOTRACE;

  // get number of nodes and number of dimensions
  DataHandle< Node*,GLOBAL > h_nodes  = s_nodes.getDataHandle();
  const CFuint nbNodes = h_nodes.size();
  const CFuint nbDim = (nbNodes? h_nodes[0]->size():0);
  cf_always_assert_desc("number of nodes is 0?",nbNodes);
  cf_always_assert_desc("number of dimensions is invalid?",nbDim==2 || nbDim==3);

  // get 'boundary' TRSs and 'inner' TRSs (only one, that's what PLaS supports)
  std::vector< SafePtr< TopologicalRegionSet > > btrs = MeshDataStack::getActive()->getFilteredTrsList("boundary");
  std::vector< SafePtr< TopologicalRegionSet > > itrs = MeshDataStack::getActive()->getFilteredTrsList("inner");
  cf_always_assert_desc("no appropriate 'inner' TRS found",itrs.size());
  CFuint nbElem = itrs[0]->getLocalNbGeoEnts();

  // access GeometricEntities and get element types
  GeometricEntityPool< StdTrsGeoBuilder >& gb = *getMethodData().getStdTrsGeoBuilder();
  std::vector< ElementTypeData >& vetypes = *MeshDataStack::getActive()->getElementTypeData();


  CFLogInfo("SetupPLaS: reset node-wise void fraction...\n");
  DataHandle< CFreal > h_nvfraction = s_nvfraction.getDataHandle();
  h_nvfraction.resize(nbNodes);
  h_nvfraction = 0.;
  CFLogInfo("SetupPLaS: reset node-wise void fraction.\n");


  CFLogInfo("SetupPLaS: reset node and element-wise volume...\n");
  DataHandle< CFreal > h_nvolume = s_nvolume.getDataHandle();
  h_nvolume.resize(nbNodes);
  h_nvolume = 0.;
  DataHandle< CFreal > h_evolume = s_evolume.getDataHandle();
  h_evolume.resize(nbElem);
  h_evolume = 0.;
  CFLogInfo("SetupPLaS: reset node and element-wise volume.\n");


  CFLogInfo("SetupPLaS: reset 'inner' and 'boundary' elements normals...\n");
  DataHandle< std::vector< RealVector > > h_ienormals = s_ienormals.getDataHandle();
  h_ienormals.resize(nbElem);
  for (CFuint i=0; i<nbElem; ++i)
    h_ienormals[i].resize(nbDim);
  DataHandle< std::vector< RealVector > > h_benormals = s_benormals.getDataHandle();
  h_benormals.resize(btrs.size());
  for (CFuint i=0; i<btrs.size(); ++i)
    h_benormals[i].assign(nbDim,RealVector(0.,btrs[i]->getLocalNbGeoEnts()));
  CFLogInfo("SetupPLaS: reset 'inner' and 'boundary' elements normals.\n");


  CFLogInfo("SetupPLaS: set 'inner' volumes and normals...\n");
  for (CFuint i=0; i<vetypes.size(); ++i) {
    //ElemAssembler *ea = getMethodData().getelem_assemblers()[ ElemID(itrs[0]->getName(),i) ];

    for (CFuint idx=0; idx<vetypes[i].getNbElems(); ++idx) {

      // build entity and prepare elemprops
      gb.getDataGE().trs = itrs[0];
      gb.getDataGE().idx = idx;
      GeometricEntity& cell = *gb.buildGE();

      // get element/nodal volume and set DataHandle(s)
      const CFreal evolume = cell.computeVolume();
      const CFreal nvolume = evolume/static_cast<CFreal>(cell.nbNodes());
      h_evolume[ cell.getID() ] = evolume;
      for (CFuint i=0; i<cell.nbNodes(); ++i)
        h_nvolume[ cell.getNode(i)->getLocalID() ] += nvolume;

      // get element normals and allocate/set DataHandle
      const CFuint nbFaces = cell.getNbFacets();
      std::vector<RealVector> avg_face_normals = cell.computeAvgFaceNormals();
      std::vector< RealVector >& normals = h_ienormals[idx];

      normals[XX].resize(nbFaces,0.);
      normals[YY].resize(nbFaces,0.);

      for(CFuint i = 0; i != nbFaces; ++i) {
        normals[XX][i] = avg_face_normals[i][XX];
        normals[YY][i] = avg_face_normals[i][YY];
      }

      if(nbDim == 3) {
        normals[ZZ].resize(nbFaces,0.);
        for(CFuint i = 0; i != nbFaces; ++i) {
          normals[ZZ][i] = avg_face_normals[i][ZZ];
      }

      // release entity
      gb.releaseGE();
    }

  } // for i
  CFLogInfo("SetupPLaS: set 'inner' volume and normals.\n");


  CFLogInfo("SetupPLaS: set 'boundary' elements normals...\n");
  // set faces normals based on element normals, converting face number from
  // coolfluid/gambit to PLaS-style (node-opposite-face)
  DataHandle< std::pair< CFuint,CFuint > > h_faceneighcell = s_faceneighcell.getDataHandle();
  for (CFuint ib=0; ib<btrs.size(); ++ib) {
    for (CFuint idx=0; idx<btrs[ib]->getLocalNbGeoEnts(); ++idx) {
      const CFuint id = btrs[ib]->getLocalGeoID(idx);
      const CFuint ic = h_faceneighcell[id].first;
      int f = (int) h_faceneighcell[id].second;
      f = (nbDim==2?
        (f==0? 2 : (f==1? 0 : (f==2? 1 : -1 ))) :
        (f==0? 3 : (f==1? 2 : (f==2? 0 : (f==3? 1 : -1 )))) );
      cf_assert_desc("problem converting face numbering!",f>=0);
      for (CFuint i=0; i<nbDim; ++i)
        h_benormals[ib][i][idx] = h_ienormals[ic][i][f];
    }
  }
  CFLogInfo("SetupPLaS: set 'boundary' elements normals.\n");
}

//////////////////////////////////////////////////////////////////////////////

void SetupPLaS::setFaceNeighCell()
{
  CFAUTOTRACE;

  vector< SafePtr<TopologicalRegionSet> > trs = MeshDataStack::getActive()->getTrsList();

  DataHandle< Node*,GLOBAL > nodes = socket_nodes.getDataHandle();

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
  DataHandle<pair<CFuint, CFuint> > faceNeighCell =
    socket_faceNeighCell.getDataHandle();
  faceNeighCell.resize(nbBFaces);

  SafePtr<vector<ElementTypeData> > elementType =
    MeshDataStack::getActive()->getElementTypeData();

  SafePtr<TopologicalRegionSet> cells =
    MeshDataStack::getActive()->getTrs("InnerCells");

    const CFuint nbElemTypes = elementType->size();

  // local connectivity face-node for each element type
  vector<Table<CFuint>*> faceNodeElement(nbElemTypes);
  for (CFuint iType = 0; iType < nbElemTypes; ++iType) {

    faceNodeElement[iType] = LocalConnectionData::getInstance().getFaceDofLocal
      ((*elementType)[iType].getGeoShape(),
       static_cast<CFPolyOrder::Type>((*elementType)[iType].getGeoOrder()),
       NODE, CFPolyForm::LAGRANGE);
  }

  // set the face shapes per element type
  vector< vector<CFGeoShape::Type> > faceShapesPerElemType(nbElemTypes);
  LocalConnectionData::getInstance().setFaceShapesPerElemType(faceShapesPerElemType);

  // state list
  DataHandle< State*,GLOBAL > states = socket_states.getDataHandle();

  // atomic number to indicate the maximum possible number
  // of nodes in a face (this allow to avoid frequent reallocations
  // of the vector nodesInFace)
  const CFuint maxNbNodesInFace = 100;
  vector<CFuint> nodesInFace(maxNbNodesInFace);

  // loop over the elements and construct faceIDs
  CFuint elemID = 0;
  typedef CFMultiMap<CFuint, CFuint>::MapIterator MapIterator;

  CFuint countF = 0;

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
          pair<MapIterator, MapIterator> faces = mapBNode2BFace.find(nodeID);

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
      }
    }
  }

  if (countF != nbBFaces) {
    throw Common::NoSuchValueException
      (FromHere(), "Not all the face neighbor cell have been detected correctly!!");
  }
}

//////////////////////////////////////////////////////////////////////////////

  }  // namespace UFEM
}  // namespace COOLFluiD

