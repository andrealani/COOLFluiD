#include <numeric>

#include "StegerWarmingCorrFluxMPI.hh"
#include "Framework/MethodStrategyProvider.hh"
#include "Framework/GeometricEntity.hh"
#include "FiniteVolumeNavierStokes/FiniteVolumeNavierStokes.hh"

#include "Common/PE.hh"
#include "Common/MPI/MPIStructDef.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

MethodStrategyProvider<StegerWarmingCorrFluxMPI,
		       CellCenterFVMData,
                       FluxSplitter<CellCenterFVMData>,
		       FiniteVolumeNavierStokesModule>
stegerWarmingCorrFluxMPIProvider("StegerWarmingCorrMPI");

//////////////////////////////////////////////////////////////////////////////

StegerWarmingCorrFluxMPI::StegerWarmingCorrFluxMPI(const std::string& name) :
  StegerWarmingCorrFlux(name),
  _faceBuilder(),
  _cellBuilder()
{
}

//////////////////////////////////////////////////////////////////////////////

StegerWarmingCorrFluxMPI::~StegerWarmingCorrFluxMPI()
{
}

//////////////////////////////////////////////////////////////////////////////

void StegerWarmingCorrFluxMPI::buildFaceBCData()
{
  setStencil();

  SafePtr<TopologicalRegionSet> pFaces =
    MeshDataStack::getActive()->getTrs(std::string("PartitionFaces"));
  SafePtr<TopologicalRegionSet> inFaces =
    MeshDataStack::getActive()->getTrs(std::string("InnerFaces"));

  const CFuint totNbFaces =
    MeshDataStack::getActive()->Statistics().getNbFaces();

  // identify the partition faces by ID
  vector<bool> isPartitionFace(totNbFaces);
  isPartitionFace.assign(totNbFaces, false);
  const CFuint nbLocalPFaces = pFaces->getLocalNbGeoEnts();
  
  DataHandle<Node*, GLOBAL> nodes = socket_nodes.getDataHandle();
  CFMultiMap<CFuint, CFuint> mapGlobalIDToLocalFaceID;

  // flag the partition faces
  for (CFuint iFace = 0; iFace < nbLocalPFaces; ++iFace) {
    const CFuint localGeoID = pFaces->getLocalGeoID(iFace);
    isPartitionFace[localGeoID] = true;
  }

  // build a mapping between global node IDs and local IDs of
  // corresponding inner faces
  const CFuint nbLocalInFaces = inFaces->getLocalNbGeoEnts();
  for (CFuint iFace = 0; iFace < nbLocalInFaces; ++iFace) {
    const CFuint nbNodesInFace = inFaces->getNbNodesInGeo(iFace);
    for (CFuint nf = 0; nf < nbNodesInFace; ++nf) {
      const CFuint globalNodeID = nodes[inFaces->getNodeID(iFace, nf)]->getGlobalID();
      mapGlobalIDToLocalFaceID.insert(globalNodeID,iFace);
    }
  }
  mapGlobalIDToLocalFaceID.sortKeys();

  // initialize the flags for the faces in direction normal to the wall
  // I.e. parallel to the wall
  _flagNormalFace.resize(totNbFaces);
  _flagNormalFace.assign(_flagNormalFace.size(), false);

  // locally built geo builders. At this stage there could be
  // something missing that doesn't let use the ones owned by
  // the method data it is safer to use local ones
  _faceBuilder.setup();
  _faceBuilder.getGeoBuilder()->setDataSockets(socket_states, socket_gstates, socket_nodes);

  SafePtr<TopologicalRegionSet> cells = MeshDataStack::getActive()->getTrs("InnerCells");

  _cellBuilder.setup();
  _cellBuilder.getGeoBuilder()->setDataSockets(socket_states, socket_gstates, socket_nodes);
  CellTrsGeoBuilder::GeoData& cellData = _cellBuilder.getDataGE();
  cellData.trs = cells;

  vector<bool> faceToConsider;
  vector<CFuint> nbLayersToAdd;
  vector<CFint> globalCellIDs;

  // at first, wall TRS are processed
  processBFaces(_wallTrsNames,
    mapGlobalIDToLocalFaceID,
    isPartitionFace,
    faceToConsider,
    nbLayersToAdd,
    globalCellIDs);

  vector<std::string> pName(1, std::string("InnerFaces"));
  CFuint localNbLayersToAdd = accumulate(nbLayersToAdd.begin(), nbLayersToAdd.end(),0);
  CFuint globalNbLayersToAdd = 0;

  CFLogDebugMin("P" << PE::GetPE().GetRank()
    << " has globalNbLayersToAdd = " << globalNbLayersToAdd << "\n");

  MPI_Allreduce(&localNbLayersToAdd, &globalNbLayersToAdd,
    1, MPIStructDef::getMPIType(&localNbLayersToAdd),MPI_MAX,PE::GetPE().GetCommunicator());

  while (globalNbLayersToAdd > 0) {
    processBFaces(pName,
      mapGlobalIDToLocalFaceID,
      isPartitionFace,
      faceToConsider,
      nbLayersToAdd,
      globalCellIDs);

    localNbLayersToAdd = accumulate(nbLayersToAdd.begin(), nbLayersToAdd.end(),0);

    MPI_Allreduce(&localNbLayersToAdd, &globalNbLayersToAdd,
      1, MPIStructDef::getMPIType(&localNbLayersToAdd),MPI_MAX,PE::GetPE().GetCommunicator());
  }
}

//////////////////////////////////////////////////////////////////////////////

void StegerWarmingCorrFluxMPI::processBFaces
(const vector<std::string>& trsNames,
 CFMultiMap<CFuint,CFuint>& mapGlobalIDToLocalFaceID,
 const vector<bool>& isPartitionFace,
 vector<bool>& faceToConsider,
 vector<CFuint>& nbLayersToAdd,
 vector<CFint>& globalCellIDs)
{
  typedef CFMultiMap<CFuint, CFuint>::MapIterator MapIterator;
  DataHandle<Node*, GLOBAL> nodes = socket_nodes.getDataHandle();
  
  PartitionFaceData sendData;
  // in the first loop only wall faces are processed
  const bool isFirstLoop = (trsNames[0] != "InnerFaces");

  // in the first loop, start from the wall faces and build up to
  // the BL limit (defined by the number of layers chosen by the user)
  // If you end up in a partition face, store the following:
  // -global nodeIDs forming that partition face,
  // -the number of face nodes(in 3D could be 3 or 4),
  // -the number of layers to add starting from such a partition face
  // -the global IDs of the cells on the internal side of the
  //  partition face
  for (CFuint j = 0; j < trsNames.size(); ++j) {
    SafePtr<TopologicalRegionSet> wallTRS = MeshDataStack::getActive()->getTrs(trsNames[j]);
    const CFuint nbTrsFaces = wallTRS->getLocalNbGeoEnts();
    
    if (isFirstLoop) {
      faceToConsider.resize(nbTrsFaces);
      nbLayersToAdd.resize(nbTrsFaces);
      globalCellIDs.resize(nbTrsFaces);
      for (CFuint f = 0; f < nbTrsFaces; ++f) {
	faceToConsider[f] = true;
	nbLayersToAdd[f] = _maxNbNormalFaces;
	globalCellIDs[f] = -1;
      }
    }
    
    // skip this completely if the number of Wall TRS faces == 0
    if (nbTrsFaces > 0) {
      scanLayers(wallTRS, isPartitionFace, faceToConsider,
		 nbLayersToAdd, globalCellIDs, sendData);
    }
  }
  
  SafePtr<TopologicalRegionSet> inFaces =
    MeshDataStack::getActive()->getTrs(std::string("InnerFaces"));
  const CFuint nbLocalInFaces = inFaces->getLocalNbGeoEnts();

  PartitionFaceData recvData;
  MPIStruct ms;
  vector<CFuint> sizeData(2);
  if (isFirstLoop) {
    // after the first loop, the size of the following storages
    // remain unchanged
    faceToConsider.resize(nbLocalInFaces);
    nbLayersToAdd.resize(nbLocalInFaces);
    globalCellIDs.resize(nbLocalInFaces);
  }
  
  // reset all the flags for faces to consider
  for (CFuint f = 0; f < nbLocalInFaces; ++f) {
    faceToConsider[f] = false;
    nbLayersToAdd[f] = 0;
    globalCellIDs[f] = -1;
  }
  
  // in this second phase, we communicate the data of the partition faces
  // from which the layers will "grow" again
  const CFuint nbProc = PE::GetPE().GetProcessorCount();
  const CFuint myRank = PE::GetPE().GetRank();
  cf_assert(nbLocalInFaces > 0);
  
  for (CFuint root = 0; root < nbProc; ++root) {
    if (root == myRank) {
      // copy send data to
      recvData = sendData;
      sizeData[0] = sendData.nodeIDs.size();
      sizeData[1] = sendData.nbFaceNodes.size();
    }
    
    // make the other processors know the size of the data to distribute
    // from the current processor
    MPI_Bcast(&sizeData[0], 2, MPIStructDef::getMPIType(&sizeData[0]), root, PE::GetPE().GetCommunicator());
    
    if (sizeData[0] > 0) {
      if (root != myRank) {
	recvData.nodeIDs.resize(sizeData[0]);
	recvData.nbFaceNodes.resize(sizeData[1]);
	recvData.globalCellIDs.resize(sizeData[1]);
	recvData.nbLayersToAdd.resize(sizeData[1]);
      }
      
      // broadcast the partition faces data
      int lnn[4];
      lnn[0] = sizeData[0];
      lnn[1] = lnn[2] = lnn[3] = sizeData[1];
      MPIStructDef::buildMPIStruct<CFuint,CFuint,CFuint,CFuint>
	(&recvData.nodeIDs[0], &recvData.nbFaceNodes[0],
	 &recvData.globalCellIDs[0], &recvData.nbLayersToAdd[0],
	 lnn, ms);
      MPI_Bcast(ms.start, 1, ms.type, root, PE::GetPE().GetCommunicator());
      
      // discard the data sent by root to itself
      if (root != myRank) {
	// if the data correspond to local partition faces build the layer data
	const CFuint nbFacesToCheck = recvData.nbFaceNodes.size();
	CFLogDebugMin("P" << PE::GetPE().GetRank() << ", nbFacesToCheck = " << nbFacesToCheck << "\n");
	CFLogDebugMin("P" << PE::GetPE().GetRank() << ", nbLocalInFaces = " << nbLocalInFaces << "\n");
	
	CFuint countNodes = 0;
	for (CFuint iFace = 0; iFace < nbFacesToCheck; ++iFace) {
	  const CFuint nbPNodes = recvData.nbFaceNodes[iFace];
	  const CFuint firstNodeID = recvData.nodeIDs[countNodes];
	  bool isFound = false;
	  pair<MapIterator, MapIterator> faces =
	    mapGlobalIDToLocalFaceID.find(firstNodeID, isFound);
	  if (isFound) {
	    for (MapIterator fPtr = faces.first; fPtr != faces.second; ++fPtr) {
	      CFuint countMatchingNodes = 1;
	      const CFuint trsFaceID = fPtr->second;
	      const CFuint nbTrsFaceNodes = inFaces->getNbNodesInGeo(trsFaceID);
	      
	      if (nbTrsFaceNodes == nbPNodes) {
		for (CFuint nn = 1; nn < nbPNodes; ++nn) {
		  const CFuint currNodeID = recvData.nodeIDs[countNodes+nn];
		  for (CFuint in = 0; in < nbPNodes; ++in) {
		    if (nodes[inFaces->getNodeID(trsFaceID,in)]->getGlobalID() == currNodeID) {
		      countMatchingNodes++;
		    }
		  }
		}
	      }
	      
	      if (countMatchingNodes == nbPNodes) {
		CFLogDebugMin("P" << PE::GetPE().GetRank() << ", matching face " << trsFaceID << "\n");
		// the partition face received corresponds to an internal face in this processor
		// store the ID of the matching face and the number of layers to add
		faceToConsider[trsFaceID] = true;
		nbLayersToAdd[trsFaceID] = recvData.nbLayersToAdd[iFace];
		globalCellIDs[trsFaceID] = recvData.globalCellIDs[iFace];
		break;
	      }
	    }
	  }
	  
	  countNodes += nbPNodes;
	}
      }
    }
  }
  
  cf_assert(nbLayersToAdd.size() == nbLocalInFaces);
}

//////////////////////////////////////////////////////////////////////////////

void StegerWarmingCorrFluxMPI::scanLayers(SafePtr<TopologicalRegionSet> wallTRS,
            const vector<bool>& isPartitionFace,
            const vector<bool>& faceToConsider,
            const vector<CFuint>& nbLayersToAdd,
            const vector<CFint>& globalCellIDs,
            PartitionFaceData& sendData)
{
  FaceTrsGeoBuilder::GeoData& faceData = _faceBuilder.getDataGE();
  faceData.trs = wallTRS;
  CellTrsGeoBuilder::GeoData& cellData = _cellBuilder.getDataGE();

  std::string fileName = std::string("blFaces-") + StringOps::to_str(PE::GetPE().GetRank()) + std::string(".dat");
  ofstream fout(fileName.c_str(),ios_base::app);

  const bool isInnerFace = (wallTRS->getName() == "InnerFaces");

  CFuint cellID = 0;
  const CFuint nbTrsFaces = wallTRS->getLocalNbGeoEnts();
  for (CFuint iFace = 0; iFace < nbTrsFaces; ++iFace) {
    if (faceToConsider[iFace]) {
      faceData.idx = iFace;
      faceData.isBFace = !isInnerFace;

      GeometricEntity *const face = _faceBuilder.buildGE();

      if (!isInnerFace) {
	cellID = face->getState(0)->getLocalID();
      }
      else {
	cf_assert(isInnerFace);
	const CFint otherCellID = static_cast<CFuint>(globalCellIDs[iFace]);
	cf_assert(otherCellID != -1);
	// consider the cell which lies on the other side of the
	// partition face coming from the donor process
	cellID = (face->getState(0)->getGlobalID() !=
		  static_cast<CFuint>(otherCellID)) ?
	  face->getState(0)->getLocalID() :
	  face->getState(1)->getLocalID();
	
	cf_assert(!face->getState(0)->isGhost());
	cf_assert(!face->getState(1)->isGhost());
      }
      
      CFuint faceID = face->getID();
      _flagNormalFace[faceID] = true;
      
      bool reachedPartitionFace = false;
      CFuint countFaces = _maxNbNormalFaces - nbLayersToAdd[iFace]+1;
      while (countFaces < _maxNbNormalFaces && (!reachedPartitionFace)) {
	cellData.idx = cellID;
	CFuint oppositeIFace = 0;
	GeometricEntity *const cell = _cellBuilder.buildGE();
	const vector<GeometricEntity*>& cellFaces = *cell->getNeighborGeos();
	const CFuint nbCellFaces =  cellFaces.size();
	
	for (CFuint f = 0; f < nbCellFaces; ++f) {
	  if (cellFaces[f]->getID() == faceID) {
	    oppositeIFace = getMethodData().getOppositeIFace
	      (f, PhysicalModelStack::getActive()->getDim(), cell->nbNodes());
	    const CFuint leftID = cellFaces[oppositeIFace]->getState(0)->getLocalID();
	    const CFuint rightID = cellFaces[oppositeIFace]->getState(1)->getLocalID();
	    cellID = (leftID == cellID) ? rightID : leftID;
	    countFaces++;
	    faceID = cellFaces[oppositeIFace]->getID();
	    _flagNormalFace[faceID] = true;
	    
	    // print cell center coordinates
	    fout << cell->getState(0)->getCoordinates() << endl;
	    
	    // if you reach the partition face without having reached the
	    // number of layers _maxNbNormalFaces, store additional data
	    if (isPartitionFace[faceID] && countFaces < _maxNbNormalFaces) {
	      const CFuint nbPFNodes = cellFaces[oppositeIFace]->nbNodes();
	      for (CFuint in = 0; in < nbPFNodes; ++in) {
		sendData.nodeIDs.push_back
		  (cellFaces[oppositeIFace]->getNode(in)->getGlobalID());
	      }
	      sendData.nbFaceNodes.push_back(nbPFNodes);
	      sendData.globalCellIDs.push_back(cell->getState(0)->getGlobalID());
	      // check the following
	      sendData.nbLayersToAdd.push_back(_maxNbNormalFaces - countFaces);
	      reachedPartitionFace = true;
	    }
	    
	    break;
	  }
	}
	_cellBuilder.releaseGE();
      }
      
      fout<< endl;
      
      _faceBuilder.releaseGE();
      faceData.isBFace = false;
    }
  }
  
  fout.close();
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
