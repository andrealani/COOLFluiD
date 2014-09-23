#include <numeric>

#include "Common/BadValueException.hh"
#include "Common/ParserException.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/SubSystemStatus.hh"
#include "Framework/PathAppender.hh"
#include "Framework/RadiationLibrary.hh"
#include "Framework/DofDataHandleIterator.hh"
#include "Environment/FileHandlerOutput.hh"
#include "Environment/DirPaths.hh"

#include "RadiativeTransferSanna/RadiativeTransferSanna.hh"
#include "RadiativeTransferSanna/Slab1DFVMCCMPI.hh"
#include "FiniteVolume/CellCenterFVM.hh"

#include "Common/PE.hh"
#include "Common/MPI/MPIStructDef.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Numerics::FiniteVolume;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace RadiativeTransfer {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<Slab1DFVMCCMPI, DataProcessingData, RadiativeTransfer> 
slab1DFVMCCMPIProvider("Slab1DFVMCCMPI");

//////////////////////////////////////////////////////////////////////////////

void Slab1DFVMCCMPI::defineConfigOptions(Config::OptionList& options)
{ 
  options.addConfigOption< CFuint >
    ("nbCellsOnLine", "Maximum number of cells in the direction orthogonal to the wall");
}
      
//////////////////////////////////////////////////////////////////////////////
      
Slab1DFVMCCMPI::Slab1DFVMCCMPI(const std::string& name) :
  Slab1DFVMCC(name),
  m_myRank(),
  m_localNbLines(),
  m_nbStatesInProc(), 
  m_sendGlobalIDs(),
  m_midFace(),
  m_mapTrsNameToNbFaces(),
  m_mapGlobalToMeshLineID(),
  m_mapGlobalToLocalStateIDs(),
  m_innerFaceIDToGlobalWallFaceID(),
  m_statesNodes(),
  m_sendArray(),
  m_qRadByLine(),
  m_qRadToSend(),
  m_qRadToRecv(),
  m_destIDsPerRank(),
  m_donorIDsPerRank(),
  m_destIDsToSend(),
  m_destIDsToRecv(),
  m_sendCount(),
  m_sendDispl(),
  m_recvCount(),
  m_recvDispl(),
  m_maxNbNormalFaces(0)
{
  addConfigOptionsTo(this);
  
  m_nbCellsOnLine = 0;
  setParameter("nbCellsOnLine",&m_nbCellsOnLine); 
}

//////////////////////////////////////////////////////////////////////////////

Slab1DFVMCCMPI::~Slab1DFVMCCMPI()
{
}

//////////////////////////////////////////////////////////////////////////////

void Slab1DFVMCCMPI::unsetup()
{
  Slab1DFVMCC::unsetup();
}
 
//////////////////////////////////////////////////////////////////////////////
      
void Slab1DFVMCCMPI::setup()
{ 
  if (m_nbCellsOnLine == 0) {
    cout << "ERROR: Slab1DFVMCCMPI::setup() nbCellsOnLine == 0 " << endl;
    abort();
  }
  
  m_myRank = PE::GetPE().GetRank();
  
  m_maxNbNormalFaces = m_nbCellsOnLine + 1;
  
  const CFuint dim = PhysicalModelStack::getActive()->getDim();
  
  // default stagnation point is NOT set
  if (m_stagnationPointXYZ.size() > 0) {
    cf_always_assert(m_stagnationPointXYZ.size() == dim);
    
    m_stagPoint.resize(dim);
    for (CFuint i = 0; i < dim; ++i) {
      m_stagPoint[i] = m_stagnationPointXYZ[i];
    }
  }
  else {
    m_stagPoint.resize(dim);
    m_stagPoint = 0.0;
  }
  
  m_midFace.resize(dim);
  
  // suppose that just one space method is available
  SafePtr<SpaceMethod> spaceMethod = getMethodData().getCollaborator<SpaceMethod>();
  SafePtr<CellCenterFVM> fvmcc = spaceMethod.d_castTo<CellCenterFVM>();
  cf_assert(fvmcc.isNotNull());
  
  m_fvmccData = fvmcc->getData();
  
  Common::SafePtr<TopologicalRegionSet> cells = MeshDataStack::getActive()->getTrs("InnerCells");
  const vector<vector<CFuint> >&  trsInfo = MeshDataStack::getActive()->getTotalTRSInfo();
  const vector<std::string>& trsNames = MeshDataStack::getActive()->getTotalTRSNames();
  vector< Common::SafePtr<TopologicalRegionSet> >& trsList = getTrsList();
  
  // compute the total number of wall faces in the whole distributed mesh
  m_nbWallFaces = 0;
  const CFuint nbWallTRS = trsList.size();
  for(CFuint iTRS = 0; iTRS < nbWallTRS; ++iTRS) {
    SafePtr<TopologicalRegionSet> currTrs = trsList[iTRS];
    for (CFuint i = 0; i < trsNames.size(); ++i) {
      if (currTrs->getName() == trsNames[i]) {
	CFuint nbTrsFaces = 0;
	const CFuint nbTRsInTRS = trsInfo[iTRS].size();
	for(CFuint tr = 0; tr < nbTRsInTRS; ++tr) {
	  m_nbWallFaces += trsInfo[iTRS][tr];
	  nbTrsFaces += trsInfo[iTRS][tr];
	}
	
	// keep accountancy of the number of faces in each wall TRS  
	m_mapTrsNameToNbFaces.insert(currTrs->getName(), nbTrsFaces);
	break;
      }
    }
  }
  m_mapTrsNameToNbFaces.sortKeys();
  
  Common::SafePtr<TopologicalRegionSet> inFaces = MeshDataStack::getActive()->getTrs("InnerFaces");
  const CFuint nbInnerFaces = inFaces->getLocalNbGeoEnts();
  m_innerFaceIDToGlobalWallFaceID.resize(nbInnerFaces, -1);
  
  buildMeshLines();
}

//////////////////////////////////////////////////////////////////////////////
           
void Slab1DFVMCCMPI::buildMeshLines()
{
  // face builder 
  Common::SafePtr<GeometricEntityPool<FaceTrsGeoBuilder> > faceBuilder = m_fvmccData->getFaceTrsGeoBuilder();
  SafePtr<FaceTrsGeoBuilder> faceBuilderPtr = faceBuilder->getGeoBuilder();
  faceBuilderPtr->setDataSockets(socket_states, socket_gstates, socket_nodes);
  
  // cell builder
  SafePtr<GeometricEntityPool<CellTrsGeoBuilder> > cellBuilder = m_fvmccData->getCellTrsGeoBuilder();
  SafePtr<CellTrsGeoBuilder> cellBuilderPtr = cellBuilder->getGeoBuilder();
  cellBuilderPtr->setDataSockets(socket_states, socket_gstates, socket_nodes);
  
  SafePtr<TopologicalRegionSet> cells = MeshDataStack::getActive()->getTrs(std::string("InnerCells"));
  CellTrsGeoBuilder::GeoData& cellData = cellBuilder->getDataGE();
  cellData.trs = cells;
  
  SafePtr<TopologicalRegionSet> pFaces = MeshDataStack::getActive()->getTrs(std::string("PartitionFaces"));
  SafePtr<TopologicalRegionSet> inFaces = MeshDataStack::getActive()->getTrs(std::string("InnerFaces"));
  
  const CFuint totNbFaces = MeshDataStack::getActive()->Statistics().getNbFaces();
  
  // identify the partition faces by ID
  vector<bool> isPartitionFace(totNbFaces);
  isPartitionFace.assign(totNbFaces, false);
  const CFuint nbLocalPFaces = pFaces->getLocalNbGeoEnts();
  
  DataHandle<Node*, GLOBAL> nodes = socket_nodes.getDataHandle();
  const CFuint nbLocalInFaces = inFaces->getLocalNbGeoEnts();
  const CFuint dim = PhysicalModelStack::getActive()->getDim();
  
  // here the size of the global map is a bit overestimated
  CFMultiMap<CFuint, CFuint> mapGlobalIDToLocalFaceID(nbLocalInFaces*dim);
  
  // flag the partition faces
  for (CFuint iFace = 0; iFace < nbLocalPFaces; ++iFace) {
    const CFuint localGeoID = pFaces->getLocalGeoID(iFace);
    isPartitionFace[localGeoID] = true;
  }
  
  // build a mapping between global node IDs and local IDs of
  // corresponding inner faces
 
  for (CFuint iFace = 0; iFace < nbLocalInFaces; ++iFace) {
    const CFuint nbNodesInFace = inFaces->getNbNodesInGeo(iFace);
    for (CFuint nf = 0; nf < nbNodesInFace; ++nf) {
      const CFuint globalNodeID = nodes[inFaces->getNodeID(iFace, nf)]->getGlobalID();
      mapGlobalIDToLocalFaceID.insert(globalNodeID,iFace);
    }
  }
  mapGlobalIDToLocalFaceID.sortKeys();
    
  vector<bool> faceToConsider;
  vector<CFuint> nbLayersToAdd;
  vector<CFint> globalCellIDs;
  vector<CFint> globalFaceIDs;
  
  // the array with all the global cellIDs starting from the boundary
  // line by line one after the other
  const CFuint nbWallTRS = getTrsNames().size();
  vector< vector<CFint> > lineGlobalCellIDs(nbWallTRS);
  for (CFuint iTRS = 0; iTRS < nbWallTRS; ++iTRS) {
    // total number of wall faces in current TRS
    const CFuint nbWallFacesInTRS = m_mapTrsNameToNbFaces.find(getTrsNames()[iTRS]);
    const CFuint sizeLineMesh = nbWallFacesInTRS*m_nbCellsOnLine;
    lineGlobalCellIDs[iTRS].resize(sizeLineMesh);
    lineGlobalCellIDs[iTRS].assign(sizeLineMesh, -1);
  }
    
  // array for storing the distances between each wall face center and the stagnation point 
  vector< RealVector > distanceFromStagPoint(nbWallTRS);
  for (CFuint iTRS = 0; iTRS < nbWallTRS; ++iTRS) {
    // total number of wall faces in current TRS
    const CFuint nbWallFacesInTRS = m_mapTrsNameToNbFaces.find(getTrsNames()[iTRS]);
    distanceFromStagPoint[iTRS].resize(nbWallFacesInTRS);
  }  
  
  // at first, wall TRS are processed
  processBFaces(getTrsNames(),
		mapGlobalIDToLocalFaceID,
		isPartitionFace,
		faceToConsider,
		nbLayersToAdd,
		globalCellIDs,
		globalFaceIDs,
		lineGlobalCellIDs,
		distanceFromStagPoint);
    
  vector<std::string> pName(1, std::string("InnerFaces"));
  CFuint localNbLayersToAdd = accumulate(nbLayersToAdd.begin(), nbLayersToAdd.end(),0);
  CFuint globalNbLayersToAdd = 0;
  
  MPI_Allreduce(&localNbLayersToAdd, &globalNbLayersToAdd,
		1, MPI_UNSIGNED,MPI_MAX,PE::GetPE().GetCommunicator());
  
  cout << "#### processor A " << m_myRank << endl;
  
  
  // check partition data at tghis point
  // string filen = "file.dat-" + StringOps::to_str(m_myRank);
//   ofstream file(filen.c_str());
  
//   file << "localNbLayersToAdd  = " << localNbLayersToAdd << endl;
//   file << "globalNbLayersToAdd = " << globalNbLayersToAdd << endl;
  
//   file << "globalCellIDs" << endl;
//   for (CFuint i = 0; i < globalCellIDs.size(); ++i) {
//     file << i << " => " << globalCellIDs[i] << endl;
//   }
  
//   file << "faceToConsider" << endl;
//   for (CFuint i = 0; i < faceToConsider.size(); ++i) {
//     if (faceToConsider[i]) {
//       file << i << " => " << true << endl;
//     }
//     else {
//       file << i << " => " << false << endl;
//     }
//   }
    
//   abort();
  
  while (globalNbLayersToAdd > 0) {
    processBFaces(pName,
		  mapGlobalIDToLocalFaceID,
		  isPartitionFace,
		  faceToConsider,
		  nbLayersToAdd,
		  globalCellIDs,
		  globalFaceIDs,
		  lineGlobalCellIDs,
		  distanceFromStagPoint);
    
    cout << "#### processor AA " << m_myRank << endl;
    
    localNbLayersToAdd = accumulate(nbLayersToAdd.begin(), nbLayersToAdd.end(),0);
    
    cout << "#### processor AB " << m_myRank << endl;
    for (;;) {}
    
    MPI_Allreduce(&localNbLayersToAdd, &globalNbLayersToAdd,
		  1, MPI_UNSIGNED,MPI_MAX,PE::GetPE().GetCommunicator());
  }
  
  cout << "#### processor " << m_myRank << endl;
  for (;;) {}
  
  vector<vector<CFint> > globalCellLines(nbWallTRS);
  for (CFuint iTRS = 0; iTRS < nbWallTRS; ++iTRS) {
    // total number of wall faces in current TRS
    const CFuint sizeLineMesh = lineGlobalCellIDs[iTRS].size();
    globalCellLines[iTRS].resize(sizeLineMesh);
    
    MPI_Allreduce(&lineGlobalCellIDs[iTRS][0], &globalCellLines[iTRS][0],
		  sizeLineMesh, MPI_INT, MPI_MAX, PE::GetPE().GetCommunicator());
  }
  
  vector<RealVector> distanceFromStagPointOut(nbWallTRS);
  for (CFuint iTRS = 0; iTRS < nbWallTRS; ++iTRS) {
    // total number of wall faces in current TRS
    const CFuint nbWallFacesInTRS = m_mapTrsNameToNbFaces.find(getTrsNames()[iTRS]);
    distanceFromStagPointOut[iTRS].resize(nbWallFacesInTRS);
    
    MPI_Allreduce(&distanceFromStagPoint[iTRS][0], &distanceFromStagPointOut[iTRS][0],
		  nbWallFacesInTRS, MPI_DOUBLE, MPI_MAX, PE::GetPE().GetCommunicator());
  }
  
  // at this point each processor has knowledge of the whole line mesh, but we have to 
  // reorder the lines following the streamwise direction
  CFMap<CFreal, CFint*> mapDistanceWithLine(m_nbWallFaces);
  for (CFuint iTRS = 0; iTRS < nbWallTRS; ++iTRS) {
    RealVector& distanceCurrLine = distanceFromStagPointOut[iTRS];
    for (CFuint i = 0; i < distanceCurrLine.size(); ++i) {
      mapDistanceWithLine.insert(distanceCurrLine[i],&globalCellLines[iTRS][i*m_nbCellsOnLine]);
    }
  }
  mapDistanceWithLine.sortKeys();
  
  // each processor can now create the corresponding part of the structured mesh 
  // which will be fed to the radiation library
  // the first processor gets more lines, if the total number of lines cannot be equally divided
  const CFuint nbProcs = PE::GetPE().GetProcessorCount();
  m_localNbLines.resize(nbProcs);
  
  const CFuint ne = m_nbWallFaces/nbProcs;
  for (CFuint i = 0; i < nbProcs; ++i) {
    m_localNbLines[i] = (i != 0) ?  ne : ne + m_nbWallFaces%nbProcs;
  }
  
  const CFuint localNbLines = m_localNbLines[m_myRank];
  const CFuint localMeshSize = localNbLines*m_nbCellsOnLine;
  m_meshByLines.reserve(localNbLines);
  
  vector<CFuint> startLine(nbProcs);
  CFuint start = 0;
  for (CFuint i = 0; i < nbProcs; ++i) {
    startLine[i] = start;
    start += m_localNbLines[i];
  }
  
  // check this !!!!!!!!!!!
  cf_assert(start == m_nbWallFaces); 
  
  m_mapGlobalToMeshLineID.reserve(localMeshSize);
  
  // allocate the local portion of the structured mesh by line 
  CFuint counter = startLine[m_myRank];
  CFuint meshByLineID = 0;
  
  for (CFuint i = 0; i < m_meshByLines.size(); ++i) {
    m_meshByLines[i] = new vector<CFuint>(m_nbCellsOnLine);
    CFint* cellIDs = mapDistanceWithLine[counter++];
    for (CFuint j = 0; j < m_nbCellsOnLine; ++j, ++meshByLineID) {
      (*m_meshByLines[i])[j] = static_cast<CFuint>(cellIDs[j]);
      m_mapGlobalToMeshLineID.insert(cellIDs[j], meshByLineID);
    }
  }
  m_mapGlobalToMeshLineID.sortKeys();
  
  // count the number of states in each processor
  const CFuint nbProc = PE::GetPE().GetProcessorCount();
  vector<CFuint> nbStatesInProc(nbProc, 0);
  nbStatesInProc[m_myRank] = socket_states.getDataHandle().size();
  m_nbStatesInProc.resize(nbProc);
  
  MPI_Allreduce(&nbStatesInProc[0], &m_nbStatesInProc[0],
		nbProc, MPI_UNSIGNED,MPI_MAX,PE::GetPE().GetCommunicator());
    
  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  const CFuint stride = nbEqs + dim;
  m_statesNodes.resize(localMeshSize*stride);
  
  // the qrad DataHandle array of qrad is per local cell
  DataHandle<CFreal> qrad = socket_qrad.getDataHandle();  
  qrad.resize(socket_states.getDataHandle().size()); 
  
  // the local (in this class) qrad array stores the oucome from the 
  // integration by line 
  m_qRadByLine.resize(localMeshSize);
  
  
  DataHandle<State*, GLOBAL> states = socket_states.getDataHandle();
  // map local to global state IDs in this process
  m_mapGlobalToLocalStateIDs.reserve(states.size());
  for (CFuint i = 0; i < states.size(); ++i) {
    m_mapGlobalToLocalStateIDs.insert(states[i]->getGlobalID(), 
				      states[i]->getLocalID());
  }
  m_mapGlobalToLocalStateIDs.sortKeys();
  
  
  // you need to take into account the processor IDs to which each global cell ID belong  
  // otherwise you can just pass around the whole qrad storage and each processor takes 
  // what it needs (first easy solution)
  
  m_destIDsPerRank.resize(nbProc);
  m_donorIDsPerRank.resize(nbProc);
  
  CFuint countDestIDs = 0;
  for (CFuint r = 0; r < nbProc; ++r) {
    const CFuint nbBcastStates = m_nbStatesInProc[r];
    m_sendGlobalIDs.resize(nbBcastStates);
    
    if (r == m_myRank) {
      cf_assert(nbBcastStates == states.size());
      for (CFuint i = 0; i < states.size(); ++i) {
	m_sendGlobalIDs[i] = states[i]->getGlobalID();
      }
    }
    
    MPI_Bcast(&m_sendGlobalIDs[0], nbBcastStates, MPI_UNSIGNED, r, PE::GetPE().GetCommunicator());
    
    for (CFuint i = 0; i < nbBcastStates; ++i) {
      bool isFound = false;
      const CFuint currGlobalID = m_sendGlobalIDs[i];
      
      pair<MapItr, MapItr> ids = m_mapGlobalToMeshLineID.find(currGlobalID, isFound);
      if (isFound) {
	const CFuint meshByLineID = ids.first->second;
	// keep track of the local state IDs in the destination processor
	m_destIDsPerRank[r].push_back(i);
	countDestIDs++;
	
	// keep track of the local meshByLineID corresponding to the m_qRadByLine value to send
	m_donorIDsPerRank[r].push_back(meshByLineID);
      }
    }
  }
  
  m_destIDsToSend.resize(countDestIDs);
  
  m_sendCount.resize(nbProc); m_sendCount.assign(nbProc,static_cast<CFuint>(0));
  m_sendDispl.resize(nbProc); m_sendDispl.assign(nbProc,static_cast<CFuint>(0));
  m_recvCount.resize(nbProc); m_recvCount.assign(nbProc,static_cast<CFuint>(0));
  m_recvDispl.resize(nbProc); m_recvDispl.assign(nbProc,static_cast<CFuint>(0));
  
  m_sendDispl[0] = 0;
  for (CFuint r = 0; r < nbProc; ++r) {
    m_sendCount[r] = m_destIDsPerRank[r].size();
    
    if (r > 0) {
      m_sendDispl[r] += m_sendCount[r-1];
    }
  }
  
  MPI_Alltoall(&m_sendCount[0], 1, MPI_UNSIGNED,
	       &m_recvCount[0], 1, MPI_UNSIGNED, 
	       PE::GetPE().GetCommunicator());
  
  m_recvDispl[0] = 0;
  CFuint count = m_recvCount[0];
  for (CFuint r = 1; r < nbProc; ++r) {
    if (m_recvCount[r] > 0) {
      m_recvDispl[r] = count;
    }
    count += m_recvCount[r];  
  }
  
  m_qRadToSend.resize(localMeshSize);
  
  const CFuint totRecvCount = std::accumulate
    (m_recvCount.begin(), m_recvCount.end(), 0);
  m_qRadToRecv.resize(totRecvCount);
  m_destIDsToRecv.resize(totRecvCount);
  m_destIDsToRecv.assign(totRecvCount, -1);
  
  CFuint countSend = 0;
  for (CFuint r = 0; r < nbProc; ++r) {
    const CFuint nr = m_destIDsPerRank[r].size();
    for (CFuint i = 0; i < nr ; ++i, ++countSend) {
      m_destIDsToSend[countSend] = m_destIDsToSend[m_donorIDsPerRank[r][i]]; 
    }
  }  
  
  MPI_Alltoallv(&m_destIDsToSend[0], &m_sendCount[0],
		&m_sendDispl[0], MPI_INT, &m_destIDsToRecv[0],
		&m_recvCount[0], &m_recvDispl[0],
		MPI_INT, PE::GetPE().GetCommunicator());
  
  cf_assert(totRecvCount == states.size()); 
  
  // sanity check
  for (CFuint i = 0; i < m_destIDsToRecv.size(); ++i) {
    if (m_destIDsToRecv[i] == -1) {
      cout << "Slab1DFVMCCMPI::buildMeshLines() => m_destIDsToRecv has -1" << endl;
      abort();
    }
  }
  
  // m_destIDsPerRank and m_donorIDsPerRank don't change anymore 
}
      
//////////////////////////////////////////////////////////////////////////////

void Slab1DFVMCCMPI::processBFaces(const vector<std::string>& trsNames,
				   CFMultiMap<CFuint,CFuint>& mapGlobalIDToLocalFaceID,
				   const vector<bool>& isPartitionFace,
				   vector<bool>& faceToConsider,
				   vector<CFuint>& nbLayersToAdd,
				   vector<CFint>& globalCellIDs,
				   vector<CFint>& globalFaceIDs,
				   vector<vector<CFint> >& lineGlobalCellIDs,
				   vector<RealVector>& distanceFromStagPoint)
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
  // -the global IDs of the cells on the internal side of the partition face
  // -the global IDs of the starting wall face (this info should survive all future data exchanges)
  // -the array of global cellIDs beloging to the line starting from the 
  //  input wall or partition boundary 
  
  for (CFuint j = 0; j < trsNames.size(); ++j) {
    SafePtr<TopologicalRegionSet> wallTRS = MeshDataStack::getActive()->getTrs(trsNames[j]);
    const CFuint nbTrsFaces = wallTRS->getLocalNbGeoEnts();
    
    if (isFirstLoop) {
      faceToConsider.resize(nbTrsFaces);
      nbLayersToAdd.resize(nbTrsFaces);
      globalCellIDs.resize(nbTrsFaces);
      globalFaceIDs.resize(nbTrsFaces);
      for (CFuint f = 0; f < nbTrsFaces; ++f) {
	faceToConsider[f] = true;
	nbLayersToAdd[f] = m_maxNbNormalFaces; //m_nbCellsOnLine;
	globalCellIDs[f] = -1;
	globalFaceIDs[f] = -1;
      }
    }
    
    cout << m_myRank << " before scan" <<  endl;
    // skip this completely if the number of Wall TRS faces == 0
    if (nbTrsFaces > 0) {
      scanLayers(wallTRS, isPartitionFace, faceToConsider,
		 nbLayersToAdd, globalCellIDs, globalFaceIDs, 
		 lineGlobalCellIDs[j], 
		 distanceFromStagPoint[j],sendData);
    } 
    cout << m_myRank << " after scan" <<  endl;
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
    globalFaceIDs.resize(nbLocalInFaces);
  }
  
  // reset all the flags for faces to consider
  for (CFuint f = 0; f < nbLocalInFaces; ++f) {
    faceToConsider[f] = false;
    nbLayersToAdd[f] = 0;
    globalCellIDs[f] = -1;
    globalFaceIDs[f] = -1;
  }
    
  // in this second phase, we communicate the data of the partition faces
  // from which the layers will "grow" again
  const CFuint nbProc = PE::GetPE().GetProcessorCount();
  const CFuint myRank = PE::GetPE().GetRank();
  cf_assert(nbLocalInFaces > 0);
  
  // string file = "par.dat-" + StringOps::to_str(m_myRank);
  // ofstream fout(file.c_str());
  
  for (CFuint root = 0; root < nbProc; ++root) {
    if (root == myRank) {
      // copy send data to
      recvData = sendData;
      
      // this size corresponds to the number of nodes on the partition boundary 
      // for the sending process
      sizeData[0] = sendData.nodeIDs.size(); 
      // this size corresponds to the number of faces on the partition boundary
      // for the sending process
      sizeData[1] = sendData.nbFaceNodes.size();
    }
        
    // make the other processors know the size of the data to distribute
    // from the current processor
    MPI_Bcast(&sizeData[0], 2, MPI_UNSIGNED, root, PE::GetPE().GetCommunicator());
        
    if (sizeData[0] > 0) {
      cf_assert(sizeData[1] > 0.);
      if (root != myRank) {
	recvData.nodeIDs.resize(sizeData[0]);
	recvData.nbFaceNodes.resize(sizeData[1]);
	recvData.globalCellIDs.resize(sizeData[1]);
	recvData.globalFaceIDs.resize(sizeData[1]);
	recvData.nbLayersToAdd.resize(sizeData[1]);
      }
          
      // broadcast the partition faces data
      int lnn[5]; // new entry here
      lnn[0] = sizeData[0];
      lnn[1] = lnn[2] = lnn[3] = lnn[4] = sizeData[1];
      MPIStructDef::buildMPIStruct<CFuint,CFuint,CFuint,int, CFuint>
	(&recvData.nodeIDs[0], &recvData.nbFaceNodes[0],
	 &recvData.globalCellIDs[0],&recvData.globalFaceIDs[0], 
	 &recvData.nbLayersToAdd[0], lnn, ms);
      
      MPI_Bcast(ms.start, 1, ms.type, root, PE::GetPE().GetCommunicator());
      
      // discard the data sent by root to itself
      if (root != myRank) {
	// if the data correspond to local partition faces build the layer data
	const CFuint nbFacesToCheck = recvData.nbFaceNodes.size();
	CFLogDebugMin("P" << PE::GetPE().GetRank() << ", nbFacesToCheck = " << nbFacesToCheck << "\n");
	CFLogDebugMin("P" << PE::GetPE().GetRank() << ", nbLocalInFaces = " << nbLocalInFaces << "\n");
	
	CFuint countNodes = 0;
	for (CFuint iFace = 0; iFace < nbFacesToCheck; ++iFace) {
	  cf_assert(iFace < nbFacesToCheck);
	  const CFuint nbPNodes = recvData.nbFaceNodes[iFace];
	  cf_assert(countNodes < recvData.nodeIDs.size());
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
		cf_assert(trsFaceID < faceToConsider.size());
		cf_assert(trsFaceID < nbLayersToAdd.size());
		cf_assert(trsFaceID < globalCellIDs.size());
		cf_assert(trsFaceID < globalFaceIDs.size());
		cf_assert(iFace < recvData.nbLayersToAdd.size());
		cf_assert(iFace < recvData.globalCellIDs.size());
		cf_assert(iFace < recvData.globalFaceIDs.size());
		
		faceToConsider[trsFaceID] = true;
		nbLayersToAdd[trsFaceID] = recvData.nbLayersToAdd[iFace];
		globalCellIDs[trsFaceID] = recvData.globalCellIDs[iFace];
		globalFaceIDs[trsFaceID] = recvData.globalFaceIDs[iFace];
		
		if (m_innerFaceIDToGlobalWallFaceID[trsFaceID] == -1) {
		  cout << "trsFaceID = " << trsFaceID << " has ID = " << recvData.globalFaceIDs[iFace] << endl;
		  // each trsFaceID should be processed only once
		  m_innerFaceIDToGlobalWallFaceID[trsFaceID] = recvData.globalFaceIDs[iFace];
		  cf_assert(recvData.globalFaceIDs[iFace] >= 0);
		}
		
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

void Slab1DFVMCCMPI::scanLayers(SafePtr<TopologicalRegionSet> wallTRS,
				const vector<bool>& isPartitionFace,
				const vector<bool>& faceToConsider,
				const vector<CFuint>& nbLayersToAdd,
				const vector<CFint>& globalCellIDs,
				const vector<CFint>& globalFaceIDs,
				vector<CFint>& lineGlobalCellIDsInTRS,
 				RealVector& distanceFromStagPoint,
				PartitionFaceData& sendData)
{
  // face builder 
  Common::SafePtr<GeometricEntityPool<FaceTrsGeoBuilder> > faceBuilder = m_fvmccData->getFaceTrsGeoBuilder();
  SafePtr<FaceTrsGeoBuilder> faceBuilderPtr = faceBuilder->getGeoBuilder();
  FaceTrsGeoBuilder::GeoData& faceData = faceBuilder->getDataGE();
  faceData.trs = wallTRS;
  
  // cell builder
  SafePtr<GeometricEntityPool<CellTrsGeoBuilder> > cellBuilder = m_fvmccData->getCellTrsGeoBuilder();
  SafePtr<CellTrsGeoBuilder> cellBuilderPtr = cellBuilder->getGeoBuilder();
  CellTrsGeoBuilder::GeoData& cellData = cellBuilder->getDataGE();
  
  //std::string fileName = std::string("blFaces-") + StringOps::to_str(PE::GetPE().GetRank()) + std::string(".dat");
  //ofstream fout(fileName.c_str(),ios_base::app);
  
  static CFuint maxID = 0;
  
  const bool isInnerFace = (wallTRS->getName() == "InnerFaces");
  const bool isWallTRS   = (wallTRS->getName() != "InnerFaces" && 
			    wallTRS->getName() != "PartitionFaces");
  
  CFuint cellID = 0;
  const CFuint nbTrsFaces = wallTRS->getLocalNbGeoEnts();
  for (CFuint iFace = 0; iFace < nbTrsFaces; ++iFace) {
    // a value must be assigned to this
    int globalWallFaceID = -1;
    
    if (faceToConsider[iFace]) {
      faceData.idx = iFace;
      faceData.isBFace = !isInnerFace;
      
      GeometricEntity *const face = faceBuilder->buildGE();
      
      CFuint cellsOnLine = 0;
      if (!isInnerFace) {
	cellID = face->getState(0)->getLocalID();
	
	// the first internal cell starting from the (wall or partition) boundary face
	// belongs to the cells line
	const CFuint globalCellID = face->getState(0)->getGlobalID();
	cf_assert(!face->getState(0)->isGhost());
	
	// set the current global wall face ID
	globalWallFaceID = wallTRS->getGlobalGeoID(iFace);
	const CFuint cID = globalWallFaceID*m_nbCellsOnLine + cellsOnLine;
	cf_assert(cID < lineGlobalCellIDsInTRS.size());
	lineGlobalCellIDsInTRS[cID] = globalCellID;
	
	// if we are considering a wall face, store the distance from the 
	// stagnation point to the corresponding face mid-point 
	// it will be useful later for ordering the lines
	if (!isPartitionFace[face->getID()]) {
	  m_midFace = 0.0;
	  const CFuint nbNodesInFace =  face->nbNodes();
	  for (CFuint n = 0; n < nbNodesInFace; ++n) {
	    m_midFace += *face->getNode(n);
	  }
	  m_midFace /= nbNodesInFace;
	  cf_assert(globalWallFaceID < distanceFromStagPoint.size());
	  cf_assert(globalWallFaceID != -1);
	  
	  distanceFromStagPoint[globalWallFaceID] = MathFunctions::getDistance(m_midFace, m_stagPoint); 
	  
	  // keep track of the global face ID of this boundary face
	  sendData.globalFaceIDs.push_back(globalWallFaceID);
	}
	
	cellsOnLine++; 
      }
      else {
	cf_assert(isInnerFace);
	cf_assert(iFace < globalCellIDs.size());
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
      
      // those are for sanity check
      bool reachedPartitionFace = false;
      CFuint countFaces = m_maxNbNormalFaces - nbLayersToAdd[iFace]+1;
      //m_nbCellsOnLine - nbLayersToAdd[iFace]+1;
      
      while (countFaces < m_maxNbNormalFaces && (!reachedPartitionFace)) {
	cellData.idx = cellID;
	CFuint oppositeIFace = 0;
	GeometricEntity *const cell = cellBuilder->buildGE();
	const vector<GeometricEntity*>& cellFaces = *cell->getNeighborGeos();
	const CFuint nbCellFaces =  cellFaces.size();
	
	bool foundFace = false;
	CFuint countMatches = 0;
	for (CFuint f = 0; f < nbCellFaces; ++f) {
	  // consider the face opposite to the current one inside the given cell
	  if (cellFaces[f]->getID() == faceID) {
	    foundFace = true;
	    countMatches++;
	    
	    oppositeIFace = m_fvmccData->getOppositeIFace
	      (f, PhysicalModelStack::getActive()->getDim(), cell->nbNodes());
	    
	    GeometricEntity *const oppositeFace = cellFaces[oppositeIFace];
	    faceID = oppositeFace->getID(); // update the faceID
	    countFaces++;
	    
	    // the leftID should always be defined because the state is internal
	    cf_assert(!oppositeFace->getState(0)->isGhost());
	    const CFuint leftID = oppositeFace->getState(0)->getLocalID();
	    const bool isBoundaryFace = oppositeFace->getState(1)->isGhost();
	    CFuint globalCellID = 0;
	    
	    if (!isPartitionFace[faceID] && (!isBoundaryFace)) {
	      if (leftID == cellID) {
		cf_assert(!oppositeFace->getState(1)->isGhost());
		const CFuint rightID = oppositeFace->getState(1)->getLocalID();
		cellID = rightID;
		globalCellID = oppositeFace->getState(1)->getGlobalID();
	      }
	      else {
		cellID = leftID;
		globalCellID = oppositeFace->getState(0)->getGlobalID();
	      }
	      
	      // print cell center coordinates
	      // fout << cell->getState(0)->getCoordinates() << endl;
	      // store the newly found cell ID
	      
	      //// FIX THIS !!!! global face ID DOESN'T exist for inner faces !!! you cannot ask it !!!
	      cf_assert(iFace < wallTRS->getLocalNbGeoEnts());
	      
	      CFuint globalFaceID = -1;
	      if (!isWallTRS) {
		if (m_innerFaceIDToGlobalWallFaceID[iFace] != -1) {
		  globalFaceID = m_innerFaceIDToGlobalWallFaceID[iFace];
		}
		else if (m_innerFaceIDToGlobalWallFaceID[iFace] == -1 && globalWallFaceID != -1) {
		  m_innerFaceIDToGlobalWallFaceID[iFace] = globalWallFaceID;
		}
		else {
		  cout << "m_innerFaceIDToGlobalWallFaceID[iFace] == -1 && globalWallFaceID == -1" << endl; 
		  abort();
		}
	      }
	      else {
		cf_assert(globalWallFaceID != -1);
		globalFaceID = globalWallFaceID;
	      }
	      
	      const CFuint globalLineCellID = globalFaceID*m_nbCellsOnLine + cellsOnLine;
	      // cout << "globalLineCellID = " << globalLineCellID << endl;
	      cf_assert(globalLineCellID < lineGlobalCellIDsInTRS.size());
	      lineGlobalCellIDsInTRS[globalLineCellID] = globalCellID;
	      	      
	      maxID = globalLineCellID;
	      
	      // check this increment
	      cellsOnLine++;
	    }
	    
	    //  if (isPartitionFace[faceID]) {
	    // 	      if (leftID != cellID) {
	    // 		cellID = leftID;
// 		globalCellID = cellFaces[oppositeIFace]->getState(0)->getGlobalID();
// 	      }
	    
// 	      // print cell center coordinates
// 	      // fout << cell->getState(0)->getCoordinates() << endl;
// 	      // store the newly found cell ID
// 	      const CFuint globalFaceID = wallTRS->getGlobalGeoID(iFace);
// 	      const CFuint globalLineCellID = globalFaceID*m_nbCellsOnLine + cellsOnLine;
// 	      //cout << "globalLineCellID = " << globalLineCellID << endl;
// 	      cf_assert(globalLineCellID < lineGlobalCellIDsInTRS.size());
// 	      lineGlobalCellIDsInTRS[globalLineCellID] = globalCellID;
	      
// 	      maxID = globalLineCellID;
	      
// 	      // check this increment
// 	      cellsOnLine++;
// 	    }
	    
	    if (isBoundaryFace) {
	      reachedPartitionFace = true;
	    }
	    
	    // if you reach the partition face without having reached the
	    // number of layers m_maxNbNormalFaces, store additional data
	    if (isPartitionFace[faceID] && countFaces < m_maxNbNormalFaces) {
	      const CFuint nbPFNodes = oppositeFace->nbNodes();
	      for (CFuint in = 0; in < nbPFNodes; ++in) {
		sendData.nodeIDs.push_back
		  (oppositeFace->getNode(in)->getGlobalID());
	      }
	      sendData.nbFaceNodes.push_back(nbPFNodes);
	      sendData.globalCellIDs.push_back(cell->getState(0)->getGlobalID());
	      //sendData.globalFaceID.push_back(globalFaceID);
	      
	      // check the following
	      sendData.nbLayersToAdd.push_back(m_maxNbNormalFaces - countFaces);
	      reachedPartitionFace = true;
	    }
	    
	    break;
	  }
	}
	
	cf_assert(foundFace);
	cf_assert(countMatches == 1);
	
	cellBuilder->releaseGE();
      }
      
      // fout<< endl;
      
      faceBuilder->releaseGE();
      faceData.isBFace = false;
    }
  }
  
  
  cout << m_myRank << " has maxID =  " << maxID << endl;
  
  // fout.close();
}

//////////////////////////////////////////////////////////////////////////////

void Slab1DFVMCCMPI::executeOnTrs()
{
  // careful here !!! one TRS at a time !!!!
  cf_assert(getTrsList().size() == 1);
  
  const CFuint dim = PhysicalModelStack::getActive()->getDim();
  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  const CFuint stride = nbEqs + dim;
  DataHandle<CFreal> qrad = socket_qrad.getDataHandle();
  MPI_Comm comm = PE::GetPE().GetCommunicator();
  
  // set the array statesNodes and the mapping
  auto_ptr<ProxyDofIterator<CFreal> > nstatesProxy
    (new DofDataHandleIterator<CFreal, CFreal>
     (&m_statesNodes, &m_mapGlobalToMeshLineID, nbEqs, dim));
  
  DataHandle<State*, GLOBAL> states = socket_states.getDataHandle();
  const CFuint nbProc = PE::GetPE().GetProcessorCount();
  
  for (CFuint r = 0; r < nbProc; ++r) {
    const CFuint rSize = m_localNbLines[r]*m_nbCellsOnLine;
    m_sendArray.resize(rSize);
    
    if (r == m_myRank) {
      // fill the send array
      for (CFuint s = 0; s < states.size(); ++s) {
	CFuint start = s*stride;
	State& currState = *states[s];  
	for (CFuint iEq = 0; iEq < nbEqs; ++iEq, start++) {
	  m_sendArray[start] = currState[iEq];
	}
	RealVector& node = currState.getCoordinates();
	for (CFuint iDim = 0; iDim < dim; ++iDim, start++) {
	  m_sendArray[start] = node[iDim];
	}
      } 
      
      m_sendGlobalIDs.resize(states.size());
      for (CFuint i = 0; i < states.size(); ++i) {
	m_sendGlobalIDs[i] = states[i]->getGlobalID();
      }
    }
    
    const CFuint nbStatesInProc = m_nbStatesInProc[r];
    
    MPIStruct ms;
    // broadcast the states + coordinates of the cell centers
    int lnn[2];
    lnn[0] = m_sendArray.size();
    lnn[1] = nbStatesInProc;
    MPIStructDef::buildMPIStruct<CFreal,CFuint>
      (&m_sendArray[0], &m_sendGlobalIDs[0], lnn, ms);
    
    // make the other processors know the size of the data to distribute
    // from the current processor
    MPI_Bcast(ms.start, 1, ms.type, r, comm);
    
    for (CFuint i = 0; i < nbStatesInProc; ++i) {
      bool isFound = false;
      const CFuint currGlobalID = m_sendGlobalIDs[i];
      
      pair<MapItr, MapItr> ids = m_mapGlobalToMeshLineID.find(currGlobalID, isFound);
      if (isFound) {
	const CFuint start = i*stride;
	// get the data (state + its position) corresponding to the current global ID
	const CFuint meshByLineID = ids.first->second;
	const CFuint startLocal = meshByLineID*stride;
	for (CFuint is = 0; is < stride; ++is) {
	  m_statesNodes[startLocal + is] = m_sendArray[start + is]; 
	}
      }
    }
  }
  
  // if (computeOnStagnationLine()) {
  //   m_sendArrayradLibrary->runOnStagnationLine(m_stagnationLineCells,
  //                                              nstatesProxy.get(),
  //                                              &m_qRadByLine[0]); 
  //}
  //else {
  
  // for the moment we use barriers, one processor at a time launches PARADE
  if (PE::GetPE().IsParallel()) {
    
    PE::GetPE().setBarrier();
    
    for (CFuint i = 0; i < PE::GetPE().GetProcessorCount(); ++i) {
      
      if (i == PE::GetPE().GetRank ()) {
	m_radLibrary->runOnStructuredMesh(m_meshByLines,
					  nstatesProxy.get(),
					  &m_qRadByLine[0]); 
      }
      PE::GetPE().setBarrier();
    }
  }
  else {
    m_radLibrary->runOnStructuredMesh(m_meshByLines,
				      nstatesProxy.get(),
				      &m_qRadByLine[0]); 
  }
  
  CFuint countSend = 0;
  for (CFuint r = 0; r < nbProc; ++r) {
    const CFuint nr = m_donorIDsPerRank[r].size();
    for (CFuint i = 0; i < nr ; ++i, ++countSend) {
      m_qRadToSend[countSend] = m_qRadByLine[m_donorIDsPerRank[r][i]]; 
    }
  }  
  
  MPI_Alltoallv(&m_qRadToSend[0], &m_sendCount[0],
		&m_sendDispl[0], MPI_DOUBLE, &m_qRadToRecv[0],
		&m_recvCount[0], &m_recvDispl[0],
		MPI_DOUBLE, comm);
  
  cf_assert(m_qRadToRecv.size() == qrad.size());
  for (CFuint i = 0; i < qrad.size(); ++i) {
    qrad[m_destIDsToRecv[i]] = m_qRadToRecv[i];
  }
}
  
//////////////////////////////////////////////////////////////////////////////

  } // namespace RadiativeTransfer

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
