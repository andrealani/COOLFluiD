#include "MathTools/MathChecks.hh"
#include "Framework/MethodCommandProvider.hh"
#include "FiniteVolume/FiniteVolume.hh"
#include "FiniteVolume/FVMCC_BCPeriodic.hh"
#include "Common/MPI/MPIStructDef.hh"

#include <iostream>
#include <map>
#include <algorithm>

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<BCPeriodic, CellCenterFVMData, FiniteVolumeModule> 
BCPeriodicFVMCCProvider("BCPeriodicFVMCC");

//////////////////////////////////////////////////////////////////////////////

void BCPeriodic::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< CFreal >("Threshold","Tolerance for considering coordinates matching.");
  options.addConfigOption< std::vector<CFreal> >("TranslationVector","TranslationVector between planes");  
}

//////////////////////////////////////////////////////////////////////////////

BCPeriodic::BCPeriodic(const std::string& name):
  FVMCC_BC(name),
  socket_uX("uX", false),
  socket_uY("uY", false),
  socket_uZ("uZ", false), 
  socket_limiter("limiter"),
  _uxExists(false),
  _uyExists(false),
  _uzExists(false),
  _faceBuilder(),
  _translationVector(),
  _localConnectivityMap()
{
  addConfigOptionsTo(this);
  _threshold = 1e-5;
  setParameter("Threshold",&_threshold);
  setParameter("TranslationVector",&_translationVector);
}

//////////////////////////////////////////////////////////////////////////////

BCPeriodic::~BCPeriodic()
{
  delete _periodicState;
}

//////////////////////////////////////////////////////////////////////////////

void BCPeriodic::setup()
{
  CFLog(VERBOSE, "BCPeriodic::setup() => start\n");
  
  CFLog(INFO, "Setup ["<<getClassName()<<"] \n");
  
  /* Setup parent class */
  FVMCC_BC::setup();
  
  const string namespaceName = MeshDataStack::getActive()->getPrimaryNamespace();
  const string uxName = namespaceName + "_uX";
  _uxExists = MeshDataStack::getActive()->getDataStorage()->checkData(uxName);
  const string uyName = namespaceName + "_uY";
  _uyExists = MeshDataStack::getActive()->getDataStorage()->checkData(uyName);
  const string uzName = namespaceName + "_uZ";
  _uzExists = MeshDataStack::getActive()->getDataStorage()->checkData(uzName);
  
  /* Check if translation vector has right dimension */
  const CFuint dim = PhysicalModelStack::getActive()->getDim();
  cf_assert(_translationVector.size()==dim);
  
  /* Convert std::vector<CFreal> to RealVector version */
  RealVector translationVector(dim), backTranslationVector(dim);
  for(CFuint iDim=0; iDim<dim; ++iDim) {
    translationVector[iDim] = _translationVector[iDim];
    backTranslationVector[iDim] = -translationVector[iDim];
  }
  CFLog(INFO, " ==> translationVector = " << translationVector << " \n");

  /* Preparing setup for faceData */
  SafePtr<TopologicalRegionSet> trs = getCurrentTRS();
  _faceBuilder.setup();
  _faceBuilder.getGeoBuilder()->setDataSockets(socket_states, socket_gstates, socket_nodes);
  FaceTrsGeoBuilder::GeoData& faceData = _faceBuilder.getDataGE();
  faceData.isBFace = true;
  faceData.trs = trs;
  _nbTrsFaces = trs->getLocalNbGeoEnts();
  _nE = PhysicalModelStack::getActive()->getNbEq();
  
  /* Setup MPI-related parameters */
  setupMPI();
  
  
  /* Check for every local face if it is an east or a west face */ 
  std::vector<FaceStruct> westFaces;
  westFaces.reserve(_nbTrsFaces); // oversized to avoid frequent memory reallocation
  
  std::vector<FaceStruct> eastFaces;
  eastFaces.reserve(_nbTrsFaces); // oversized to avoid frequent memory reallocation
  
  _localWestEastMap.reserve(_nbTrsFaces);
  
  DataHandle<CFreal> normals = socket_normals.getDataHandle();
  RealVector faceNormal(dim);
  RealVector centreNode(dim);
  FaceStruct thisFace;

  Stopwatch<WallTime> stp;
  stp.start();
  
  //std::cout<<"_nbTrsFaces = "<<_nbTrsFaces<<"\n"; 
  for (CFuint iFace = 0; iFace<_nbTrsFaces; iFace++){
    
    /* Face setup */
    faceData.idx = iFace;
    GeometricEntity *const face = _faceBuilder.buildGE();
    const CFuint faceGlobalID = face->getID();
    _globalToLocalTRSFaceID.insert(faceGlobalID,iFace);
    
    
    /* Calculate coordinate */
    const std::vector<Node*>& nodes = *face->getNodes();
    const CFuint nbNodesInFace = nodes.size();
    centreNode = 0.;
    for(CFuint iNode=0; iNode< nbNodesInFace; ++iNode) {
      centreNode += *(nodes[iNode]);
    }
    centreNode /= nbNodesInFace;
    
    /* Load this face in a structure */
    thisFace.setCentreCoordinates(centreNode);
    thisFace.setGlobalFaceID(faceGlobalID);
    thisFace.setLocalFaceID(iFace);
    

    /* Check if the projection of the face-normal on the translation vector
     * is positive (westFace) or negative (eastFace)
     */
    const CFuint startID = faceGlobalID*dim;
    DataHandle<CFreal> normals = socket_normals.getDataHandle();
    RealVector faceNormal(dim);
    for(CFuint iDim=0; iDim<dim; ++iDim) {
      faceNormal[iDim]=normals[startID+iDim];
    }
    
    CFreal proj = MathTools::MathFunctions::innerProd(translationVector, faceNormal);
    if (proj > 0) {
      eastFaces.push_back(thisFace);
      _localWestEastMap.insert(iFace,1);
      //std::cout<<"EAST \t rank = "<<_rank<<" faceGlobalID["<<iFace<<"] = "<< faceGlobalID<<"\n";
    }
    else if (proj < 0){
      westFaces.push_back(thisFace);
      _localWestEastMap.insert(iFace,0);
      //std::cout<<"WEST \t rank = "<<_rank<<" faceGlobalID["<<iFace<<"] = "<< faceGlobalID<<"\n";
    }
    
    /* release geometry */
    _faceBuilder.releaseGE();
    
  }
  _localWestEastMap.sortKeys();
  
  CFLog(INFO,"BCPeriodic::setup() => step 1 took " << stp.read() << "s\n");
  stp.start();
   
  /* Gather the number of Faces for each process */
  CFuint nbWestFaces = westFaces.size();
  cf_assert(nbWestFaces <= westFaces.capacity());
  CFuint nbEastFaces = eastFaces.size();
  cf_assert(nbEastFaces <= eastFaces.capacity());
  
  std::vector<CFuint> nbWestFacesPerProcess(_nbProcesses,0);
  std::vector<CFuint> nbEastFacesPerProcess(_nbProcesses,0);
  MPI_Allgather(&nbWestFaces, 1, MPIStructDef::getMPIType(&nbWestFaces), 
	&nbWestFacesPerProcess[0], 1, MPIStructDef::getMPIType(&nbWestFaces), _comm);
  MPI_Allgather(&nbEastFaces, 1, MPIStructDef::getMPIType(&nbEastFaces), 
	&nbEastFacesPerProcess[0], 1, MPIStructDef::getMPIType(&nbEastFaces), _comm);
  
  CFuint totalNbWestFaces = 0;
  CFuint totalNbEastFaces = 0;
  
  if (_nbProcesses > 1) {
    CFLog(INFO, " ==> number of west faces on each processor: \n");
    for(CFuint iP=0; iP<_nbProcesses; ++iP) {
      CFLog(INFO, "     "<<iP<<" --> " << nbWestFacesPerProcess[iP] << " \n");
      totalNbWestFaces += nbWestFacesPerProcess[iP];
    }
    CFLog(INFO, "     sum = " << totalNbWestFaces << " \n");
    CFLog(INFO, " ==> number of east faces on each processor: \n");
    for(CFuint iP=0; iP<_nbProcesses; ++iP) {
      CFLog(INFO, "     "<<iP<<" --> " << nbEastFacesPerProcess[iP] << " \n");
      totalNbEastFaces += nbEastFacesPerProcess[iP];
    }
    CFLog(INFO, "     sum = " << totalNbEastFaces << " \n");
  }
  else{
    CFLog(INFO, " ==> number of west faces: " << nbWestFaces << " \n");
    CFLog(INFO, " ==> number of east faces: " << nbEastFaces << " \n");
  }
  
  CFLog(INFO,"BCPeriodic::setup() => step 2 took " << stp.read() << "s\n");
  stp.start();
  
  /* Definition of MPI struct */
  FaceMPIStruct faceMPIStruct;

  /* 
   * Assemble local connectivity map using MPI. It stores for each local face
   0* on which processor the periodic face is, and its local and global ID
   */
  _localConnectivityMap.resize(_nbTrsFaces);

  //Sending information of all the faces
  CFuint matches = 0;
  if(_nbProcesses == 1){
    for(CFuint iP=0; iP<_nbProcesses; ++iP) {

      /* Find WestFace connectivity */
      const CFuint nbW = nbWestFacesPerProcess[iP];
      for(CFuint iFace=0; iFace<nbW; iFace++) {

	/* fill the mpi_struct */
	if(iP == _rank) faceMPIStruct.copy(westFaces[iFace]);

	/* MPI Broadcast struct */
	faceMPIStruct.broadcast(iP);

	/* Unload the buffer */
	FaceStruct westFace = static_cast<FaceStruct>(faceMPIStruct);

	/* Find eastFace corresponding to this westFace */
	std::vector<FaceStruct>::iterator found = std::find_if(eastFaces.begin(), eastFaces.end(), findTranslated(westFace,translationVector,_threshold));
	if (found != eastFaces.end()) {
	  if (_localConnectivityMap[found->getLocalFaceID()].size() == 0)   matches++;
	  _localConnectivityMap[found->getLocalFaceID()].push_back(PairStruct(iP,westFace.getLocalFaceID()));
	}

      } // end westFace connectivity


      /* Find eastFace connectivity */
      const CFuint nbE = nbEastFacesPerProcess[iP];
      for(CFuint iFace=0; iFace<nbE; iFace++) {

	/* fill the mpi_struct */
	if(iP == _rank) faceMPIStruct.copy(eastFaces[iFace]);

	/* MPI Broadcast struct */
	faceMPIStruct.broadcast(iP);

	/* Unload the buffer */
	FaceStruct eastFace = static_cast<FaceStruct>(faceMPIStruct);

	/* Find westFace corresponding to this eastFace */
	std::vector<FaceStruct>::iterator found = std::find_if(westFaces.begin(), westFaces.end(), findTranslated(eastFace,backTranslationVector,_threshold));
	if (found != westFaces.end()) {
	  if (_localConnectivityMap[found->getLocalFaceID()].size() == 0)   matches++;
	  _localConnectivityMap[found->getLocalFaceID()].push_back(PairStruct(iP,eastFace.getLocalFaceID()));
	}
      } // end eastFace connectivity

    } // end Assembly of localConnectivityMap

  }
  else { 
    std::vector<CFuint> sendGlobalIDW(nbWestFaces, 0);
    std::vector<CFuint> sendLocalIDW(nbWestFaces, 0);
    std::vector<CFreal> sendCoordsW(nbWestFaces*dim, 0);
    std::vector<CFuint> recvGlobalIDW(totalNbWestFaces, 0);
    std::vector<CFuint> recvLocalIDW(totalNbWestFaces, 0);
    std::vector<CFreal> recvCoordsW(totalNbWestFaces*dim, 0);

    std::vector<int> displsW(_nbProcesses, 0);
    std::vector<int> rcountsW(_nbProcesses, 0);
    std::vector<int> displsWDim(_nbProcesses, 0);
    std::vector<int> rcountsWDim(_nbProcesses, 0);

    int localLengthW = nbWestFaces;
    int localLengthWDim = nbWestFaces*dim;

    std::vector<CFuint> sendGlobalIDE(nbEastFaces, 0);
    std::vector<CFuint> sendLocalIDE(nbEastFaces, 0);
    std::vector<CFreal> sendCoordsE(nbEastFaces*dim, 0);
    std::vector<CFuint> recvGlobalIDE(totalNbEastFaces, 0);
    std::vector<CFuint> recvLocalIDE(totalNbEastFaces, 0);
    std::vector<CFreal> recvCoordsE(totalNbEastFaces*dim, 0);

    std::vector<int> displsE(_nbProcesses, 0);
    std::vector<int> rcountsE(_nbProcesses, 0);
    std::vector<int> displsEDim(_nbProcesses, 0);
    std::vector<int> rcountsEDim(_nbProcesses, 0);

    int localLengthE = nbEastFaces;
    int localLengthEDim = nbEastFaces*dim;
    
    displsW[0] = 0;
    displsE[0] = 0;
    for (CFuint iP=1;iP<_nbProcesses;iP++) {
      displsW[iP] = displsW[iP-1] + nbWestFacesPerProcess[iP-1];
      displsWDim[iP] = displsW[iP]*dim;
      displsE[iP] = displsE[iP-1] + nbEastFacesPerProcess[iP-1];
      displsEDim[iP] = displsE[iP]*dim;
    }

    for(CFuint iP=0; iP<_nbProcesses; ++iP) {
      rcountsW[iP] = nbWestFacesPerProcess[iP];
      rcountsWDim[iP] = rcountsW[iP]*dim;
      rcountsE[iP] = nbEastFacesPerProcess[iP];
      rcountsEDim[iP] = rcountsE[iP]*dim;
    }
      
    /* Fill the send vector */
    for(CFuint iFace=0; iFace<nbWestFaces; iFace++) {

      /* fill the sentd info */
      sendGlobalIDW[iFace] = westFaces[iFace].getGlobalFaceID();
      sendLocalIDW[iFace] = westFaces[iFace].getLocalFaceID();
      for (CFuint iDim = 0; iDim<dim; iDim++){
	sendCoordsW[iFace*dim + iDim] = westFaces[iFace].getCentreCoordinates()[iDim];
      }
    }
    for(CFuint iFace=0; iFace<nbEastFaces; iFace++) {

      /* fill the sentd info */
      sendGlobalIDE[iFace] = eastFaces[iFace].getGlobalFaceID();
      sendLocalIDE[iFace] = eastFaces[iFace].getLocalFaceID();
      for (CFuint iDim = 0; iDim<dim; iDim++){
	sendCoordsE[iFace*dim + iDim] = eastFaces[iFace].getCentreCoordinates()[iDim];
      }
    }
    MPI_Barrier(_comm);

    MPI_Datatype MPI_CFUINT = Common::MPIStructDef::getMPIType(&sendGlobalIDW[0]);
    MPI_Datatype MPI_CFREAL = Common::MPIStructDef::getMPIType(&sendCoordsW[0]);
    MPI_Allgatherv(&sendGlobalIDW[0], localLengthW, MPI_CFUINT,
		    &recvGlobalIDW[0], &rcountsW[0], &displsW[0], MPI_CFUINT, _comm);
    MPI_Allgatherv(&sendLocalIDW[0], localLengthW, MPI_CFUINT,
		    &recvLocalIDW[0], &rcountsW[0], &displsW[0], MPI_CFUINT, _comm);
    MPI_Allgatherv(&sendCoordsW[0], localLengthWDim, MPI_CFREAL,
		    &recvCoordsW[0], &rcountsWDim[0], &displsWDim[0], MPI_CFREAL, _comm);
    MPI_Allgatherv(&sendGlobalIDE[0], localLengthE, MPI_CFUINT,
		    &recvGlobalIDE[0], &rcountsE[0], &displsE[0], MPI_CFUINT, _comm);
    MPI_Allgatherv(&sendLocalIDE[0], localLengthE, MPI_CFUINT,
		    &recvLocalIDE[0], &rcountsE[0], &displsE[0], MPI_CFUINT, _comm);
    MPI_Allgatherv(&sendCoordsE[0], localLengthEDim, MPI_CFREAL,
		    &recvCoordsE[0], &rcountsEDim[0], &displsEDim[0], MPI_CFREAL, _comm);
    MPI_Barrier(_comm);

    // Constructing all the faces
    /// WEST
    std::vector<FaceStruct> totalWestFaces;
    RealVector coordinates(dim);
    totalWestFaces.resize(totalNbWestFaces);
    for (CFuint iFace=0; iFace<totalNbWestFaces; iFace++) {
      for (CFuint iDim=0; iDim<dim; iDim++) {
	coordinates[iDim] = recvCoordsW[iFace*dim + iDim];
      }
      /* Load this face in a structure */
      totalWestFaces[iFace].setCentreCoordinates(coordinates);
      totalWestFaces[iFace].setGlobalFaceID(recvGlobalIDW[iFace]);
      totalWestFaces[iFace].setLocalFaceID(recvLocalIDW[iFace]);
    }
    ///EAST
    std::vector<FaceStruct> totalEastFaces(totalNbEastFaces);
    for (CFuint iFace=0; iFace<totalNbEastFaces; iFace++) {
      for (CFuint iDim=0; iDim<dim; iDim++) {
	coordinates[iDim] = recvCoordsE[iFace*dim + iDim]; 
      }
      /* Load this face in a structure */
      totalEastFaces[iFace].setCentreCoordinates(coordinates);
      totalEastFaces[iFace].setGlobalFaceID(recvGlobalIDE[iFace]);
      totalEastFaces[iFace].setLocalFaceID(recvLocalIDE[iFace]);
    }

    /* Find WestFace connectivity */
    for(CFuint iP=0; iP<_nbProcesses; iP++){
      const CFuint nbW = nbWestFacesPerProcess[iP]; 
      for(CFuint iFace=0; iFace<nbW; iFace++) {
	FaceStruct westFace = totalWestFaces[displsW[iP] + iFace];
	/* Find eastFace corresponding to this westFace */
	std::vector<FaceStruct>::iterator found = std::find_if(eastFaces.begin(), eastFaces.end(), findTranslated(westFace,translationVector,_threshold));
	if (found != eastFaces.end()) {
	  if (_localConnectivityMap[found->getLocalFaceID()].size() == 0)   matches++;
	  _localConnectivityMap[found->getLocalFaceID()].push_back(PairStruct(iP,westFace.getLocalFaceID()));
	}
      } 
      
      /* Find eastFace connectivity */
      const CFuint nbE = nbEastFacesPerProcess[iP];
      for(CFuint iFace=0; iFace<nbE; iFace++) {
	//CFLog(INFO,"iD = "<<displsE[iP] + iFace <<"\n");
	FaceStruct eastFace = totalEastFaces[displsE[iP] + iFace]; 
	
	/* Find westFace corresponding to this eastFace */
	std::vector<FaceStruct>::iterator found = std::find_if(westFaces.begin(), westFaces.end(), findTranslated(eastFace,backTranslationVector,_threshold));
	if (found != westFaces.end()) {
	  if (_localConnectivityMap[found->getLocalFaceID()].size() == 0)   matches++;
	  _localConnectivityMap[found->getLocalFaceID()].push_back(PairStruct(iP,eastFace.getLocalFaceID()));
	}
      } // end eastFace connectivity
    }
    
    MPI_Barrier(_comm);    
  }
  
  CFLog(INFO,"BCPeriodic::setup() => step 3 took " << stp.read() << "s\n");
  
  if(matches!=_nbTrsFaces)
    CFLog(INFO, "Only " << matches << "/"<< _nbTrsFaces << " matches found. Wrong TranslationVector or try increasing Threshold \n");  
  
  CFLog(VERBOSE, "BCPeriodic::setup() => end\n");
}

//////////////////////////////////////////////////////////////////////////////

void BCPeriodic::setupMPI() 
{
  const std::string nsp = MeshDataStack::getActive()->getPrimaryNamespace();
  _nbProcesses = PE::GetPE().GetProcessorCount(nsp);
  _rank = PE::GetPE().GetRank(nsp);
  _comm = PE::GetPE().GetCommunicator(nsp);
  _nbFacesPerProcess.resize(_nbProcesses,0);
  _periodicState = new Framework::State;
  _periodicState->resize(_nE);
  
  MPI_Barrier(_comm);
  MPI_Allgather(&_nbTrsFaces, 1, MPIStructDef::getMPIType(&_nbTrsFaces), 
	&_nbFacesPerProcess[0], 1, MPIStructDef::getMPIType(&_nbTrsFaces), _comm);

  
  // build counts and displacements for preProcess step
  
  _sendcounts.resize(_nbProcesses);
  _recvcounts.resize(_nbProcesses);
  _recvdispls.resize(_nbProcesses);
  _senddispls.resize(_nbProcesses);
  
  _recvdispls[0] = 0;
  _senddispls[0] = 0;
  _LastDisplacement = 0;
  for(CFuint iP=0; iP<_nbProcesses; iP++){
    _sendcounts[iP] = _nE*_nbFacesPerProcess[_rank];
    _recvcounts[iP] = _nE*_nbFacesPerProcess[iP];
  }
  if(_nbProcesses>1){
    for(CFuint iP=1; iP<_nbProcesses; iP++){
      _recvdispls[iP] = _recvdispls[iP-1] + _recvcounts[iP-1];
      _senddispls[iP] = _senddispls[iP-1] + _sendcounts[iP-1];
    }
  }
  _LastDisplacement = _recvdispls[_nbProcesses-1] + _recvcounts[_nbProcesses-1];
  
  _recvbuf.resize(_LastDisplacement, 0.);
  _sendbuf.reserve(_nbTrsFaces*_nE*_nbProcesses);
  
  _recvbufLimiter.resize(_LastDisplacement, 0.);
  _sendbufLimiter.reserve(_nbTrsFaces*_nE*_nbProcesses);
  
  // allocate data for the transfer of cell-based gradients
  const CFuint dim = PhysicalModelStack::getActive()->getDim();
  _recvbufGrad.resize(dim);
  for (CFuint i = 0; i < dim; ++i) {
    _recvbufGrad[i].resize(_LastDisplacement, 0.);
  }
  _sendbufGrad.resize(dim);
  for (CFuint i = 0; i < dim; ++i) {
    _sendbufGrad[i].reserve(_nbTrsFaces*_nE*_nbProcesses);
  }
}


//////////////////////////////////////////////////////////////////////////////

void BCPeriodic::preProcess()
{
  CFLog(VERBOSE, "BCPeriodic::preProcess() => start\n");
  
  if (_nbProcesses > 1) {
    
    _sendbuf.clear();
    
    FaceTrsGeoBuilder::GeoData& faceData = _faceBuilder.getDataGE();
    faceData.isBFace = true;

    std::vector<CFreal> boundaryState;
    if (_nbTrsFaces > 0) {
      boundaryState.reserve(_nbTrsFaces*_nE); 
    }
    
    // Loop for every face
    for(CFuint iFace=0; iFace<_nbTrsFaces; iFace++) {
      faceData.idx = iFace;
      GeometricEntity *const face = _faceBuilder.buildGE();
      State *const innerState = face->getState(0);
      for(CFuint e=0; e<_nE; e++){
        boundaryState.push_back( (*innerState)[e] );
      }
      _faceBuilder.releaseGE();
    }

    const CFuint nbTrsNE = _nbTrsFaces*_nE;
    for(CFuint iP=0; iP<_nbProcesses; iP++){
      for(CFuint i=0; i< nbTrsNE; i++){
        _sendbuf.push_back(boundaryState[i]);
      }
    }

    MPI_Barrier(_comm);    
    
    MPI_Datatype MPI_CFREAL = Common::MPIStructDef::getMPIType(&_sendbuf[0]);
    
    MPI_Alltoallv(&_sendbuf[0], &_sendcounts[0], &_senddispls[0], MPI_CFREAL,
		  &_recvbuf[0], &_recvcounts[0], &_recvdispls[0], MPI_CFREAL, _comm);
    MPI_Barrier(_comm);
  }
  
  CFLog(VERBOSE, "BCPeriodic::preProcess() => end\n");
}

//////////////////////////////////////////////////////////////////////////////

Framework::State* BCPeriodic::computePeriodicState(GeometricEntity *const face)
{
  CFLog(DEBUG_MIN, "BCPeriodic::computePeriodicState() => start\n");
  
  const CFuint faceGlobalID = face->getID();
  const CFuint faceLocalID = _globalToLocalTRSFaceID.find(faceGlobalID);
  const CFuint j = _localConnectivityMap[faceLocalID].size() - 1;
  const CFuint periodicFaceID = _localConnectivityMap[faceLocalID][j].getFaceID();
  
  if (_nbProcesses == 1) {
    FaceTrsGeoBuilder::GeoData& faceData = _faceBuilder.getDataGE();
    faceData.isBFace = true;
    faceData.idx = periodicFaceID;
    GeometricEntity *const periodicFace = _faceBuilder.buildGE();
    State *const periodicState = periodicFace->getState(0);
    _faceBuilder.releaseGE();
    CFLog(DEBUG_MIN, "BCPeriodic::computePeriodicState() => end\n");
    return periodicState;
  }
  else {
    cf_assert(faceLocalID < _localConnectivityMap.size());
    cf_assert(j < _localConnectivityMap[faceLocalID].size()); 
    CFLog(DEBUG_MIN, "BCPeriodic::computePeriodicState() => 1 end\n");
    const CFuint periodicProcess = _localConnectivityMap[faceLocalID][j].getProcess();
    for(CFuint h=0; h<_nE; h++){
      cf_assert(periodicProcess < _recvdispls.size());
      CFLog(DEBUG_MIN, "BCPeriodic::computePeriodicState() => 2 start\n");
      const CFuint idx = _recvdispls[periodicProcess] + _nE*periodicFaceID + h;
      cf_assert(idx < _recvbuf.size());
      (*_periodicState)[h] = _recvbuf[idx];
      CFLog(DEBUG_MIN, "BCPeriodic::computePeriodicState() => 2 end\n");
    }
    CFLog(DEBUG_MIN, "BCPeriodic::computePeriodicState() => end\n");
    return _periodicState;
  }
}

//////////////////////////////////////////////////////////////////////////////

void BCPeriodic::setGhostState(GeometricEntity *const face)
{
  CFLog(DEBUG_MIN, "BCPeriodic::setGhostState() => start\n");
  State *const ghostState = face->getState(1);
  State *const periodicState = computePeriodicState(face);
  for(CFuint h=0; h<_nE; h++){
    (*ghostState)[h] = (*periodicState)[h];
  }
  CFLog(DEBUG_MIN, "BCPeriodic::setGhostState() => end\n");
}
      
//////////////////////////////////////////////////////////////////////////////

void BCPeriodic::transferGradientsData()
{
  const CFuint dim = PhysicalModelStack::getActive()->getDim();
  if (_uxExists) {transferArray(socket_uX.getDataHandle(), _sendbufGrad[XX], _recvbufGrad[XX]);} 
  if (_uyExists && dim >= DIM_2D) {transferArray(socket_uY.getDataHandle(), _sendbufGrad[YY], _recvbufGrad[YY]);}
  if (_uzExists && dim == DIM_3D) {transferArray(socket_uZ.getDataHandle(), _sendbufGrad[ZZ], _recvbufGrad[ZZ]);} 
  
  DataHandle<CFreal> limiter = socket_limiter.getDataHandle();
  if (limiter.size() > 0) {transferArray(limiter, _sendbufLimiter, _recvbufLimiter);}
}
      
//////////////////////////////////////////////////////////////////////////////
      
void BCPeriodic::transferArray(const DataHandle<CFreal>& gradients, 
			       vector<CFreal>& sendbuf, 
			       vector<CFreal>& recvbuf)
{
  if (_nbProcesses > 1) {
    sendbuf.clear();
    
    FaceTrsGeoBuilder::GeoData& faceData = _faceBuilder.getDataGE();
    faceData.isBFace = true;
    
    std::vector<CFreal> boundaryState;
    if (_nbTrsFaces > 0) {
      boundaryState.reserve(_nbTrsFaces*_nE); 
    }
    
    // Loop for every face
    for(CFuint iFace=0; iFace<_nbTrsFaces; iFace++) {
      faceData.idx = iFace;
      GeometricEntity *const face = _faceBuilder.buildGE();
      const CFuint stateID = face->getState(0)->getLocalID();
      const CFuint start = stateID*_nE;
      for(CFuint e=0; e<_nE; e++){
        boundaryState.push_back( gradients[start+e]);
      }
      
      _faceBuilder.releaseGE();
    }
    
    const CFuint nbTrsNE = _nbTrsFaces*_nE;
    for(CFuint iP=0; iP<_nbProcesses; iP++){
      for(CFuint i=0; i< nbTrsNE; i++){
        sendbuf.push_back(boundaryState[i]);
      }
    }

    MPI_Barrier(_comm);    
    
    MPI_Datatype MPI_CFREAL = Common::MPIStructDef::getMPIType(&sendbuf[0]);
    
    MPI_Alltoallv(&sendbuf[0], &_sendcounts[0], &_senddispls[0], MPI_CFREAL,
		  &recvbuf[0], &_recvcounts[0], &_recvdispls[0], MPI_CFREAL, _comm);
    MPI_Barrier(_comm);
  }
}

//////////////////////////////////////////////////////////////////////////////

void BCPeriodic::computePeriodicGradient(GeometricEntity *const face, 
					 vector<CFreal*>& cellGradient,
					 CFreal*& cellLimiter)
{
  CFLog(DEBUG_MIN, "BCPeriodic::computePeriodicGradient() => start\n");
  
  const CFuint dim = PhysicalModelStack::getActive()->getDim();
  CFreal* uX = (_uxExists) ?                  &socket_uX.getDataHandle()[0] : CFNULL;
  CFreal* uY = (_uyExists && dim >= DIM_2D) ? &socket_uY.getDataHandle()[0] : CFNULL;
  CFreal* uZ = (_uzExists && dim == DIM_3D) ? &socket_uZ.getDataHandle()[0] : CFNULL;
  DataHandle<CFreal> limiter = socket_limiter.getDataHandle();
  CFreal* lim = (limiter.size() > 0.) ? &limiter[0] : CFNULL;
  
  const CFuint faceGlobalID = face->getID();
  const CFuint faceLocalID = _globalToLocalTRSFaceID.find(faceGlobalID);
  const CFuint j = _localConnectivityMap[faceLocalID].size() - 1;
  const CFuint periodicFaceID = _localConnectivityMap[faceLocalID][j].getFaceID();
  
  if (_nbProcesses == 1) {
    FaceTrsGeoBuilder::GeoData& faceData = _faceBuilder.getDataGE();
    faceData.isBFace = true;
    faceData.idx = periodicFaceID;
    GeometricEntity *const periodicFace = _faceBuilder.buildGE();
    const CFuint stateID = periodicFace->getState(0)->getLocalID();
    const CFuint start = stateID*_nE;
    if (uX != CFNULL) {cellGradient[XX] = &uX[start];}
    if (uY != CFNULL) {cellGradient[YY] = &uY[start];}
    if (uZ != CFNULL) {cellGradient[ZZ] = &uZ[start];}
    if (limiter.size() > 0.) {cellLimiter = &lim[start];}
    _faceBuilder.releaseGE();
  }
  else {
    cf_assert(faceLocalID < _localConnectivityMap.size());
    cf_assert(j < _localConnectivityMap[faceLocalID].size()); 
    CFLog(DEBUG_MIN, "BCPeriodic::computePeriodicState() => 1 end\n");
    const CFuint periodicProcess = _localConnectivityMap[faceLocalID][j].getProcess();
    cf_assert(periodicProcess < _recvdispls.size());
    const CFuint idx = _recvdispls[periodicProcess] + _nE*periodicFaceID;
    if (uX != CFNULL) {
      cf_assert(idx < _recvbufGrad[XX].size());
      cellGradient[XX] = &_recvbufGrad[XX][idx];
    }
    if (uY != CFNULL) {
      cf_assert(idx < _recvbufGrad[YY].size());
      cellGradient[YY] = &_recvbufGrad[YY][idx];
    }
    if (uZ != CFNULL) {
      cf_assert(idx < _recvbufGrad[ZZ].size());
      cellGradient[ZZ] = &_recvbufGrad[ZZ][idx];
    }    
    if (limiter.size() > 0.) {cellLimiter = &_recvbufLimiter[idx];}
  }
  
  CFLog(DEBUG_MIN, "BCPeriodic::computePeriodicState() => end\n");
}
      
//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> > BCPeriodic::needsSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSink> > result =
    FVMCC_BC::needsSockets();
  
  // Add the needed DataSocketSinks
  result.push_back(&socket_uX);
  result.push_back(&socket_uY);
  result.push_back(&socket_uZ);
  result.push_back(&socket_limiter);
  
  return result;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

