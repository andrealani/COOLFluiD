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

MethodCommandProvider<BCPeriodic, CellCenterFVMData, FiniteVolumeModule> BCPeriodicFVMCCProvider("BCPeriodicFVMCC");

//////////////////////////////////////////////////////////////////////////////

void BCPeriodic::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< CFreal >("Threshold","Tolerance for considering coordinates matching.");
  options.addConfigOption< std::vector<CFreal> >("TranslationVector","TranslationVector between planes");  
}

//////////////////////////////////////////////////////////////////////////////

BCPeriodic::BCPeriodic(const std::string& name):
  FVMCC_BC(name),
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
  CFLog(INFO, "Setup ["<<getClassName()<<"] \n");
  
  /* Setup parent class */
  FVMCC_BC::setup();
  
  /* Check if translation vector has right dimension */
  CFuint dim = PhysicalModelStack::getActive()->getDim();
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
  std::vector<FaceStruct> eastFaces;
  _localWestEastMap.reserve(_nbTrsFaces);
  for (CFuint iFace = 0; iFace<_nbTrsFaces; iFace++){

    /* Face setup */
    faceData.idx = iFace;
    GeometricEntity *const face = _faceBuilder.buildGE();
    const CFuint faceGlobalID = face->getID();
    _globalToLocalTRSFaceID.insert(faceGlobalID,iFace);


    /* Calculate coordinate */
    const std::vector<Node*>& nodes = *face->getNodes();
    RealVector centreNode(dim);
    for(CFuint iNode=0; iNode<nodes.size(); ++iNode) {
      for(CFuint iDim=0; iDim<dim; ++iDim) {
        centreNode[iDim] += (*(nodes[iNode]))[iDim]/nodes.size();
      }
    }
    
    /* Load this face in a structure */
    FaceStruct thisFace;
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
    }
    else {
      westFaces.push_back(thisFace);
      _localWestEastMap.insert(iFace,0);
    }

    /* release geometry */
    _faceBuilder.releaseGE();

  }
  _localWestEastMap.sortKeys();
    
    
  /* Gather the number of Faces for each process */
  CFuint nbWestFaces = westFaces.size();
  CFuint nbEastFaces = eastFaces.size();
  
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
    
  /* Definition of MPI struct */
  FaceMPIStruct faceMPIStruct;

  /* 
   * Assemble local connectivity map using MPI. It stores for each local face
   * on which processor the periodic face is, and its local and global ID
   */
  _localConnectivityMap.resize(_nbTrsFaces);
CFuint matches(0);
  for(CFuint iP=0; iP<_nbProcesses; ++iP) {

    /* Find WestFace connectivity */
    for(CFuint iFace=0; iFace<nbWestFacesPerProcess[iP]; iFace++) {
    
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
    for(CFuint iFace=0; iFace<nbEastFacesPerProcess[iP]; iFace++) {

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
  
  MPI_Barrier(_comm);    
  
  if(matches!=_nbTrsFaces)
    CFLog(INFO, "Only " << matches << "/"<< _nbTrsFaces << " matches found. Wrong TranslationVector or try increasing Threshold \n");  
}

//////////////////////////////////////////////////////////////////////////////

void BCPeriodic::setupMPI() 
{
  _nbProcesses = PE::GetPE().GetProcessorCount();
  _rank = PE::GetPE().GetRank();
  _comm = PE::GetPE().GetCommunicator();
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
  _recvbuf.resize(_LastDisplacement);
  for(CFuint i=0; i<_LastDisplacement; i++){
    _recvbuf[i] = 0.0;
  }
  _sendbuf.reserve(_nbTrsFaces*_nE*_nbProcesses);
  
}


//////////////////////////////////////////////////////////////////////////////


void BCPeriodic::preProcess()
{
  if (_nbProcesses > 1) {
    
    _sendbuf.clear();
    
    FaceTrsGeoBuilder::GeoData& faceData = _faceBuilder.getDataGE();
    faceData.isBFace = true;

    std::vector<CFreal> boundaryState;
    boundaryState.reserve(_nbTrsFaces*_nE);  
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

    for(CFuint iP=0; iP<_nbProcesses; iP++){
      for(CFuint i=0; i<_nbTrsFaces*_nE; i++){
        _sendbuf.push_back(boundaryState[i]);
      }
    }

    MPI_Barrier(_comm);    
    
    MPI_Datatype MPI_CFREAL = Common::MPIStructDef::getMPIType(&_sendbuf[0]);
    
    MPI_Alltoallv(&_sendbuf[0], &_sendcounts[0], &_senddispls[0], MPI_CFREAL,
		  &_recvbuf[0], &_recvcounts[0], &_recvdispls[0], MPI_CFREAL, _comm);
    MPI_Barrier(_comm);
  }
}

//////////////////////////////////////////////////////////////////////////////

Framework::State* BCPeriodic::computePeriodicState(GeometricEntity *const face)
{
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
    return periodicState;
  }
  else {
    const CFuint periodicProcess = _localConnectivityMap[faceLocalID][j].getProcess();
    for(CFuint h=0; h<_nE; h++){
      (*_periodicState)[h] = _recvbuf[_recvdispls[periodicProcess] + _nE*periodicFaceID + h];
    }
    return _periodicState;
  }
}

//////////////////////////////////////////////////////////////////////////////

bool BCPeriodic::isOutletFace(GeometricEntity *const face)
{
  const CFuint faceGlobalID = face->getID();
  const CFuint faceLocalID = _globalToLocalTRSFaceID.find(faceGlobalID);
  return _localWestEastMap[faceLocalID];
}

//////////////////////////////////////////////////////////////////////////////

bool BCPeriodic::isInletFace(GeometricEntity *const face)
{
  const CFuint faceGlobalID = face->getID();
  const CFuint faceLocalID = _globalToLocalTRSFaceID.find(faceGlobalID);
  return !_localWestEastMap[faceLocalID];
}

//////////////////////////////////////////////////////////////////////////////

bool BCPeriodic::isOutletFace(const CFuint& faceLocalID)
{
  return _localWestEastMap[faceLocalID];
}

//////////////////////////////////////////////////////////////////////////////

bool BCPeriodic::isInletFace(const CFuint& faceLocalID)
{
  return !_localWestEastMap[faceLocalID];
}

//////////////////////////////////////////////////////////////////////////////

void BCPeriodic::setGhostState(GeometricEntity *const face)
{
  State *const ghostState = face->getState(1);
  State *const periodicState = computePeriodicState(face);
  for(CFuint h=0; h<_nE; h++){
    (*ghostState)[h] = (*periodicState)[h];
  }
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

