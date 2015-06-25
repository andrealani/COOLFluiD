#include "MathTools/MathChecks.hh"
#include "Framework/MethodCommandProvider.hh"
#include "FiniteVolume/FiniteVolume.hh"
#include "FiniteVolume/PeriodicturboMPI.hh"
#include "mpi.h"
#include "Common/PE.hh"
#include "Common/MPI/MPIStructDef.hh"
#include <iostream>
#include <map>

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<PeriodicturboMPI, CellCenterFVMData, FiniteVolumeModule> periodicturboMPIFVMCCProvider("PeriodicturboMPIFVMCC");

//////////////////////////////////////////////////////////////////////////////

void PeriodicturboMPI::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< CFreal >("Threshold","Tolerance for considering coordinates matching.");
}

//////////////////////////////////////////////////////////////////////////////

PeriodicturboMPI::PeriodicturboMPI(const std::string& name):
  FVMCC_BC(name),
  _comm(),
  _my_nP(),
  _n_P(),
  _nE(), 
  _LastDisplacement(),
  _nbTrsFaces(),
  _nbDim(), 
  _countFpP(),
  _faceBuilder(),
  _globalToLocalTRSFaceID(),
  _ConnectionFacePeriodic(),
  _ConnectionProcessPeriodic(),
  _RecvDis2(),
  _sendcounts2(),
  _recvcounts2(),
  _rdispls2(),
  _sdispls2(),
  _rbuf2(),
  _BoundaryState()
{
  addConfigOptionsTo(this);
  _threshold = 10e-3;
  setParameter("Threshold",&_threshold);
}

//////////////////////////////////////////////////////////////////////////////

PeriodicturboMPI::~PeriodicturboMPI()
{
}

//////////////////////////////////////////////////////////////////////////////
void PeriodicturboMPI::setup()
{
  const std::string nsp = this->getMethodData().getNamespace();
  
  // MPI parameters
  CFuint nP = PE::GetPE().GetProcessorCount(nsp);
  CFuint myP = PE::GetPE().GetRank(nsp);
  _comm = PE::GetPE().GetCommunicator(nsp);
  _my_nP = myP;
  _n_P = nP;

  // TRS
  FVMCC_BC::setup();
  SafePtr<TopologicalRegionSet> trs = getCurrentTRS();
  _faceBuilder.setup();
  _faceBuilder.getGeoBuilder()->setDataSockets(socket_states, socket_gstates, socket_nodes);
  FaceTrsGeoBuilder::GeoData& faceData = _faceBuilder.getDataGE();
  faceData.isBFace = true;
  faceData.trs = trs;
  CFuint nbTrsFaces = trs->getLocalNbGeoEnts();
  _nbTrsFaces = nbTrsFaces;
  std::vector<CFint> nFpP(nP,-1);
  MPI_Allgather(&nbTrsFaces, 1, MPI_UNSIGNED, &nFpP[0], 1, MPI_UNSIGNED, _comm);
  _countFpP = nFpP;
  const CFuint nEquation = PhysicalModelStack::getActive()->getNbEq();
  _nE = nEquation;
  _nbDim = PhysicalModelStack::getActive()->getDim();
  _globalToLocalTRSFaceID.reserve(nbTrsFaces);
  _ConnectionFacePeriodic.reserve(nbTrsFaces);
  _ConnectionProcessPeriodic.reserve(nbTrsFaces);

  // Find which surface (west or east) a face belong to, moreover calculate the coordinate
  Common::CFMap<CFuint,CFuint> WestEastSurfaceMap;
  WestEastSurfaceMap.reserve(_nbTrsFaces);
  CFuint nWf = 0;
  CFuint nEf = 0;
  for (CFuint iFace = 0; iFace<nbTrsFaces; iFace++){

    // Face setup
    CFLogDebugMed( "iFace = " << iFace << "\n");
    faceData.idx = iFace;
    GeometricEntity *const face = _faceBuilder.buildGE();
    const CFuint faceGlobalID = face->getID();
    _globalToLocalTRSFaceID.insert(faceGlobalID,iFace);

    // Normals
    const CFuint startID = faceGlobalID*PhysicalModelStack::getActive()->getDim();
    DataHandle<CFreal> normals = socket_normals.getDataHandle();
    CFreal nx = normals[startID];
    CFreal ny = normals[startID + 1];
    CFreal nz = 0;
    if(_nbDim == 3){
      nz = normals[startID + 2];
    }
    CFreal invFaceLength = 1./sqrt(nx*nx + ny*ny + nz*nz);
    nx *= invFaceLength;
    ny *= invFaceLength;
    nz *= invFaceLength;

    CFint flagSurface = (ny == 0) ? 0 : (ny<0 ? -1 : 1);

    // Calculate the number of face belonging to east and west surface
    if (flagSurface == -1){
      nWf = nWf + 1;
    }
    else {
      nEf = nEf + 1;
    }

    // Build map of west and east surface
    WestEastSurfaceMap.insert(iFace, flagSurface);

    // release geometry
    _faceBuilder.releaseGE();

  }
  _globalToLocalTRSFaceID.sortKeys();
  WestEastSurfaceMap.sortKeys();
  std::vector<CFint> nFWpP(nP,0);
  MPI_Allgather(&nWf, 1, MPI_UNSIGNED, &nFWpP[0], 1, MPI_UNSIGNED, _comm);
  std::vector<CFint> nFEpP(nP,0);
  MPI_Allgather(&nEf, 1, MPI_UNSIGNED, &nFEpP[0], 1, MPI_UNSIGNED, _comm);

  vector<CFreal> FaceHeightCwest;
  FaceHeightCwest.reserve(nWf);
  vector<CFreal> FaceAxialCwest;
  FaceAxialCwest.reserve(nWf);
  vector<CFuint> FaceFacewest;
  FaceFacewest.reserve(nWf);
  vector<CFreal> FaceHeightCeast;
  FaceHeightCeast.reserve(nEf);
  vector<CFreal> FaceAxialCeast;
  FaceAxialCeast.reserve(nEf);
  vector<CFuint> FaceFaceeast;
  FaceFaceeast.reserve(nEf);
  for (CFuint iFace = 0; iFace<nbTrsFaces; iFace++){

    // Face setup
    CFLogDebugMed( "iFace = " << iFace << "\n");
    faceData.idx = iFace;
    GeometricEntity *const face = _faceBuilder.buildGE();
    const CFuint faceGlobalID = face->getID();
    _globalToLocalTRSFaceID.insert(faceGlobalID,iFace);

    // Normals
    const CFuint startID = faceGlobalID*PhysicalModelStack::getActive()->getDim();
    DataHandle<CFreal> normals = socket_normals.getDataHandle();
    CFreal nx = normals[startID];
    CFreal ny = normals[startID + 1];
    CFreal nz = 0;
    if(_nbDim == 3){
      nz = normals[startID + 2];
    }
    CFreal invFaceLength = 1./sqrt(nx*nx + ny*ny + nz*nz);
    nx *= invFaceLength;
    ny *= invFaceLength;
    nz *= invFaceLength;

    // Calculate coordinate
    const vector<Node*>& nodes = *face->getNodes();
    CFreal zCentreCoordinate = 0.;
    CFreal xCentreCoordinate = 0.;
    if(_nbDim == 2){
      zCentreCoordinate = 0;
      xCentreCoordinate = 0.5*( (*(nodes[0]))[XX] + (*(nodes[1]))[XX] );
    }
    if(_nbDim == 3){
      zCentreCoordinate = 0.25*( (*(nodes[0]))[ZZ] + (*(nodes[1]))[ZZ] + (*(nodes[2]))[ZZ] + (*(nodes[3]))[ZZ] );
      xCentreCoordinate = 0.25*( (*(nodes[0]))[XX] + (*(nodes[1]))[XX] + (*(nodes[2]))[XX] + (*(nodes[3]))[XX] );
    }
    CFreal angleCoeff = sin(acos(nz));
    CFreal FaceHeightCoordinate = zCentreCoordinate*angleCoeff;

    CFint flagSurface = WestEastSurfaceMap.find(iFace);
    if (flagSurface == -1){
      FaceHeightCwest.push_back(FaceHeightCoordinate);
      FaceAxialCwest.push_back(xCentreCoordinate);
      FaceFacewest.push_back(iFace);
    }
    else {
      FaceHeightCeast.push_back(FaceHeightCoordinate);
      FaceAxialCeast.push_back(xCentreCoordinate);
      FaceFaceeast.push_back(iFace);
    }

    // release geometry
    _faceBuilder.releaseGE();

  }

  CFuint iterEast = 0;
  CFuint iterWest = 0;
  for(CFuint iP=0; iP<_n_P; iP++){

    MPI_Barrier(_comm);

    // root's parameters for the broadcast
    int root = iP;
    int rootnWf = nFWpP[root];
    int rootnEf = nFEpP[root];

    // build vectors to be sent
    vector<CFreal> sendFHCwest(rootnWf,0);
    vector<CFreal> sendFACwest(rootnWf,0);
    vector<CFuint> sendFFwest(rootnWf,0);

    vector<CFreal> sendFHCeast(rootnEf,0);
    vector<CFreal> sendFACeast(rootnEf,0);
    vector<CFuint> sendFFeast(rootnEf,0);

    if(_my_nP == root){
      for(CFuint i=0; i<rootnWf; i++){
        sendFHCwest[i] = FaceHeightCwest[i];
        sendFACwest[i] = FaceAxialCwest[i];
        sendFFwest[i] = FaceFacewest[i];
      }
      for(CFuint i=0; i<rootnEf; i++){
        sendFHCeast[i] = FaceHeightCeast[i];
        sendFACeast[i] = FaceAxialCeast[i];
        sendFFeast[i] = FaceFaceeast[i];
      }
    }
    
    MPI_Datatype MPI_CFREAL = Common::MPIStructDef::getMPIType(&sendFHCwest[0]);
    
    // Broadcast
    MPI_Bcast(&sendFHCwest[0], rootnWf, MPI_CFREAL, root, _comm);
    MPI_Bcast(&sendFACwest[0], rootnWf, MPI_CFREAL, root, _comm);
    MPI_Bcast(&sendFFwest[0], rootnWf, MPI_UNSIGNED, root, _comm);
    MPI_Bcast(&sendFHCeast[0], rootnEf, MPI_CFREAL, root, _comm);
    MPI_Bcast(&sendFACeast[0], rootnEf, MPI_CFREAL, root, _comm);
    MPI_Bcast(&sendFFeast[0], rootnEf, MPI_UNSIGNED, root, _comm);

    // Build the received vectors
    vector<CFreal> RecvFHCwest(sendFHCwest);
    vector<CFreal> RecvFACwest(sendFACwest);
    vector<CFuint> RecvFFwest(sendFFwest);
    vector<CFreal> RecvFHCeast(sendFHCeast);
    vector<CFreal> RecvFACeast(sendFACeast);
    vector<CFuint> RecvFFeast(sendFFeast);
    
    MPI_Barrier(_comm);

    // Check the connectivity of west faces
    if(nWf>0){
      for(CFuint i=0; i<nWf; i++){
        // data of the face
        CFuint iFaceFFwest = FaceFacewest[i];
        CFreal iFaceFHCwest = FaceHeightCwest[i];
        CFreal iFaceFACwest = FaceAxialCwest[i];
        if(rootnEf>0){
          for(CFuint ii=0; ii<rootnEf; ii++){
            //data of the periodic face in the east surface
            CFuint iPeriodicFaceFFwest = RecvFFeast[ii];
            CFreal iPeriodicFaceFHCwest = RecvFHCeast[ii];
            CFreal iPeriodicFaceFACwest = RecvFACeast[ii];
            if(MathTools::MathChecks::isEqualWithError(iFaceFHCwest,iPeriodicFaceFHCwest, _threshold)){
              if(MathTools::MathChecks::isEqualWithError(iFaceFACwest,iPeriodicFaceFACwest, _threshold)){
                _ConnectionFacePeriodic.insert(iFaceFFwest,iPeriodicFaceFFwest);
                _ConnectionProcessPeriodic.insert(iFaceFFwest,root);
                iterWest++;
                break;
              }
            }
          }
        }
      }
    }

    // Check the connectivity of west faces
    if(nEf>0){
      for(CFuint i=0; i<nEf; i++){
        // data of the face
        CFuint iFaceFFeast = FaceFaceeast[i];
        CFreal iFaceFHCeast = FaceHeightCeast[i];
        CFreal iFaceFACeast = FaceAxialCeast[i];
        if(rootnWf>0){
          for(CFuint ii=0; ii<rootnWf; ii++){
            //data of the periodic face in the east surface
            CFuint iPeriodicFaceFFeast = RecvFFwest[ii];
            CFreal iPeriodicFaceFHCeast = RecvFHCwest[ii];
            CFreal iPeriodicFaceFACeast = RecvFACwest[ii];
            if(MathTools::MathChecks::isEqualWithError(iFaceFHCeast,iPeriodicFaceFHCeast, _threshold)){
              if(MathTools::MathChecks::isEqualWithError(iFaceFACeast,iPeriodicFaceFACeast, _threshold)){
                _ConnectionFacePeriodic.insert(iFaceFFeast,iPeriodicFaceFFeast);
                _ConnectionProcessPeriodic.insert(iFaceFFeast,root);
                iterEast++;
                break;
              }
            }
          }
        }
      }
    }

  }
  _ConnectionFacePeriodic.sortKeys();
  _ConnectionProcessPeriodic.sortKeys();

  // output to verify if, for each face, a periodic face has been found
  // pay attention that it is verified if:
  // 1) iterWest >= nWf (not necessary only equal due to the overposition of processes's mesh partitions)
  // 2) iterEast >= nEf
  // (not necessary only equal due to the overposition of processes's mesh partitions)
  // if it is not verify try to decrease the tolerance
  cout<<" _my_nP = "<<_my_nP<<" iterWest = "<<iterWest<<" nWf "<<nWf<<endl;
  cout<<" _my_nP = "<<_my_nP<<" iterEast = "<<iterEast<<" nEf "<<nEf<<endl;

  // build counts and displacements for preProcess step
  _sendcounts2.resize(_n_P);
  _recvcounts2.resize(_n_P);
  _rdispls2.resize(_n_P);
  _sdispls2.resize(_n_P);
  _rdispls2[0] = 0;
  _sdispls2[0] = 0;
  CFuint _LastDisplacement = 0;
  for(CFuint p=0; p<_n_P; p++){
    _sendcounts2[p] = _nE*_countFpP[_my_nP];
    _recvcounts2[p] = _nE*_countFpP[p];
  }
  if(_n_P>1){
    for(CFuint p=1; p<_n_P; p++){
      _rdispls2[p] = _rdispls2[p-1] + _recvcounts2[p-1];
      _sdispls2[p] = _sdispls2[p-1] + _sendcounts2[p-1];
    }
  }
  _LastDisplacement = _rdispls2[_n_P-1] + _recvcounts2[_n_P-1];
  _rbuf2.resize(_LastDisplacement);
  _BoundaryState.resize(_LastDisplacement);
  for(CFuint i=0; i<_LastDisplacement; i++){
    _rbuf2[i] = 0.0;
    _BoundaryState[i] = 0.0;
  }
  for(CFuint i=0; i<_n_P; i++){
    _RecvDis2.push_back(_rdispls2[i]);
  }
  cout<<" _my_nP = "<<_my_nP<<" adiossssssssssssssssssssssssssssssssssssss 3"<<endl;

}


//////////////////////////////////////////////////////////////////////////////

void PeriodicturboMPI::preProcess()
{

  FaceTrsGeoBuilder::GeoData& faceData = _faceBuilder.getDataGE();
  faceData.isBFace = true;

  const CFuint dimq = _nbTrsFaces*_nE*_n_P;
  const CFuint dimq2 = _nbTrsFaces*_nE;// <<-----------------------------------with allgatherv
  std::vector<CFreal> boundaryState(dimq2);
  for(CFuint iFace=0; iFace<_nbTrsFaces; iFace++){
    faceData.idx = iFace;
    GeometricEntity *const face = _faceBuilder.buildGE();
    State *const innerState = face->getState(0);
    for(CFuint e=0; e<_nE; e++){
      boundaryState[_nE*iFace + e] = (*innerState)[e];
    }
    _faceBuilder.releaseGE();
  }
  std::vector<CFreal> SendBoundaryState;
  SendBoundaryState.reserve(dimq);
  for(CFuint p=0; p<_n_P; p++){
    for(CFuint pp=0; pp<_nbTrsFaces*_nE; pp++){
      SendBoundaryState.push_back(boundaryState[pp]);
    }
  }

  vector<CFreal> sbuf2(SendBoundaryState);
  vector<CFreal> sbuf2AG(boundaryState);// <<-----------------------------------with allgatherv
 
  MPI_Datatype MPI_CFREAL = Common::MPIStructDef::getMPIType(&sbuf2[0]);
  
  MPI_Alltoallv(&sbuf2[0], &_sendcounts2[0], &_sdispls2[0], MPI_CFREAL, &_rbuf2[0], &_recvcounts2[0], &_rdispls2[0], MPI_CFREAL, _comm);
  //MPI_Allgatherv(&sbuf2AG[0], dimq2, MPI_CFREAL, &_rbuf2[0], &_recvcounts2[0], &_rdispls2[0], MPI_CFREAL, _comm);

  _BoundaryState.clear();
  _BoundaryState.reserve(_LastDisplacement);
  for(CFuint i=0; i<_LastDisplacement; i++){
    //_BoundaryState[i] = _rbuf2[i];
    _BoundaryState.push_back(_rbuf2[i]);
  }
}

//////////////////////////////////////////////////////////////////////////////

void PeriodicturboMPI::setGhostState(GeometricEntity *const face)
{

   State *const ghostState = face->getState(1);
   const CFuint faceGlobalID = face->getID();
   const CFuint iFace = _globalToLocalTRSFaceID.find(faceGlobalID);
   const CFuint ProcessPeriodic = _ConnectionProcessPeriodic[iFace];
   const CFuint FacePeriodic = _ConnectionFacePeriodic[iFace];
   for(CFuint h=0; h<_nE; h++){
     //(*ghostState)[h] = _BoundaryState[_RecvDis2[ProcessPeriodic] + _nE*FacePeriodic + h];
     (*ghostState)[h] = _rbuf2[_RecvDis2[ProcessPeriodic] + _nE*FacePeriodic + h];
   }

}

//////////////////////////////////////////////////////////////////////////////


void PeriodicturboMPI::Barrier()
{
  MPI_Barrier(_comm);
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD
















