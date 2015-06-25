#include "MathTools/MathChecks.hh"
#include "Framework/MethodCommandProvider.hh"
#include "FiniteVolume/FiniteVolume.hh"
#include "FiniteVolume/PeriodicX2DMPI.hh"
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

MethodCommandProvider<PeriodicX2DMPI, CellCenterFVMData, FiniteVolumeModule> periodicX2DMPIFVMCCProvider("PeriodicX2DMPIFVMCC");

//////////////////////////////////////////////////////////////////////////////

void PeriodicX2DMPI::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< CFreal >("Threshold","Tolerance for considering coordinates matching.");
}

//////////////////////////////////////////////////////////////////////////////

PeriodicX2DMPI::PeriodicX2DMPI(const std::string& name):
  FVMCC_BC(name),
  _globalToLocalTRSFaceID(),
  _ConnectionFacePeriodic(),
  _ConnectionProcessPeriodic(),
  _RecvDis2(),
  _sendcounts2(),
  _recvcounts2(),
  _rdispls2(),
  _sdispls2(),
  _LastDisplacement(),
  _rbuf2(),
  _BoundaryState(),
  _countFpP(),
  _faceBuilder(),
  _nE(), 
  _nbTrsFaces(),
  _my_nP(),
  _n_P(),
  _comm()
{
  addConfigOptionsTo(this);
  _threshold = 10e-2;
  setParameter("Threshold",&_threshold);
}

//////////////////////////////////////////////////////////////////////////////

PeriodicX2DMPI::~PeriodicX2DMPI()
{
}

//////////////////////////////////////////////////////////////////////////////

void PeriodicX2DMPI::setup()
{
  const std::string nsp = this->getMethodData().getNamespace();
  
  CFuint nP = PE::GetPE().GetProcessorCount(nsp);
  CFuint myP = PE::GetPE().GetRank(nsp);
  _comm = PE::GetPE().GetCommunicator(nsp);
  _my_nP = myP;
  _n_P = nP;
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
  std::vector<CFreal> nodeXCoordinate;
  nodeXCoordinate.reserve(2*nbTrsFaces);
  std::vector<CFreal> nodeYCoordinate;
  nodeYCoordinate.reserve(2*nbTrsFaces);
  std::vector<CFreal> NodeCoordinate;
  NodeCoordinate.reserve(4*nbTrsFaces);
  std::vector<CFuint> FaceGlobalID;
  FaceGlobalID.reserve(nbTrsFaces);
  CFreal nodeBiggerXCoordinate = 0.0;
  CFreal nodeSmallerXCoordinate = 0.0;
  CFreal nodeYCoordinateBiggerX = 0.0;
  CFreal nodeYCoordinateSmallerX = 0.0;
  _globalToLocalTRSFaceID.clear();
  for (CFuint iFace = 0; iFace<nbTrsFaces; iFace++){
     CFLogDebugMed( "iFace = " << iFace << "\n");
     faceData.idx = iFace;
     GeometricEntity *const face = _faceBuilder.buildGE();
     const CFuint faceGlobalID = face->getID();
     _globalToLocalTRSFaceID.insert(faceGlobalID,iFace);
     const vector<Node*>& nodes = *face->getNodes();
     cf_assert(nodes.size() == 2);
     if((*(nodes[0]))[XX] > (*(nodes[1]))[XX]){
       nodeBiggerXCoordinate = (*(nodes[0]))[XX];
       nodeSmallerXCoordinate = (*(nodes[1]))[XX];
       nodeYCoordinateBiggerX = (*(nodes[0]))[YY];
       nodeYCoordinateSmallerX = (*(nodes[1]))[YY];
     }
     else if((*(nodes[1]))[XX] > (*(nodes[0]))[XX]){
	    nodeBiggerXCoordinate = (*(nodes[1]))[XX];
	    nodeSmallerXCoordinate = (*(nodes[0]))[XX];
            nodeYCoordinateBiggerX = (*(nodes[1]))[YY];
            nodeYCoordinateSmallerX = (*(nodes[0]))[YY];
     }
     nodeXCoordinate.push_back(nodeBiggerXCoordinate);
     nodeXCoordinate.push_back(nodeSmallerXCoordinate);
     nodeYCoordinate.push_back(nodeYCoordinateBiggerX);
     nodeYCoordinate.push_back(nodeYCoordinateSmallerX);
     NodeCoordinate.push_back(nodeBiggerXCoordinate);
     NodeCoordinate.push_back(nodeSmallerXCoordinate);
     NodeCoordinate.push_back(nodeYCoordinateBiggerX);
     NodeCoordinate.push_back(nodeYCoordinateSmallerX);
     FaceGlobalID.push_back(faceGlobalID);
     _faceBuilder.releaseGE();
  }
  _globalToLocalTRSFaceID.sortKeys();
  CFuint nC = 4; //total number of coordinate of the nodes belonging to a face
  const CFuint dimq = _n_P*nbTrsFaces*nC;
  const CFuint dimq2 = nbTrsFaces*nC;
  std::vector<CFreal> SendNodeCoordinate;
  SendNodeCoordinate.reserve(dimq);
  for(CFuint p=0; p<_n_P; p++){
    for(CFuint pp=0; pp<dimq2; pp++){
      SendNodeCoordinate.push_back(NodeCoordinate[pp]);
    }
  }
  vector<int> sendcounts(_n_P);
  vector<int> recvcounts(_n_P);
  vector<int> rdispls(_n_P);
  vector<int> sdispls(_n_P);
  rdispls[0] = 0;
  sdispls[0] = 0;
  CFuint LastDisplacement;
  for(CFuint p=0; p<_n_P; p++){
    sendcounts[p] = nC*_countFpP[_my_nP];
    recvcounts[p] = nC*_countFpP[p];
  }
  if(_n_P>1){
    for(CFuint p=1; p<_n_P; p++){
      rdispls[p] = rdispls[p-1] + recvcounts[p-1];
      sdispls[p] = sdispls[p-1] + sendcounts[p-1];
    }
  }
  _ConnectionFacePeriodic.reserve(nbTrsFaces);
  _ConnectionProcessPeriodic.reserve(nbTrsFaces);
  vector<int> RecvDis(rdispls);
  LastDisplacement = rdispls[_n_P-1] + recvcounts[_n_P-1];
  vector<CFreal> rbuf(LastDisplacement,0.0);
  vector<CFreal> sbuf(SendNodeCoordinate);
  
  MPI_Datatype MPI_CFREAL = Common::MPIStructDef::getMPIType(&rbuf[0]);
  
  MPI_Alltoallv(&sbuf[0], &sendcounts[0], &sdispls[0], MPI_CFREAL, &rbuf[0], &recvcounts[0], &rdispls[0], MPI_CFREAL, _comm);
  std::vector<CFreal> RecvNodeCoordinate(rbuf);
  for(CFuint p=0; p<nP; ++p){
     CFuint count = _countFpP[p];
     for(CFuint iFace = 0; iFace<nbTrsFaces; iFace++){
        CFreal nodeBiggerXCoordinate = nodeXCoordinate[2*iFace];
        CFreal nodeSmallerXCoordinate = nodeXCoordinate[2*iFace + 1];
        CFreal nodeYCoordinateBiggerX = nodeYCoordinate[2*iFace];
        CFreal nodeYCoordinateSmallerX = nodeYCoordinate[2*iFace + 1];
        for (CFuint iPeriodicFace = 0; iPeriodicFace<count; ++iPeriodicFace){
          CFreal nodePeriodicBiggerXCoordinate = RecvNodeCoordinate[RecvDis[p] + nC*iPeriodicFace + 0];
          CFreal nodePeriodicSmallerXCoordinate = RecvNodeCoordinate[RecvDis[p] + nC*iPeriodicFace + 1];
          CFreal nodePeriodicYCoordinateBiggerX = RecvNodeCoordinate[RecvDis[p] + nC*iPeriodicFace + 2];
          CFreal nodePeriodicYCoordinateSmallerX = RecvNodeCoordinate[RecvDis[p] + nC*iPeriodicFace + 3];
          if(MathTools::MathChecks::isEqualWithError(nodeBiggerXCoordinate,nodePeriodicBiggerXCoordinate, _threshold) &&
	  MathTools::MathChecks::isEqualWithError(nodeSmallerXCoordinate,nodePeriodicSmallerXCoordinate, _threshold)){
            if(MathTools::MathChecks::isNotEqualWithError(nodeYCoordinateBiggerX,nodePeriodicYCoordinateBiggerX, _threshold) &&     MathTools::MathChecks::isNotEqualWithError(nodeYCoordinateSmallerX,nodePeriodicYCoordinateSmallerX, _threshold)){
              _ConnectionFacePeriodic.insert(iFace,iPeriodicFace);
              _ConnectionProcessPeriodic.insert(iFace,p);
            break;
          }
        }
      }
    }
 }
  _ConnectionFacePeriodic.sortKeys();
  _ConnectionProcessPeriodic.sortKeys();

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
}

//////////////////////////////////////////////////////////////////////////////

void PeriodicX2DMPI::preProcess()
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

void PeriodicX2DMPI::setGhostState(GeometricEntity *const face)
{

   State *const ghostState = face->getState(1);
   const CFuint faceGlobalID = face->getID();
   const CFuint iFace = _globalToLocalTRSFaceID.find(faceGlobalID);
   const CFuint ProcessPeriodic = _ConnectionProcessPeriodic.find(iFace);
   const CFuint FacePeriodic = _ConnectionFacePeriodic.find(iFace);
   for(CFuint h=0; h<_nE; h++){
     //(*ghostState)[h] = _BoundaryState[_RecvDis2[ProcessPeriodic] + _nE*FacePeriodic + h];
     (*ghostState)[h] = _rbuf2[_RecvDis2[ProcessPeriodic] + _nE*FacePeriodic + h];
   }

}

//////////////////////////////////////////////////////////////////////////////


void PeriodicX2DMPI::Barrier()
{
  MPI_Barrier(_comm);
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD
















