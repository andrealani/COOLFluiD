#include "MathTools/MathChecks.hh"
#include "Framework/MethodCommandProvider.hh"
#include "FiniteVolume/FiniteVolume.hh"
#include "FiniteVolume/Periodic3DMPI.hh"
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

MethodCommandProvider<Periodic3DMPI, CellCenterFVMData, FiniteVolumeModule> periodic3DMPIFVMCCProvider("Periodic3DMPIFVMCC");

//////////////////////////////////////////////////////////////////////////////

void Periodic3DMPI::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< CFreal >("Threshold","Tolerance for considering coordinates matching.");
}

//////////////////////////////////////////////////////////////////////////////

Periodic3DMPI::Periodic3DMPI(const std::string& name):
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
  _countNWf(),
  _countNEf(),
  _countFpP(),
  _faceBuilder(),
  _nE(), 
  _nbTrsFaces(),
  _my_nP(),
  _n_P(),
  _comm()
{
  addConfigOptionsTo(this);
  _threshold = 0.0;
  setParameter("Threshold",&_threshold);
}

//////////////////////////////////////////////////////////////////////////////

Periodic3DMPI::~Periodic3DMPI()
{
}

//////////////////////////////////////////////////////////////////////////////

void Periodic3DMPI::setup()
{
  const std::string nsp = getMethodData().getNamespace();
  
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
  _globalToLocalTRSFaceID.clear();

  // count the number of face of west and east surface
  CFuint nWf = 0;
  CFuint nEf = 0;

  for (CFuint iFace = 0; iFace<nbTrsFaces; iFace++){

     // Face setup
     CFLogDebugMed( "iFace = " << iFace << "\n");
     faceData.idx = iFace;
     GeometricEntity *const face = _faceBuilder.buildGE();
     CFuint faceGlobalID = face->getID();
     _globalToLocalTRSFaceID.insert(faceGlobalID,iFace);

     // Normals
     CFuint startID = faceGlobalID*PhysicalModelStack::getActive()->getDim();
     DataHandle<CFreal> normals = socket_normals.getDataHandle();
     CFreal nx = normals[startID];
     CFreal ny = normals[startID + 1];
     CFreal nz = normals[startID + 2];
     CFreal invFaceLength = 1./sqrt(nx*nx + ny*ny + nz*nz);
     nx *= invFaceLength;
     ny *= invFaceLength;
     nz *= invFaceLength;
     CFint flagSurface = (ny == 0) ? 0 : (ny<0 ? -1 : 1);
     if (flagSurface == -1){
       nWf = nWf + 1;
     }
     else {
       nEf = nEf + 1;
     }
     // release geometry
     _faceBuilder.releaseGE();

  }

  std::vector<CFint> NWf(nP,-1);
  MPI_Allgather(&nWf, 1, MPI_UNSIGNED, &NWf[0], 1, MPI_UNSIGNED, _comm);
  _countNWf = NWf;
  std::vector<CFint> NEf(nP,-1);
  MPI_Allgather(&nEf, 1, MPI_UNSIGNED, &NEf[0], 1, MPI_UNSIGNED, _comm);
  _countNEf = NEf;

  // define vector
  vector< pair<CFreal,CFuint> > VecXwest;
  vector< pair<CFreal,CFuint> > VecYwest;
  vector< pair<CFreal,CFuint> > VecZwest;

  vector< pair<CFreal,CFuint> > VecXeast;
  vector< pair<CFreal,CFuint> > VecYeast;
  vector< pair<CFreal,CFuint> > VecZeast;

  VecXwest.reserve(nWf);
  VecYwest.reserve(nWf);
  VecZwest.reserve(nWf);

  VecXeast.reserve(nEf);
  VecYeast.reserve(nEf);
  VecZeast.reserve(nEf);

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
     CFreal nz = normals[startID + 2];
     const CFreal invFaceLength = 1./sqrt(nx*nx + ny*ny + nz*nz);
     nx *= invFaceLength;
     ny *= invFaceLength;
     nz *= invFaceLength;
     CFint flagSurface = (ny == 0) ? 0 : (ny<0 ? -1 : 1);

     // Nodes
     const vector<Node*>& nodes = *face->getNodes();
     const CFuint nbNodesInFace = nodes.size();
     // cout<<" face = "<<iFace<<" nbNodesInFace = "<<nbNodesInFace<<endl;

     CFreal TempNodeXcoordinate = 0.0;
     CFreal TempNodeYcoordinate = 0.0;
     CFreal TempNodeZcoordinate = 0.0;
     // cout<<" Face = "<<iFace<<endl;
     for(CFuint i=0; i<nbNodesInFace; i++){
       CFreal nodeXcoordinate = (*(nodes[i]))[XX];
       CFreal nodeYcoordinate = (*(nodes[i]))[YY];
       CFreal nodeZcoordinate = (*(nodes[i]))[ZZ];
       TempNodeXcoordinate = TempNodeXcoordinate + nodeXcoordinate;
       TempNodeYcoordinate = TempNodeYcoordinate + nodeYcoordinate;
       TempNodeZcoordinate = TempNodeZcoordinate + nodeZcoordinate;
       // cout<<" node = "<<i<<endl;
       // cout<<" nodeXcoordinate = "<<nodeXcoordinate<<endl;
       // cout<<" nodeYcoordinate = "<<nodeYcoordinate<<endl;
       // cout<<" nodeZcoordinate = "<<nodeZcoordinate<<endl;
     }
     CFreal centralPointXcoordinate = TempNodeXcoordinate/nbNodesInFace;
     CFreal centralPointYcoordinate = TempNodeYcoordinate/nbNodesInFace;
     CFreal centralPointZcoordinate = TempNodeZcoordinate/nbNodesInFace;

     /*
     CFreal TempNodeXcoordinate = 10e5;
     CFreal TempNodeYcoordinate = 10e5;
     CFreal TempNodeZcoordinate = 10e5;
     cout<<" Face = "<<iFace<<endl;
     for(CFuint i=0; i<nbNodesInFace; i++){
       CFreal nodeXcoordinate = (*(nodes[i]))[XX];
       CFreal nodeYcoordinate = (*(nodes[i]))[YY];
       CFreal nodeZcoordinate = (*(nodes[i]))[ZZ];
       if(TempNodeXcoordinate > nodeXcoordinate){TempNodeXcoordinate = nodeXcoordinate;}
       if(TempNodeYcoordinate > nodeYcoordinate){TempNodeYcoordinate = nodeYcoordinate;}
       if(TempNodeZcoordinate > nodeZcoordinate){TempNodeZcoordinate = nodeZcoordinate;}
     }
     cout<<" node = "<<i<<endl;
     cout<<" TempNodeXcoordinate = "<<TempNodeXcoordinate<<endl;
     cout<<" TempNodeYcoordinate = "<<TempNodeYcoordinate<<endl;
     cout<<" TempNodeZcoordinate = "<<TempNodeZcoordinate<<endl;
     CFreal centralPointXcoordinate = TempNodeXcoordinate;
     CFreal centralPointYcoordinate = TempNodeYcoordinate;
     CFreal centralPointZcoordinate = TempNodeZcoordinate;*/
     //////////////////////////////////////////////////////////////////////////////////////////////////////////

     // build vectors
     pair<CFreal,CFuint> PairX;
     pair<CFreal,CFuint> PairY;
     pair<CFreal,CFuint> PairZ;
     pair<CFint,CFuint> PairFlagSurface;

     PairX = make_pair(centralPointXcoordinate,iFace);
     PairY = make_pair(centralPointYcoordinate,iFace);
     PairZ = make_pair(centralPointZcoordinate,iFace);
     PairFlagSurface = make_pair(flagSurface,iFace);

     if (flagSurface == -1){
       VecXwest.push_back(PairX);
       VecYwest.push_back(PairY);
       VecZwest.push_back(PairZ);
     }
     else {
       VecXeast.push_back(PairX);
       VecYeast.push_back(PairY);
       VecZeast.push_back(PairZ);
     }
     // release geometry
     _faceBuilder.releaseGE();

  }

  // sort
  _globalToLocalTRSFaceID.sortKeys();

  // build vector of west surface to be sent
  const CFuint dimqWest = _n_P*nWf;
  const CFuint dimq2West = nWf;

  vector<CFreal> SendVecXwestfirst;
  vector<CFreal> SendVecXwestsecond;
  SendVecXwestfirst.reserve(dimqWest);
  SendVecXwestsecond.reserve(dimqWest);
  for(CFuint p=0; p<_n_P; p++){
    for(CFuint pp=0; pp<dimq2West; pp++){
      SendVecXwestfirst.push_back(VecXwest[pp].first);
      SendVecXwestsecond.push_back(VecXwest[pp].second);
    }
  }

  vector<CFreal> SendVecYwestfirst;
  vector<CFreal> SendVecYwestsecond;
  SendVecYwestfirst.reserve(dimqWest);
  SendVecYwestsecond.reserve(dimqWest);
  for(CFuint p=0; p<_n_P; p++){
    for(CFuint pp=0; pp<dimq2West; pp++){
      SendVecYwestfirst.push_back(VecYwest[pp].first);
      SendVecYwestsecond.push_back(VecYwest[pp].second);
    }
  }

  vector<CFreal> SendVecZwestfirst;
  vector<CFreal> SendVecZwestsecond;
  SendVecZwestfirst.reserve(dimqWest);
  SendVecZwestsecond.reserve(dimqWest);
  for(CFuint p=0; p<_n_P; p++){
    for(CFuint pp=0; pp<dimq2West; pp++){
      SendVecZwestfirst.push_back(VecZwest[pp].first);
      SendVecZwestsecond.push_back(VecZwest[pp].second);
    }
  }

  // build vector of east surface to be sent
  const CFuint dimqEast = _n_P*nEf;
  const CFuint dimq2East = nEf;

  vector<CFreal> SendVecXeastfirst;
  vector<CFreal> SendVecXeastsecond;
  SendVecXeastfirst.reserve(dimqEast);
  SendVecXeastsecond.reserve(dimqEast);
  for(CFuint p=0; p<_n_P; p++){
    for(CFuint pp=0; pp<dimq2East; pp++){
      SendVecXeastfirst.push_back(VecXeast[pp].first);
      SendVecXeastsecond.push_back(VecXeast[pp].second);
    }
  }

  vector<CFreal> SendVecYeastfirst;
  vector<CFreal> SendVecYeastsecond;
  SendVecYeastfirst.reserve(dimqEast);
  SendVecYeastsecond.reserve(dimqEast);
  for(CFuint p=0; p<_n_P; p++){
    for(CFuint pp=0; pp<dimq2East; pp++){
      SendVecYeastfirst.push_back(VecYeast[pp].first);
      SendVecYeastsecond.push_back(VecYeast[pp].second);
    }
  }

  vector<CFreal> SendVecZeastfirst;
  vector<CFreal> SendVecZeastsecond;
  SendVecZeastfirst.reserve(dimqEast);
  SendVecZeastsecond.reserve(dimqEast);
  for(CFuint p=0; p<_n_P; p++){
    for(CFuint pp=0; pp<dimq2East; pp++){
      SendVecZeastfirst.push_back(VecZeast[pp].first);
      SendVecZeastsecond.push_back(VecZeast[pp].second);
    }
  }

  // build counts and displacements for west vectors
  vector<int> sendcountsWest(_n_P);
  vector<int> recvcountsWest(_n_P);
  vector<int> rdisplsWest(_n_P);
  vector<int> sdisplsWest(_n_P);
  rdisplsWest[0] = 0;
  sdisplsWest[0] = 0;
  CFuint LastDisplacementWest;
  for(CFuint p=0; p<_n_P; p++){
    sendcountsWest[p] = _countNWf[_my_nP];
    recvcountsWest[p] = _countNWf[p];
  }
  if(_n_P>1){
    for(CFuint p=1; p<_n_P; p++){
      rdisplsWest[p] = rdisplsWest[p-1] + recvcountsWest[p-1];
      sdisplsWest[p] = sdisplsWest[p-1] + sendcountsWest[p-1];
    }
  }
  vector<int> RecvDisWest(rdisplsWest);
  LastDisplacementWest = rdisplsWest[_n_P-1] + recvcountsWest[_n_P-1];

  // buil vector of processes for west face
  vector<CFreal> VecProcWest;
  VecProcWest.reserve(LastDisplacementWest);
  for(CFuint p=0; p<_n_P; p++){
    CFuint countPwest = _countNWf[p];
    if(countPwest != 0){
      for(CFuint pp=0; pp<countPwest; pp++){
        VecProcWest.push_back(p);
      }
    }
  }

  // build buffer for west vectors
  vector<CFreal> sbufVecXwestfirst(SendVecXwestfirst);
  vector<CFreal> sbufVecYwestfirst(SendVecYwestfirst);
  vector<CFreal> sbufVecZwestfirst(SendVecZwestfirst);
  vector<CFreal> sbufVecXwestsecond(SendVecXwestsecond);
  vector<CFreal> sbufVecYwestsecond(SendVecYwestsecond);
  vector<CFreal> sbufVecZwestsecond(SendVecZwestsecond);

  vector<CFreal> rbufVecXwestfirst(LastDisplacementWest,0.0);
  vector<CFreal> rbufVecYwestfirst(LastDisplacementWest,0.0);
  vector<CFreal> rbufVecZwestfirst(LastDisplacementWest,0.0);
  vector<CFreal> rbufVecXwestsecond(LastDisplacementWest,0.0);
  vector<CFreal> rbufVecYwestsecond(LastDisplacementWest,0.0);
  vector<CFreal> rbufVecZwestsecond(LastDisplacementWest,0.0);

  // build counts and displacements for east vectors
  vector<int> sendcountsEast(_n_P);
  vector<int> recvcountsEast(_n_P);
  vector<int> rdisplsEast(_n_P);
  vector<int> sdisplsEast(_n_P);
  rdisplsEast[0] = 0;
  sdisplsEast[0] = 0;
  CFuint LastDisplacementEast;
  for(CFuint p=0; p<_n_P; p++){
    sendcountsEast[p] = _countNEf[_my_nP];
    recvcountsEast[p] = _countNEf[p];
  }
  if(_n_P>1){
    for(CFuint p=1; p<_n_P; p++){
      rdisplsEast[p] = rdisplsEast[p-1] + recvcountsEast[p-1];
      sdisplsEast[p] = sdisplsEast[p-1] + sendcountsEast[p-1];
    }
  }
  vector<int> RecvDisEast(rdisplsEast);
  LastDisplacementEast = rdisplsEast[_n_P-1] + recvcountsEast[_n_P-1];

  // buil vector of processes for east face
  vector<CFreal> VecProcEast;
  VecProcEast.reserve(LastDisplacementEast);
  for(CFuint p=0; p<_n_P; p++){
    CFuint countPeast = _countNEf[p];
    if(countPeast != 0){
      for(CFuint pp=0; pp<countPeast; pp++){
        VecProcEast.push_back(p);
      }
    }
  }

  // build buffer for east vectors
  vector<CFreal> sbufVecXeastfirst(SendVecXeastfirst);
  vector<CFreal> sbufVecYeastfirst(SendVecYeastfirst);
  vector<CFreal> sbufVecZeastfirst(SendVecZeastfirst);
  vector<CFreal> sbufVecXeastsecond(SendVecXeastsecond);
  vector<CFreal> sbufVecYeastsecond(SendVecYeastsecond);
  vector<CFreal> sbufVecZeastsecond(SendVecZeastsecond);

  vector<CFreal> rbufVecXeastfirst(LastDisplacementEast,0.0);
  vector<CFreal> rbufVecYeastfirst(LastDisplacementEast,0.0);
  vector<CFreal> rbufVecZeastfirst(LastDisplacementEast,0.0);
  vector<CFreal> rbufVecXeastsecond(LastDisplacementEast,0.0);
  vector<CFreal> rbufVecYeastsecond(LastDisplacementEast,0.0);
  vector<CFreal> rbufVecZeastsecond(LastDisplacementEast,0.0);
 
  MPI_Datatype MPI_CFREAL = Common::MPIStructDef::getMPIType(&rbufVecXeastfirst[0]);
  
  // MPI_Alltoallv for west surface
  MPI_Alltoallv(&sbufVecXwestfirst[0], &sendcountsWest[0], &sdisplsWest[0], MPI_CFREAL, &rbufVecXwestfirst[0], &recvcountsWest[0], &rdisplsWest[0], MPI_CFREAL, _comm);
  std::vector<CFreal> RecvVecXwestfirst(rbufVecXwestfirst);
  MPI_Barrier(_comm);
  MPI_Alltoallv(&sbufVecXwestsecond[0], &sendcountsWest[0], &sdisplsWest[0], MPI_CFREAL, &rbufVecXwestsecond[0], &recvcountsWest[0], &rdisplsWest[0], MPI_CFREAL, _comm);
  std::vector<CFreal> RecvVecXwestsecond(rbufVecXwestsecond);
  MPI_Barrier(_comm);

  MPI_Alltoallv(&sbufVecYwestfirst[0], &sendcountsWest[0], &sdisplsWest[0], MPI_CFREAL, &rbufVecYwestfirst[0], &recvcountsWest[0], &rdisplsWest[0], MPI_CFREAL, _comm);
  std::vector<CFreal> RecvVecYwestfirst(rbufVecYwestfirst);
  MPI_Barrier(_comm);
  MPI_Alltoallv(&sbufVecYwestsecond[0], &sendcountsWest[0], &sdisplsWest[0], MPI_CFREAL, &rbufVecYwestsecond[0], &recvcountsWest[0], &rdisplsWest[0], MPI_CFREAL, _comm);
  std::vector<CFreal> RecvVecYwestsecond(rbufVecYwestsecond);
  MPI_Barrier(_comm);

  MPI_Alltoallv(&sbufVecZwestfirst[0], &sendcountsWest[0], &sdisplsWest[0], MPI_CFREAL, &rbufVecZwestfirst[0], &recvcountsWest[0], &rdisplsWest[0], MPI_CFREAL, _comm);
  std::vector<CFreal> RecvVecZwestfirst(rbufVecZwestfirst);
  MPI_Barrier(_comm);
  MPI_Alltoallv(&sbufVecZwestsecond[0], &sendcountsWest[0], &sdisplsWest[0], MPI_CFREAL, &rbufVecZwestsecond[0], &recvcountsWest[0], &rdisplsWest[0], MPI_CFREAL, _comm);
  std::vector<CFreal> RecvVecZwestsecond(rbufVecZwestsecond);
  MPI_Barrier(_comm);

  // MPI_Alltoallv for east surface
  MPI_Alltoallv(&sbufVecXeastfirst[0], &sendcountsEast[0], &sdisplsEast[0], MPI_CFREAL, &rbufVecXeastfirst[0], &recvcountsEast[0], &rdisplsEast[0], MPI_CFREAL, _comm);
  std::vector<CFreal> RecvVecXeastfirst(rbufVecXeastfirst);
  MPI_Barrier(_comm);
  MPI_Alltoallv(&sbufVecXeastsecond[0], &sendcountsEast[0], &sdisplsEast[0], MPI_CFREAL, &rbufVecXeastsecond[0], &recvcountsEast[0], &rdisplsEast[0], MPI_CFREAL, _comm);
  std::vector<CFreal> RecvVecXeastsecond(rbufVecXeastsecond);
  MPI_Barrier(_comm);

  MPI_Alltoallv(&sbufVecYeastfirst[0], &sendcountsEast[0], &sdisplsEast[0], MPI_CFREAL, &rbufVecYeastfirst[0], &recvcountsEast[0], &rdisplsEast[0], MPI_CFREAL, _comm);
  std::vector<CFreal> RecvVecYeastfirst(rbufVecYeastfirst);
  MPI_Barrier(_comm);
  MPI_Alltoallv(&sbufVecYeastsecond[0], &sendcountsEast[0], &sdisplsEast[0], MPI_CFREAL, &rbufVecYeastsecond[0], &recvcountsEast[0], &rdisplsEast[0], MPI_CFREAL, _comm);
  std::vector<CFreal> RecvVecYeastsecond(rbufVecYeastsecond);
  MPI_Barrier(_comm);

  MPI_Alltoallv(&sbufVecZeastfirst[0], &sendcountsEast[0], &sdisplsEast[0], MPI_CFREAL, &rbufVecZeastfirst[0], &recvcountsEast[0], &rdisplsEast[0], MPI_CFREAL, _comm);
  std::vector<CFreal> RecvVecZeastfirst(rbufVecZeastfirst);
  MPI_Barrier(_comm);
  MPI_Alltoallv(&sbufVecZeastsecond[0], &sendcountsEast[0], &sdisplsEast[0], MPI_CFREAL, &rbufVecZeastsecond[0], &recvcountsEast[0], &rdisplsEast[0], MPI_CFREAL, _comm);
  std::vector<CFreal> RecvVecZeastsecond(rbufVecZeastsecond);
  MPI_Barrier(_comm);

  /*
if(_my_nP == 0){
  cout<<" rbufVecZeastsecond[0] = "<<rbufVecZwestsecond[0]<<endl;
  cout<<" RecvVecZeastsecond[0] = "<<RecvVecZwestsecond[0]<<endl;
  cout<<" rbufVecZeastsecond[1] = "<<rbufVecZwestsecond[1]<<endl;
  cout<<" RecvVecZeastsecond[1] = "<<RecvVecZwestsecond[1]<<endl;
  cout<<" rbufVecZeastsecond[3] = "<<rbufVecZwestsecond[3]<<endl;
  cout<<" RecvVecZeastsecond[3] = "<<RecvVecZwestsecond[3]<<endl;
  cout<<" rbufVecZeastsecond[5] = "<<rbufVecZwestsecond[5]<<endl;
  cout<<" RecvVecZeastsecond[5] = "<<RecvVecZwestsecond[5]<<endl;
  cout<<" rbufVecZeastsecond[6] = "<<rbufVecZwestsecond[6]<<endl;
  cout<<" RecvVecZeastsecond[6] = "<<RecvVecZwestsecond[6]<<endl;
  cout<<" rbufVecZeastsecond[9] = "<<rbufVecZwestsecond[9]<<endl;
  cout<<" RecvVecZeastsecond[9] = "<<RecvVecZwestsecond[9]<<endl;

  cout<<" rbufVecZeastfirst[0] = "<<rbufVecZwestfirst[0]<<endl;
  cout<<" rbufVecZeastfirst[1] = "<<rbufVecZwestfirst[1]<<endl;
  cout<<" rbufVecZeastfirst[3] = "<<rbufVecZwestfirst[3]<<endl;
  cout<<" rbufVecZeastfirst[5] = "<<rbufVecZwestfirst[5]<<endl;
  cout<<" rbufVecZeastfirst[6] = "<<rbufVecZwestfirst[6]<<endl;
  cout<<" rbufVecZeastfirst[9] = "<<rbufVecZwestfirst[9]<<endl;

  cout<<" rbufVecYeastsecond[0] = "<<rbufVecYwestsecond[0]<<endl;
  cout<<" rbufVecYeastsecond[1] = "<<rbufVecYwestsecond[1]<<endl;
  cout<<" rbufVecYeastsecond[3] = "<<rbufVecYwestsecond[3]<<endl;
  cout<<" rbufVecYeastsecond[5] = "<<rbufVecYwestsecond[5]<<endl;
  cout<<" rbufVecYeastsecond[6] = "<<rbufVecYwestsecond[6]<<endl;
  cout<<" rbufVecYeastsecond[9] = "<<rbufVecYwestsecond[9]<<endl;

  cout<<" rbufVecYeastfirst[0] = "<<rbufVecYwestfirst[0]<<endl;
  cout<<" rbufVecYeastfirst[1] = "<<rbufVecYwestfirst[1]<<endl;
  cout<<" rbufVecYeastfirst[3] = "<<rbufVecYwestfirst[3]<<endl;
  cout<<" rbufVecYeastfirst[5] = "<<rbufVecYwestfirst[5]<<endl;
  cout<<" rbufVecYeastfirst[6] = "<<rbufVecYwestfirst[6]<<endl;
  cout<<" rbufVecYeastfirst[9] = "<<rbufVecYwestfirst[9]<<endl;

  cout<<" rbufVecXeastsecond[0] = "<<rbufVecXwestsecond[0]<<endl;
  cout<<" rbufVecXeastsecond[1] = "<<rbufVecXwestsecond[1]<<endl;
  cout<<" rbufVecXeastsecond[3] = "<<rbufVecXwestsecond[3]<<endl;
  cout<<" rbufVecXeastsecond[5] = "<<rbufVecXwestsecond[5]<<endl;
  cout<<" rbufVecXeastsecond[6] = "<<rbufVecXwestsecond[6]<<endl;
  cout<<" rbufVecXeastsecond[9] = "<<rbufVecXwestsecond[9]<<endl;

  cout<<" rbufVecXeastfirst[0] = "<<rbufVecXwestfirst[0]<<endl;
  cout<<" rbufVecXeastfirst[1] = "<<rbufVecXwestfirst[1]<<endl;
  cout<<" rbufVecXeastfirst[3] = "<<rbufVecXwestfirst[3]<<endl;
  cout<<" rbufVecXeastfirst[5] = "<<rbufVecXwestfirst[5]<<endl;
  cout<<" rbufVecXeastfirst[6] = "<<rbufVecXwestfirst[6]<<endl;
  cout<<" rbufVecXeastfirst[9] = "<<rbufVecXwestfirst[9]<<endl;
}

if(_my_nP == 1){
  cout<<"1 rbufVecZeastsecond[0] = "<<rbufVecZwestsecond[0]<<endl;
  cout<<"1 RecvVecZeastsecond[0] = "<<RecvVecZwestsecond[0]<<endl;
  cout<<"1 rbufVecZeastsecond[1] = "<<rbufVecZwestsecond[1]<<endl;
  cout<<"1 RecvVecZeastsecond[1] = "<<RecvVecZwestsecond[1]<<endl;
  cout<<"1 rbufVecZeastsecond[3] = "<<rbufVecZwestsecond[3]<<endl;
  cout<<"1 RecvVecZeastsecond[3] = "<<RecvVecZwestsecond[3]<<endl;
  cout<<"1 rbufVecZeastsecond[5] = "<<rbufVecZwestsecond[5]<<endl;
  cout<<"1 RecvVecZeastsecond[5] = "<<RecvVecZwestsecond[5]<<endl;
  cout<<"1 rbufVecZeastsecond[6] = "<<rbufVecZwestsecond[6]<<endl;
  cout<<"1 RecvVecZeastsecond[6] = "<<RecvVecZwestsecond[6]<<endl;
  cout<<"1 rbufVecZeastsecond[9] = "<<rbufVecZwestsecond[9]<<endl;
  cout<<"1 RecvVecZeastsecond[9] = "<<RecvVecZwestsecond[9]<<endl;

  cout<<"1 rbufVecZeastfirst[0] = "<<rbufVecZwestfirst[0]<<endl;
  cout<<"1 rbufVecZeastfirst[1] = "<<rbufVecZwestfirst[1]<<endl;
  cout<<"1 rbufVecZeastfirst[3] = "<<rbufVecZwestfirst[3]<<endl;
  cout<<"1 rbufVecZeastfirst[5] = "<<rbufVecZwestfirst[5]<<endl;
  cout<<"1 rbufVecZeastfirst[6] = "<<rbufVecZwestfirst[6]<<endl;
  cout<<"1 rbufVecZeastfirst[9] = "<<rbufVecZwestfirst[9]<<endl;

  cout<<"1 rbufVecYeastsecond[0] = "<<rbufVecYwestsecond[0]<<endl;
  cout<<"1 rbufVecYeastsecond[1] = "<<rbufVecYwestsecond[1]<<endl;
  cout<<"1 rbufVecYeastsecond[3] = "<<rbufVecYwestsecond[3]<<endl;
  cout<<"1 rbufVecYeastsecond[5] = "<<rbufVecYwestsecond[5]<<endl;
  cout<<"1 rbufVecYeastsecond[6] = "<<rbufVecYwestsecond[6]<<endl;
  cout<<"1 rbufVecYeastsecond[9] = "<<rbufVecYwestsecond[9]<<endl;

  cout<<"1 rbufVecYeastfirst[0] = "<<rbufVecYwestfirst[0]<<endl;
  cout<<"1 rbufVecYeastfirst[1] = "<<rbufVecYwestfirst[1]<<endl;
  cout<<"1 rbufVecYeastfirst[3] = "<<rbufVecYwestfirst[3]<<endl;
  cout<<"1 rbufVecYeastfirst[5] = "<<rbufVecYwestfirst[5]<<endl;
  cout<<"1 rbufVecYeastfirst[6] = "<<rbufVecYwestfirst[6]<<endl;
  cout<<"1 rbufVecYeastfirst[9] = "<<rbufVecYwestfirst[9]<<endl;

  cout<<"1 rbufVecXeastsecond[0] = "<<rbufVecXwestsecond[0]<<endl;
  cout<<"1 rbufVecXeastsecond[1] = "<<rbufVecXwestsecond[1]<<endl;
  cout<<"1 rbufVecXeastsecond[3] = "<<rbufVecXwestsecond[3]<<endl;
  cout<<"1 rbufVecXeastsecond[5] = "<<rbufVecXwestsecond[5]<<endl;
  cout<<"1 rbufVecXeastsecond[6] = "<<rbufVecXwestsecond[6]<<endl;
  cout<<"1 rbufVecXeastsecond[9] = "<<rbufVecXwestsecond[9]<<endl;

  cout<<"1 rbufVecXeastfirst[0] = "<<rbufVecXwestfirst[0]<<endl;
  cout<<"1 rbufVecXeastfirst[1] = "<<rbufVecXwestfirst[1]<<endl;
  cout<<"1 rbufVecXeastfirst[3] = "<<rbufVecXwestfirst[3]<<endl;
  cout<<"1 rbufVecXeastfirst[5] = "<<rbufVecXwestfirst[5]<<endl;
  cout<<"1 rbufVecXeastfirst[6] = "<<rbufVecXwestfirst[6]<<endl;
  cout<<"1 rbufVecXeastfirst[9] = "<<rbufVecXwestfirst[9]<<endl;
}
MPI_Barrier(_comm);
*/
  // build MAPS for west surface
  vector< pair<CFreal,CFuint> > MapXfacesWest;
  MapXfacesWest.reserve(LastDisplacementWest);
  vector< pair<CFreal,CFuint> > MapXprocessesWest;
  MapXprocessesWest.reserve(LastDisplacementWest);
  vector< pair<CFreal,CFuint> > MapYfacesWest;
  MapYfacesWest.reserve(LastDisplacementWest);
  vector< pair<CFreal,CFuint> > MapYprocessesWest;
  MapYprocessesWest.reserve(LastDisplacementWest);
  vector< pair<CFreal,CFuint> > MapZfacesWest;
  MapZfacesWest.reserve(LastDisplacementWest);
  vector< pair<CFreal,CFuint> > MapZprocessesWest;
  MapZprocessesWest.reserve(LastDisplacementWest);

  for(CFuint i=0; i<LastDisplacementWest; i++){

    CFreal xCoordinateWest = RecvVecXwestfirst[i];
    CFuint xFaceWest = RecvVecXwestsecond[i];
    CFuint xProcessWest = VecProcWest[i];
    MapXfacesWest.push_back(make_pair(xCoordinateWest, xFaceWest));
    MapXprocessesWest.push_back(make_pair(xProcessWest, xFaceWest));

    CFreal yCoordinateWest = RecvVecYwestfirst[i];
    CFuint yFaceWest = RecvVecYwestsecond[i];
    CFuint yProcessWest = VecProcWest[i];
    MapYfacesWest.push_back(make_pair(yCoordinateWest, yFaceWest));
    MapYprocessesWest.push_back(make_pair(yProcessWest, yFaceWest));

    CFreal zCoordinateWest = RecvVecZwestfirst[i];
    CFuint zFaceWest = RecvVecZwestsecond[i];
    CFuint zProcessWest = VecProcWest[i];
    MapZfacesWest.push_back(make_pair(zCoordinateWest, zFaceWest));
    MapZprocessesWest.push_back(make_pair(zProcessWest, zFaceWest));

  }

  // build MAPS for east surface
  vector< pair<CFreal,CFuint> > MapXfacesEast;
  MapXfacesEast.reserve(LastDisplacementEast);
  vector< pair<CFreal,CFuint> > MapXprocessesEast;
  MapXprocessesEast.reserve(LastDisplacementEast);
  vector< pair<CFreal,CFuint> > MapYfacesEast;
  MapYfacesEast.reserve(LastDisplacementEast);
  vector< pair<CFreal,CFuint> > MapYprocessesEast;
  MapYprocessesEast.reserve(LastDisplacementEast);
  vector< pair<CFreal,CFuint> > MapZfacesEast;
  MapZfacesEast.reserve(LastDisplacementEast);
  vector< pair<CFreal,CFuint> > MapZprocessesEast;
  MapZprocessesEast.reserve(LastDisplacementEast);

  for(CFuint i=0; i<LastDisplacementEast; i++){

    CFreal xCoordinateEast = RecvVecXeastfirst[i];
    CFuint xFaceEast = RecvVecXeastsecond[i];
    CFuint xProcessEast = VecProcEast[i];
    MapXfacesEast.push_back(make_pair(xCoordinateEast, xFaceEast));
    MapXprocessesEast.push_back(make_pair(xProcessEast, xFaceEast));//(xCoordinateEast, xProcessEast));

    CFreal yCoordinateEast = RecvVecYeastfirst[i];
    CFuint yFaceEast = RecvVecYeastsecond[i];
    CFuint yProcessEast = VecProcEast[i];
    MapYfacesEast.push_back(make_pair(yCoordinateEast, yFaceEast));
    MapYprocessesEast.push_back(make_pair(yProcessEast, yFaceEast));//(yCoordinateEast, yProcessEast));

    CFreal zCoordinateEast = RecvVecZeastfirst[i];
    CFuint zFaceEast = RecvVecZeastsecond[i];
    CFuint zProcessEast = VecProcEast[i];
    MapZfacesEast.push_back(make_pair(zCoordinateEast, zFaceEast));
    MapZprocessesEast.push_back(make_pair(zProcessEast, zFaceEast));//(zCoordinateEast, zProcessEast));

  }

  // sort west surface
  typedef pair<CFreal,CFuint> ThePair1;
  typedef pair<CFuint,CFuint> ThePair2;
  typedef pair<ThePair1,ThePair2> TheCFrealPair;
  typedef pair<TheCFrealPair,TheCFrealPair> TheTriplePair;
  typedef pair<TheCFrealPair,TheTriplePair> TheQuadPair;

  vector<TheCFrealPair> MapXwest;
  MapXwest.reserve(LastDisplacementWest);
  vector<TheCFrealPair> MapYwest;
  MapYwest.reserve(LastDisplacementWest);
  vector<TheCFrealPair> MapZwest;
  MapZwest.reserve(LastDisplacementWest);
  for(CFuint i=0; i<LastDisplacementWest; i++){
    TheCFrealPair xCFrealCorewest;
    xCFrealCorewest = make_pair(MapXfacesWest[i], MapXprocessesWest[i]);
    MapXwest.push_back(xCFrealCorewest);
    TheCFrealPair yCFrealCorewest;
    yCFrealCorewest = make_pair(MapYfacesWest[i], MapYprocessesWest[i]);
    MapYwest.push_back(yCFrealCorewest);
    TheCFrealPair zCFrealCorewest;
    zCFrealCorewest = make_pair(MapZfacesWest[i], MapZprocessesWest[i]);
    MapZwest.push_back(zCFrealCorewest);
  }

  vector<TheQuadPair> xQuadCorewest;
  xQuadCorewest.reserve(LastDisplacementWest);
  for(CFuint y=0; y<LastDisplacementWest; y++){
    TheTriplePair xTriplePairwest;
    TheQuadPair xQuadPairwest;
    xTriplePairwest = make_pair(MapYwest[y], MapZwest[y]);
    xQuadPairwest = make_pair(MapXwest[y],xTriplePairwest);
    xQuadCorewest.push_back(xQuadPairwest);
  }

  sort(xQuadCorewest.begin(), xQuadCorewest.end());

  vector<TheQuadPair> yQuadCorewest;
  yQuadCorewest.reserve(LastDisplacementWest);
  for(CFuint y=0; y<LastDisplacementWest; y++){
    TheTriplePair yTriplePairwest;
    TheQuadPair yQuadPairwest;
    yTriplePairwest = make_pair((xQuadCorewest[y].second).second, xQuadCorewest[y].first);
    yQuadPairwest = make_pair((xQuadCorewest[y].second).first, yTriplePairwest);
    yQuadCorewest.push_back(yQuadPairwest);
  }

  int indexArrayYwest[2];
  int indexWhileYwest = 0;
  for(CFuint i=1; i<LastDisplacementWest; i++){
    while(MathTools::MathChecks::isEqualWithError( (((yQuadCorewest[i].second).second).first).first, (((yQuadCorewest[i-1].second).second).first).first, _threshold)){
      indexArrayYwest[indexWhileYwest] = i - 1;
      indexArrayYwest[1] = i;
      indexWhileYwest = 1;
      i++;
    }
    if(indexWhileYwest != 0){
      sort(yQuadCorewest.begin() + indexArrayYwest[0], yQuadCorewest.begin() + indexArrayYwest[indexWhileYwest] + 1);
    }
    indexWhileYwest = 0;
  }

  vector<TheQuadPair> zQuadCorewest;
  zQuadCorewest.reserve(LastDisplacementWest);
  for(CFuint y=0; y<LastDisplacementWest; y++){
    TheTriplePair zTriplePairwest;
    TheQuadPair zQuadPairwest;
    zTriplePairwest = make_pair((yQuadCorewest[y].second).second, yQuadCorewest[y].first);
    zQuadPairwest = make_pair((yQuadCorewest[y].second).first, zTriplePairwest);
    zQuadCorewest.push_back(zQuadPairwest);
  }

  int indexArrayZwest[2];
  int indexWhileZwest = 0;
  for(CFuint i=1; i<LastDisplacementWest; i++){
    while((MathTools::MathChecks::isEqualWithError( (((zQuadCorewest[i].second).first).first).first, (((zQuadCorewest[i-1].second).first).first).first, _threshold)) && (MathTools::MathChecks::isEqualWithError( (((zQuadCorewest[i].second).second).first).first,  (((zQuadCorewest[i-1].second).second).first).first, _threshold))){
      indexArrayZwest[indexWhileZwest] = i - 1;
      indexArrayZwest[1] = i;
      indexWhileZwest = 1;
      i++;
    }
    if(indexWhileZwest != 0){
      sort(zQuadCorewest.begin() + indexArrayZwest[0], zQuadCorewest.begin() + indexArrayZwest[indexWhileZwest] + 1);
    }
    indexWhileZwest = 0;
  }

  // sort east surface
  vector<TheCFrealPair> MapXeast;
  MapXeast.reserve(LastDisplacementEast);
  vector<TheCFrealPair> MapYeast;
  MapYeast.reserve(LastDisplacementEast);
  vector<TheCFrealPair> MapZeast;
  MapZeast.reserve(LastDisplacementEast);
  for(CFuint i=0; i<LastDisplacementEast; i++){
    TheCFrealPair xCFrealCoreeast;
    xCFrealCoreeast = make_pair(MapXfacesEast[i], MapXprocessesEast[i]);
    MapXeast.push_back(xCFrealCoreeast);
    TheCFrealPair yCFrealCoreeast;
    yCFrealCoreeast = make_pair(MapYfacesEast[i], MapYprocessesEast[i]);
    MapYeast.push_back(yCFrealCoreeast);
    TheCFrealPair zCFrealCoreeast;
    zCFrealCoreeast = make_pair(MapZfacesEast[i], MapZprocessesEast[i]);
    MapZeast.push_back(zCFrealCoreeast);
  }

  vector<TheQuadPair> xQuadCoreeast;
  xQuadCoreeast.reserve(LastDisplacementEast);
  for(CFuint y=0; y<LastDisplacementEast; y++){
    TheTriplePair xTriplePaireast;
    TheQuadPair xQuadPaireast;
    xTriplePaireast = make_pair(MapYeast[y], MapZeast[y]);
    xQuadPaireast = make_pair(MapXeast[y],xTriplePaireast);
    xQuadCoreeast.push_back(xQuadPaireast);
  }

  sort(xQuadCoreeast.begin(), xQuadCoreeast.end());

  vector<TheQuadPair> yQuadCoreeast;
  yQuadCoreeast.reserve(LastDisplacementEast);
  for(CFuint y=0; y<LastDisplacementEast; y++){
    TheTriplePair yTriplePaireast;
    TheQuadPair yQuadPaireast;
    yTriplePaireast = make_pair((xQuadCoreeast[y].second).second, xQuadCoreeast[y].first);
    yQuadPaireast = make_pair((xQuadCoreeast[y].second).first, yTriplePaireast);
    yQuadCoreeast.push_back(yQuadPaireast);
  }

  int indexArrayYeast[2];
  int indexWhileYeast = 0;
  for(CFuint i=1; i<LastDisplacementEast; i++){
    while(MathTools::MathChecks::isEqualWithError( (((yQuadCoreeast[i].second).second).first).first, (((yQuadCoreeast[i-1].second).second).first).first, _threshold)){
      indexArrayYeast[indexWhileYeast] = i - 1;
      indexArrayYeast[1] = i;
      indexWhileYeast = 1;
      i++;
    }
    if(indexWhileYeast != 0){
      sort(yQuadCoreeast.begin() + indexArrayYeast[0], yQuadCoreeast.begin() + indexArrayYeast[indexWhileYeast] + 1);
    }
    indexWhileYeast = 0;
  }

  vector<TheQuadPair> zQuadCoreeast;
  zQuadCoreeast.reserve(LastDisplacementEast);
  for(CFuint y=0; y<LastDisplacementEast; y++){
    TheTriplePair zTriplePaireast;
    TheQuadPair zQuadPaireast;
    zTriplePaireast = make_pair((yQuadCoreeast[y].second).second, yQuadCoreeast[y].first);
    zQuadPaireast = make_pair((yQuadCoreeast[y].second).first, zTriplePaireast);
    zQuadCoreeast.push_back(zQuadPaireast);
  }

  int indexArrayZeast[2];
  int indexWhileZeast = 0;
  for(CFuint i=1; i<LastDisplacementEast; i++){
    while((MathTools::MathChecks::isEqualWithError( (((zQuadCoreeast[i].second).first).first).first, (((zQuadCoreeast[i-1].second).first).first).first, _threshold)) && (MathTools::MathChecks::isEqualWithError( (((zQuadCoreeast[i].second).second).first).first,  (((zQuadCoreeast[i-1].second).second).first).first, _threshold))){
      indexArrayZeast[indexWhileZeast] = i - 1;
      indexArrayZeast[1] = i;
      indexWhileZeast = 1;
      i++;
    }
    if(indexWhileZeast != 0){
      sort(zQuadCoreeast.begin() + indexArrayZeast[0], zQuadCoreeast.begin() + indexArrayZeast[indexWhileZeast] + 1);
    }
    indexWhileZeast = 0;
  }

  for(CFuint i=0; i<LastDisplacementEast; i++){
    // cout<<endl;
    // cout<<" east x["<<i<<"] = "<<(((zQuadCoreeast[i].second).first).first).first<<" west x["<<i<<"] = "<<(((zQuadCorewest[i].second).first).first).first<<endl;
    // cout<<" east y["<<i<<"] = "<<(((zQuadCoreeast[i].second).second).first).first<<" west y["<<i<<"] = "<<(((zQuadCorewest[i].second).second).first).first<<endl;
    // cout<<" east z["<<i<<"] = "<<((zQuadCoreeast[i].first).first).first<<" west z["<<i<<"] = "<<((zQuadCorewest[i].first).first).first<<endl;
    // cout<<endl;
  }

  // build connectivity vectors of west surface
  CFuint sizeWest = LastDisplacementWest;
  vector<CFuint> VectorProcessesWest;
  vector<CFuint> VectorFacesWest;
  vector<CFuint> VectorPeriodicFacesWest;
  VectorProcessesWest.reserve(LastDisplacementWest);
  VectorFacesWest.reserve(LastDisplacementWest);
  VectorPeriodicFacesWest.reserve(LastDisplacementWest);
  for(CFuint i=0; i<sizeWest; i++){
    if(MapXprocessesWest[i].first == _my_nP){
      VectorProcessesWest.push_back(((zQuadCorewest[i].first).second).first);
      VectorFacesWest.push_back(((zQuadCorewest[i].first).second).second);
      VectorPeriodicFacesWest.push_back(((zQuadCoreeast[i].first).second).second);
    }
  }

  // build connectivity vectors of east surface
  CFuint sizeEast = LastDisplacementEast;
  vector<CFuint> VectorProcessesEast;
  vector<CFuint> VectorFacesEast;
  vector<CFuint> VectorPeriodicFacesEast;
  VectorProcessesEast.reserve(LastDisplacementEast);
  VectorFacesEast.reserve(LastDisplacementEast);
  VectorPeriodicFacesEast.reserve(LastDisplacementEast);
  for(CFuint i=0; i<sizeEast; i++){
    if(MapXprocessesEast[i].first == _my_nP){
      VectorProcessesEast.push_back(((zQuadCoreeast[i].first).second).first);
      VectorFacesEast.push_back(((zQuadCoreeast[i].first).second).second);
      VectorPeriodicFacesEast.push_back(((zQuadCorewest[i].first).second).second);
    }
  }

  // build connectivity maps
  CFuint totalSize = sizeWest + sizeEast;
  _ConnectionFacePeriodic.reserve(totalSize);
  _ConnectionProcessPeriodic.reserve(totalSize);

  for(CFuint i=0; i<sizeWest; i++){
    CFuint iFaceWest = VectorFacesWest[i];
    CFuint iPeriodicFaceWest = VectorPeriodicFacesWest[i];
    CFuint pWest = VectorProcessesWest[i];
    _ConnectionFacePeriodic.insert(iFaceWest,iPeriodicFaceWest);
    _ConnectionProcessPeriodic.insert(iFaceWest,pWest);
  }
  for(CFuint i=0; i<sizeEast; i++){
    CFuint iFaceEast = VectorFacesEast[i];
    CFuint iPeriodicFaceEast = VectorPeriodicFacesEast[i];
    CFuint pEast = VectorProcessesEast[i];
    _ConnectionFacePeriodic.insert(iFaceEast,iPeriodicFaceEast);
    _ConnectionProcessPeriodic.insert(iFaceEast,pEast);
  }

  _ConnectionFacePeriodic.sortKeys();
  _ConnectionProcessPeriodic.sortKeys();

 // cout<<" _ConnectionFacePeriodic.print(): "<<endl;
  _ConnectionFacePeriodic.print();
  // cout<<" _ConnectionProcessPeriodic.print(): "<<endl;
  _ConnectionProcessPeriodic.print();

  // build counts and displacements for preProcess step
  _LastDisplacement = LastDisplacementEast + LastDisplacementWest;
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

void Periodic3DMPI::preProcess()
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

void Periodic3DMPI::setGhostState(GeometricEntity *const face)
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


void Periodic3DMPI::Barrier()
{
  MPI_Barrier(_comm);
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD
















