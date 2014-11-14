#include "Common/CFMultiMap.hh"
#include "Common/BadValueException.hh"
#include "Common/NoSuchValueException.hh"
#include "MathTools/RealMatrix.hh"
#include "Framework/MeshData.hh"
#include "Framework/MapGeoToTrsAndIdx.hh"
#include "Framework/MethodCommandProvider.hh"

#include "UFEM/UFEM.hh"
#include "UFEM/StdSetup.hh"
 
//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace UFEM {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<StdSetup, UFEMSolverData, UFEMPlugin> stdSetupProvider("StdSetup");

//////////////////////////////////////////////////////////////////////////////

StdSetup::StdSetup(const std::string& name) :
  UFEMSolverCom(name),
  _nodeIdToStateId(),
  socket_nstatesProxy("nstatesProxy"),
  socket_isUpdated("isUpdated"),
  socket_states("states"),
  socket_nodes("nodes"),
  socket_connBndFace2InnerCell("connBndFace2InnerCell"),
  socket_connBndFace2BndState("connBndFace2BndState"),
  socket_wallNearestSegment("wallNearestSegment"),
  socket_wallNearestDistance("wallNearestDistance"),
  socket_wallNearestDistanceState("wallNearestDistanceState"),
  socket_wallNearestVelocityGradient("wallNearestVelocityGradient"),
  socket_appliedStrongBC("appliedStrongBC")
{
  CFAUTOTRACE;
  wallNearestTrs=-1;
}

//////////////////////////////////////////////////////////////////////////////

StdSetup::~StdSetup()
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSource> >
StdSetup::providesSockets()
{
  CFAUTOTRACE;
  std::vector<Common::SafePtr<BaseDataSocketSource> > result;

  result.push_back(&socket_nstatesProxy);
  result.push_back(&socket_isUpdated);
  result.push_back(&socket_connBndFace2InnerCell);
  result.push_back(&socket_connBndFace2BndState);
  result.push_back(&socket_wallNearestSegment);
  result.push_back(&socket_wallNearestDistance);
  result.push_back(&socket_wallNearestDistanceState);
  result.push_back(&socket_wallNearestVelocityGradient);
  result.push_back(&socket_appliedStrongBC);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

void StdSetup::execute()
{
  CFAUTOTRACE;

  Common::SafePtr<TopologicalRegionSet> cells =
    MeshDataStack::getActive()->getTrs("InnerCells");

  DataHandle<State*,GLOBAL> states  = socket_states.getDataHandle();
  DataHandle<Node*,GLOBAL> nodes  = socket_nodes.getDataHandle();
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
  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  isUpdated.resize(nbStates*nbEqs);

  DataHandle<vector<bool> > appliedStrongBC = socket_appliedStrongBC.getDataHandle();
  appliedStrongBC.resize(nbStates);

  // this does not work in Windows so for the moment we
  // only compile it if not Win32
#if (!defined WIN32)
  setBndFace2InnerCellAndBndFace2BndStateConnectivity();
  
  if (getMethodData().getCalcWallDistance()) setWallNearestSegment();

  if (getMethodData().getReadWallDistFromFile())
  {
    setWallNearestSegment(); //this has to be called anyway to fill also the other sockets..
    readWallDistanceFromFile();
  }
    
#endif
}

//////////////////////////////////////////////////////////////////////////////

void StdSetup::setBndFace2InnerCellAndBndFace2BndStateConnectivity()
{
  CFAUTOTRACE;

  // nomenclature: everything starts with inner is the property of trs "InnerCells"
  //               then each other trs is sequentially processed, their vars starting with bnd

  vector< SafePtr<TopologicalRegionSet> > trs = MeshDataStack::getActive()->getTrsList();

  // this is what this function is made for
  DataHandle< std::vector<CFuint> > connBndFace2InnerCell = socket_connBndFace2InnerCell.getDataHandle();
  connBndFace2InnerCell.resize(trs.size());
  DataHandle< Common::ConnectivityTable<CFuint> > connBndFace2BndState = socket_connBndFace2BndState.getDataHandle();
  connBndFace2BndState.resize(trs.size());

  // finding innercells
  SafePtr<TopologicalRegionSet> innerTrs=0;
  for (CFuint i = 0; i < trs.size(); ++i) {
    if (trs[i]->getName() == "InnerCells") {
      innerTrs = trs[i];
      break;
    }
  }
  if (innerTrs==0) Common::NoSuchValueException(FromHere(),"Trs 'InnerCells' not found.");

  // local vars can be outside loop
  DataHandle< Node*,GLOBAL > innerNodes = socket_nodes.getDataHandle();
  const CFuint innerNbNodes = innerNodes.size();
  CFVec<CFint> bndNodeGlobal2Local(innerNbNodes);
  const CFuint innerNbGeos=innerTrs->getLocalNbGeoEnts();
  const Common::SafePtr< Common::ConnectivityTable<CFuint> > innerGeo2Node=innerTrs->getGeo2NodesConn();

  // looping on all other trss
  for (CFuint itrs = 0; itrs < trs.size(); ++itrs) {
    SafePtr<TopologicalRegionSet> bndTrs= trs[itrs];

    if (bndTrs->getName() != "InnerCells") {
/*
cout << "TRS: " << endl << bndTrs->getName() << endl << flush;
*/
      // getting nodes and preliminary settings
      Common::SafePtr< std::vector< CFuint > > bndNodes = bndTrs->getNodesInTrs();
      const CFuint bndNbNodes= bndTrs->getNbNodesInTrs();
      const CFuint bndNbGeos= bndTrs->getLocalNbGeoEnts();
      bndNodeGlobal2Local=-1;
      std::vector< std::vector<CFuint> > bndNod2Elm(bndNbNodes);
      connBndFace2InnerCell[itrs].resize(bndNbGeos);

      // cycle all the nodes in the TRS to set isBndNode flag
      CFuint bndCheckNbNode=0;
      for (std::vector< CFuint >::iterator itd = bndNodes->begin(); itd != bndNodes->end(); ++itd)
        bndNodeGlobal2Local[*itd]=bndCheckNbNode++;
      if (bndCheckNbNode!=bndNbNodes) Common::BadValueException(FromHere(),"bndCheckNbNode!=bndNbNodes");

      // build boundary node -> innercells elem connectivity
      CFuint bndNod2ElmMax=0;
      for (CFuint i=0; i<innerNbGeos; i++) {
        const CFuint innerGeoLocalID=innerTrs->getLocalGeoID(i);
        const CFuint innerNbGeoNodes=innerGeo2Node->nbCols(i);
        for (CFuint j=0; j<innerNbGeoNodes; j++)
          if (bndNodeGlobal2Local[(*innerGeo2Node)(i,j)]!=-1)
            bndNod2Elm[bndNodeGlobal2Local[(*innerGeo2Node)(i,j)]].push_back(innerGeoLocalID);
        if (bndNod2Elm[i].size()>bndNod2ElmMax) bndNod2ElmMax=bndNod2Elm[i].size();
      }

      // resetting trs face and trs state connectivity table
      connBndFace2BndState[itrs]=*bndTrs->getGeo2StatesConn();
      for (CFuint igeo=0; igeo<bndNbGeos; igeo++){
        const CFuint igeoNbStates=bndTrs->getNbStatesInGeo(igeo);
        for (CFuint istate=0; istate<igeoNbStates; istate++){
          connBndFace2BndState[itrs](igeo,istate)=bndNodeGlobal2Local[connBndFace2BndState[itrs](igeo,istate)];
        }
      }

//cout << "connBndFace2BndState[itrs]" << endl << connBndFace2BndState[itrs] << flush;

      // looping on each bnd element, from bndNod2Elm counting which element is connected and how many times,
      // the winner with the most shared nodes is the wall adjacent cell of that bnd face
      // also filling current trs's fac to state connectivity with trs local numbering
      std::vector<CFuint> bndGeoIdentifier;
      CFVec<CFuint> bndGeoIdentifierCounter;
      for (CFuint igeo=0; igeo<bndNbGeos; igeo++){

        const CFuint igeoNbNodes=bndTrs->getNbNodesInGeo(igeo);

        bndGeoIdentifier.resize(0);
        for (CFuint inod=0; inod<igeoNbNodes; inod++){
          const CFuint inodLocalID= bndTrs->getNodeID(igeo,inod);
          const CFuint inodNbCon=bndNod2Elm[bndNodeGlobal2Local[inodLocalID]].size();
          for (CFuint icon=0; icon<inodNbCon; icon++){
            bndGeoIdentifier.push_back(bndNod2Elm[bndNodeGlobal2Local[inodLocalID]][icon]);
          }
        }
        std::sort(bndGeoIdentifier.begin(),bndGeoIdentifier.end());
        bndGeoIdentifierCounter.resize(bndGeoIdentifier.size());
        bndGeoIdentifierCounter=1;
        const CFuint bndNbGeoIdentifier=bndGeoIdentifier.size();
        for (CFuint iid=1; iid<bndNbGeoIdentifier; iid++)
          if (bndGeoIdentifier[iid]==bndGeoIdentifier[iid-1])
            bndGeoIdentifierCounter[iid]+=bndGeoIdentifierCounter[iid-1];
        CFuint nbGeoMatch=0;
        CFuint lastGeoMatch=0;
        for (CFuint iid=1; iid<bndNbGeoIdentifier; iid++){
          if (bndGeoIdentifierCounter[iid]> igeoNbNodes)
            Common::BadValueException(FromHere(),"More nodes matching than the boundary element actually has.");
          if (bndGeoIdentifierCounter[iid]==igeoNbNodes){
            nbGeoMatch++;
            lastGeoMatch=iid;
          }
        }
        if (nbGeoMatch<1) Common::BadValueException(FromHere(),"No wall adjacent cell has been found.");
        if (nbGeoMatch>1) Common::BadValueException(FromHere(),"More than one adjacent cell has been found.");

        connBndFace2InnerCell[itrs][igeo]=bndGeoIdentifier[lastGeoMatch]; // hurray 1

//cout << igeo << " :   " << flush;
//for (CFuint i=0; i<bndGeoIdentifier.size(); i++)
//  cout << " " << bndGeoIdentifier[i] << flush;
//cout <<  " matching= " << connBndFace2InnerCell[itrs][igeo] << endl << flush;

      }
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

//#include <fstream>

void StdSetup::setWallNearestSegment()
{
  CFAUTOTRACE;
  CFLog(INFO, "Computing wall distances using simple state-node distance...\n");
/*
  // sockets...
  DataHandle< CFuint > wallNearestSegment = socket_wallNearestSegment.getDataHandle();
  DataHandle< CFreal > wallNearestDistance = socket_wallNearestDistance.getDataHandle();
  DataHandle< CFreal > wallNearestDistanceState = socket_wallNearestDistanceState.getDataHandle();
  DataHandle< CFreal > wallNearestVelocityGradient = socket_wallNearestVelocityGradient.getDataHandle();

  // finding trs named as "wall" and the trs of innercells
  std::vector<Common::SafePtr<TopologicalRegionSet> > trs = MeshDataStack::getActive()->getTrsList();
  SafePtr< TopologicalRegionSet > bndTrs=0;
  SafePtr< TopologicalRegionSet > innerTrs=0;
  for (CFuint i = 0; i < trs.size(); ++i) {
    if (trs[i]->hasTag("inner")) {
      innerTrs = trs[i];
    }
    if (trs[i]->getName() == "wall" ) {
      bndTrs=trs[i];
      wallNearestTrs=i;
    }
  }
  const CFuint nbInnerGeoEnt=innerTrs->getLocalNbGeoEnts();
  const CFuint nbInnerStates=innerTrs->getNbStatesInTrs();
  if (bndTrs==0) {
    wallNearestSegment.resize(nbInnerGeoEnt);
    wallNearestDistance.resize(nbInnerGeoEnt);
    wallNearestDistanceState.resize(nbInnerStates);
    wallNearestVelocityGradient.resize(0);
    for (CFuint i = 0; i != nbInnerGeoEnt; ++i) {
      wallNearestSegment[i] = 0;
      wallNearestDistance[i] = 0.;
      wallNearestDistanceState[i] = 0.;
    }
  } else {
    // calculating boundary segments centroids
    CFuint nbBndGeoEnt=bndTrs->getLocalNbGeoEnts();
    std::vector<CFreal> bndCentroidsX(nbBndGeoEnt,0.);
    std::vector<CFreal> bndCentroidsY(nbBndGeoEnt,0.);
    std::vector<CFreal> bndCentroidsZ(nbBndGeoEnt,0.);
    std::vector<CFreal> bndNormalsX(nbBndGeoEnt,0.);
    std::vector<CFreal> bndNormalsY(nbBndGeoEnt,0.);
    std::vector<CFreal> bndNormalsZ(nbBndGeoEnt,0.);
    SafePtr< GeometricEntityPool< StdTrsGeoBuilder > > geoBuilder = getMethodData().getStdTrsGeoBuilder();
    StdTrsGeoBuilder::GeoData& geoData = geoBuilder->getDataGE();
    geoData.trs = bndTrs;
    for (CFuint iEnt = 0; iEnt < nbBndGeoEnt; ++iEnt){
      geoData.idx = iEnt;
      GeometricEntity& cell = *geoBuilder->buildGE();
      RealVector centroid=cell.computeCentroid();
      bndCentroidsX[iEnt]=centroid[0];
      bndCentroidsY[iEnt]=centroid[1];
      if (centroid.size()==3) bndCentroidsZ[iEnt]=centroid[2];
      RealVector normal=cell.computeAvgCellNormal();
      normal.normalize();
      bndNormalsX[iEnt]=normal[0];
      bndNormalsY[iEnt]=normal[1];
      if (normal.size()==3) bndNormalsZ[iEnt]=normal[2];
      geoBuilder->releaseGE();
    }

    // brute-force loop for identifying bnd cell for all cells
    wallNearestSegment.resize(nbInnerGeoEnt);
    wallNearestDistance.resize(nbInnerGeoEnt);
    wallNearestDistanceState.resize(nbInnerStates);
    wallNearestVelocityGradient.resize(nbBndGeoEnt);
    geoData.trs = innerTrs;
    for (CFuint is=0; is<nbInnerStates; is++) wallNearestDistanceState[is]=1.e30;
    for (CFuint iEnt = 0; iEnt < nbInnerGeoEnt; ++iEnt){
      // first identifying nearest cell
      CFreal mindist=1.e30,currdist;
      CFuint minbnd=-1;
      wallNearestSegment[iEnt]=-1;
      wallNearestDistance[iEnt]=1.e30;
      geoData.idx = iEnt;
      GeometricEntity& cell = *geoBuilder->buildGE();
      RealVector centroid=cell.computeCentroid();
      for (CFuint iBnd = 0; iBnd < nbBndGeoEnt; ++iBnd){
        CFreal X=centroid[0]-bndCentroidsX[iBnd];
        CFreal Y=centroid[1]-bndCentroidsY[iBnd];
        CFreal Z=0.;
        if (centroid.size()==3) Z=centroid[2]-bndCentroidsZ[iBnd];
        currdist=X*X;
        if (currdist<mindist) {
          currdist+=Y*Y;
          if (currdist<mindist) {
            currdist+=Z*Z;
            if (currdist<mindist) {
              mindist=currdist;
              minbnd=iBnd;
            }
          }
        }
      }

      CFreal X=centroid[0]-bndCentroidsX[minbnd];
      CFreal Y=centroid[1]-bndCentroidsY[minbnd];
      CFreal Z=0.;
      if (centroid.size()==3) Z=centroid[2]-bndCentroidsZ[minbnd];
      X*=bndNormalsX[minbnd];
      Y*=bndNormalsY[minbnd];
      Z*=bndNormalsZ[minbnd];
      wallNearestDistance[iEnt]=fabs(X+Y+Z);
      wallNearestSegment[iEnt]=minbnd;
      std::vector< State* >* states = cell.getStates();
      for (CFuint istate=0; istate<states->size(); istate++){
        CFuint lid=(*states)[istate]->getLocalID();
        Node& inode=(*states)[istate]->getCoordinates();
        CFreal X=inode[0]-bndCentroidsX[minbnd];
        CFreal Y=inode[1]-bndCentroidsY[minbnd];
        CFreal Z=0.;
        if (centroid.size()==3) Z=inode[2]-bndCentroidsZ[minbnd];
        X*=bndNormalsX[minbnd];
        Y*=bndNormalsY[minbnd];
        Z*=bndNormalsZ[minbnd];
        wallNearestDistanceState[lid]=min(fabs(X+Y+Z),wallNearestDistanceState[lid]);

      }

      geoBuilder->releaseGE();
    }

    // making sure, on bndTrs that state-wise wall distances are exactly zero
    Common::SafePtr< std::vector<CFuint> > bndTrsStates=bndTrs->getStatesInTrs();
    for (CFuint i=0; i<bndTrs->getNbStatesInTrs(); i++){
      const CFuint lid=(*bndTrsStates)[i];
      wallNearestDistanceState[lid]=0.;
    }

  }


//{ CFuint count=wallNearestDistanceState.size();
//cout << "\nwallNearestDistanceState(" << count << "): \n";
//for (CFuint i=0; i<count; i++) 
//  cout << (wallNearestDistanceState)[i] << " ";
//cout << "\n" << flush;}
*/

  // this breaks standard k-epsilon


  // sockets...
  DataHandle< CFuint > wallNearestSegment = socket_wallNearestSegment.getDataHandle();
  DataHandle< CFreal > wallNearestDistance = socket_wallNearestDistance.getDataHandle();
  DataHandle< CFreal > wallNearestDistanceState = socket_wallNearestDistanceState.getDataHandle();
  DataHandle< CFreal > wallNearestVelocityGradient = socket_wallNearestVelocityGradient.getDataHandle();

  // finding trs named as "wall" and the trs of innercells
  std::vector<Common::SafePtr<TopologicalRegionSet> > trs = MeshDataStack::getActive()->getTrsList();
  SafePtr< TopologicalRegionSet > bndTrs=0;
  SafePtr< TopologicalRegionSet > innerTrs=0;
  for (CFuint i = 0; i < trs.size(); ++i) {
    if (trs[i]->hasTag("inner")) {
      innerTrs = trs[i];
    }
    if (trs[i]->getName() == "wall" ) {
      bndTrs=trs[i];
      wallNearestTrs=i;
    }
  }
  const CFuint nbInnerGeoEnt=innerTrs->getLocalNbGeoEnts();
  const CFuint nbInnerStates=innerTrs->getNbStatesInTrs();
  wallNearestSegment.resize(nbInnerGeoEnt);
  wallNearestDistance.resize(nbInnerGeoEnt);
  wallNearestDistanceState.resize(nbInnerStates);
  wallNearestVelocityGradient.resize(0);
  for (CFuint i = 0; i != nbInnerGeoEnt; ++i) {
    wallNearestSegment[i] = 0;
    wallNearestDistance[i] = 0.;
  }
  for (CFuint i = 0; i != nbInnerStates; ++i) {
    wallNearestDistanceState[i] = 0.;
  }
  if (bndTrs!=0) {
    // calculating boundary segments centroids
    CFuint nbBndGeoEnt=bndTrs->getLocalNbGeoEnts();
    CFuint nbBndGeoEntAll=0;
    MPI_Allreduce(&nbBndGeoEnt,&nbBndGeoEntAll,1,MPIStructDef::getMPIType(&nbBndGeoEnt), MPI_SUM,MPI_COMM_WORLD);
    std::vector<CFreal> bndCentroidsX(nbBndGeoEntAll,0.);
    std::vector<CFreal> bndCentroidsY(nbBndGeoEntAll,0.);
    std::vector<CFreal> bndCentroidsZ(nbBndGeoEntAll,0.);
    std::vector<CFreal> bndNormalsX(nbBndGeoEntAll,0.);
    std::vector<CFreal> bndNormalsY(nbBndGeoEntAll,0.);
    std::vector<CFreal> bndNormalsZ(nbBndGeoEntAll,0.);
    SafePtr< GeometricEntityPool< StdTrsGeoBuilder > > geoBuilder = getMethodData().getStdTrsGeoBuilder();
    StdTrsGeoBuilder::GeoData& geoData = geoBuilder->getDataGE();
    geoData.trs = bndTrs;
    for (CFuint iEnt = 0; iEnt < nbBndGeoEnt; ++iEnt){
      geoData.idx = iEnt;
      GeometricEntity& cell = *geoBuilder->buildGE();
      RealVector centroid=cell.computeCentroid();
      bndCentroidsX[iEnt]=centroid[0];
      bndCentroidsY[iEnt]=centroid[1];
      if (centroid.size()==3) bndCentroidsZ[iEnt]=centroid[2];
      RealVector normal=cell.computeAvgCellNormal();
      normal.normalize();
      bndNormalsX[iEnt]=normal[0];
      bndNormalsY[iEnt]=normal[1];
      if (normal.size()==3) bndNormalsZ[iEnt]=normal[2];
      geoBuilder->releaseGE();
    }

    int *nbBndGeos=new int[PE::GetPE().GetProcessorCount()];
    int *dispBndGeos=new int[PE::GetPE().GetProcessorCount()];
    MPI_Allgather(&nbBndGeoEnt,1,MPI_INT,nbBndGeos,1,MPI_INT,MPI_COMM_WORLD);
    dispBndGeos[0]=0;
    for (int i=1; i<PE::GetPE().GetProcessorCount(); i++) dispBndGeos[i]=dispBndGeos[i-1]+nbBndGeos[i-1];
/*
cout << "proccount=" << PE::GetPE().GetProcessorCount() << "\n" << flush;

for (int i=1; i<PE::GetPE().GetProcessorCount(); i++) cout << nbBndGeos[i] << " " <<  dispBndGeos[i] << "\n" << flush;

{ CFuint count=nbBndGeoEnt;
cout << "\nbndCentroidsX(" << count << "): \n";
for (CFuint i=0; i<count; i++) 
  cout << (bndCentroidsX)[i] << " ";
cout << "\n" << flush;}
*/
    MPI_Allgatherv(&bndCentroidsX[0],nbBndGeoEnt,MPI_DOUBLE,&bndCentroidsX[0],nbBndGeos,dispBndGeos,MPI_DOUBLE,MPI_COMM_WORLD);
    MPI_Allgatherv(&bndCentroidsY[0],nbBndGeoEnt,MPI_DOUBLE,&bndCentroidsY[0],nbBndGeos,dispBndGeos,MPI_DOUBLE,MPI_COMM_WORLD);
    MPI_Allgatherv(&bndCentroidsZ[0],nbBndGeoEnt,MPI_DOUBLE,&bndCentroidsZ[0],nbBndGeos,dispBndGeos,MPI_DOUBLE,MPI_COMM_WORLD);
    MPI_Allgatherv(&bndNormalsX[0],  nbBndGeoEnt,MPI_DOUBLE,&bndNormalsX[0],  nbBndGeos,dispBndGeos,MPI_DOUBLE,MPI_COMM_WORLD);
    MPI_Allgatherv(&bndNormalsY[0],  nbBndGeoEnt,MPI_DOUBLE,&bndNormalsY[0],  nbBndGeos,dispBndGeos,MPI_DOUBLE,MPI_COMM_WORLD);
    MPI_Allgatherv(&bndNormalsZ[0],  nbBndGeoEnt,MPI_DOUBLE,&bndNormalsZ[0],  nbBndGeos,dispBndGeos,MPI_DOUBLE,MPI_COMM_WORLD);
/*
{ CFuint count=bndCentroidsX.size();
cout << "\nbndCentroidsX(" << count << "): \n";
for (CFuint i=0; i<count; i++) 
  cout << (bndCentroidsX)[i] << " ";
cout << "\n" << flush;}
*/
    delete []nbBndGeos;
    delete []dispBndGeos;



    // brute-force loop for identifying bnd cell for all cells
    wallNearestSegment.resize(nbInnerGeoEnt);
    wallNearestDistance.resize(nbInnerGeoEnt);
    wallNearestDistanceState.resize(nbInnerStates);
    wallNearestVelocityGradient.resize(nbBndGeoEntAll);
    geoData.trs = innerTrs;
    for (CFuint is=0; is<nbInnerStates; is++) wallNearestDistanceState[is]=1.e30;
    for (CFuint iEnt = 0; iEnt < nbInnerGeoEnt; ++iEnt){
      // first identifying nearest cell
      CFreal mindist=1.e30,currdist;
      CFuint minbnd=-1;
      wallNearestSegment[iEnt]=-1;
      wallNearestDistance[iEnt]=1.e30;
      geoData.idx = iEnt;
      GeometricEntity& cell = *geoBuilder->buildGE();
      RealVector centroid=cell.computeCentroid();
      for (CFuint iBnd = 0; iBnd < nbBndGeoEntAll; ++iBnd){
        CFreal X=centroid[0]-bndCentroidsX[iBnd];
        CFreal Y=centroid[1]-bndCentroidsY[iBnd];
        CFreal Z=0.;
        if (centroid.size()==3) Z=centroid[2]-bndCentroidsZ[iBnd];
        currdist=X*X;
        if (currdist<mindist) {
          currdist+=Y*Y;
          if (currdist<mindist) {
            currdist+=Z*Z;
            if (currdist<mindist) {
              mindist=currdist;
              minbnd=iBnd;
            }
          }
        }
      }

      CFreal X=centroid[0]-bndCentroidsX[minbnd];
      CFreal Y=centroid[1]-bndCentroidsY[minbnd];
      CFreal Z=0.;
      if (centroid.size()==3) Z=centroid[2]-bndCentroidsZ[minbnd];
      X*=bndNormalsX[minbnd];
      Y*=bndNormalsY[minbnd];
      Z*=bndNormalsZ[minbnd];
      wallNearestDistance[iEnt]=fabs(X+Y+Z);
      wallNearestSegment[iEnt]=minbnd;
      std::vector< State* >* states = cell.getStates();
      for (CFuint istate=0; istate<states->size(); istate++){
        CFuint lid=(*states)[istate]->getLocalID();
        Node& inode=(*states)[istate]->getCoordinates();
        CFreal X=inode[0]-bndCentroidsX[minbnd];
        CFreal Y=inode[1]-bndCentroidsY[minbnd];
        CFreal Z=0.;
        if (centroid.size()==3) Z=inode[2]-bndCentroidsZ[minbnd];
        X*=bndNormalsX[minbnd];
        Y*=bndNormalsY[minbnd];
        Z*=bndNormalsZ[minbnd];
        wallNearestDistanceState[lid]=min(fabs(X+Y+Z),wallNearestDistanceState[lid]);

      }

      geoBuilder->releaseGE();
    }

    // making sure, on bndTrs that state-wise wall distances are exactly zero
    Common::SafePtr< std::vector<CFuint> > bndTrsStates=bndTrs->getStatesInTrs();
    for (CFuint i=0; i<bndTrs->getNbStatesInTrs(); i++){
      const CFuint lid=(*bndTrsStates)[i];
      wallNearestDistanceState[lid]=0.;
    }

  }

/*
{ CFuint count=wallNearestDistanceState.size();
cout << "\nwallNearestDistanceState(" << count << "): \n";
for (CFuint i=0; i<count; i++) 
  cout << (wallNearestDistanceState)[i] << " ";
cout << "\n" << flush;}
*/

  CFLog(INFO, "Wall distances computation finished...\n");

}

//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////
void StdSetup::readWallDistanceFromFile()
{

  // finding trs named as "wall" and the trs of innercells
  std::vector<Common::SafePtr<TopologicalRegionSet> > trs = MeshDataStack::getActive()->getTrsList();
  SafePtr< TopologicalRegionSet > innerTrs=0;
  for (CFuint i = 0; i < trs.size(); ++i) {
    if (trs[i]->hasTag("inner")) {
      innerTrs = trs[i];
    }
  }
//   const CFuint nbInnerGeoEnt=innerTrs->getLocalNbGeoEnts();
  const CFuint nbInnerStates=innerTrs->getNbStatesInTrs();
  
  
  
  
  
   int irank, nrank;
   MPI_Comm_size(MPI_COMM_WORLD,&nrank);
   MPI_Comm_rank(MPI_COMM_WORLD,&irank);
   if(!irank) cout<< "Reading Wall Distance.. " << nbInnerStates <<endl;

   
   int GlobNbNodes = 0;
   MPI_Allreduce((void*)&nbInnerStates,(void*)&GlobNbNodes,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
   vector< CFreal > WallDistVec(GlobNbNodes);
   
   //soket..
   DataHandle< CFreal > wallNearestDistanceState = socket_wallNearestDistanceState.getDataHandle();
   wallNearestDistanceState.resize(nbInnerStates);
      
   Common::SafePtr<std::vector<CFuint> > NodesID;
   NodesID = innerTrs->getNodesInTrs();


 //---------- read the wallDist file.. ------------
   ifstream ParamFile;
   ParamFile.open( getMethodData().getWallDistFileName().c_str() );
   ParamFile.clear();
   ParamFile.seekg ( 0, ios::beg ); 
   CFuint count =0;
   CFreal strR;
   
   cout << "File open!" <<endl;
 
   while ( !ParamFile.eof() )
   {     
    ParamFile >> strR;
    if(count<GlobNbNodes) WallDistVec[count] = strR;
    count++;
   }
 //---------------------------------------------
 
 //---------- search for ID's ------------------
  Common::SafePtr<GeometricEntityPool<StdTrsGeoBuilder> > GeoBuilder = getMethodData().getStdTrsGeoBuilder();
  StdTrsGeoBuilder::GeoData& GeoData = GeoBuilder->getDataGE();
  
   for (CFuint igeo=0; igeo<innerTrs->getLocalNbGeoEnts(); igeo++){ //loop over the inner cells
    GeoData.trs = innerTrs;
    GeoData.idx = igeo;
    GeometricEntity *const Geo = GeoBuilder->buildGE();
    const vector< State* >* states = Geo->getStates();
    const CFuint nbInnerGeoStates= Geo->nbStates();

    for (CFuint istate=0; istate<nbInnerGeoStates; istate++)
    { 
     wallNearestDistanceState[((*states)[istate])->getLocalID()] = WallDistVec[ ((*states)[istate])->getGlobalID() ];
    }

    GeoBuilder->releaseGE();
  }
 //---------------------------------------------
 
   if(!irank) cout<< "Wall Distance Read!!.." <<endl;          

}
//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
StdSetup::needsSockets()
{
  CFAUTOTRACE;
  std::vector<Common::SafePtr<BaseDataSocketSink> > result;

  result.push_back(&socket_states);
  result.push_back(&socket_nodes);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

 } // namespace UFEM

} // namespace COOLFluiD
