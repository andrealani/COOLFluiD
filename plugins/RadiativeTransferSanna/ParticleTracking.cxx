#include <iostream>
#include <cmath>
#include <algorithm>

#include "MathTools/MathChecks.hh"
#include "RadiativeTransferSanna/ParticleTracking.hh"
#include "FiniteVolume/CellCenterFVM.hh"
#include "Framework/MeshData.hh"
#include "Framework/PhysicalModel.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::MathTools;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace RadiativeTransfer {

//////////////////////////////////////////////////////////////////////////////
    
ParticleTracking::ParticleTracking() :
  socket_states("Null"),
  socket_gstates("Null"),
  socket_nodes("Null"),
  socket_nstates("Null"),
  socket_normals("Null"),
  socket_isOutward("Null"),
  _cellBuilder(),
  _faceBuilder(),
  _direction(),
  _endPoint(),
  _startPoint(),
  _internalPoint(),
  _intersectionPoint(),
  _barycentreActualCell(),
  _barycentreActualFace(),
  _Lx(),
  _Ly(),
  _Lz(),
  _Rx(),
  _Ry(),
  _Rz(),
  _nn(),
  _barycentreFace(), 
  _nodeXcoordinate(),
  _nodeYcoordinate(),
  _nodeZcoordinate(),
  _foundIt(),
  _dim(),
  _set(),
  _startCellID(),
  _cellIdx(),
  _FaceNb(),
  _ExitFaceID(),
  _wallFaceID(),
  _CellIDmap()
{
}

//////////////////////////////////////////////////////////////////////////////

ParticleTracking::~ParticleTracking()
{
}

//////////////////////////////////////////////////////////////////////////////

void ParticleTracking::setDataSockets(DataSocketSink<Framework::State*, Framework::GLOBAL> statesSocket,
                                      DataSocketSink<State*> gstatesSocket,
                                      DataSocketSink<Framework::Node*, Framework::GLOBAL> nodesSocket,
                                      DataSocketSink<RealVector> nstatesSocket,
                                      Framework::DataSocketSink< CFreal> normalsSocket,
                                      Framework::DataSocketSink<CFint> isOutwardsSocket
                                      )
{
  socket_states = statesSocket;
  socket_gstates = gstatesSocket;
  socket_nodes = nodesSocket;
  socket_nstates = nstatesSocket;
  socket_normals = normalsSocket;
  socket_isOutward = isOutwardsSocket;
  
  // cell builder initialization
  _cellBuilder.setup();
  _cellBuilder.getGeoBuilder()->setDataSockets(socket_states, socket_gstates, socket_nodes);
  CellTrsGeoBuilder::GeoData& cellData = _cellBuilder.getDataGE();
  cellData.trs = MeshDataStack::getActive()->getTrs("InnerCells");
  
  // face builder initialization
  _faceBuilder.setup();
  _faceBuilder.getGeoBuilder()->setDataSockets(socket_states, socket_gstates, socket_nodes);
  
  _dim = PhysicalModelStack::getActive()->getDim();
  _direction.resize(_dim);
  _endPoint.resize(_dim);
  _startPoint.resize(_dim);
  _internalPoint.resize(_dim);
  _intersectionPoint.resize(_dim);
  _barycentreActualCell.resize(_dim);
  _barycentreActualFace.resize(_dim);
  _Lx.resize(_dim);
  _Ly.resize(_dim);
  _Lz.resize(_dim);
  _Rx.resize(_dim);
  _Ry.resize(_dim);
  _Rz.resize(_dim);
  _nn.resize(_dim);
  _barycentreFace.resize(_dim);
  
  // AL: atomic number here => this corresponds to the number of faces in a hexahedron
  const CFuint maxNbNodesInCell = 8;
  _nodeXcoordinate.reserve(maxNbNodesInCell);
  _nodeYcoordinate.reserve(maxNbNodesInCell);
  _nodeZcoordinate.reserve(maxNbNodesInCell);				
  
  // build the cell IDs mapping
  buildCellIDmap();
}

//////////////////////////////////////////////////////////////////////////////

void ParticleTracking::buildCellIDmap()
{
  if (_CellIDmap.size() > 0) {
    _CellIDmap.clear();
  }
  
  CellTrsGeoBuilder::GeoData& cellData = _cellBuilder.getDataGE();
  const CFuint nCells = cellData.trs->getLocalNbGeoEnts();
  for(CFuint i=0; i<nCells; ++i){
    cellData.idx = i;
    GeometricEntity *const cell = _cellBuilder.buildGE();
    // AL: mapping is based on global IDs 
    _CellIDmap.insert(cell->getState(0)->getGlobalID(), i);
    _cellBuilder.releaseGE();
  }
  _CellIDmap.sortKeys();
}
    
//////////////////////////////////////////////////////////////////////////////
    
void ParticleTracking::setupStartPoint(CFuint startCellID)
{
  _startCellID = startCellID;
  // compute the coordinate of the start point: it is assumed that the start point is the baricentre of the startCell
  _cellBuilder.getDataGE().idx = _CellIDmap.find(_startCellID);
  GeometricEntity *const cell = _cellBuilder.buildGE();
  computeAverage(*cell->getNodes(), cell->nbNodes(), _startPoint);
  _cellBuilder.releaseGE();
}
    
//////////////////////////////////////////////////////////////////////////////

void ParticleTracking::setupStartPoint(CFuint startFaceLocalTrsID, const string& trs)
{  
  // compute the coordinate of the start point: it is assumed that the start point is the baricentre of the startFace
  FaceTrsGeoBuilder::GeoData& faceData = _faceBuilder.getDataGE();
  SafePtr<TopologicalRegionSet> wallFaces = MeshDataStack::getActive()->getTrs(trs);
  
  // cf_assert(trs != "InnerFaces");
  faceData.isBFace = true;
  faceData.trs = wallFaces;
  faceData.idx = startFaceLocalTrsID;
  
  GeometricEntity *const face = _faceBuilder.buildGE();
  // local face ID
  _wallFaceID = face->getID();                  
  // internal cell ID
  _startCellID = face->getState(0)->getGlobalID();
  computeAverage(*face->getNodes(), face->nbNodes(), _startPoint);
  _faceBuilder.releaseGE();
  
  CFLog(DEBUG_MAX, "ParticleTracking::setupStartPoint() => _wallFaceID = " << _wallFaceID << ", _startCellID = " << _startCellID << "\n");
}
    
//////////////////////////////////////////////////////////////////////////////
    
void ParticleTracking::setupPath(CFreal* pathInformation, 
				 InfoType InformationType)
{ 
  if(InformationType == ENDPOINT){
    _set = 1;
    for (CFuint d = 0; d < _dim; ++d) {
      _endPoint[d] = pathInformation[d];
    }
    _foundIt = false;
    _direction = _endPoint - _startPoint;
    const CFreal invNorm = 1./_direction.norm2();
    _direction *= invNorm;
    _internalPoint=_startPoint;
  }
  
  if(InformationType == DIRECTION){
    _set = 2;
    CFreal norm = 0.0;
    for (CFuint d = 0; d < _dim; ++d) {
      norm += pathInformation[d]*pathInformation[d];
    }
    norm = sqrt(norm);
    for (CFuint d = 0; d < _dim; ++d) {
      _direction[d] = pathInformation[d]/norm;
    }
    _foundIt = false;
  }
}
    
/////////////////////////////////////////////////////////////////////////////

//Maybe use the cell->getNeighborGeos() for caching.
//TODO: USE 2D algorithm with last ray tangent for first aproximation
//TODO: Don't need to calculate the intersections on an entry face(same as the last cell's exit face)
GeoEntOut ParticleTracking::myAxiRayTracing(Ray& ray, CFint currentCellID){
  //cout<<"MY_AXY_RAYTRACER"<<endl;
  CFreal a,b,c,x0,y0,z0,invA,D,D1,D2,B1,B2,z1,z2,r1,r2,dr,dz,zz,dr_2,dz_2,t[2],s[2];
  bool isS1,isS2;
  CFreal D22, cx0, cy0,a1,a2,ca2,c_2;
  GeoEntOut out;
  vector<CFreal>tCantidates;
  vector<CFuint>fCandidates;
  vector<CFreal>tdebug;
  RealVector faceOutNormal(2), rayTangent(2);
  tCantidates.reserve(5);
  fCandidates.reserve(5);
  //cout<<"on cell with ID: "<<currentCellID<<endl;
  DataHandle<CFint> faceIsOutwards= socket_isOutward.getDataHandle();

  a =ray.direction[0];   b=ray.direction[1];   c=ray.direction[2];
  x0=ray.startPoint[0]; y0=ray.startPoint[1]; z0=ray.startPoint[2];

    //Only need to calculate it once
    D22=pow(a*y0-b*x0,2);
    cx0=c*x0;
    cy0=c*y0;
    a1=pow(a,2)+pow(b,2);
    a2=(x0*a+y0*b);
    ca2=c*a2;
    c_2=pow(c,2);
  CellTrsGeoBuilder::GeoData& cellData = _cellBuilder.getDataGE();
  _cellIdx = _CellIDmap.find(currentCellID);
  cellData.idx = _cellIdx;

  GeometricEntity *const cell = _cellBuilder.buildGE();
  CFuint nFaces = cell->nbNeighborGeos();
  //cout<<"***********************************************"<<endl;
  for(CFuint f=0; f<nFaces; ++f){
    GeometricEntity* const face = cell->getNeighborGeo(f);
    const vector<Node*>& nodes = *face->getNodes();
    CFint faceID=face->getID();
    z1 = (*(nodes[0]))[XX]; z2 = (*(nodes[1]))[XX];
    r1 = (*(nodes[0]))[YY]; r2 = (*(nodes[1]))[YY];
    //cout<<"line: ( "<< z1 <<" , "<< r1 <<" ) to ( "<< z2 <<" , "<< r2 <<" )"<<endl;

    dr=r2-r1; dz=z2-z1; zz=z0-z1;
    dr_2=pow(dr,2); dz_2=pow(dz,2);

    faceOutNormal[0] =  dr;
    faceOutNormal[1] = -dz;
    // this way we assume that the parametrization runs clockwise over the cell's facets.
    // so that the face normal points outwards
    // CF doesn't impose this condition ( bug? ), so we have to correct it
    faceOutNormal*=(faceIsOutwards[faceID]==_cellIdx) ? 1.:-1.;

    D=pow(dr*(a*zz-cx0)+a*dz*r1,2)+pow(b*dz*r1+dr*(b*zz-cy0),2)-dz_2*D22;

    if (D>=0){
      invA=1/(dz_2*a1-c_2*dr_2);
      D=sqrt(D)*invA;
      D1=c*D;
      B1=(dz*(zz*a1-ca2)+c_2*dr*r1)*invA;

      s[0]=B1+D1; s[1]=B1-D1;

      isS1=s[0]<=1. && s[0]>=0.; // if the intersection point
      isS2=s[1]<=1. && s[1]>=0.; // lays between the vertices

      if(isS1|| isS2){
        //RealVector faceNormal(2);
        
        //cout<<"found something"<<endl;
        t[0]=-1.;t[1]=-1.;
        D2=dz*D;
        B2=(c*(zz*dr_2+r1*dz*dr)-dz_2*a2)*invA;
        CFuint kk=0;
        if (isS1){
          t[kk]=B2+D2;
          kk=kk+1;
        }
        if (isS2){
          t[kk]=B2-D2;
          kk=kk+1;
        }
        for (CFuint k=0;k<kk;++k){
          rayTangent[0]=c;
          rayTangent[1]=( a*( x0+a*t[k] ) + b*( y0+b*t[k]) )/ sqrt( pow(x0+a*t[k],2) + pow(y0+b*t[k],2) );

          //cout<<"fount t= "<<t[k]<<endl;
          tdebug.push_back(t[k]);
          if (t[k]>ray.tt && MathFunctions::innerProd(rayTangent,faceOutNormal)>=0){ //found exit point
            //cout<<"current t:"<<t[k]<<endl;
            tCantidates.push_back(t[k]);
            fCandidates.push_back(f);
          }
        }
      }
    }
    //cout<<currentCellID<<" got the nodes"<<endl;
  }
  //cout<<"number of candidates: "<<tCantidates.size()<<endl;
  if (tCantidates.size()<=0){
    CFLog(INFO,"Can't find an exit point! \n");
    out.exitCellID=-1;
    out.exitFaceID=-1;
    _cellBuilder.releaseGE();
    return out;
  }

  CFreal myT=tCantidates[0];
  CFuint index=0;
  for(CFuint i=1;i<tCantidates.size();++i){
    if(tCantidates[i]<myT){
      myT=tCantidates[i];
      index=i;
    }
  }

  ray.tt=myT;
  //cout<<"ray.t: "<<ray.tt<<endl;
  //cout<<"localGE face ID: "<<fCandidates[index]<<endl;
  //cout<<currentCellID<<" END2"<<endl;
  GeometricEntity* exitFace = cell->getNeighborGeo(fCandidates[index]);
  out.exitFaceID = exitFace->getID();

  out.exitCellID = exitFace->getState(0)->getGlobalID();
  if (out.exitCellID == currentCellID && !exitFace->getState(1)->isGhost()){
    out.exitCellID=exitFace->getState(1)->getGlobalID();
  }

  _cellBuilder.releaseGE();
  return out;

}



/////////////////////////////////////////////////////////////////////////////

void ParticleTracking::tracking(vector<CFuint>& CellList, CFuint actualCellID) // use this if _set = 1.
{ 
  //if(_set == 1){
    actualCellID = _startCellID;
    do{
      CellList.push_back(actualCellID);
      if(_dim == 2){ 
        actualCellID = tracking2D(actualCellID);
      }
      
      if(_dim == 3){
        actualCellID = tracking3D(actualCellID);
      }
    } while(_foundIt == false);
    //}
}

//////////////////////////////////////////////////////////////////////////////

CFuint ParticleTracking::tracking2D(CFuint cellID)
{ 
  // initialization
  _ExitFaceID = -1;
  
  CellTrsGeoBuilder::GeoData& cellData = _cellBuilder.getDataGE();
  _cellIdx = _CellIDmap.find(cellID);
  cellData.idx = _cellIdx;
  
  CFLog(DEBUG_MAX, "ParticleTracking::tracking2D() => global cellID = " << cellID << ", _cellIdx = " << _cellIdx << "\n");
  
  GeometricEntity *const cell = _cellBuilder.buildGE();
  
  // collect the data of the cell
  CFuint nFaces = cell->nbNeighborGeos();
  CFreal h = 0.0;
  CFreal t = 0.0;
  _barycentreActualCell = 0.;
  
  _nodeXcoordinate.clear();
  _nodeYcoordinate.clear();
  for(CFuint f=0; f<nFaces; ++f){
    GeometricEntity* const face = cell->getNeighborGeo(f);
    //cout<<"nb of nodes"<<face->nbNodes()<<endl;;
    const vector<Node*>& nodes = *face->getNodes();
    
    _nodeXcoordinate.push_back( (*(nodes[0]))[XX] );
    _nodeXcoordinate.push_back( (*(nodes[1]))[XX] );
    _nodeYcoordinate.push_back( (*(nodes[0]))[YY] );
    _nodeYcoordinate.push_back( (*(nodes[1]))[YY] );
    
    h = max(h, abs(MathFunctions::getDistance(*(nodes[0]),*(nodes[1]))));
    _barycentreActualCell += (*(nodes[0])) + (*(nodes[1]));
  }
  h *= nFaces; // h it is just a distance bigger than whatever segment obtained cutting the part
               // that lies into the actual cell of whatever line that cross the actual cell

  const CFreal inv2nFaces = 1./(2.*nFaces);
  _barycentreActualCell *= inv2nFaces;
  
  // compute the exit point if _set = 2
  if(_set == 2){
    if(_startCellID == cellID){
      _internalPoint = _startPoint;
    }
    else{
      // t is the distance between the startPoint and the projection of the barycentre of the actual cell on the path of the particle
      t = ((_barycentreActualCell[0] - _startPoint[0])*_direction[0]) + ((_barycentreActualCell[1] - _startPoint[1])*_direction[1]);
      
      // internalPoint is a point that lies into the actual cell and on the particle path
      _internalPoint = _startPoint + _direction*t;
    }
  
    // in case set=2, the end point has to be a whatever point beyond the actual cell on the particle path 
    _endPoint = _internalPoint + _direction*h;
  }

  // find the exit face
  sortNodes2D(_nodeXcoordinate, _nodeYcoordinate, _barycentreActualCell, nFaces);
  vector<CFuint> intersectedFaces;
  T2Ltest(_nodeXcoordinate, _nodeYcoordinate, nFaces, intersectedFaces);
  CFuint exitFace = 0;
  _foundIt  = P2Ltest(_nodeXcoordinate, _nodeYcoordinate, intersectedFaces, exitFace);
  
  // an exit face has to be found
  // cf_assert(exitFaceFound);
  _FaceNb = exitFace;
  
  // identifies the neighbour cell where is now the particle
  CFuint newCellID = cellID;
  if(_foundIt == false){
    GeometricEntity* ExitFace = cell->getNeighborGeo(exitFace);
    _ExitFaceID = ExitFace->getID();
    CFLog(DEBUG_MAX,"ParticleTracking::tracking2D() => #(foundIt == false)# _ExitFaceID = " <<_ExitFaceID <<" \n");
    const CFuint CellsFaceZero = ExitFace->getState(0)->getGlobalID();
    CFuint CellsFaceUno = 0;
    if (!ExitFace->getState(1)->isGhost()) {
      CellsFaceUno = ExitFace->getState(1)->getGlobalID();
      newCellID = (CellsFaceZero == cellID)? CellsFaceUno : CellsFaceZero;
    }
    else {
      // AL: this is a default behaviour not necessarily consistent ...
      newCellID = CellsFaceZero;
    }
  }
  
  _cellBuilder.releaseGE();
  return newCellID;

}

//////////////////////////////////////////////////////////////////////////////

CFuint ParticleTracking::tracking3D(CFuint cellID)
{ 
  // initialization
  _ExitFaceID = -1;
  
  CellTrsGeoBuilder::GeoData& cellData = _cellBuilder.getDataGE();
  _cellIdx = _CellIDmap.find(cellID);
  cellData.idx = _cellIdx;
  GeometricEntity *const cell = _cellBuilder.buildGE();
  const CFuint nFaces = cell->nbNeighborGeos();
  
  // compute the barycentre of the cell, it is needed for sorting the various face nodes
  const vector<Node*>& cellNodes = *cell->getNodes();
  computeAverage(cellNodes, cell->nbNodes(),_barycentreActualCell);
  const CFreal invDirLength =  1./_direction.norm2();
  
  // compute the internal and end point if _set = 2 
  // (it is not necessary for the algorithm the computation of the internal point, 
  // but in case of set=2 the user of ParticleTracking could need it)
  if(_set == 2){
    if(_startCellID == cellID){
      _internalPoint = _startPoint;
    }
    else{
      /*vector<CFreal> intersectionPoint = getIntersectionPoint();
      // t is the distance between the intersectionPoint and the projection of the barycentre of the actual cell on the path of the particle
      CFreal t = ((_barycentreActualCell[0] - intersectionPoint[0])*_direction[0])/(pow(pow(_direction[0],2) + pow(_direction[1],2) + pow(_direction[2],2),0.5)) +
                 ((_barycentreActualCell[1] - intersectionPoint[1])*_direction[1])/(pow(pow(_direction[0],2) + pow(_direction[1],2) + pow(_direction[2],2),0.5)) +
                 ((_barycentreActualCell[2] - intersectionPoint[2])*_direction[2])/(pow(pow(_direction[0],2) + pow(_direction[1],2) + pow(_direction[2],2),0.5));
      t = std::abs(t);
      
      // internalPoint is a point that lies into the actual cell and on the particle path
      _internalPoint[0] = intersectionPoint[0];// + _direction[0]*t;
      _internalPoint[1] = intersectionPoint[1];// + _direction[1]*t;
      _internalPoint[2] = intersectionPoint[2];// + _direction[2]*t;*/
      
      
      // t is the distance between the startPoint and the projection of the barycentre of the actual cell on the path of the particle
      const CFreal t = std::abs((((_barycentreActualCell[0] - _startPoint[0])*_direction[0]) +
				 ((_barycentreActualCell[1] - _startPoint[1])*_direction[1]) +
				 ((_barycentreActualCell[2] - _startPoint[2])*_direction[2]))*invDirLength);
      
      // internalPoint is a point that lies into the actual cell and on the particle path
      _internalPoint = _startPoint + _direction*t;
    }
    
    // in case set=2, the end point has to be a whatever point beyond the actual cell on the particle path 
    
    // AL: why 5. ????
    const CFreal h = 5.*std::abs((((_barycentreActualCell[0] - (*(cellNodes[0]))[XX])*_direction[0]) +
				  ((_barycentreActualCell[1] - (*(cellNodes[0]))[YY])*_direction[1]) +
				  ((_barycentreActualCell[2] - (*(cellNodes[0]))[ZZ])*_direction[2]))*invDirLength);
    
    _endPoint = _internalPoint + _direction*h;
  }
  
  // find the exit face
  bool exitFaceFound = false;
  CFuint exitFace = nFaces-1;
  for(CFuint f=0; f<nFaces; ++f){
    // collect the data of the face
    GeometricEntity* const face = cell->getNeighborGeo(f);
    const CFuint nFaceNodes = face->nbNodes();
    const vector<Node*>& faceNodes = *face->getNodes();
    computeAverage(faceNodes, nFaceNodes, _barycentreActualFace);

    _nodeXcoordinate.clear();
    _nodeYcoordinate.clear();
    _nodeZcoordinate.clear();
    for(CFuint n=0; n<nFaceNodes; ++n){
      _nodeXcoordinate.push_back( (*(faceNodes[n]))[XX] );
      _nodeYcoordinate.push_back( (*(faceNodes[n]))[YY] );
      _nodeZcoordinate.push_back( (*(faceNodes[n]))[ZZ] );
    }
    
    // compute the exit face
    sortNodes3D(_nodeXcoordinate, _nodeYcoordinate, _nodeZcoordinate, _barycentreActualCell, _barycentreActualFace, 0, nFaceNodes-1);
    exitFaceFound = T2Itest(_nodeXcoordinate, _nodeYcoordinate, _nodeZcoordinate);
    
    if(exitFaceFound){
      exitFace = f;
      if(_set == 1){
        _foundIt = P2Itest(_nodeXcoordinate, _nodeYcoordinate, _nodeZcoordinate);
      }
      break;
    }
  }
  
  // an exit face has to be found
  // cf_assert(exitFaceFound);
  _FaceNb = exitFace;

  //if(!exitFaceFound){
  //  cout<<"EXIT FACE NOT FOUND"<<endl;
  //}
  
  // identifies the neighbour cell where is now the particle
  CFuint newCellID = cellID;
  if(exitFaceFound){
    GeometricEntity* ExitFace = cell->getNeighborGeo(exitFace);
    _ExitFaceID = ExitFace->getID();
    const CFuint CellsFaceZero = ExitFace->getState(0)->getGlobalID(); // here global IDs
    CFuint CellsFaceUno = 0;
    if (!ExitFace->getState(1)->isGhost()) {
      CellsFaceUno = ExitFace->getState(1)->getGlobalID(); // here global IDs
      newCellID = (CellsFaceZero == cellID)? CellsFaceUno : CellsFaceZero;
    }
    else { 
      // AL: this is a default behaviour not necessarily consistent ...
      newCellID = CellsFaceZero;
    }
  }
  _cellBuilder.releaseGE();
  return newCellID;
}

//////////////////////////////////////////////////////////////////////////////

void ParticleTracking::sortNodes2D(vector<CFreal>& NodeXcoordinate, 
				   vector<CFreal>& NodeYcoordinate, 
				   const RealVector& CellBarycentre,
				   CFuint nFaces)
{
  const CFreal Xg = CellBarycentre[0];
  const CFreal Yg = CellBarycentre[1];
  for(CFuint i=0; i<nFaces; ++i){
    const CFreal G0x = NodeXcoordinate[2*i] - Xg;
    const CFreal G0y = NodeYcoordinate[2*i] - Yg;
    const CFreal G1x = NodeXcoordinate[2*i+1] - Xg;
    const CFreal G1y = NodeYcoordinate[2*i+1] - Yg;
    const CFreal V = G0x*G1y - G0y*G1x;
    if(V < 0){
      swap(NodeXcoordinate[2*i],NodeXcoordinate[2*i+1]);
      swap(NodeYcoordinate[2*i],NodeYcoordinate[2*i+1]);
    }
  }
}
    
//////////////////////////////////////////////////////////////////////////////

void ParticleTracking::sortNodes3D(vector<CFreal>& nodeXcoordinate, vector<CFreal>& nodeYcoordinate, vector<CFreal>& nodeZcoordinate, 
                                   const RealVector& barycentreActualCell, const RealVector& barycentreActualFace, CFuint p, CFuint r)
{  
  // It applies the merge-sort algorithm to sort the nodes in a recurrence way. 
  // In the merge step the nodes are sorted counterclockwise (looking at the face from the outside of the cell).
  if(p < r) { //p=index of the first element of the vectors; r=index of the last element of the vectors
    // recurrence step
    const CFuint q = (p + r)/2;
    sortNodes3D(nodeXcoordinate, nodeYcoordinate, nodeZcoordinate, barycentreActualCell, barycentreActualFace, p, q);
    sortNodes3D(nodeXcoordinate, nodeYcoordinate, nodeZcoordinate, barycentreActualCell, barycentreActualFace, q+1, r);

    // merge step
    const CFuint n1 = q - p + 1;
    const CFuint n2 = r - q;
    
    if (n1 > _Lx.size()) {
      _Lx.resize(n1);
      _Ly.resize(n1);
      _Lz.resize(n1);
    }
    
    if (n2 > _Rx.size()) {
      _Rx.resize(n2);
      _Ry.resize(n2);
      _Rz.resize(n2);
    }

    for(CFuint i=0; i<n1; ++i){
      _Lx[i] = nodeXcoordinate[p+i];
      _Ly[i] = nodeYcoordinate[p+i];
      _Lz[i] = nodeZcoordinate[p+i];
    }
    for(CFuint j=0; j<n2; ++j){
      _Rx[j] = nodeXcoordinate[q+j+1];
      _Ry[j] = nodeYcoordinate[q+j+1];
      _Rz[j] = nodeZcoordinate[q+j+1];
    }
    
    const CFreal GcGfx = barycentreActualFace[0] - barycentreActualCell[0];
    const CFreal GcGfy = barycentreActualFace[1] - barycentreActualCell[1];
    const CFreal GcGfz = barycentreActualFace[2] - barycentreActualCell[2];
    CFuint i = 0;
    CFuint j = 0;
    for(CFuint k=p; k<r+1; ++k){
      CFreal V = 0.;
      if(i<n1 && j<n2){
        const CFreal GfLx = _Lx[i] - barycentreActualFace[0];
        const CFreal GfLy = _Ly[i] - barycentreActualFace[1];
        const CFreal GfLz = _Lz[i] - barycentreActualFace[2];
        const CFreal GfRx = _Rx[j] - barycentreActualFace[0];
        const CFreal GfRy = _Ry[j] - barycentreActualFace[1];
        const CFreal GfRz = _Rz[j] - barycentreActualFace[2];
        V = GcGfx*(GfLy*GfRz - GfLz*GfRy) + GcGfy*(GfLz*GfRx - GfLx*GfRz) + GcGfz*(GfLx*GfRy - GfLy*GfRx);
      }
      if(i>=n1) V = -1.;
      if(j>=n2) V = 1.;

      if(V >= 0.){
        nodeXcoordinate[k] = _Lx[i];
        nodeYcoordinate[k] = _Ly[i];
        nodeZcoordinate[k] = _Lz[i];
        i++;
      }
      else{
        nodeXcoordinate[k] = _Rx[j];
        nodeYcoordinate[k] = _Ry[j];
        nodeZcoordinate[k] = _Rz[j];
        j++;
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

bool ParticleTracking::P2Ltest(const vector<CFreal>& nodeXcoordinate,
			       const vector<CFreal>& nodeYcoordinate, 
			       const vector<CFuint>& intersectedFaces, 
			       CFuint & exitFace)
{
  bool foundIt = true; // the endPoint is into the actualCell (supposing that the exitFace does not exist)
  CFuint nExitFaces = intersectedFaces.size();
  for(CFuint i=0; i<nExitFaces; ++i){
    const CFuint f = intersectedFaces[i];
    const CFreal Om = (nodeXcoordinate[2*f+1] - nodeXcoordinate[2*f])*(_endPoint[1] - nodeYcoordinate[2*f]) - 
      (_endPoint[0] - nodeXcoordinate[2*f])*(nodeYcoordinate[2*f+1] - nodeYcoordinate[2*f]);
    if(Om<0){
      exitFace = f;
      foundIt = false; // the exitFace exists and then the endPoint is outside the actualCell
      break;
    }
  }
  return foundIt; // always false if _set = 2
}

//////////////////////////////////////////////////////////////////////////////

bool ParticleTracking::T2Itest(const vector<CFreal>& nodeXcoordinate,
			                   const vector<CFreal>& nodeYcoordinate, 
			                   const vector<CFreal>& nodeZcoordinate)
{
  bool exitFaceFound = true;
  const CFreal SEx = _endPoint[0] - _startPoint[0];
  const CFreal SEy = _endPoint[1] - _startPoint[1];
  const CFreal SEz = _endPoint[2] - _startPoint[2];
  const CFuint nNodes = nodeXcoordinate.size();
  assert(nodeXcoordinate.size() <= 4);
  
  for(CFuint i=0; i<nNodes; ++i){
    CFuint j = i + 1;
    if(i == (nNodes-1) ) j = 0;
    const CFreal SAx = nodeXcoordinate[i] - _startPoint[0];
    const CFreal SAy = nodeYcoordinate[i] - _startPoint[1];
    const CFreal SAz = nodeZcoordinate[i] - _startPoint[2];
    const CFreal SBx = nodeXcoordinate[j] - _startPoint[0];
    const CFreal SBy = nodeYcoordinate[j] - _startPoint[1];
    const CFreal SBz = nodeZcoordinate[j] - _startPoint[2];
    CFreal V = SEx*(SAy*SBz - SAz*SBy) + SEy*(SAz*SBx - SAx*SBz) + SEz*(SAx*SBy - SAy*SBx);
    V /= std::abs(V);
        
    if (V <= 0.){
      exitFaceFound = false;
      break;
    }
  }
  return exitFaceFound;
}

//////////////////////////////////////////////////////////////////////////////

bool ParticleTracking::P2Itest(const vector<CFreal>& nodeXcoordinate,
			       const vector<CFreal>& nodeYcoordinate, 
			       const vector<CFreal>& nodeZcoordinate)
{

  bool foundIt = true;
  const CFuint nNodes = nodeXcoordinate.size();
  for(CFuint i=0; i<nNodes; ++i){
    const CFuint m = (i+2)%nNodes;  // i+2
    const CFuint n = (i+1)%nNodes;  // i+1
    const CFuint o = i;             // i
    const CFreal MNx = nodeXcoordinate[n] - nodeXcoordinate[m];
    const CFreal MNy = nodeYcoordinate[n] - nodeYcoordinate[m];
    const CFreal MNz = nodeZcoordinate[n] - nodeZcoordinate[m];
    const CFreal NOx = nodeXcoordinate[o] - nodeXcoordinate[n];
    const CFreal NOy = nodeYcoordinate[o] - nodeYcoordinate[n];
    const CFreal NOz = nodeZcoordinate[o] - nodeZcoordinate[n];
    const CFreal NEx = _endPoint[0] - nodeXcoordinate[n];
    const CFreal NEy = _endPoint[1] - nodeYcoordinate[n];
    const CFreal NEz = _endPoint[2] - nodeZcoordinate[n];
    const CFreal V = NEx*(MNy*NOz - MNz*NOy) + NEy*(MNz*NOx - MNx*NOz) + NEz*(MNx*NOy - MNy*NOx);
    if(V<0){
      foundIt = false;
      break;
    }
  }
  return foundIt;
}

//////////////////////////////////////////////////////////////////////////////

RealVector& ParticleTracking::getIntersectionPoint()
{
  CellTrsGeoBuilder::GeoData& cellData = _cellBuilder.getDataGE();
  cellData.idx = _cellIdx;
  GeometricEntity *const cell = _cellBuilder.buildGE();
  GeometricEntity* exitFace = cell->getNeighborGeo(_FaceNb);
  cf_assert(_ExitFaceID >= 0);
  const CFuint startID = _ExitFaceID*_dim;
  
  // compute face normal
  DataHandle<CFreal> normals = socket_normals.getDataHandle();
    
  _nn[0] = normals[startID];
  _nn[1] = normals[startID + 1];
  if(_dim == 3){
    _nn[2] = normals[startID + 2];
  }
  const CFreal invNLength = 1./_nn.norm2();
  _nn *= invNLength;
  
  // compute face baricentre
  computeAverage(*exitFace->getNodes(), exitFace->nbNodes(), _barycentreFace);
  _cellBuilder.releaseGE();
  
  const CFreal A = MathFunctions::innerProd(_nn,_internalPoint);
  const CFreal B = MathFunctions::innerProd(_nn,_direction);
  const CFreal D = MathFunctions::innerProd(_nn,_barycentreFace);
  const CFreal t = (D - A)/B;
  _intersectionPoint = _internalPoint + _direction*t;
  return _intersectionPoint;
}

//////////////////////////////////////////////////////////////////////////////

CFuint  ParticleTracking::getExitFaceGlobalGhostID(GeoEntOut out)
{
  //std::cout<<"Im here!"<<std::endl;
  //FaceTrsGeoBuilder::GeoData& faceData = _faceBuilder.getDataGE();
  //SafePtr<TopologicalRegionSet> ghostFaces = MeshDataStack::getActive()->getTrs("InnerFaces");
  // cf_assert(trs != "InnerFaces");
  //faceData.isBFace = false;
  //faceData.trs = ghostFaces;
  //faceData.idx = exitFaceID;
  //GeometricEntity const* exitFace = _faceBuilder.buildGE();
  CellTrsGeoBuilder::GeoData& cellData = _cellBuilder.getDataGE();
  _cellIdx = _CellIDmap.find(out.exitCellID);
  cellData.idx = _cellIdx;
  GeometricEntity *const cell = _cellBuilder.buildGE();
  CFuint nFaces = cell->nbNeighborGeos();
  //cout<<currentCellID<<" END1"<<endl;
  for(CFuint f=0; f<nFaces; ++f){
    GeometricEntity* const face = cell->getNeighborGeo(f);
    if (face->getID()==out.exitFaceID){
      CFuint ghostGlobalID;
      if(face->getState(1)->isGhost()){
        ghostGlobalID= face->getState(1)->getGlobalID();
        //std::cout<<"state 1 is ghost"<<std::endl;
      }
      else{
        ghostGlobalID= face->getState(0)->getGlobalID();
        //std::cout<<"state 0 is ghost"<<std::endl;
      }
      //cout<<"ExitFace global ID: "<<ghostGlobalID<<endl;
      _cellBuilder.releaseGE();
      return ghostGlobalID;

    }

  }
  cout<<"Error in getExitFaceGlobalGhostID !"<<endl;
  _cellBuilder.releaseGE();
  return 0;
}
/////////////////////////////////////////////////////////////////////////////
  } // namespace RadiativeTransfer

} // namespace COOLFluiD



