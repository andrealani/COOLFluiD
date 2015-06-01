#include "MathTools/RealMatrix.hh"
#include "MathTools/InverterT.hh"

#include "Framework/MeshData.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/GeometricEntity.hh"
#include "Framework/GeometricEntityPool.hh"
#include "Framework/FaceTrsGeoBuilder.hh"
#include "Framework/CellTrsGeoBuilder.hh"

#include "FiniteVolume/FiniteVolume.hh"
#include "FiniteVolume/StdSetup.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<StdSetup, CellCenterFVMData, FiniteVolumeModule> stdSetupProvider("StdSetup");

//////////////////////////////////////////////////////////////////////////////

void StdSetup::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< std::string >("stencil","stencil type");
  options.addConfigOption< std::string >
    ("InitLimiterSocket","name of the input limiter socket");
}

//////////////////////////////////////////////////////////////////////////////

StdSetup::StdSetup(const std::string& name) :
  BaseSetupFVMCC<CellCenterFVMCom>(name),
  socket_activeDiffusion("activeDiffusion"),
  _v1(DIM_3D),
  _v2(DIM_3D),
  _v3(DIM_3D),
  _xcCellApprox(DIM_3D)
{
  addConfigOptionsTo(this);
  
  _stencilType = "Face";
  setParameter("stencil",&_stencilType);
  
  _limiterSocketName = "Null";
  setParameter("InitLimiterSocket",&_limiterSocketName);
}

//////////////////////////////////////////////////////////////////////////////

StdSetup::~StdSetup()
{
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSource> > StdSetup::providesSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSource> > result = 
    BaseSetupFVMCC<CellCenterFVMCom>::providesSockets();
  
  result.push_back(&socket_activeDiffusion);
  
  return result;
}

//////////////////////////////////////////////////////////////////////////////

void StdSetup::configure ( Config::ConfigArgs& args )
{
  BaseSetupFVMCC<CellCenterFVMCom>::configure(args);
  
  // here call the mesh reader to ask if _limiterSocketName is one of the extra fields
  if (_limiterSocketName != "Null") {
    m_dynamicSockets.createSocketSink<CFreal>(_limiterSocketName, true);
  }
}      

//////////////////////////////////////////////////////////////////////////////

void StdSetup::execute()
{
  CFAUTOTRACE;
  
  CFLog(VERBOSE, "StdSetup::execute() => start\n");
  
  BaseSetupFVMCC<CellCenterFVMCom>::execute();
  
  // ---------------------------
  // AL: experiment
  // if (dim == DIM_3D) {
  //     cout << "BEFORE  compute3DGeometricData()" << endl;
  //     compute3DGeometricData();
  //     cout << "AFTER   compute3DGeometricData()" << endl;
  //   }
  // ---------------------------
  
  DataHandle<CFreal> activeDiffusion = socket_activeDiffusion.getDataHandle();
  activeDiffusion.resize(socket_states.getDataHandle().size());
  activeDiffusion = 1.;
  
  CFLog(VERBOSE, "StdSetup::execute() => end\n");
}
      
//////////////////////////////////////////////////////////////////////////////

void StdSetup::compute3DGeometricData()
{     
  // face normals
  DataHandle<CFreal> normals = socket_normals.getDataHandle();
  
  // face centroids
  DataHandle<CFreal> faceCenters = socket_faceCenters.getDataHandle();
    
  // temporary wrappers for pointers
  RealVector normalPtr(DIM_3D, static_cast<CFreal*>(NULL));
  RealVector xcFacePtr(DIM_3D, static_cast<CFreal*>(NULL));
  
  // list of TRSs
  vector<SafePtr<TopologicalRegionSet> > trs = MeshDataStack::getActive()->getTrsList();
    
  // face builder is created and setup
  SafePtr<GeometricEntityPool<FaceTrsGeoBuilder> > faceBuilder = getMethodData().getFaceTrsGeoBuilder();
  FaceTrsGeoBuilder::GeoData& geoData = faceBuilder->getDataGE();
  SafePtr<FaceTrsGeoBuilder> faceBuilderPtr = faceBuilder->getGeoBuilder();
  faceBuilderPtr->setDataSockets(socket_states, socket_gstates, socket_nodes);

  // here we loop over all (inner and boundary) faces
  const CFuint nbTRSs = trs.size();
  for (CFuint iTRS = 0; iTRS < nbTRSs; ++iTRS) {
    assert(iTRS < trs.size());
    SafePtr<TopologicalRegionSet> currTrs = trs[iTRS];
    
    if (currTrs->getName() != "InnerCells") {
      geoData.trs = currTrs;
      const CFuint nbFaces = currTrs->getLocalNbGeoEnts();
      for (CFuint iFace = 0; iFace < nbFaces; ++iFace) {
	geoData.idx = iFace;
	const GeometricEntity *const face = faceBuilder->buildGE();
	
	// we can compute the face normal directly with this approximated centroid
	const CFuint faceID  = face->getID();
	const CFuint startID = faceID*DIM_3D;
	
	// wrap the normal pointer
	assert(startID < normals.size());
	normalPtr.wrap(DIM_3D, &normals[startID]);
	normalPtr = 0.0;
	
	// now the centroid can be repositioned
	assert(startID < faceCenters.size());
	xcFacePtr.wrap(DIM_3D, &faceCenters[startID]);
	xcFacePtr = 0.0;
	
	computeFaceGeometry(face, normalPtr, xcFacePtr);
	
	faceBuilder->releaseGE();
      }
    }
  }
  
  // --------------------------
  // cell centroids and volumes
  // --------------------------
  
  DataHandle<CFreal> volumes = socket_volumes.getDataHandle();
  SafePtr<TopologicalRegionSet> cells = MeshDataStack::getActive()->getTrs("InnerCells");
  
  // create and setup the cell builder
  SafePtr<GeometricEntityPool<CellTrsGeoBuilder> > cellBuilder = getMethodData().getCellTrsGeoBuilder();
  SafePtr<CellTrsGeoBuilder> cellBuilderPtr = cellBuilder->getGeoBuilder();
  cellBuilderPtr->setDataSockets(socket_states, socket_gstates, socket_nodes);
  CellTrsGeoBuilder::GeoData& cellData = cellBuilder->getDataGE();
  cellData.trs = cells;
    
  ofstream oldv("volumes.old");
  ofstream newv("volumes.new");
  ofstream oldc("centroid.old");
  ofstream newc("centroid.new");
  
  for (CFuint iElem = 0; iElem < volumes.size(); ++iElem) {
    oldv << iElem << " => "; 
    oldv.precision(14); oldv.setf(ios::scientific,ios::floatfield); oldv << volumes[iElem] << endl;
  }
  oldv.close();
  
  CFuint count = 0;
  CFuint countNegativeVolumes = 0;
  const CFuint nbElems = cells->getLocalNbGeoEnts();
  for (CFuint iElem = 0; iElem < nbElems; ++iElem) {
    // build the GeometricEntity
    cellData.idx = iElem;
    GeometricEntity *const cell = cellBuilder->buildGE();
    
    // this is just a test
    // assignFaceNodes(cell);
    
    computeCellGeometry3D(iElem, cell, xcFacePtr, countNegativeVolumes, count);
    
    cellBuilder->releaseGE();
  }
  
  cout << "TEST 0: negative tet-volumes detected "  << countNegativeVolumes << endl;
  cout << "TEST 1: centroid is not in Cell for "  << count << " / " << volumes.size() << " cells" << endl;
  
  isCentroidInside();
  
  // for (;;) {}
  
  newv.close();
  oldc.close();
  newc.close();
  
  cout << "StdSetup::compute3DGeometricData()" << endl;
}

//////////////////////////////////////////////////////////////////////////////

bool StdSetup::isCentroidInside(CFuint cellID, const vector<GeometricEntity*>& facesInCell, const RealVector& xc)
{ 
  // face normals
  DataHandle<CFreal> normals = socket_normals.getDataHandle();
  
  // face centroids
  DataHandle<CFreal> faceCenters = socket_faceCenters.getDataHandle();
  
  // index of the element for which normal is outward
  DataHandle<CFint> isOutward = socket_isOutward.getDataHandle();
  
  // for the moment this works only on hexahedra
  RealVector vfc(DIM_3D);
  RealVector vn(DIM_3D);
  RealVector normalPtr(DIM_3D, static_cast<CFreal*>(NULL));
  RealVector midFacelPtr(DIM_3D, static_cast<CFreal*>(NULL));
    
  const CFuint nbFaces = facesInCell.size();
  for (CFuint iFace = 0; iFace < nbFaces; ++iFace) {
    // assume nodes are numbered in a counter-clockwise manner when viewed from inside the cell
    const GeometricEntity& face = *facesInCell[iFace];
    const CFuint faceID = face.getID();
    const CFuint startID = faceID*DIM_3D;
    normalPtr.wrap(DIM_3D, &normals[startID]);
    midFacelPtr.wrap(DIM_3D, &faceCenters[startID]);
    
    vfc = xc - midFacelPtr;
    
    for (CFuint d = 0; d < 3; ++d) {
      vn[d] = vfc[d];
    }
    
    if (isOutward[faceID] != static_cast<CFint>(cellID)) vn *= -1.;
    
    if (MathFunctions::innerProd(vn,normalPtr) >= 0.) {
      return false;
    }
  }
  return true;
}
      
      
//////////////////////////////////////////////////////////////////////////////
      
// bool StdSetup::isCentroidInside(CFuint iElem, const vector<GeometricEntity*>& facesInCell, const RealVector& xc)
// {
//   // for the moment this works only on hexahedra
//   static RealVector xa(DIM_3D);
//   static RealVector b(DIM_3D);
//   static RealVector c(DIM_3D);  
//   static RealVector d(DIM_3D);
//   static RealVector p(DIM_3D);
//   static RealVector xc2(DIM_3D);
//   static RealVector temp(DIM_3D); 
  
//   static RealMatrix matA(DIM_3D, DIM_3D, 0.);
//   static RealMatrix matB(DIM_3D, DIM_3D, 0.);
//   static RealMatrix matBi(DIM_3D, DIM_3D, 0.);
//   matA(0,0) = matA(0,2) = matA(1,1) = matA(1,2) = 1.;
//   static RealMatrix mm(DIM_3D, DIM_3D, 0.);
  
//   static InverterT<3> mInv;
    
//   const CFuint nbFaces = facesInCell.size();
//   for (CFuint iFace = 0; iFace < nbFaces; ++iFace) {
//     // assume nodes are numbered in a counter-clockwise manner when viewed from inside the cell
//     const GeometricEntity& face = *facesInCell[iFace];
//     const vector<Node*>& nodesInFace = face.getNodes();
//     const CFuint nbNodesInFace = nodesInFace.size();
    
    
    
    
//     // pick first node xa 
//     xa = *nodesInFace[0];
    
//     // translate all nodes to new coordinates system
//     // a = xa - xa = {0,0,0}  by construction 
//     b = *nodesInFace[1] - xa;
//     c = *nodesInFace[2] - xa;
//     d = *nodesInFace[3] - xa;
//     p = xc - xa;
    
//     const CFreal h = MathFunctions::mixedProd(b,d,c,temp);
//     matA(2,2) = h;
    
//     if (MathChecks::isZero(h)) {
//       // AL: we assume here that if the face is planar the centroid is internal
//       //     but a test is needed !!!
//       continue;
//     } 
//     else {      
      
//       for (CFuint i = 0; i < 3; ++i) {
// 	matB(i,0) =  b[i];
// 	matB(i,1) =  d[i];
// 	matB(i,2) =  c[i];
//       }
      
//       mInv.invert(matB, matBi);
//       mm = matA*matBi;
      
//       //   xc2 = mm*b; cout << "b = " << xc2  << endl;
//       //     xc2 = mm*c; cout << "c = " << xc2  << endl;
//       //     xc2 = mm*d; cout << "d = " << xc2  << endl;
      
//       xc2 = mm*p;
      
//       if (!(xc2[ZZ] > h*xc2[XX]*xc2[YY])) {
// 	// cout << "Face N. "<< iFace << " => " << h << ", determ3 = " << mm.determ3() << endl;
//         return false;
//       }
//     }
//   }
//   return true;
// }

//////////////////////////////////////////////////////////////////////////////

void StdSetup::isCentroidInside()
{
  // face centroids
  DataHandle<CFreal> faceCenters = socket_faceCenters.getDataHandle();
  
  // temporary wrappers for pointers
  RealVector xcFacePtr(DIM_3D, static_cast<CFreal*>(NULL));
  
  // list of TRSs
  vector<SafePtr<TopologicalRegionSet> > trs = MeshDataStack::getActive()->getTrsList();
  
  // face builder is created and setup
  SafePtr<GeometricEntityPool<FaceTrsGeoBuilder> > faceBuilder = getMethodData().getFaceTrsGeoBuilder();
  FaceTrsGeoBuilder::GeoData& geoData = faceBuilder->getDataGE();
  SafePtr<FaceTrsGeoBuilder> faceBuilderPtr = faceBuilder->getGeoBuilder();
  faceBuilderPtr->setDataSockets(socket_states, socket_gstates, socket_nodes);

  // here we loop over all (inner and boundary) faces
  const CFuint nbTRSs = trs.size();
  
  CFuint nbFaces = 0;
  CFuint count = 0;
  for (CFuint iTRS = 0; iTRS < nbTRSs; ++iTRS) {
    SafePtr<TopologicalRegionSet> currTrs = trs[iTRS];
    
    if (currTrs->getName() == "InnerFaces") {
      geoData.trs = currTrs;
      nbFaces = currTrs->getLocalNbGeoEnts();
      for (CFuint iFace = 0; iFace < nbFaces; ++iFace) {
	geoData.idx = iFace;
	const GeometricEntity *const face = faceBuilder->buildGE();
	const CFuint faceID  = face->getID();
	const CFuint startID = faceID*DIM_3D;
	
	// wrap the normal pointer
	xcFacePtr.wrap(DIM_3D,&faceCenters[startID]);
	
	_v1 = face->getState(0)->getCoordinates() - xcFacePtr;
	_v2 = face->getState(1)->getCoordinates() - xcFacePtr;
	
	if (MathFunctions::innerProd(_v1,_v2) >= 0.) {
	  count++;
	}
	
	faceBuilder->releaseGE();
      }
    }
  }
  
  if (count == 0) {
    cout << "StdSetup::isCentroidInside() => all centroids are internal to both cells for each face\n";
  }
  else {  
    cout << "StdSetup::isCentroidInside() => centroid is not in Cell for "  << count << " / " << nbFaces << " internal faces" << endl;
  }
}
      
//////////////////////////////////////////////////////////////////////////////

void StdSetup::assignFaceNodes(GeometricEntity *const cell)
{
  vector<GeometricEntity*>&  facesInCell = const_cast<vector<GeometricEntity*>&>(*cell->getNeighborGeos());
  static RealMatrix tnodes(8,3);
  const CFreal ah = 1e-6;
  
  tnodes(0,XX) = 0.; tnodes(0,YY) = 0.; tnodes(0,ZZ) = 0.;
  tnodes(1,XX) = 0.; tnodes(1,YY) = ah; tnodes(1,ZZ) = 0.;
  tnodes(2,XX) = ah; tnodes(2,YY) = ah; tnodes(2,ZZ) = 0.;
  tnodes(3,XX) = ah; tnodes(3,YY) = 0.; tnodes(3,ZZ) = 0.;
  tnodes(4,XX) = 0.; tnodes(4,YY) = 0.; tnodes(4,ZZ) = ah;
  tnodes(5,XX) = 0.; tnodes(5,YY) = ah; tnodes(5,ZZ) = ah;
  tnodes(6,XX) = ah; tnodes(6,YY) = ah; tnodes(6,ZZ) = ah/2.;
  tnodes(7,XX) = ah; tnodes(7,YY) = 0.; tnodes(7,ZZ) = ah;
  
  *facesInCell[0]->getNode(0) = tnodes.getRow<RealVector>(0);
  *facesInCell[0]->getNode(1) = tnodes.getRow<RealVector>(3);
  *facesInCell[0]->getNode(2) = tnodes.getRow<RealVector>(2);
  *facesInCell[0]->getNode(3) = tnodes.getRow<RealVector>(1);
  
  *facesInCell[1]->getNode(0) = tnodes.getRow<RealVector>(4);
  *facesInCell[1]->getNode(1) = tnodes.getRow<RealVector>(5);
  *facesInCell[1]->getNode(2) = tnodes.getRow<RealVector>(6);
  *facesInCell[1]->getNode(3) = tnodes.getRow<RealVector>(7);
  
  *facesInCell[2]->getNode(0) = tnodes.getRow<RealVector>(0);
  *facesInCell[2]->getNode(1) = tnodes.getRow<RealVector>(1);
  *facesInCell[2]->getNode(2) = tnodes.getRow<RealVector>(5);
  *facesInCell[2]->getNode(3) = tnodes.getRow<RealVector>(4);
  
  *facesInCell[3]->getNode(0) = tnodes.getRow<RealVector>(1);
  *facesInCell[3]->getNode(1) = tnodes.getRow<RealVector>(2);
  *facesInCell[3]->getNode(2) = tnodes.getRow<RealVector>(6);
  *facesInCell[3]->getNode(3) = tnodes.getRow<RealVector>(5);
  
  *facesInCell[4]->getNode(0) = tnodes.getRow<RealVector>(2);
  *facesInCell[4]->getNode(1) = tnodes.getRow<RealVector>(3);
  *facesInCell[4]->getNode(2) = tnodes.getRow<RealVector>(7);
  *facesInCell[4]->getNode(3) = tnodes.getRow<RealVector>(6);
  
  *facesInCell[5]->getNode(0) = tnodes.getRow<RealVector>(0);
  *facesInCell[5]->getNode(1) = tnodes.getRow<RealVector>(4);
  *facesInCell[5]->getNode(2) = tnodes.getRow<RealVector>(7);
  *facesInCell[5]->getNode(3) = tnodes.getRow<RealVector>(3);
}
    
//////////////////////////////////////////////////////////////////////////////

void StdSetup::computeFaceGeometry(const GeometricEntity *const face, 
				   RealVector& normalPtr, 
				   RealVector& xcFacePtr)				   
{
  // loop over the nodes of the current face and get a centroid approximation
  static RealVector xcFaceApprox(DIM_3D);
  const vector<Node*>& nodesInFace = face->getNodes();
  const CFuint nbNodesInFace = nodesInFace.size();
  assert(nbNodesInFace == 4);
  const CFreal nbNodesInFaceInv = 1./static_cast<CFreal>(nbNodesInFace);

  xcFaceApprox = 0.; 
  for (CFuint iNode = 0; iNode < nbNodesInFace; ++iNode) {
    xcFaceApprox += (*nodesInFace[iNode])*nbNodesInFaceInv;
  }
  
  CFuint ino2 = nbNodesInFace-1; // last node in the list
  for (CFuint iNode = 0; iNode < nbNodesInFace; ++iNode) {
    const CFuint ino1 = ino2;
    ino2 = iNode;
    
    assert(ino1 < nodesInFace.size());
    assert(ino2 < nodesInFace.size());
    
    _v1 = (*nodesInFace[ino1]) - xcFaceApprox;
    _v2 = (*nodesInFace[ino2]) - xcFaceApprox;
    
    normalPtr[XX] += 0.5*(_v1[YY]*_v2[ZZ] - _v1[ZZ]*_v2[YY]);
    normalPtr[YY] += 0.5*(_v1[ZZ]*_v2[XX] - _v1[XX]*_v2[ZZ]);
    normalPtr[ZZ] += 0.5*(_v1[XX]*_v2[YY] - _v1[YY]*_v2[XX]);
  }
  
  // loop through the edges again to get the normal pieces...
  CFreal sumArea = 0.0;
  ino2 = nbNodesInFace-1; // last node in the list
  for (CFuint iNode = 0; iNode < nbNodesInFace; ++iNode) {
    const CFuint ino1 = ino2;
    ino2 = iNode;
    assert(ino1 < nodesInFace.size());
    assert(ino2 < nodesInFace.size());
    
    _v1 = (*nodesInFace[ino1]) - xcFaceApprox;
    _v2 = (*nodesInFace[ino2]) - xcFaceApprox;
    
    const CFreal currArea = normalPtr[XX] * (_v1[YY] * _v2[ZZ] - _v1[ZZ] * _v2[YY]) +
                            normalPtr[YY] * (_v1[ZZ] * _v2[XX] - _v1[XX] * _v2[ZZ]) +
                            normalPtr[ZZ] * (_v1[XX] * _v2[YY] - _v1[YY] * _v2[XX]);
    
    if (currArea <= 0.0) {cout << " Warning: area negative -> " << currArea << endl;}
    sumArea += currArea;
    xcFacePtr += currArea*(xcFaceApprox + (*nodesInFace[ino1]) + (*nodesInFace[ino2]));
  }
  
  //	CF_DEBUG_OBJ(sumArea);
  assert(MathChecks::isNotZero(sumArea));	
  xcFacePtr /= (3.0*sumArea);
  
  //  if (xcFacePtr != xcFaceApprox) {
  //     cout.precision(14); cout.setf(ios::scientific,ios::floatfield); cout << xcFacePtr << " vs ";
  //     cout.precision(14); cout.setf(ios::scientific,ios::floatfield); cout << xcFaceApprox << "\n";
  //   }
  
  // remove this
  //  xcFacePtr = xcFaceApprox;
}
      
//////////////////////////////////////////////////////////////////////////////

void StdSetup::computeFaceGeometryIter(const GeometricEntity *const face, 
				       RealVector& normalPtr, 
				       RealVector& xcFacePtr)
{
  // loop over the nodes of the current face and get a centroid approximation
  static RealVector xcFaceApprox(DIM_3D);
  static RealVector dxc(DIM_3D);
  
  // sanity test: must converge in 1 iteration
  // vector<Node*>& nodesInFace = const_cast< vector<Node*>& > (face->getNodes());
  //   (*nodesInFace[0])[XX] = 0.; (*nodesInFace[0])[YY] = 0.; (*nodesInFace[0])[ZZ] = 0.;
  //   (*nodesInFace[1])[XX] = 1.; (*nodesInFace[1])[YY] = 0.; (*nodesInFace[1])[ZZ] = 0.;
  //   (*nodesInFace[2])[XX] = 1.; (*nodesInFace[2])[YY] = 1.; (*nodesInFace[2])[ZZ] = 0.;
  //   (*nodesInFace[3])[XX] = 0.; (*nodesInFace[3])[YY] = 1.; (*nodesInFace[3])[ZZ] = 0.;
  
  const vector<Node*>& nodesInFace = face->getNodes();
  const CFuint nbNodesInFace = nodesInFace.size();
  const CFreal nbNodesInFaceInv = 1./static_cast<CFreal>(nbNodesInFace);
  const CFreal ov3 = 1./3.;
  const CFreal eps = 1e-12;
  
  static RealVector dAdXterm(DIM_3D);
  static RealVector dAdX(DIM_3D);
  static RealVector dFdX(DIM_3D);
  
  // initial position of the centroid
  xcFaceApprox = 0.; 
  for (CFuint iNode = 0; iNode < nbNodesInFace; ++iNode) {
    xcFaceApprox += (*nodesInFace[iNode])*nbNodesInFaceInv;
  }
  
  // cout.precision(14); cout.setf(ios::scientific,ios::floatfield); cout << xcFaceApprox << ", after ";
  
  dxc = 0.;
  xcFacePtr = xcFaceApprox;
  CFuint iter = 0;
  
  do {
    xcFacePtr += dxc;
    dxc = 0.;
    
    // derivatives:  dAtdX*(2*x - xt1 - xt2)
    dAdXterm = 0.; 
    dAdX = 0.; 
    normalPtr = 0.;
        
    CFuint ino2 = nbNodesInFace-1; // last node in the list
    for (CFuint iNode = 0; iNode < nbNodesInFace; ++iNode) {
      const CFuint ino1 = ino2;
      ino2 = iNode;
      
      assert(ino1 < nodesInFace.size());
      assert(ino2 < nodesInFace.size());
      
      _v1 = (*nodesInFace[ino1]) - xcFacePtr;
      _v2 = (*nodesInFace[ino2]) - xcFacePtr;
      
      normalPtr[XX] += 0.5*(_v1[YY]*_v2[ZZ] - _v1[ZZ]*_v2[YY]);
      normalPtr[YY] += 0.5*(_v1[ZZ]*_v2[XX] - _v1[XX]*_v2[ZZ]);
      normalPtr[ZZ] += 0.5*(_v1[XX]*_v2[YY] - _v1[YY]*_v2[XX]);
    }
    
    // loop through the edges again to get the normal pieces...
    CFreal sumArea = 0.0;
    ino2 = nbNodesInFace-1; // last node in the list
    for (CFuint iNode = 0; iNode < nbNodesInFace; ++iNode) {
      const CFuint ino1 = ino2;
      ino2 = iNode;
      assert(ino1 < nodesInFace.size());
      assert(ino2 < nodesInFace.size());
      
      _v1 = (*nodesInFace[ino1]) - xcFacePtr;
      _v2 = (*nodesInFace[ino2]) - xcFacePtr;
      
      //  const CFreal currArea = normalPtr[XX] * (_v1[YY] * _v2[ZZ] - _v1[ZZ] * _v2[YY]) +
      //                             normalPtr[YY] * (_v1[ZZ] * _v2[XX] - _v1[XX] * _v2[ZZ]) +
      //                             normalPtr[ZZ] * (_v1[XX] * _v2[YY] - _v1[YY] * _v2[XX]);
      
      const CFreal ax = _v1[YY] * _v2[ZZ] - _v1[ZZ] * _v2[YY];
      const CFreal ay = _v1[ZZ] * _v2[XX] - _v1[XX] * _v2[ZZ];
      const CFreal az = _v1[XX] * _v2[YY] - _v1[YY] * _v2[XX];
      const CFreal currArea = 0.5*(ax*ax + ay*ay + az*az);
      
      if (currArea <= 0.0) {cout << " Warning: area negative -> " << currArea << endl;}
      sumArea += currArea;
      dxc += ov3*currArea*(xcFacePtr + (*nodesInFace[ino1]) + (*nodesInFace[ino2]));
      
      // dAdX = 0.5*(2*ax*daxdX + 2*ay*daydX + 2*az*dazdX) = ax*daxdX + ay*daydX + az*dazdX;
      dAdX[XX] = ay*(-_v1[ZZ] + _v2[ZZ]) + az*(-_v2[YY] + _v1[YY]);
      dAdX[YY] = ax*(-_v2[ZZ] + _v1[ZZ]) + az*(-_v1[XX] + _v2[XX]);
      dAdX[ZZ] = ax*(-_v1[YY] + _v2[YY]) + ay*(-_v2[XX] + _v1[XX]);
      
      dAdXterm += dAdX*(2.*xcFacePtr - (*nodesInFace[ino1]) - (*nodesInFace[ino2]));
    }
    
    assert(MathChecks::isNotZero(sumArea));	
    
    dFdX = ov3*(2.*sumArea + dAdXterm);
    dxc -= (xcFacePtr*sumArea);
    
    for (CFuint i = 0; i < 3; ++i) {
      assert(MathChecks::isNotZero(dFdX[i]));
      dxc[i] /= dFdX[i];
    }
    
    iter++;
  } while (dxc.norm2() > eps && iter < 30);
  
  // cout.precision(14); cout.setf(ios::scientific,ios::floatfield); cout << dxc.norm2() << endl;
  
  // cout.precision(14); cout.setf(ios::scientific,ios::floatfield); cout << " Iter[" << iter << "] => " << xcFacePtr << "\n";
}
      
//////////////////////////////////////////////////////////////////////////////

void StdSetup::computeCellGeometry3D(CFuint iElem, 
				     const GeometricEntity *const cell, 
				     RealVector& xcFacePtr,
				     CFuint& countNegativeVolumes, 
				     CFuint& count)
{
  // face centroids
  DataHandle<CFreal> faceCenters = socket_faceCenters.getDataHandle();
  
  // compute an approximate center as the mean of the face xcFace's...
  const vector<GeometricEntity*>& facesInCell = *cell->getNeighborGeos();
  const CFuint nbFacesInCell = facesInCell.size();
  const CFreal nbFacesInCellInv = 1./static_cast<CFreal>(facesInCell.size());
  
  _xcCellApprox = 0.; 
  for (CFuint iFace = 0; iFace < nbFacesInCell; ++iFace) {
    const CFuint faceID = facesInCell[iFace]->getID();
    xcFacePtr.wrap(DIM_3D, &faceCenters[faceID*DIM_3D]);
    _xcCellApprox += xcFacePtr*nbFacesInCellInv;
  }
  
  // oldc << iElem << " => "; 
  // oldc.precision(14); oldc.setf(ios::scientific,ios::floatfield); oldc << _xcCellApprox << endl;

  DataHandle<CFreal> volumes = socket_volumes.getDataHandle();
  
  // add tetrahedra volumes
  volumes[iElem] = 0.0;
  
  Node& xcCell = cell->getState(0)->getCoordinates();
  xcCell = 0.0;
  
  // loop on the faces of this cell
  for (CFuint iFace = 0; iFace < nbFacesInCell; ++iFace) {
    const GeometricEntity *const face = facesInCell[iFace];
    const CFuint faceID = facesInCell[iFace]->getID();
    xcFacePtr.wrap(DIM_3D, &faceCenters[faceID*DIM_3D]);
    _v1 = xcFacePtr - _xcCellApprox;
    
    const vector<Node*>& nodesInFace = face->getNodes();
    const CFuint nbNodesInFace = nodesInFace.size();
    
    // check if the normal to the face is outward with respect to this cell
    const CFuint cellHavingNormalOutward = static_cast<CFuint>(socket_isOutward.getDataHandle()[faceID]);
    
    if ((cellHavingNormalOutward == iElem && nbFacesInCell == 6) || 
	(cellHavingNormalOutward != iElem && nbFacesInCell == 4)) {
      // face is outward - loop through edges in forward direction...
      CFuint ino2 = nbNodesInFace - 1;
      for (CFuint iNode = 0; iNode < nbNodesInFace; ++iNode) {
	const CFuint ino1 = ino2;
	ino2 = iNode;
	
	_v2 = (*nodesInFace[ino1]) - _xcCellApprox;
	_v3 = (*nodesInFace[ino2]) - _xcCellApprox;
	
	// 2 nodes, the face, and the approx cv form a tet...
	const CFreal currVolume = _v1[XX] * (_v2[YY] * _v3[ZZ] - _v2[ZZ] * _v3[YY]) + 
	  _v1[YY] * (_v2[ZZ] * _v3[XX] - _v2[XX] * _v3[ZZ]) + 
	  _v1[ZZ] * (_v2[XX] * _v3[YY] - _v2[YY] * _v3[XX]);
	
	if (currVolume <= 0.0) {
	  countNegativeVolumes++;
	}
	
	// update volume
	volumes[iElem] += currVolume;
        
	// update cell centroid
	xcCell += currVolume*(_xcCellApprox + xcFacePtr + (*nodesInFace[ino1]) + (*nodesInFace[ino2]));
      }
    }
    else {
      // face is outward, loop through edges in backward direction...
      CFuint ino2 = 0;
      for (CFint iNode = nbNodesInFace - 1; iNode >= 0; iNode--) {
	const CFuint ino1 = ino2;
	ino2 = iNode;
	
	_v2 = (*nodesInFace[ino1]) - _xcCellApprox;
	_v3 = (*nodesInFace[ino2]) - _xcCellApprox;
	
	const CFreal currVolume = _v1[XX] * (_v2[YY] * _v3[ZZ] - _v2[ZZ] * _v3[YY]) + 
	  _v1[YY] * (_v2[ZZ] * _v3[XX] - _v2[XX] * _v3[ZZ]) + 
	  _v1[ZZ] * (_v2[XX] * _v3[YY] - _v2[YY] * _v3[XX]);
	
	if (currVolume <= 0.0) {
	  countNegativeVolumes++;
	}
	
	// update volume
	volumes[iElem] += currVolume;
	
	// update cell centroid
	xcCell += currVolume*(_xcCellApprox + xcFacePtr + (*nodesInFace[ino1]) + (*nodesInFace[ino2]));
      }
    }
  }
  
  // normalize both ...
  assert(volumes[iElem] > 0.);
  xcCell /= 4.0*volumes[iElem];
  
//   xcCell = 0.;
//   for (int i = 0; i < 8; ++i) {
//     xcCell += *cell->getNode(i);
//   }
//   xcCell /= 8.;
  
  if (!isCentroidInside(iElem, facesInCell, xcCell)) count++;
  
  //newc << iElem << " => "; 
  //newc.precision(14); newc.setf(ios::scientific,ios::floatfield); newc << xcCell << endl;
  
  volumes[iElem] /= 6.0;
  
  //newv << iElem << " => "; 
  //newv.precision(14); newv.setf(ios::scientific,ios::floatfield); newv << volumes[iElem] << endl;
  
  assert(volumes[iElem] > 0.);
}

//////////////////////////////////////////////////////////////////////////////

void StdSetup::computeCellGeometry3DIter(CFuint iElem, 
					 const GeometricEntity *const cell, 
					 RealVector& xcFacePtr,
					 CFuint& countNegativeVolumes, 
					 CFuint& count)
{
  // face centroids
  DataHandle<CFreal> faceCenters = socket_faceCenters.getDataHandle();
  
  // compute an approximate center as the mean of the face xcFace's...
  const vector<GeometricEntity*>& facesInCell = *cell->getNeighborGeos();
  const CFuint nbFacesInCell = facesInCell.size();
  const CFreal nbFacesInCellInv = 1./static_cast<CFreal>(facesInCell.size());
  
  _xcCellApprox = 0.; 
  for (CFuint iFace = 0; iFace < nbFacesInCell; ++iFace) {
    const CFuint faceID = facesInCell[iFace]->getID();
    xcFacePtr.wrap(DIM_3D, &faceCenters[faceID*DIM_3D]);
    _xcCellApprox += xcFacePtr*nbFacesInCellInv;
  }
  
  DataHandle<CFreal> volumes = socket_volumes.getDataHandle();
  // add tetrahedra volumes
  volumes[iElem] = 0.0;
  
  Node& xcCell = cell->getState(0)->getCoordinates();
  xcCell = 0.0;
  
  static RealVector dxc(DIM_3D);
  static RealVector dAdXterm(DIM_3D);
  static RealVector dAdX(DIM_3D);
  static RealVector dFdX(DIM_3D);
  
  const CFuint countNegativeVolumesBkp = countNegativeVolumes;
  const CFreal ov6 = 1./6.;
  const CFreal eps = 1e-12;
   
  dxc = 0.;
  xcCell = _xcCellApprox;
  CFuint iter = 0;
  
  do {
    xcCell += dxc;
    dxc = 0.;
    
    // derivatives:  dAtdX*(3*x - xf - xt1 - xt2)
    dAdXterm = 0.; 
    dAdX = 0.; 
    volumes[iElem] = 0.;
    countNegativeVolumes = countNegativeVolumesBkp;
    
    // loop on the faces of this cell
    for (CFuint iFace = 0; iFace < nbFacesInCell; ++iFace) {
      const GeometricEntity *const face = facesInCell[iFace];
      const CFuint faceID = facesInCell[iFace]->getID();
      xcFacePtr.wrap(DIM_3D, &faceCenters[faceID*DIM_3D]);
      _v1 = xcFacePtr - xcCell;
      
      const vector<Node*>& nodesInFace = face->getNodes();
      const CFuint nbNodesInFace = nodesInFace.size();
      
      // check if the normal to the face is outward with respect to this cell
      const CFuint cellHavingNormalOutward = static_cast<CFuint>(socket_isOutward.getDataHandle()[faceID]);
      
      if ((cellHavingNormalOutward == iElem && nbFacesInCell == 6) || 
	  (cellHavingNormalOutward != iElem && nbFacesInCell == 4)) {
	// face is outward - loop through edges in forward direction...
	CFuint ino2 = nbNodesInFace - 1;
	for (CFuint iNode = 0; iNode < nbNodesInFace; ++iNode) {
	  const CFuint ino1 = ino2;
	  ino2 = iNode;
	  
	  _v2 = (*nodesInFace[ino1]) - xcCell;
	  _v3 = (*nodesInFace[ino2]) - xcCell;
	  
	  // 2 nodes, the face, and the approx cv form a tet...
	  const CFreal ax = _v2[YY] * _v3[ZZ] - _v2[ZZ] * _v3[YY];
	  const CFreal ay = _v2[ZZ] * _v3[XX] - _v2[XX] * _v3[ZZ];
	  const CFreal az = _v2[XX] * _v3[YY] - _v2[YY] * _v3[XX];
	  const CFreal currVolume = ov6*(_v1[XX]*ax + _v1[YY]*ay + _v1[ZZ]*az);
	  
	  if (currVolume <= 0.0) {
	    countNegativeVolumes++;
	  }
	  
	  // update volume
	  volumes[iElem] += currVolume;
	  
	  // update cell centroid
	  dxc += 0.25*currVolume*(xcCell + xcFacePtr + (*nodesInFace[ino1]) + (*nodesInFace[ino2]));
	  dAdX[XX] = ov6*(-ax + _v1[YY]*(-_v2[ZZ] + _v3[ZZ]) + _v1[ZZ]*(-_v3[YY] + _v2[YY]));
	  dAdX[YY] = ov6*(_v1[XX]*(-_v3[ZZ] + _v2[ZZ]) - ay + _v1[ZZ]*(-_v2[XX] + _v3[XX]));
	  dAdX[ZZ] = ov6*(_v1[XX]*(-_v2[YY] + _v3[YY]) + _v1[YY]*(-_v3[XX] + _v2[XX]) - az);
	  dAdXterm += dAdX*(3.*xcCell - xcFacePtr - (*nodesInFace[ino1]) - (*nodesInFace[ino2]));
	}
      }
      else {
	// face is outward, loop through edges in backward direction...
	CFuint ino2 = 0;
	for (CFint iNode = nbNodesInFace - 1; iNode >= 0; iNode--) {
	  const CFuint ino1 = ino2;
	  ino2 = iNode;
	  
	  _v2 = (*nodesInFace[ino1]) - xcCell;
	  _v3 = (*nodesInFace[ino2]) - xcCell;
	  
	  // 2 nodes, the face, and the approx cv form a tet...
	  const CFreal ax = _v2[YY] * _v3[ZZ] - _v2[ZZ] * _v3[YY];
	  const CFreal ay = _v2[ZZ] * _v3[XX] - _v2[XX] * _v3[ZZ];
	  const CFreal az = _v2[XX] * _v3[YY] - _v2[YY] * _v3[XX];
	  const CFreal currVolume = ov6*(_v1[XX]*ax + _v1[YY]*ay + _v1[ZZ]*az);
	  
	  if (currVolume <= 0.0) {
	    countNegativeVolumes++;
	  }
	  
	  // update volume
	  volumes[iElem] += currVolume;
	  
	  // update cell centroid
	  dxc += 0.25*currVolume*(xcCell + xcFacePtr + (*nodesInFace[ino1]) + (*nodesInFace[ino2]));
	  
	  dAdX[XX] = ov6*(-ax + _v1[YY]*(-_v2[ZZ] + _v3[ZZ]) + _v1[ZZ]*(-_v3[YY] + _v2[YY]));
	  dAdX[YY] = ov6*(_v1[XX]*(-_v3[ZZ] + _v2[ZZ]) - ay + _v1[ZZ]*(-_v2[XX] + _v3[XX]));
	  dAdX[ZZ] = ov6*(_v1[XX]*(-_v2[YY] + _v3[YY]) + _v1[YY]*(-_v3[XX] + _v2[XX]) - az);
	  dAdXterm += dAdX*(3.*xcCell - xcFacePtr - (*nodesInFace[ino1]) - (*nodesInFace[ino2]));
	}
      }
    }
    
    assert(volumes[iElem] > 0.);
    
    dFdX = 0.25*(3.*volumes[iElem] + dAdXterm);
    dxc -= (xcCell*volumes[iElem]);
    
    for (CFuint i = 0; i < 3; ++i) {
      assert(MathChecks::isNotZero(dFdX[i]));
      dxc[i] /= dFdX[i];
    }
    
    iter++;
  } while (dxc.norm2() > eps && iter < 30);
  
  cout << "Iter = " << iter << endl;
  
  if (!isCentroidInside(iElem, facesInCell, xcCell)) count++;
  assert(volumes[iElem] > 0.);
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD
