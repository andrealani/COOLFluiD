#include "DiamondVolume3DDerivative.hh"
#include "Framework/MethodStrategyProvider.hh"
#include "Framework/GeometricEntity.hh"
#include "Common/NotImplementedException.hh"
#include "FiniteVolume/FiniteVolume.hh"

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

MethodStrategyProvider<DiamondVolume3DDerivative,
                       CellCenterFVMData,
		       DerivativeComputer,
                       FiniteVolumeModule>
diamondVolume3DDerivativeProvider("DiamondVolume3D");

//////////////////////////////////////////////////////////////////////////////

DiamondVolume3DDerivative::DiamondVolume3DDerivative(const std::string& name) :
  DerivativeComputer(name),
  _third(1./3.),
  socket_states("states"),
  socket_nodes("nodes"),
  socket_gstates("gstates"), 
  socket_cellFlag("cellFlag"),
  _normalCalculator(),
  _n031(),
  _n132(),
  _n023(),
  _n041(),
  _n142(),
  _n024(),
  _midNode(),
  _tmpNormal(),
  _v1(),
  _v2(),
  _v3(),
  _v4(),
  _xproj(),
  _ncoord(),
  _lcoord(),
  _rcoord(),
  _tetraCoord(),
  _pyramCoord(),
  _pyramNormals1(),
  _pyramNormals2(),
  _cVolume(),
  _cFace(),
  _ncF(),
  _cVcF(),
  _n4cF(),
  _n5cF(),
  _cvolumes(),
  _changeNormalSign(),
  _nonPlanarFaceNodeID()
{
}

//////////////////////////////////////////////////////////////////////////////

DiamondVolume3DDerivative::~DiamondVolume3DDerivative()
{
}

//////////////////////////////////////////////////////////////////////////////

void DiamondVolume3DDerivative::setup()
{
  DerivativeComputer::setup();

  _n031.resize(PhysicalModelStack::getActive()->getDim());
  _n132.resize(PhysicalModelStack::getActive()->getDim());
  _n023.resize(PhysicalModelStack::getActive()->getDim());
  _n041.resize(PhysicalModelStack::getActive()->getDim());
  _n142.resize(PhysicalModelStack::getActive()->getDim());
  _n024.resize(PhysicalModelStack::getActive()->getDim());
  _midNode.resize(PhysicalModelStack::getActive()->getDim());
  _tmpNormal.resize(PhysicalModelStack::getActive()->getDim());
  _v1.resize(PhysicalModelStack::getActive()->getDim());
  _v2.resize(PhysicalModelStack::getActive()->getDim());
  _v3.resize(PhysicalModelStack::getActive()->getDim());
  _v4.resize(PhysicalModelStack::getActive()->getDim());
  _xproj.resize(PhysicalModelStack::getActive()->getDim());
  _ncoord.resize(PhysicalModelStack::getActive()->getDim());
  _lcoord.resize(PhysicalModelStack::getActive()->getDim());
  _rcoord.resize(PhysicalModelStack::getActive()->getDim());
  _tetraCoord.resize(4, PhysicalModelStack::getActive()->getDim());
  _pyramCoord.resize(5, PhysicalModelStack::getActive()->getDim());
  _pyramNormals1.resize(5, PhysicalModelStack::getActive()->getDim());
  _pyramNormals2.resize(5, PhysicalModelStack::getActive()->getDim());
  _cVolume.resize(3);
  _cFace.resize(3);
  _ncF.resize(3);
  _cVcF.resize(3);
  _n4cF.resize(3);
  _n5cF.resize(3);
  
  // control volumes around faces are here precomputed to enhance efficiency
  // acording to profiler-driven indications
  const CFuint nbFaces = MeshDataStack::getActive()->Statistics().getNbFaces();
  _cvolumes.resize(nbFaces);
  _changeNormalSign.resize(nbFaces);
  _changeNormalSign.assign(nbFaces, false);
  
  _nonPlanarFaceNodeID.resize(nbFaces);
  _nonPlanarFaceNodeID.assign(nbFaces, -1);
  
  // set the list of faces
  vector<Common::SafePtr<TopologicalRegionSet> > trs = MeshDataStack::getActive()->getTrsList();
  const CFuint nbTRSs = trs.size();
  
  // prepare the building of the faces
  GeometricEntityPool<FaceCellTrsGeoBuilder> geoBuilder;
  geoBuilder.getGeoBuilder()->setCellFlagSocket(socket_cellFlag);
  geoBuilder.getGeoBuilder()->setDataSockets(socket_states, socket_gstates, socket_nodes);
  geoBuilder.setup();
  FaceCellTrsGeoBuilder::GeoData& geoData = geoBuilder.getDataGE();
  geoData.allCells = true;
  geoData.cells = MeshDataStack::getActive()->getTrs("InnerCells");
  
  // fill the _cvolumes and _changeNormalSign arrays
  for (CFuint iTRS = 0; iTRS < nbTRSs; ++iTRS) {
    SafePtr<TopologicalRegionSet> currTrs = trs[iTRS];
    geoData.faces = currTrs;
    
    // the faces on the boundary of the partition don't have to
    // be processed (their fluxes could give NaN)
    if (currTrs->getName() != "PartitionFaces" && currTrs->getName() != "InnerCells") {
      if (currTrs->hasTag("writable")) {
	geoData.isBFace = true;
      }
      else {
        geoData.isBFace = false;
      }
      
      const CFuint nbTrsFaces = currTrs->getLocalNbGeoEnts();
      for (CFuint iFace = 0; iFace < nbTrsFaces; ++iFace) {
	
	// build the GeometricEntity
        geoData.idx = iFace;
        GeometricEntity *const geo = geoBuilder.buildGE();
	
	if (geo->nbNodes() == 3) {
	  compute2TetraVolume(geo);
	}
	else if (geo->nbNodes() == 4) {
	  compute2PyramidVolume(geo);
	}
	
	_cvolumes[geo->getID()] = _volume;
	cf_always_assert(_volume > 0.0);
	
	geoBuilder.releaseGE();
      }
    }
  }
}
      
//////////////////////////////////////////////////////////////////////////////

void DiamondVolume3DDerivative::computeGradients(const RealMatrix& values,
						 vector<RealVector*>& gradients)
{
  if (getMethodData().getCurrentFace()->nbNodes() == 3) {
    compute2TetraGradient(values, gradients);
  }
  else {
    cf_assert(getMethodData().getCurrentFace()->nbNodes() == 4);
    compute2PyramidGradient(values, gradients);
  }
}

//////////////////////////////////////////////////////////////////////////////

void DiamondVolume3DDerivative::compute2TetraGradient
(const RealMatrix& values, vector<RealVector*>& gradients)
{
  GeometricEntity& geo = *getMethodData().getCurrentFace();
  
  // construct the 2 tetras volume and the normals around the face
  Node& node0 = *geo.getNode(0);
  Node& node1 = *geo.getNode(1);
  Node& node2 = *geo.getNode(2);
  Node& node3 = geo.getState(0)->getCoordinates();
  Node& node4 = geo.getState(1)->getCoordinates();
  
  // set the coordinates of the tetra 0132 in a matrix
  _tetraCoord.setRow(node0,0);
  _tetraCoord.setRow(node1,1);
  _tetraCoord.setRow(node2,2);
  _tetraCoord.setRow(node3,3);
  
  CFreal factor = (!_changeNormalSign[geo.getID()]) ? -0.5 : 0.5;
  
  _v1 = node3 - node0;
  _v2 = node1 - node0;
  MathFunctions::crossProd(_v1, _v2, _v3);
  _n031 = factor*_v3;
  
  _v1 = node3 - node1;
  _v2 = node2 - node1;
  MathFunctions::crossProd(_v1, _v2, _v3);
  _n132 = factor*_v3;
  
  _v1 = node2 - node0;
  _v2 = node3 - node0;
  MathFunctions::crossProd(_v1, _v2, _v3);
  _n023 = factor*_v3;
  
  // set the coordinates of the tetra 0124 in a matrix
  _tetraCoord.setRow(node4,3);
  
  // the factor is surely the opposite of the previous one
  factor *= -1;
  
  _v1 = node4 - node0;
  _v2 = node1 - node0;
  MathFunctions::crossProd(_v1, _v2, _v3);
  _n041 = factor*_v3;
  
  _v1 = node4 - node1;
  _v2 = node2 - node1;
  MathFunctions::crossProd(_v1, _v2, _v3);
  _n142 = factor*_v3;
  
  _v1 = node2 - node0;
  _v2 = node4 - node0;
  MathFunctions::crossProd(_v1, _v2, _v3);
  _n024 = factor*_v3;
  
  // Green Gauss is applied in the diamond volume
  // to compute the gradients
  const CFreal invVolume = _third/_volume;
  const CFuint nbVars = values.nbRows();
  cf_assert(nbVars == PhysicalModelStack::getActive()->getNbEq());
  cf_assert(gradients.size() == nbVars);
  
  for (CFuint i = 0; i < nbVars; ++i) {
    RealVector& grad = *gradients[i];
    
    const CFreal v031 = values(i,0) + values(i,3) + values(i,1);
    const CFreal v132 = values(i,1) + values(i,3) + values(i,2);
    const CFreal v023 = values(i,0) + values(i,2) + values(i,3);
    const CFreal v041 = values(i,0) + values(i,4) + values(i,1);
    const CFreal v142 = values(i,1) + values(i,4) + values(i,2);
    const CFreal v024 = values(i,0) + values(i,2) + values(i,4);
    
    for (CFuint ix = 0; ix < DIM_3D; ++ix) {
      grad[ix] = invVolume*(_n031[ix]*v031 + _n132[ix]*v132 + _n023[ix]*v023 +
			    _n041[ix]*v041 + _n142[ix]*v142 + _n024[ix]*v024);
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void DiamondVolume3DDerivative::compute2PyramidGradient
(const RealMatrix& values, vector<RealVector*>& gradients)
{
  // Green Gauss is applied in the diamond volume
  // to compute the gradients
  const CFreal invVolume = _third/_volume;
  const CFuint nbVars = values.nbRows();
  cf_assert(nbVars == PhysicalModelStack::getActive()->getNbEq());
  cf_assert(gradients.size() == nbVars);
  
  for (CFuint i = 0; i < nbVars; ++i) {
    RealVector& grad = *gradients[i];

    // pyramid 01234
    const CFreal v014 = values(i,0) + values(i,1) + values(i,4);
    const CFreal v124 = values(i,1) + values(i,2) + values(i,4);
    const CFreal v234 = values(i,2) + values(i,3) + values(i,4);
    const CFreal v043 = values(i,0) + values(i,4) + values(i,3);

    // pyramid 01235
    const CFreal v015 = values(i,0) + values(i,1) + values(i,5);
    const CFreal v125 = values(i,1) + values(i,2) + values(i,5);
    const CFreal v235 = values(i,2) + values(i,3) + values(i,5);
    const CFreal v053 = values(i,0) + values(i,5) + values(i,3);

    for (CFuint ix = 0; ix < DIM_3D; ++ix) {
      grad[ix] = invVolume*(_pyramNormals1(1,ix)*v014 +
			    _pyramNormals1(2,ix)*v124 +
			    _pyramNormals1(3,ix)*v234 +
			    _pyramNormals1(4,ix)*v043 +
			    _pyramNormals2(1,ix)*v015 +
			    _pyramNormals2(2,ix)*v125 +
			    _pyramNormals2(3,ix)*v235 +
			    _pyramNormals2(4,ix)*v053);
    }
  }

 //  /// just for testing .....
/*  static RealVector n1(3);
  CFreal maxArea = 0.;
  CFreal sumNormalsU0 = 0.;
  CFreal sumNormalsU1 = 0.;
  
  for (CFuint i = 1; i < 5; ++i) {
    _pyramNormals1.putRow(i,n1); 
    const CFreal normN = n1.norm2();
    maxArea = max(normN, maxArea);  
    sumNormalsU0 += normN;
  }

  for (CFuint i = 1; i < 5; ++i) {
    _pyramNormals2.putRow(i,n1); 
    const CFreal normN = n1.norm2();
    maxArea = max(normN, maxArea);
    sumNormalsU1 += normN;
  }

 
  _maxCVFaceArea = maxArea;
  
  _refCVArea[0] = sumNormalsU0/3.;
  _refCVArea[1] = sumNormalsU1/3.;*/
}


// //////////////////////////////////////////////////////////////////////////////

void DiamondVolume3DDerivative::compute2PyramidNormals()
{
  GeometricEntity& geo = *getMethodData().getCurrentFace();
  
  // construct the 2 pyramid volume and the normals around the face
  Node& node0 = *geo.getNode(0);
  Node& node1 = *geo.getNode(1);
  Node& node2 = *geo.getNode(2);
  Node& node3 = *geo.getNode(3);
  Node& node4 = geo.getState(0)->getCoordinates();
  Node& node5 = geo.getState(1)->getCoordinates();

  // set the coordinates of the pyramid 01234 in a matrix
  _pyramCoord.setRow(node0,0);
  _pyramCoord.setRow(node1,1);
  _pyramCoord.setRow(node2,2);
  _pyramCoord.setRow(node3,3);
  _pyramCoord.setRow(node4,4);
  
  CFreal factor = (!_changeNormalSign[geo.getID()]) ? 1.0 : -1.0;
  
  // compute the face normals (here avoid to compute
  // the normal on the quad face !!!)
  _normalCalculator.computePyramNormals(_pyramCoord,_pyramNormals1);
  _pyramNormals1 *= factor;
  
  // set the coordinates of the pyramid 01235 in a matrix
  _pyramCoord.setRow(node5,4);
  
  // compute the face normals (here avoid to compute
  // the normal on the quad face !!!)
  _normalCalculator.computePyramNormals(_pyramCoord,_pyramNormals2);
  
  // the factor is surely the opposite of the previous one
  factor *= -1.;
  _pyramNormals2 *= factor;

  // on the boundary it should never happen that the right (ghost) state is inside the domain
  //  if (geo.getState(1)->isGhost()) {
    // check if the left (4) and ghost (5) states lie on two sides of the base face center (F)
    // by checking the dot product of the two vectors (4-F, 5-F)
  //  _cFace = 0.25*(node0 + node1 + node2 + node3);
  //  _n4cF = node4 - _cFace; 
  //  _n5cF = node5 - _cFace;
  //  if(MathFunctions::innerProd(_n4cF,_n5cF) > 0.0) {
  //    cout << "DiamondVolume3DDerivative::compute2PyramidNormals() => L and R states are on the same side of the face"<< endl;
      // abort();
   // }
 // }
}
      
//////////////////////////////////////////////////////////////////////////////

void DiamondVolume3DDerivative::computeControlVolume
(vector<RealVector*>& states, GeometricEntity *const geo)
{
  DataHandle<RealVector> nodalStates = this->socket_nstates.getDataHandle();
  
  states[0] = &nodalStates[geo->getNode(0)->getLocalID()];
  states[1] = &nodalStates[geo->getNode(1)->getLocalID()];
  states[2] = &nodalStates[geo->getNode(2)->getLocalID()];
  
  if (geo->nbNodes() == 3) {
    states[3] = geo->getState(0);
    states[4] = geo->getState(1);
    
    //    compute2TetraVolume(geo);
  }
  else if (geo->nbNodes() == 4) {
    states[3] = &nodalStates[geo->getNode(3)->getLocalID()];
    states[4] = geo->getState(0);
    states[5] = geo->getState(1);
    
    //     compute2PyramidVolume(geo);
    compute2PyramidNormals();
  }
  
  cf_assert(geo->getID() < _cvolumes.size());
  _volume = _cvolumes[geo->getID()];
}

//////////////////////////////////////////////////////////////////////////////
      
void DiamondVolume3DDerivative::compute2TetraVolume
(GeometricEntity *const geo)
{
  // construct the 2 tetras volume and the normals around the face
  Node& node0 = *geo->getNode(0);
  Node& node1 = *geo->getNode(1);
  Node& node2 = *geo->getNode(2);
  Node& node3 = geo->getState(0)->getCoordinates();
  Node& node4 = geo->getState(1)->getCoordinates();
  
  // set the coordinates of the tetra 0132 in a matrix
  _tetraCoord.setRow(node0,0);
  _tetraCoord.setRow(node1,1);
  _tetraCoord.setRow(node2,2);
  _tetraCoord.setRow(node3,3);
  
  CFreal volume0123 = _volumeCalculator.calculateTetraVolume(_tetraCoord);
  if (volume0123 < 0.0) {
    _changeNormalSign[geo->getID()] = true;
    volume0123 = std::abs(volume0123);
  }
  
  // set the coordinates of the tetra 0124 in a matrix
  _tetraCoord.setRow(node0,0);
  _tetraCoord.setRow(node1,1);
  _tetraCoord.setRow(node2,2);
  _tetraCoord.setRow(node4,3);
  
  CFreal volume0124 = _volumeCalculator.calculateTetraVolume(_tetraCoord);
  if (volume0124 < 0.0) {
    volume0124 = std::abs(volume0124);
  }
  
  _volume = volume0123 + volume0124;
  cf_assert(_volume > 0.0);
}

//////////////////////////////////////////////////////////////////////////////

void DiamondVolume3DDerivative::compute2PyramidVolume
(GeometricEntity *const geo)
{
  // planarize face (some nodes might be moved)
 // planarizeFaceStart(*geo);
  
  // construct the 2 tetras volume and the normals around the face
  Node& node0 = *geo->getNode(0);
  Node& node1 = *geo->getNode(1);
  Node& node2 = *geo->getNode(2);
  Node& node3 = *geo->getNode(3);
  Node& node4 = geo->getState(0)->getCoordinates();
  Node& node5 = geo->getState(1)->getCoordinates();
  
  //cout << "=> " << node4 << ", " << node5 <<  endl;
  
  // set the coordinates of the pyramid 01234 in a matrix
  _pyramCoord.setRow(node0,0);
  _pyramCoord.setRow(node1,1);
  _pyramCoord.setRow(node2,2);
  _pyramCoord.setRow(node3,3);
  _pyramCoord.setRow(node4,4);
  
  CFreal volume01234 = _volumeCalculator.calculatePyramVolume(_pyramCoord);
  
  CFreal prodVol = volume01234;
  if (volume01234 < 0.0) {
    _changeNormalSign[geo->getID()] = true;
    volume01234 = std::abs(volume01234);
  }
  // set the coordinates of the pyramid 01235 in a matrix
  _pyramCoord.setRow(node5,4);
  
  CFreal volume01235 = _volumeCalculator.calculatePyramVolume(_pyramCoord);
  prodVol *=  volume01235;
  if (volume01235 < 0.0) {
    volume01235 = std::abs(volume01235);
  }
  
  ofstream file("centroid.err", ios::app);
  
  if (prodVol > 0.0) {
    static CFuint count = 0;
    count++;
    
    //    cout << "DiamondVolume3DDerivative::compute2PyramidVolume() => pyramids product volumes (should be < 0) = " << prodVol  << ", " << count << endl;
    if (!geo->getState(1)->isGhost()) {
      file << count <<", INNER => " << geo->getState(0)->getCoordinates() << " ; " << geo->getState(1)->getCoordinates() << endl;
      
      //       RealVector n1(3);
      //       RealVector n2(3);
      //       RealVector nF(3);
      //       nF =  0.25*(node0 + node1 + node2 + node3);
      //       n1 = node4 - nF;
      //       n2 = node5 - nF;
      
      //       cout << "node0 = " << node0 << endl;
      //       cout << "node1 = " << node1 << endl;
      //       cout << "node2 = " << node2 << endl;
      //       cout << "node3 = " << node3 << endl;
      //       cout << "node4 = " << node4 << endl;
      //       cout << "node5 = " << node5 << endl;
      
      //       assert(MathFunctions::innerProd(n1,n2) < 0.0);
      //     abort();   
    } 
    
    if (geo->getState(1)->isGhost()) {
      file << count <<", BOUND => " << geo->getState(0)->getCoordinates() << " ; " << geo->getState(1)->getCoordinates() << endl;
    } 
    
    const CFuint MAXCOUNT = 100;
    if (count < MAXCOUNT) {
      const string filename = "trouble_cell.plt-" + StringOps::to_str(PE::GetPE().GetRank("Default")); 
      // write a TECPLOT file with : face nodes, L/R nodes, face normal (translated to L), left cell nodes  
      ofstream fout(filename.c_str(), ios::app);
      
      if (count == 1) {
	fout << "TITLE = Unstructured grid data" << endl;
	fout << "VARIABLES  =  \"x0\" \"x1\" \"x2\" " << endl;
      }
      
      for (CFuint iCell = 0; iCell < 2; ++iCell) {
	const string pref = (!geo->getState(1)->isGhost()) ? "INNER-" : "BOUND-";
	const string zone = pref + "ZONE" + StringOps::to_str(count-1) + "-" + StringOps::to_str(iCell);
	fout << "ZONE   T=\" "<< zone  <<" Cell\", N=8, E=1, F=FEPOINT, ET=BRICK" << endl;
	GeometricEntity *const cell = geo->getNeighborGeo(iCell);
	const CFuint nbNodes = cell->nbNodes();
	for (CFuint iNode = 0; iNode < nbNodes; ++iNode) {
	  fout.precision(14);  fout.setf(ios::scientific,ios::floatfield); fout << *cell->getNode(iNode) << endl;
	}
	
	if (nbNodes < 8) {
	  for (CFuint iNode = nbNodes; iNode < 8; ++iNode) {
	    fout.precision(14);  fout.setf(ios::scientific,ios::floatfield); fout << *cell->getNode(nbNodes-1) << endl;
	  }
	}
	fout << "1 2 3 4 5 6 7 8" << endl;
      }
      
      if (count == MAXCOUNT-1) {
	fout.close();
      }
    }
  }
  
  // set the volume
  _volume = volume01234 + volume01235; 
  
  // restore face as before planarization
  //planarizeFaceEnd(*geo);
}

//////////////////////////////////////////////////////////////////////////////

void DiamondVolume3DDerivative::computeAverageValues
(GeometricEntity *const geo,
 const vector<RealVector*>& values,
 RealVector& avValues)
{
  const CFuint nbValues = avValues.size();
  if (geo->nbNodes() == 3) {
    for (CFuint i = 0; i < nbValues; ++i) {
      avValues[i] = _third*((*values[0])[i] +
			   (*values[1])[i] +
			   (*values[2])[i]);
    }
  }
  else {
    for (CFuint i = 0; i < nbValues; ++i) {
      avValues[i] = 0.25*((*values[0])[i] +
			  (*values[1])[i] +
			  (*values[2])[i] +
			  (*values[3])[i]);
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

CFuint DiamondVolume3DDerivative::getNbVerticesInControlVolume
(GeometricEntity *const geo) const
{
  return geo->nbNodes() + 2;
}

//////////////////////////////////////////////////////////////////////////////

SafePtr<vector<RealVector> > DiamondVolume3DDerivative::getGradientsJacob()
{
  const CFreal invVolume = _third/_volume;
  
  if (getMethodData().getCurrentFace()->nbNodes() == 3) {
    _gradientsJacob[0] = invVolume*(_n031 + _n132 + _n023);
    _gradientsJacob[1] = invVolume*(_n041 + _n142 + _n024);
  }
  else if (getMethodData().getCurrentFace()->nbNodes() == 4) {
    for (CFuint ix = 0; ix < DIM_3D; ++ix) {
      _gradientsJacob[0][ix] = invVolume*(_pyramNormals1(1,ix) + _pyramNormals1(2,ix) +
					  _pyramNormals1(3,ix) + _pyramNormals1(4,ix));
      _gradientsJacob[1][ix] = invVolume*(_pyramNormals2(1,ix) + _pyramNormals2(2,ix) +
					  _pyramNormals2(3,ix) + _pyramNormals2(4,ix));
    }
  }
  
  return &_gradientsJacob;
}
      
//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > 
DiamondVolume3DDerivative::needsSockets()
{
  std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > result =
   DerivativeComputer::needsSockets();
  result.push_back(&socket_states);
  result.push_back(&socket_nodes);
  result.push_back(&socket_gstates);
  result.push_back(&socket_cellFlag);
  return result;
}
      
//////////////////////////////////////////////////////////////////////////////

void DiamondVolume3DDerivative::planarizeFaceStart(GeometricEntity& geo) 
{ 
  Node& node0 = *geo.getNode(0);
  Node& node1 = *geo.getNode(1);
  Node& node2 = *geo.getNode(2);
  Node& node3 = *geo.getNode(3);
  
  // the mixed product tells us how much the face is non-planar
  const CFreal h1 = MathFunctions::checkPlanarity(node0, node1, node3, node2);
  const CFreal h2 = MathFunctions::checkPlanarity(node1, node2, node0, node3);
  const CFreal h3 = MathFunctions::checkPlanarity(node2, node3, node1, node0);
  const CFreal h4 = MathFunctions::checkPlanarity(node3, node0, node2, node1);
  
  const CFreal maxH = max(h1,max(h2,max(h3,h4)));
  const CFreal minH = min(h1,min(h2,min(h3,h4)));
  
  if (maxH > 1e-14) {
    // choose the plane for which the forth point is less off-plane
    if (h1 <= h2 && h1 <= h3 && h1 <= h4) {
      // project point 2 on the plane defined by 013
      MathFunctions::computeProjectedPoint(node0, node1, node3, node2, _xproj);
      correctNonPlanarFace(geo, 2);  
      cout << "(" << minH << ", " << maxH << "), 1 From " << h1  << " to " << MathFunctions::checkPlanarity(node0, node1, node3, node2) << endl;
    } 
    else if (h2 <= h1 && h2 <= h3 && h2 <= h4) {
      // project point 3 on the plane defined by 120
      MathFunctions::computeProjectedPoint(node1, node2, node0, node3, _xproj);
      correctNonPlanarFace(geo, 3); 
      cout << "(" << minH << ", " << maxH << "), 2 From " << h2  << " to " << MathFunctions::checkPlanarity(node1, node2, node0, node3) << endl;
    } 
    else if (h3 <= h1 && h3 <= h2 && h3 <= h4) {
      // project point 0 on the plane defined by 231
      MathFunctions::computeProjectedPoint(node2, node3, node1, node0, _xproj);
      correctNonPlanarFace(geo, 0); 
      cout << "(" << minH << ", " << maxH << "), 3 From " << h3  << " to " << MathFunctions::checkPlanarity(node2, node3, node1, node0) << endl;
    } 
    else if (h4 <= h1 && h4 <= h2 && h4 <= h3) {
      // project point 1 on the plane defined by 302
      MathFunctions::computeProjectedPoint(node3, node0, node2, node1, _xproj);
      correctNonPlanarFace(geo, 1);
      cout << "(" << minH << ", " << maxH << "), 4 From " << h4  << " to " << MathFunctions::checkPlanarity(node3, node0, node2, node1) << endl;
    }
    
    // recompute the coordinates of the cell centers on the left and on the right of the face (ghost cells are excluded)
    //cout << geo.getState(0)->getCoordinates() << " BEFORE " << geo.getState(1)->getCoordinates() << endl;
    
    //  cout << "HEXA L [0,1,2,3,4,5,6,7]" << endl;
    recomputeCellCenter(geo.getNeighborGeo(0), _lcoord);
    if (!geo.getState(1)->isGhost()) {
      //  cout << "HEXA R [0,1,2,3,4,5,6,7]" << endl;
      recomputeCellCenter(geo.getNeighborGeo(1), _rcoord);
      //cout << geo.getState(0)->getCoordinates() << " AFTER " << geo.getState(1)->getCoordinates() << endl;
      //   cout << MathFunctions::getDistance(geo.getNeighborGeo(0)->getState(0)->getCoordinates(),
      //   		 geo.getNeighborGeo(1)->getState(0)->getCoordinates()) << endl;
    }
  }
}

//////////////////////////////////////////////////////////////////////////////
      
} // namespace FiniteVolume

} // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
