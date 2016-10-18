#include "Framework/ComputeDummyStates.hh"
#include "Framework/MeshData.hh"
#include "Framework/PhysicalModel.hh"
#include "Framework/State.hh"
#include "Framework/NormalsCalculator.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

ComputeDummyStates::ComputeDummyStates() :
  socket_normals("Null"),
  socket_gstates("Null"),
  socket_states("Null"),
  socket_nodes("Null"),
  _coord(), 
  _faceMidPoint(),
  _updating(false),
  _isGhostOnFace(false)
{
}

//////////////////////////////////////////////////////////////////////////////

ComputeDummyStates::~ComputeDummyStates()
{
}

//////////////////////////////////////////////////////////////////////////////

void ComputeDummyStates::setup()
{
  _coord.resize(PhysicalModelStack::getActive()->getDim());
  _faceMidPoint.resize(PhysicalModelStack::getActive()->getDim());
}

//////////////////////////////////////////////////////////////////////////////

void ComputeDummyStates::setDataSockets(DataSocketSink< CFreal> normals,
                                        DataSocketSink<State*> gstates,
                                        DataSocketSink<State*, GLOBAL> states,
                                        DataSocketSink<Node*, GLOBAL> nodes)
{
  socket_normals = normals;
  socket_gstates = gstates;
  socket_states = states;
  socket_nodes = nodes;
}

//////////////////////////////////////////////////////////////////////////////

void ComputeDummyStates::operator() (const std::vector<std::string>& trssWithGhostOnFace)
{
  CFLog(VERBOSE, "ComputeDummyStates::operator() => START\n");
  
  _updating = false;
  
  // normals to the boundary are ALWAYS computed outward, by construction
  DataHandle<CFreal> normals = socket_normals.getDataHandle();
  DataHandle<State*> gstates = socket_gstates.getDataHandle();
  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();
  DataHandle < Framework::Node*, Framework::GLOBAL > nodes = socket_nodes.getDataHandle();
  vector< Common::SafePtr<TopologicalRegionSet> > alltrs = MeshDataStack::getActive()->getTrsList();
  RealVector newState(PhysicalModelStack::getActive()->getNbEq());
  const CFuint dim = PhysicalModelStack::getActive()->getDim();
  RealVector faceNormal(dim);
  vector<Node*> fnodes;
  
  const bool isGhost = true;
  CFuint nbGhostStates = 0;
  // loop over all the BOUNDARY  topological region sets
  vector< Common::SafePtr<TopologicalRegionSet> >::iterator itrs;
  for (itrs = alltrs.begin(); itrs != alltrs.end(); ++itrs) {

    Common::SafePtr<TopologicalRegionSet> trs = (*itrs);
    if (trs->getName() != "InnerCells" &&
        trs->getName() != "InnerFaces") {
      
      _isGhostOnFace = isGhostOnFace(trssWithGhostOnFace, trs->getName());
      
      const CFuint nbGeos = trs->getLocalNbGeoEnts();
      // loop over all the boundary faces
      for (CFuint iGeo = 0; iGeo < nbGeos; ++iGeo, ++nbGhostStates) {
        // create a ghost state
        State* ghostState = new State(newState, isGhost);
	const State* const state = states[trs->getStateID(iGeo, 0)];
        const CFuint faceID = trs->getLocalGeoID(iGeo);
        const CFuint startID = faceID*dim;
	
        // set the current normal
        for (CFuint i = 0; i < dim; ++i) {
          faceNormal[i] = normals[startID + i];
        }
	
	// get all nodes of the current face
       	const CFuint nbFaceNodes = trs->getNbNodesInGeo(iGeo);
	fnodes.resize(nbFaceNodes);
	for (CFuint f = 0; f < nbFaceNodes; ++f) {
	  fnodes[f] = nodes[trs->getNodeID(iGeo, f)];
	}
	
        // set the coordinates in the ghost node
        setCoordinatesInState(fnodes,
                              *state,
                              faceNormal,
                              *ghostState);
	
	// override the current normal 
        for (CFuint i = 0; i < dim; ++i) { 
	  normals[startID + i] = faceNormal[i]; 
        } 
	
        // set the local ID corresponding to this ghost state
        trs->setStateID(iGeo, 1, nbGhostStates);
        ghostState->setLocalID(nbGhostStates);
        // set the new state in the state list that will keep ownership on them
        gstates[nbGhostStates] = ghostState;
      }
    }
  } 
  
  CFLog(VERBOSE, "ComputeDummyStates::operator() => END\n");
}

//////////////////////////////////////////////////////////////////////////////

void ComputeDummyStates::updateAllDummyStates()
{
  CFLog(VERBOSE, "ComputeDummyStates::updateAllDummyStates() => START\n");
  
  _updating = true;

  // normals to the boundary are ALWAYS computed outward, by construction
  DataHandle<CFreal> normals = socket_normals.getDataHandle();
  DataHandle<State*> gstates = socket_gstates.getDataHandle();
  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();
  DataHandle < Framework::Node*, Framework::GLOBAL > nodes = socket_nodes.getDataHandle();

  RealVector newState(PhysicalModelStack::getActive()->getNbEq());
  const CFuint dim = PhysicalModelStack::getActive()->getDim();
  RealVector faceNormal(dim);

  // unused //  const bool isGhost = true;
  CFuint nbGhostStates = 0;
  // loop over all the BOUNDARY  topological region sets
  vector< Common::SafePtr<TopologicalRegionSet> > alltrs = MeshDataStack::getActive()->getTrsList();
  vector<Node*> fnodes;
  
  vector< Common::SafePtr<TopologicalRegionSet> >::iterator itrs;
  for (itrs = alltrs.begin(); itrs != alltrs.end(); ++itrs) {

    SafePtr<TopologicalRegionSet> trs = (*itrs);
    if (trs->getName() != "InnerCells" &&
        trs->getName() != "InnerFaces") {

      const CFuint nbGeos = trs->getLocalNbGeoEnts();
      // loop over all the boundary faces
      for (CFuint iGeo = 0; iGeo < nbGeos; ++iGeo, ++nbGhostStates) {
        // create a ghost state
        State* ghostState = gstates[nbGhostStates];
	const State* const state = states[trs->getStateID(iGeo, 0)];
        const CFuint faceID = trs->getLocalGeoID(iGeo);
        const CFuint startID = faceID*dim;
        // set the current normal
        for (CFuint i = 0; i < dim; ++i) {
          faceNormal[i] = normals[startID + i];
        }
	
	// get all nodes of the current face
       	const CFuint nbFaceNodes = trs->getNbNodesInGeo(iGeo);
	fnodes.resize(nbFaceNodes);
	for (CFuint f = 0; f < nbFaceNodes; ++f) {
	  fnodes[f] = nodes[trs->getNodeID(iGeo, f)];
	}
	
        // update the coordinates in the ghost node
        setCoordinatesInState(fnodes,
                              *state,
                              faceNormal,
                              *ghostState);
	
	// override the current normal 
        for (CFuint i = 0; i < dim; ++i) {  
          normals[startID + i] = faceNormal[i];  
        }  
      }
    }
  }
  
  CFLog(VERBOSE, "ComputeDummyStates::updateAllDummyStates() => END\n");
}

//////////////////////////////////////////////////////////////////////////////

void ComputeDummyStates::setCoordinatesInState(const std::vector<Node*>& nodes,
					       const State& state,
					       RealVector& normal,
					       State& ghostState)
{
  if (nodes.size() == 4) {    
    // consider the average of all normals to the 2 or 4 subtriangles in which the 
    // quadrilateral face can be decomposed
    
    static NormalsCalculator nc;
    static RealVector tnormal0(3);
    static RealVector tnormal1(3);
    static RealVector tnormal2(3);
    static RealVector tnormal3(3);
    static RealMatrix tcoord(3,3);
    static RealVector newNormal(3);
    CFreal a0 =0.0;
    CFreal a1 =0.0;
    CFreal a2 =0.0;
    CFreal a3 =0.0;
    
    // here unit normals are computed
    
    tcoord.setRow(*nodes[0], 0);
    tcoord.setRow(*nodes[1], 1);
    tcoord.setRow(*nodes[3], 2);
    nc.computeTriagNormal(tnormal0,tcoord);
    a0 = tnormal0.norm2();
    tnormal0 /= a0;

    tcoord.setRow(*nodes[1], 0);
    tcoord.setRow(*nodes[2], 1);
    tcoord.setRow(*nodes[3], 2);
    nc.computeTriagNormal(tnormal1,tcoord);
    a1 = tnormal1.norm2();
    tnormal1 /= a1;

    tcoord.setRow(*nodes[0], 0);
    tcoord.setRow(*nodes[1], 1);
    tcoord.setRow(*nodes[2], 2);
    nc.computeTriagNormal(tnormal2,tcoord);
    a2 = tnormal2.norm2();
    tnormal2 /= a2;

    tcoord.setRow(*nodes[0], 0);
    tcoord.setRow(*nodes[2], 1);
    tcoord.setRow(*nodes[3], 2);
    nc.computeTriagNormal(tnormal3,tcoord);
    a3 = tnormal3.norm2();
    tnormal3 /= a3;
    
    // average face unit normal
    newNormal = 0.25*(tnormal0 + tnormal1 + tnormal2 + tnormal3);
    // newNormal = 0.5*(tnormal0 + tnormal1);
  
    // here we assign the original area of the face
    newNormal *= normal.norm2();
    
    // average face normal 
    //    newNormal = (2.0*0.25)*(tnormal0 + tnormal1 + tnormal2 + tnormal3); 
    //    newNormal = (2.0*0.5)*(tnormal0 + tnormal1); 
    
    // New and old normal must be both outward with respect to the computational domain
    if (MathFunctions::innerProd(normal,newNormal) < 0.) {
      normal = (-1.)*newNormal;
    }
    else {
      normal = newNormal;
    }
  }
  
  // The equation of the plane containing the boundary face
  // and the given node (xp, yp, zp) (first node of the face),
  // with normal (a,b,c) is
  // a*x + b*y + c*z + k = 0
  // with k = -a*xp - b*yp - c*zp
  
  CFLogDebugMax( "normal = " << normal << "\n");
  CFLogDebugMax( "state coord = " << state.getCoordinates() << "\n");
  
  _faceMidPoint = 0.;
  for (CFuint i = 0; i < nodes.size(); ++i) {
    _faceMidPoint += *nodes[i];
  }
  _faceMidPoint /= (CFreal) nodes.size();
  
  CFreal k = 0.0;
  // if (nodes.size() < 4) {
  //k = - MathFunctions::innerProd(normal, *nodes[0]);
  k = - MathFunctions::innerProd(normal, _faceMidPoint);  
  
  // }
  // else if (nodes.size() == 4) {
  // if the face has 4 nodes, they could not lie all on the same plane
  // so we take an average
  //const CFreal p0 = - MathFunctions::innerProd(normal, *nodes[0]);
  // const CFreal p1 = - MathFunctions::innerProd(normal, *nodes[1]);
  //  const CFreal p2 = - MathFunctions::innerProd(normal, *nodes[2]);
  //  const CFreal p3 = - MathFunctions::innerProd(normal, *nodes[3]);
  //  k = 0.25*(p0 + p1 + p2 + p3);
  //}
  
  // t is parameter for vectorial representation of a straight line
  // in space
  // t = (a*xM + b*yM + c*zM + k)/(a*a + b*b + c*c)
  
  const CFreal n2 = MathFunctions::innerProd(normal, normal);
  cf_assert(std::abs(n2) > 0.0);
  
  const RealVector& stateCoord = state.getCoordinates();
  const CFreal t = (MathFunctions::innerProd(normal,stateCoord) + k)/n2;
    
  // The point N symmetric to M (position of the other state,
  // internal neighbor of the face) with respect to the
  // given plane is given by
  // (xN, yN, zN) = (xM, yM, zM) - 2*t*(a,b,c)
  const CFreal fac = (_isGhostOnFace) ? 1. : 2.;
  _coord =  stateCoord - fac*t*normal;
  /*  if (nodes.size() == 4) {
    // we check if the candidate ghost point G is on the other side with respect to the face mid point F
    // and the inner state I 
    // consider the two vectors (G-F) and (I-F): if their dot product is > 0 they are on the same 
    // side of the face 

    static RealVector ghostF(3);
    static RealVector innerF(3);
    static RealVector midF(3);
  
    midF = 0.25*((*nodes[0]) + (*nodes[1]) + (*nodes[2]) + (*nodes[3]));
    ghostF = _coord - midF ;
    innerF = stateCoord - midF;
    
    static CFuint count = 0; 
    if (MathFunctions::innerProd(ghostF,innerF) > 0.) {
      cout << "ComputeDummyStates => ghost wrong !! "<< endl; cout << count++ << endl;
      
      // compute ghost point as the symmetric point of the internal state with respect to the face centroid 
      // xG = xI + 2 * (xF - xI) = xI - 2 * (xI - xF) 
      _coord = stateCoord - 2.0*innerF;
    }
    }*/
  // the ghost node is not on the mesh
  bool isOnMesh = false;
  if(!_updating){
    Node* coord = new Node(_coord, isOnMesh);
    ghostState.setSpaceCoordinates(coord);
  }
  else{
    Node* tmpNode = ghostState.getNodePtr();
    (*tmpNode) = _coord;
    ghostState.setSpaceCoordinates(tmpNode);
  }
  
  CFLogDebugMax( "G coord = " << _coord << "\n" << "\n");
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
