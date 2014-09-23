#include "CorrectedDerivative2D.hh"
#include "Framework/MethodStrategyProvider.hh"
#include "Framework/GeometricEntity.hh"
#include "FiniteVolume/FiniteVolume.hh"


#include "ComputeFaceVertexNeighborsPlusGhost.hh"

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

MethodStrategyProvider<CorrectedDerivative2D,
                       CellCenterFVMData,
		       DerivativeComputer,
                       FiniteVolumeModule>
corrrectedDerivative2DProvider("Corrected2D");

//////////////////////////////////////////////////////////////////////////////

CorrectedDerivative2D::CorrectedDerivative2D(const std::string& name) :
  DerivativeComputer(name),
  socket_gstates("gstates"),
  socket_volumes("volumes"),
  socket_states("states"),
  socket_nodes("nodes"),
  socket_stencil("stencil"),
  socket_uX("uX"),
  socket_uY("uY"),
  _cellsTRS(CFNULL),
  _dr(0.0),
  _eRL(),
  _eRLdotNN(),
  _oEminusNN(),
  _gradLR(),
  _lf1(),
  _lf2(),
  _cellBuilder(),
  _stencil()
{
}

//////////////////////////////////////////////////////////////////////////////

CorrectedDerivative2D::~CorrectedDerivative2D()
{
}

//////////////////////////////////////////////////////////////////////////////

void CorrectedDerivative2D::setup()
{
  DerivativeComputer::setup();

  const CFuint dim = PhysicalModelStack::getActive()->getDim();
  _eRL.resize(dim);
  _eRLdotNN.resize(dim);
  _oEminusNN.resize(dim,dim);
  _gradLR.resize(dim);

  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  _lf1.resize(nbEqs);
  _lf2.resize(nbEqs);

  // set up the cell builder
  _cellBuilder.setup();

  _cellsTRS = MeshDataStack::getActive()->getTrs("InnerCells");
  TrsGeoWithNodesBuilder::GeoData& geoData = _cellBuilder.getDataGE();
  geoData.trs = _cellsTRS;
}

//////////////////////////////////////////////////////////////////////////////

void CorrectedDerivative2D::computeGradients(const RealMatrix& values,
					     vector<RealVector*>& gradients)
{
  if (!getMethodData().getCurrentFace()->getState(1)->isGhost()) {
    computeInnerGradients(values,gradients);
  }
  else {
    computeBoundaryGradients(values,gradients);
  }
}

//////////////////////////////////////////////////////////////////////////////

void CorrectedDerivative2D::computeInnerGradients(const RealMatrix& values,
						  vector<RealVector*>& gradients)
{
  const CFuint nbVars = values.nbRows();
  cf_assert(gradients.size() == nbVars);

  State* const stateL = getMethodData().getCurrentFace()->getState(0);
  State* const stateR = getMethodData().getCurrentFace()->getState(1);
  
//   string file = "outgrad-" + StringOps::to_str(PE::GetPE().GetRank());
//   ofstream fout(file.c_str(), ios::app);
  
 //  if (!stateL->isParUpdatable() && !stateR->isParUpdatable()) {
  //     fout << "stateL and stateR on INNER are not parallel updatable "<< endl;
  //     abort();
  //   }
  
//   if (!stateL->isParUpdatable() && socket_stencil.getDataHandle()[stateL->getLocalID()].size() < 7) {
//     fout << "P" << PE::GetPE().GetRank() <<  " stateL on INNER is not parallel updatable for [" << stateL->getGlobalID() 
// 	 << "], stencil = "  << socket_stencil.getDataHandle()[stateL->getLocalID()].size() << endl;
//     // cf_assert(stateL->isParUpdatable());
//   } 
//   if (!stateR->isParUpdatable() && socket_stencil.getDataHandle()[stateR->getLocalID()].size() < 7) {
//     fout << "P" << PE::GetPE().GetRank() <<  " stateR on INNER is not parallel updatable for [" << stateR->getGlobalID() 
// 	 << "], stencil = "  << socket_stencil.getDataHandle()[stateR->getLocalID()].size() << endl;
//     // cf_assert(stateR->isParUpdatable());
//   }
  
  DataHandle<CFreal> uX = socket_uX.getDataHandle();
  DataHandle<CFreal> uY = socket_uY.getDataHandle();
  
  // if the face is a boundary face, consider the stencil made by
  // the ghost state and the inner cells sharing at least a node
  // with the face*
  for (CFuint i = 0; i < nbVars; ++i) {
    RealVector& grad = *gradients[i];
    cf_assert(grad.size() == 2);
   
    const CFreal dv = (values(i,1) - values(i,0))/_dr;
    const CFuint stateIDR = stateR->getLocalID();
    const CFuint stateIDL = stateL->getLocalID();
    
    _gradLR[XX] = 0.5*(uX(stateIDL,i,nbVars) + uX(stateIDR,i,nbVars));
    _gradLR[YY] = 0.5*(uY(stateIDL,i,nbVars) + uY(stateIDR,i,nbVars));
    correctGradient(dv,_gradLR, grad);
  }
}

//////////////////////////////////////////////////////////////////////////////

void CorrectedDerivative2D::computeBoundaryGradients(const RealMatrix& values,
						     vector<RealVector*>& gradients)
{
  // if the face is a boundary face, consider the stencil made by
  // the ghost state and the inner cells sharing at least a node
  // with the face and apply a LS (least square) algorithm
  
  // left + right + nodal states + neighbor states
  GeometricEntity *const face = getMethodData().getCurrentFace();
  State* const stateL = face->getState(0);
  State* const stateR = face->getState(1);
  
 //  string file = "outgrad-" + StringOps::to_str(PE::GetPE().GetRank());
//   ofstream fout(file.c_str(), ios::app);
//   if (!stateL->isParUpdatable() && socket_stencil.getDataHandle()[stateL->getLocalID()].size() < 7) {
//     fout << "P" << PE::GetPE().GetRank() 
// 	 <<  " stateL (" << stateL->getCoordinates()[XX] << "," << stateL->getCoordinates()[YY] 
// 	 << ") on BOUNDARY is not parallel updatable for [" << stateL->getGlobalID() 
// 	 << "], stencil = "  << socket_stencil.getDataHandle()[stateL->getLocalID()].size() << endl;
//     // cf_assert(stateL->isParUpdatable());
//   } 
  
  
  cf_assert(stateR->isGhost());
  const CFuint nbVars = values.nbRows();
  cf_assert(gradients.size() == nbVars);
  
  // storage of states
  const CFuint startID = 2 + face->nbNodes();
  const CFuint stencilSize = _stencil.size();
  const RealVector& nodeR = stateR->getCoordinates();
  
  //  if (PE::GetPE().GetRank() == 3 &&  stateL->getGlobalID() == 19946) {
  //     cout << "P3, state 19946, stencil L,R = " << socket_stencil.getDataHandle()[stateL->getLocalID()].size() << ", " << _stencil.size() << endl;
  //   }
  
    
//   if (stateL->getGlobalID() == ComputeFaceVertexNeighborsPlusGhost::GLOBAL_STATE_ID.second) {
//     cout << "P" << ComputeFaceVertexNeighborsPlusGhost::GLOBAL_STATE_ID.first  << 
//       " =>=>=> " << ComputeFaceVertexNeighborsPlusGhost::GLOBAL_STATE_ID.second << ", "<< stateL->getLocalID() << " => [ ";
//     for (int i = 0; i < stencilSize; ++i) {
//       cout << _stencil[i]->getLocalID() << " ";
  
//     }
//     cout << "]\n";
//   }
    
  CFreal l11 = 0.0;
  CFreal l12 = 0.0;
  CFreal l22 = 0.0;
  _lf1 = 0.0;
  _lf2 = 0.0;
  
 // static CFuint minStencil = 10000;
 // static CFuint maxStencil = 0;
 // static CFuint minGlobalStateID = 0;
  
  // loop over the neighbor cells belonging to the chosen stencil
  for(CFuint in = 0; in < stencilSize; ++in) {
   // if (stencilSize < minStencil) {
   //   minStencil = min(minStencil, stencilSize);
    //  minGlobalStateID = stateL->getGlobalID();
    //}
    //maxStencil = max(maxStencil, stencilSize);
    
    const RealVector& nodeLast = _stencil[in]->getCoordinates();
    const CFreal deltaR = MathFunctions::getDistance(nodeR,nodeLast);
    const CFreal weight = 1.0/deltaR;
    const CFreal dx = weight*(nodeLast[XX] - nodeR[XX]);
    const CFreal dy = weight*(nodeLast[YY] - nodeR[YY]);
    l11 += dx*dx;
    l12 += dx*dy;
    l22 += dy*dy;
    
    const CFuint nID = startID + in;
    for (CFuint iVar = 0; iVar < nbVars; ++iVar) {
      const CFreal du = weight*(values(iVar,nID) - values(iVar,1));
      _lf1[iVar] += dx*du;
      _lf2[iVar] += dy*du;
    }
  }

  DataHandle<CFreal> uX = socket_uX.getDataHandle();
  DataHandle<CFreal> uY = socket_uY.getDataHandle();
  const CFreal invDet = 1./(l11*l22 - l12*l12);

  for (CFuint i = 0; i < nbVars; ++i) {
    RealVector& grad = *gradients[i];
    const CFreal dv = (values(i,1) - values(i,0))/_dr;
    const CFreal uXR = (l22*_lf1[i] - l12*_lf2[i])*invDet;
    const CFreal uYR = (l11*_lf2[i] - l12*_lf1[i])*invDet;
    const CFuint stateID = stateL->getLocalID();
    _gradLR[XX] = 0.5*(uX(stateID,i,nbVars) + uXR);
    _gradLR[YY] = 0.5*(uY(stateID,i,nbVars) + uYR);
    correctGradient(dv,_gradLR, grad);
  }
  
  //  cout << "Left stateID  = " << stateL->getGlobalID() << ", (min,max) = ( " << minStencil << ", "  << maxStencil << " ) " << endl;
}
      
//////////////////////////////////////////////////////////////////////////////

void CorrectedDerivative2D::computeControlVolume(vector<RealVector*>& states,
						 GeometricEntity *const geo)
{ 
  // set all states needed in the reconstruction stencil
  setControlVolumeStates(states);
  
  State* const stateL = getMethodData().getCurrentFace()->getState(0);
  State* const stateR = getMethodData().getCurrentFace()->getState(1);
  Node& nodeL = stateL->getCoordinates();
  Node& nodeR = stateR->getCoordinates();
  _dr = MathTools::MathFunctions::getDistance(nodeL, nodeR);
  cf_assert(_dr > 0.);
  
  // unit vector joining the two cell centers
  _eRL = (nodeR - nodeL)/_dr;
  cf_assert(MathChecks::isEqualWithError(_eRL.norm2(),1.0,1e-12));
  
  RealVector& n = getMethodData().getUnitNormal();
  cf_assert (socket_isOutward.getDataHandle()[geo->getID()] == (CFint)stateL->getLocalID());  
  
  const CFreal eDotN = MathFunctions::innerProd(_eRL, n);
  if (eDotN <= 0.) {	
    if (PhysicalModelStack::getActive()->getDim() == 2) {
      // only quadrilateral cells can be invalid
      cf_assert(geo->getNeighborGeo(0)->nbNodes() == 4);
      
      // write a TECPLOT file with : face nodes, L/R nodes, face normal (translated to L), left cell nodes  
      ofstream fout("invalid_cell.plt");
      fout << "TITLE = Unstructured grid data" << endl;
      fout << "VARIABLES  =  \"x0\" \"x1\"" << endl;
      fout << "ZONE   T=\"ZONE0 Face\", N=2, E=1, F=FEPOINT, ET=LINESEG" << endl;
      fout.precision(14); fout << *geo->getNode(0) << endl;
      fout.precision(14); fout << *geo->getNode(1) << endl;
      fout << "1 2" << endl;
      fout << "ZONE   T=\"ZONE1 L/R\", N=2, E=1, F=FEPOINT, ET=LINESEG" << endl;
      fout.precision(14); fout << nodeL << endl;
      fout.precision(14); fout << nodeR << endl;
      fout << "1 2" << endl; 
      fout << "ZONE   T=\"ZONE2 Normal\", N=2, E=1, F=FEPOINT, ET=LINESEG" << endl;
      fout.precision(14); fout << nodeL[XX] << " " << nodeL[YY] << endl;
      fout.precision(14); fout << n[XX]+nodeL[XX] << " " << n[YY]+nodeL[YY] << endl;
      fout << "1 2" << endl;
      fout << "ZONE   T=\"ZONE3 Cell\", N=4, E=1, F=FEPOINT, ET=QUADRILATERAL" << endl;
      fout.precision(14); fout << *geo->getNeighborGeo(0)->getNode(0) << endl;
      fout.precision(14); fout << *geo->getNeighborGeo(0)->getNode(1) << endl;
      fout.precision(14); fout << *geo->getNeighborGeo(0)->getNode(2) << endl;
      fout.precision(14); fout << *geo->getNeighborGeo(0)->getNode(3) << endl;
      fout << "1 2 3 4" << endl;
    }
  }
  
//  cf_assert(eDotN > 0.); // the angle between face normal and LR should be < 90 degrees 
  _eRLdotNN = eDotN*n;
  
  computeOEminusNN(n);
  
  DataHandle<CFreal> volumes = socket_volumes.getDataHandle();
  if (!stateR->isGhost()) {
    _volume = 0.5*(volumes[stateL->getLocalID()] + volumes[stateR->getLocalID()]);
  }
  else {
    _volume = volumes[stateL->getLocalID()];
  }
}

//////////////////////////////////////////////////////////////////////////////

void CorrectedDerivative2D::computeAverageValues(GeometricEntity *const geo,
						 const vector<RealVector*>& values,
						 RealVector& avValues)
{
  // average state on the face is based on the average of the
  // nodal values
  const CFuint nbValues = avValues.size();
  const CFuint nbNodesInFace = geo->nbNodes();
  const CFreal invNbNodes = 1./static_cast<CFreal>(nbNodesInFace);
  
  for (CFuint i = 0; i < nbValues; ++i) {
    avValues[i] = 0.0;
    for (CFuint iNode = 0; iNode < nbNodesInFace; ++iNode) {
      avValues[i] += invNbNodes*(*values[iNode + 2])[i];
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void CorrectedDerivative2D::setControlVolumeStates(vector<RealVector*>& states)
{ 
  // left + right + nodal states + neighbor states used in the least 
  // square reconstruction
  GeometricEntity *const face = getMethodData().getCurrentFace();
  State* const stateL = face->getState(0);
  State* const stateR = face->getState(1);
  const CFuint nbNodesInFace = face->nbNodes();
  
  // first state is the left state
  states[0] = stateL;
  // second state is the right state
  states[1] = stateR;
  
  // nodal states for the given face (needed for computing the average state) 
  // start in position 2
  CFuint counter = 2;
  DataHandle<RealVector> nstates = socket_nstates.getDataHandle();
  for (CFuint iNode = 0; iNode < nbNodesInFace; ++iNode, counter++) {
    states[counter] = &nstates[face->getNode(iNode)->getLocalID()];
  } 
  
  cf_assert((counter == 4 && PhysicalModelStack::getActive()->getDim() == 2) ||
	    ((counter == 5 || counter == 6) && PhysicalModelStack::getActive()->getDim() == 3));
  
  // if the face is on the boundary, build a local stencil to estimate the gradients
  if (stateR->isGhost()) {
    _stencil.clear();
    cf_assert(_stencil.size() == 0);
    
    //    // alternative algorithm: get all the vertex neighbors of the face nodes (on the boundary)
    //     static set<State*> setStates;
    //     const CFuint nbBNodes = face->nbNodes();
    //     for (CFuint iNode = 0; iNode < nbBNodes; ++iNode) {
    //       const CFuint nodeID = face->getNode(iNode)->getLocalID();
    //       const vector<State*>& s = getMethodData().
    // 	getNodalStatesExtrapolator()->getNodalStateNeighbors(nodeID);
    //       for (CFuint is = 0; is < s.size(); ++is) {
    // 	setStates.insert(s[is]);
    //       }
    //     }
    
    //     set<State*>::iterator itr;
    //     for (itr = setStates.begin(); itr != setStates.end(); ++itr) {
    //       if ((*itr != stateL) && (*itr != stateR)) {
    // 	_stencil.push_back(*itr);
    // 	states[counter++] = *itr;
    //       }
    //     }
    //     setStates.clear();
    
    DataHandle<vector<State*> > stencilSocket = socket_stencil.getDataHandle();
    const CFuint stateIDL = stateL->getLocalID();
    const CFuint stencilSize = stencilSocket[stateIDL].size();
    _stencil.push_back(stateL);
    
    for(CFuint in = 0; in < stencilSize; ++in) {
      State* const last = stencilSocket[stateIDL][in];
      
      // if the neighbor state of the current left state is a ghost state 
      // and is different from the right state, include it in the 
      // reconstruction stencil
      if (last->isGhost()) { 
	if (last != stateR) {
	  states[counter++] = last;
	  _stencil.push_back(last);
	}
      }
      else {
	// check if the corresponding internal neighbor cell shares at least 
	// one node with current face here a simple TrsGeoWithNodesBuilder 
	// can be used to save time
	_cellBuilder.getDataGE().idx = last->getLocalID();
	GeometricEntity *const cell = _cellBuilder.buildGE();
	const CFuint nbFaceNodes = face->nbNodes();
	for (CFuint iNode = 0; iNode < nbFaceNodes; ++iNode) {
	  if (cell->containNode(face->getNode(iNode))) {
	    states[counter++] = last;
	    _stencil.push_back(last);
	    break;
	  }
	}
	
	_cellBuilder.releaseGE();
      }
    }   
  }
  
//   string fname = "outfile-" + StringOps::to_str(ComputeFaceVertexNeighborsPlusGhost::GLOBAL_STATE_ID.first);
//   ofstream fout(fname.c_str(), ios::app);
//   DataHandle<vector<State*> > stencilSocket = socket_stencil.getDataHandle();
//   const CFuint stateIDL = stateL->getLocalID();
//   const CFuint stencilSizeL = stencilSocket[stateIDL].size();
//   if (stateL->getGlobalID() == ComputeFaceVertexNeighborsPlusGhost::GLOBAL_STATE_ID.second) {
//     fout << "P" << ComputeFaceVertexNeighborsPlusGhost::GLOBAL_STATE_ID.first  << 
//       " =>=>=> " << ComputeFaceVertexNeighborsPlusGhost::GLOBAL_STATE_ID.second << ", "<< stateL->getLocalID() << " => [ ";
//     for (int i = 0; i < stencilSizeL; ++i) {
//       fout << stencilSocket(stateIDL,i,nbVars)->getLocalID() << " ";
//    }
//     fout << "]\n";
//   } 
//   else if (!stateR->isGhost()) {
//     if (stateR->getGlobalID() == ComputeFaceVertexNeighborsPlusGhost::GLOBAL_STATE_ID.second) {
//       const CFuint stateIDR = stateR->getLocalID();
//       const CFuint stencilSizeR = stencilSocket[stateIDR].size();
//       fout << "RIGHT!!!!!!" << endl;
//       fout << "P" << ComputeFaceVertexNeighborsPlusGhost::GLOBAL_STATE_ID.first  << 
// 	" =>=>=> " << ComputeFaceVertexNeighborsPlusGhost::GLOBAL_STATE_ID.second << ", "<< stateR->getLocalID() << " => [ ";
//       for (int i = 0; i < stencilSizeR; ++i) {
// 	fout << stencilSocket(stateIDR,i,nbVars)->getLocalID() << " ";
//       }
//       fout << "]\n";
//     }
//   }
    
  _nbDofsInCV = counter;
  assert(_nbDofsInCV >= 4);
}

//////////////////////////////////////////////////////////////////////////////

SafePtr<vector<RealVector> > CorrectedDerivative2D::getGradientsJacob()
{
  /// @TODO to be implemented properly
  const CFreal ovDr = 1./_dr;
  _gradientsJacob[0] = (-ovDr)*_eRLdotNN;
  _gradientsJacob[1] = ovDr*_eRLdotNN;
  
  return &_gradientsJacob;
}

//////////////////////////////////////////////////////////////////////////////

vector<SafePtr<BaseDataSocketSink> > CorrectedDerivative2D::needsSockets()
{
  vector<SafePtr<BaseDataSocketSink> > result =
    DerivativeComputer::needsSockets();

  result.push_back(&socket_gstates);
  result.push_back(&socket_volumes);
  result.push_back(&socket_states);
  result.push_back(&socket_nodes);
  result.push_back(&socket_stencil);
  result.push_back(&socket_uX);
  result.push_back(&socket_uY);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

void CorrectedDerivative2D::computeOEminusNN(const RealVector& n)
{
  const CFreal nx = n[XX];
  const CFreal ny = n[YY];
  _oEminusNN(0,0) = 1.0 - nx*nx;
  _oEminusNN(0,1) = _oEminusNN(1,0) = - nx*ny;
  _oEminusNN(1,1) = 1.0 - ny*ny;
}

//////////////////////////////////////////////////////////////////////////////

} // namespace FiniteVolume

} // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
