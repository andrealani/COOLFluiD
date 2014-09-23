#include "CorrectedDerivativeGG2D.hh"
#include "Framework/MethodStrategyProvider.hh"
#include "Framework/GeometricEntity.hh"
#include "Framework/LocalConnectionData.hh"
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

MethodStrategyProvider<CorrectedDerivativeGG2D,
                       CellCenterFVMData,
		       DerivativeComputer,
                       FiniteVolumeModule>
corrrectedDerivativeGG2DProvider("CorrectedGG2D");

//////////////////////////////////////////////////////////////////////////////

CorrectedDerivativeGG2D::CorrectedDerivativeGG2D(const std::string& name) :
  DerivativeComputer(name),
  socket_normals("normals"),
  socket_gstates("gstates"),
  socket_volumes("volumes"),
  socket_states("states"),
  socket_nodes("nodes"),
  _cellsTRS(CFNULL),
  _dr(0.0),
  _eRL(),
  _eRLdotNN(),
  _oEminusNN(),
  _gradLR(),
  _lf1(),
  _lf2(),
  _uxL(),
  _uxR(),
  _fnormal(),
  _cellBuilder(),
  _stencil(),
  _faceNodesL(CFNULL),
  _faceNodesR(CFNULL),
  _cellFaceIDsL(),
  _cellFaceIDsR(),
  _nbNodesInCellL(0),
  _nbNodesInCellR(0)
{
  addConfigOptionsTo(this);
  _useWeights = false;
  setParameter("useWeights",&_useWeights);
}

//////////////////////////////////////////////////////////////////////////////

CorrectedDerivativeGG2D::~CorrectedDerivativeGG2D()
{
}

//////////////////////////////////////////////////////////////////////////////

void CorrectedDerivativeGG2D::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< bool >
    ("useWeights", "Tell if to use weights in the LS algo");
}

//////////////////////////////////////////////////////////////////////////////

void CorrectedDerivativeGG2D::setup()
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
  _uxL.resize(dim,nbEqs);
  _uxR.resize(dim,nbEqs);
  _fnormal.resize(dim);

  // set up the cell builder
  _cellBuilder.getGeoBuilder()->setDataSockets
    (socket_states, socket_gstates, socket_nodes);
  _cellBuilder.setup();

  _cellsTRS = MeshDataStack::getActive()->getTrs("InnerCells");
  CellTrsGeoBuilder::GeoData& geoData = _cellBuilder.getDataGE();
  geoData.trs = _cellsTRS;
}

//////////////////////////////////////////////////////////////////////////////

void CorrectedDerivativeGG2D::computeGradients(const RealMatrix& values,
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

void CorrectedDerivativeGG2D::computeInnerGradients(const RealMatrix& values,
						    vector<RealVector*>& gradients)
{
  const CFuint nbVars = values.nbRows();
  cf_assert(gradients.size() == nbVars);
  State* const stateL = getMethodData().getCurrentFace()->getState(0);
  State* const stateR = getMethodData().getCurrentFace()->getState(1);
  
  CFuint counter = 2;
  _uxL = 0.0;
  computeGradientsGG(_faceNodesL, _cellFaceIDsL, stateL, values, counter, _uxL);

  // increment the counter
  counter += _nbNodesInCellL;

  _uxR = 0.0;
  computeGradientsGG(_faceNodesR, _cellFaceIDsR, stateR, values, counter, _uxR);

  const CFuint dim = PhysicalModelStack::getActive()->getDim();
  for (CFuint i = 0; i < nbVars; ++i) {
    RealVector& grad = *gradients[i];
    const CFreal dv = (values(i,1) - values(i,0))/_dr;
    for (CFuint iDim = 0; iDim < dim; ++iDim) {
      _gradLR[iDim] = 0.5*(_uxL(iDim,i) + _uxR(iDim,i));
    }
    correctGradient(dv,_gradLR, grad);
  }
}

//////////////////////////////////////////////////////////////////////////////

void CorrectedDerivativeGG2D::computeBoundaryGradients(const RealMatrix& values,
						     vector<RealVector*>& gradients)
{
  // if the face is a boundary face, consider the stencil made by
  // the ghost state and the inner cells sharing at least a node
  // with the face and apply a LS (least square) algorithm

  // left + right + nodal states + neighbor states
  const CFuint nbVars = values.nbRows();
  cf_assert(gradients.size() == nbVars);
//   const CFuint dim = PhysicalModelStack::getActive()->getDim();
  State* const stateL = getMethodData().getCurrentFace()->getState(0);
  State* const stateR = getMethodData().getCurrentFace()->getState(1);
  DataHandle<CFreal> volumes = socket_volumes.getDataHandle();
  DataHandle<CFreal> normals = socket_normals.getDataHandle();
  DataHandle<CFint> isOutward = socket_isOutward.getDataHandle();

  // storage of states
  const CFuint startID = 2 + _nbNodesInCellL;
  const CFuint stencilSize = _stencil.size();
  const RealVector& nodeR = stateR->getCoordinates();

  CFreal l11 = 0.0;
  CFreal l12 = 0.0;
  CFreal l22 = 0.0;
  _lf1 = 0.0;
  _lf2 = 0.0;

  // loop over the neighbor cells belonging to the chosen stencil
  for (CFuint in = 0; in < stencilSize; ++in) {
    const RealVector& nodeLast = _stencil[in]->getCoordinates();
    const CFreal weight = (!_useWeights) ? 1.0 :
      1.0/MathFunctions::getDistance(nodeR,nodeLast);
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

  const CFreal invDet = 1./(l11*l22 - l12*l12);

  CFuint counter = 2;
  _uxL = 0.0;
  computeGradientsGG(_faceNodesL, _cellFaceIDsL, stateL, values, counter, _uxL);

  for (CFuint i = 0; i < nbVars; ++i) {
    RealVector& grad = *gradients[i];
    const CFreal dv = (values(i,1) - values(i,0))/_dr;
    const CFreal uXR = (l22*_lf1[i] - l12*_lf2[i])*invDet;
    const CFreal uYR = (l11*_lf2[i] - l12*_lf1[i])*invDet;

    _gradLR[XX] = 0.5*(_uxL(XX,i) + uXR);
    _gradLR[YY] = 0.5*(_uxL(YY,i) + uYR);

    correctGradient(dv,_gradLR, grad);
  }
}

//////////////////////////////////////////////////////////////////////////////

void CorrectedDerivativeGG2D::computeControlVolume(vector<RealVector*>& states,
						   GeometricEntity *const geo)
{
  setControlVolumeStates(states);
  
  State* const stateL = getMethodData().getCurrentFace()->getState(0);
  State* const stateR = getMethodData().getCurrentFace()->getState(1);
  Node& node0 = stateL->getCoordinates();
  Node& node1 = stateR->getCoordinates();
  _dr = MathTools::MathFunctions::getDistance(node0, node1);

  const CFuint dim = PhysicalModelStack::getActive()->getDim();
  for (CFuint ix = 0; ix < dim; ++ix) {
    _eRL[ix] =  node1[ix] - node0[ix];
  }
  _eRL /= _dr;
  
  RealVector& n = getMethodData().getUnitNormal();
  const CFreal eDotN = MathFunctions::innerProd(_eRL, n);
  _eRLdotNN = eDotN*n;

  computeOEminusNN(n);

  // check this !!!!!!!!!

  //  // handle to ID of the cell for which the normal is outward
  //   DataHandle<CFint> isOutward = socket_isOutward.getDataHandle();
  //   const CFuint faceID = geo->getID();
  //   if (static_cast<CFuint>(isOutward[faceID]) != geo->getState(0)->getLocalID()) {
  //     _eRL *= -1.;
  //   }

  DataHandle<CFreal> volumes = socket_volumes.getDataHandle();
  if (!stateR->isGhost()) {
    _volume = 0.5*(volumes[stateL->getLocalID()] + volumes[stateR->getLocalID()]);
  }
  else {
    _volume = volumes[stateL->getLocalID()];
  }
}

//////////////////////////////////////////////////////////////////////////////

void CorrectedDerivativeGG2D::computeAverageValues
(GeometricEntity *const geo,
 const vector<RealVector*>& values,
 RealVector& avValues)
{
  // average state on the face is based on the average of the
  // nodal values
  const CFuint nbValues = avValues.size();
  const CFuint nbNodesInFace = geo->nbNodes();
  const CFreal invNbNodes = 1./static_cast<CFreal>(nbNodesInFace);
  cf_assert(_currFaceNodeIDs.size() == nbNodesInFace);


  for (CFuint i = 0; i < nbValues; ++i) {
    avValues[i] = 0.0;
    for (CFuint iNode = 0; iNode < _currFaceNodeIDs.size(); ++iNode) {
      avValues[i] += (*values[_currFaceNodeIDs[iNode]])[i];
    }
    avValues[i] *= invNbNodes;
  }
}

//////////////////////////////////////////////////////////////////////////////

void CorrectedDerivativeGG2D::setControlVolumeStates(vector<RealVector*>& states)
{  
  // left + right + face nodal nstates + other nstates in the cell 
  GeometricEntity *const face = getMethodData().getCurrentFace();
  State* const stateL = face->getState(0);
  State* const stateR = face->getState(1);
  const CFuint nbNodesInFace = face->nbNodes();
  
  // last state is the ghost state of the current face
  states[0] = stateL;
  states[1] = stateR;
  
  _currFaceNodeIDs.clear(); 
  _cellFaceIDsL.clear();
  _cellFaceIDsR.clear(); 
  _stencil.clear();
  cf_assert(_stencil.size() == 0);
  
  CFuint counter = 2;
  // left state computed always with Green-Gauss applied to the cell itself
  _cellBuilder.getDataGE().idx = stateL->getLocalID();
  GeometricEntity *const cell = _cellBuilder.buildGE();
  
  DataHandle<RealVector> nstates = socket_nstates.getDataHandle();
  const vector<Node*>& nodesL = *cell->getNodes();
  _nbNodesInCellL = nodesL.size();
  
  for (CFuint iNode = 0; iNode < _nbNodesInCellL; ++iNode) {
    const CFuint nodeID = nodesL[iNode]->getLocalID();
    states[counter++] = &nstates[nodeID];

    if (face->containNode(nodesL[iNode])) {
      _currFaceNodeIDs.push_back(counter);
    }
  }
  
  const vector<GeometricEntity*>& cellFacesL = *cell->getNeighborGeos();
  const CFuint nbCellFacesL =  cellFacesL.size();
  for (CFuint iFace = 0; iFace < nbCellFacesL; ++iFace) {
    _cellFaceIDsL.push_back(cellFacesL[iFace]->getID());
  }

  // connectivity table for the local face-node onnectivities
  _faceNodesL = LocalConnectionData::getInstance().
    getFaceDofLocal(cell->getShape(), CFPolyOrder::ORDER1, NODE, CFPolyForm::LAGRANGE);

  if (!stateR->isGhost()) {
    // release the left cell
    _cellBuilder.releaseGE();

    // here you should check that nodes are not doubled stored (this creates other complications ...)
    _cellBuilder.getDataGE().idx = stateR->getLocalID();

    // left state computed always with Green-Gauss applied to the cell itself
    GeometricEntity *const cell = _cellBuilder.buildGE();
    DataHandle<RealVector> nstates = socket_nstates.getDataHandle();
    const vector<Node*>& nodesR = *cell->getNodes();

    _nbNodesInCellR = nodesR.size();

    for (CFuint iNode = 0; iNode < _nbNodesInCellR; ++iNode) {
      const CFuint nodeID = nodesR[iNode]->getLocalID();
      states[counter++] = &nstates[nodeID];
    }

    const vector<GeometricEntity*>& cellFacesR = *cell->getNeighborGeos();
    const CFuint nbCellFacesR =  cellFacesR.size();
    for (CFuint iFace = 0; iFace < nbCellFacesR; ++iFace) {
      _cellFaceIDsR.push_back(cellFacesR[iFace]->getID());
    }

    // connectivity table for the local face-node onnectivities
    _faceNodesR = LocalConnectionData::getInstance().
      getFaceDofLocal(cell->getShape(), CFPolyOrder::ORDER1, NODE, CFPolyForm::LAGRANGE);
    _cellBuilder.releaseGE();
  }
  else {
    // if the face is a boundary face, consider the stencil made by
    // the ghost state and the inner cells sharing at least a node
    // with the face
    const vector<GeometricEntity*>& cellFaces = *cell->getNeighborGeos();
    const CFuint nbCellFaces =  cellFaces.size();

    for (CFuint iFace = 0; iFace < nbCellFaces; ++iFace) {
      GeometricEntity* currFace = cellFaces[iFace];

      if (currFace->getID() != face->getID()) {
	// if the face shares at least one node, consider its state
	bool sharesNode = false;
	for (CFuint iNode = 0; iNode < nbNodesInFace; ++iNode) {
	  if (currFace->containNode(face->getNode(iNode))) {
	    sharesNode = true;
	    break;
	  }
	}

	if (sharesNode) {
	  State* neighState = (currFace->getState(0) != stateL) ?
	    currFace->getState(0) : currFace->getState(1);
	  states[counter++] = neighState;
	  _stencil.push_back(neighState);
	}
      }
    }

    _cellBuilder.releaseGE();
  }

  _nbDofsInCV = counter;  
}

//////////////////////////////////////////////////////////////////////////////

SafePtr<vector<RealVector> > CorrectedDerivativeGG2D::getGradientsJacob()
{
  /// @TODO to be implemented properly
  return &_gradientsJacob;
}

//////////////////////////////////////////////////////////////////////////////

vector<SafePtr<BaseDataSocketSink> > CorrectedDerivativeGG2D::needsSockets()
{
  vector<SafePtr<BaseDataSocketSink> > result =
    DerivativeComputer::needsSockets();

  result.push_back(&socket_normals);
  result.push_back(&socket_gstates);
  result.push_back(&socket_volumes);
  result.push_back(&socket_states);
  result.push_back(&socket_nodes);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

void CorrectedDerivativeGG2D::computeOEminusNN(const RealVector& n)
{
  const CFreal nx = n[XX];
  const CFreal ny = n[YY];
  _oEminusNN(0,0) = 1.0 - nx*nx;
  _oEminusNN(0,1) = - nx*ny;
  _oEminusNN(1,0) = - nx*ny;
  _oEminusNN(1,1) = 1.0 - ny*ny;
}

//////////////////////////////////////////////////////////////////////////////

void CorrectedDerivativeGG2D::computeGradientsGG(SafePtr<Table<CFuint> > faceNodes,
						 vector<CFuint> cellFaceIDs,
						 State *const state,
						 const RealMatrix& values,
 						 CFuint counter,
						 RealMatrix& ux)
{
  DataHandle<CFreal> volumes = socket_volumes.getDataHandle();
  DataHandle<CFreal> normals = socket_normals.getDataHandle();
  DataHandle<CFint> isOutward = socket_isOutward.getDataHandle();

  const CFuint nbVars = values.nbRows();
  const CFuint nbFacesIDs = cellFaceIDs.size();
  const CFuint dim = PhysicalModelStack::getActive()->getDim();

  for (CFuint iFace = 0; iFace < nbFacesIDs; ++iFace) {
    const CFuint faceID = cellFaceIDs[iFace];
    const CFuint startNormalID = faceID*dim;

    // backup normal
    for (CFuint iDim = 0; iDim < dim; ++iDim) {
      _fnormal[iDim] = normals[startNormalID + iDim];
      if ( (CFuint) isOutward[faceID] != state->getLocalID()) {
	_fnormal[iDim] *= -1.;
      }
    }

    const CFuint nbNodesInFace = faceNodes->nbCols(iFace);
    const CFreal invNodesInFace = 1./nbNodesInFace;
    for (CFuint i = 0; i < nbVars; ++i) {
      for (CFuint iNode = 0; iNode < nbNodesInFace; ++iNode) {
	const CFuint colID = counter + (*faceNodes)(iFace,iNode);
        const CFreal vOvN = values(i,colID)*invNodesInFace;
	for (CFuint iDim = 0; iDim < dim; ++iDim) {
	  ux(iDim,i) += vOvN*_fnormal[iDim];
	}
      }
    }
  }
  ux /= volumes[state->getLocalID()];
}

//////////////////////////////////////////////////////////////////////////////

} // namespace FiniteVolume

} // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
