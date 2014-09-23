#include "WMUSCLPolyRec.hh"
#include "Framework/VolumeIntegrator.hh"
#include "Framework/MethodStrategyProvider.hh"
#include "Common/CFLog.hh"
#include "Framework/MeshData.hh"
#include "Framework/SubSystemStatus.hh"
#include "Framework/PhysicalModel.hh"

#include "FiniteVolume/FiniteVolume.hh"
#include "FiniteVolume/CellCenterFVMData.hh"

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

MethodStrategyProvider<WMUSCLPolyRec, 
		       CellCenterFVMData, 
		       PolyReconstructor<CellCenterFVMData>, 
		       FiniteVolumeModule> 
wmusclPolyRecProvider("WMUSCL");

//////////////////////////////////////////////////////////////////////////////

WMUSCLPolyRec::WMUSCLPolyRec(const std::string& name) :
  FVMCC_PolyRec(name),
  socket_stencil("stencil"),
  _upStates(4),
  _lStateBkp(),
  _rStateBkp(),
  _vLimiter(),
  _drUlUll(),
  _drUrUl(),
  _drUrrUr(),
  _r1(),
  _r2(),
  _faceNormal(),
  _midNode(),
  _dx(),
  _stencilPattern()
{
}

//////////////////////////////////////////////////////////////////////////////

WMUSCLPolyRec::~WMUSCLPolyRec()
{
}
      
//////////////////////////////////////////////////////////////////////////////

void WMUSCLPolyRec::configure ( Config::ConfigArgs& args )
{
  FVMCC_PolyRec::configure(args);
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
WMUSCLPolyRec::needsSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSink> > result =
    FVMCC_PolyRec::needsSockets();

  // Add the needed DataSocketSinks
  result.push_back(&socket_stencil);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

void WMUSCLPolyRec::computeGradients()
{
}

//////////////////////////////////////////////////////////////////////////////

void WMUSCLPolyRec::extrapolateImpl(GeometricEntity* const face)
{ 
  FVMCC_PolyRec::baseExtrapolateImpl(face);
  
  getMethodData().getUpdateVar()->setExtraData(true);
  (!isBoundaryFace()) ? extrapolateInner(face) : extrapolateBoundary(face);
  
  getBackupValues(LEFT) = getValues(LEFT);
  getBackupValues(RIGHT) =  getValues(RIGHT);
  
  getMethodData().getUpdateVar()->setExtraData(false);
}

//////////////////////////////////////////////////////////////////////////////

void WMUSCLPolyRec::extrapolateInner(GeometricEntity* const face)
{
  // U_i-1  U_i   U_i+1  U_i+2
  // c[2]   c[0]  c[1]   c[3]
  
  DataHandle<CFreal> newLimiter = socket_limiter.getDataHandle();
  DataHandle<vector<State*> > stencil = socket_stencil.getDataHandle();
  const CFuint dim = PhysicalModelStack::getActive()->getDim();
  SafePtr<Limiter<CellCenterFVMData> > limiter = getMethodData().getLimiter();
  
  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  const CFreal residual = SubSystemStatusStack::getActive()->getResidual();
  const CFuint faceID = face->getID();
  cf_assert (face->getID() < stencil.size());
  const vector<State*>& cs = stencil[faceID];
  
  cf_assert (_zeroGradient != CFNULL);
  
  cf_assert (cs.size() == 4);
  _drUlUll = (*cs[0]) - (*cs[2]);
  _drUrUl  = (*cs[1]) - (*cs[0]);
  _drUrrUr = (*cs[3]) - (*cs[1]);
  
  // fix for division by zero
  for  (CFuint iVar = 0; iVar < nbEqs; ++iVar) { 
    _r1[iVar] = _drUlUll[iVar]/max(_drUrUl[iVar],MathTools::MathConsts::CFrealEps());
    _r2[iVar] = _drUrUl[iVar]/max(_drUrrUr[iVar],MathTools::MathConsts::CFrealEps());
  }
  
  const CFuint startID = faceID*nbEqs;
  const CFuint twoNbEqs = 2*nbEqs;
  
  if (residual > _limitRes) {
    for (CFuint iVar = 0; iVar < twoNbEqs; ++iVar) {
      newLimiter[startID + iVar] = 1.0;
    }
    limiter->limitOnFace(_r1, _r2, &newLimiter[startID]);
  }
  else {
    if (!_freezeLimiter) { 
      _vLimiter = 1.0;
      
      // historical modification of the limiter
      limiter->limitOnFace(_r1, _r2, &_vLimiter[0]);
      CFuint currID = startID;
      for (CFuint iVar = 0; iVar < twoNbEqs; ++iVar, ++currID) {
	newLimiter[currID] = min(_vLimiter[iVar],newLimiter[currID]);
      }
    }
  }

  if (dim > DIM_1D) { 
    _midNode = 0.5*(*face->getNode(0) + *face->getNode(1));
  }
  else {
    _midNode =  (*face->getNode(0));
  }

  _dx(0,0) = MathFunctions::getDistance(cs[0]->getCoordinates(), cs[2]->getCoordinates());
  _dx(1,0) = MathFunctions::getDistance(cs[1]->getCoordinates(), cs[3]->getCoordinates());
  _dx(0,1) = MathFunctions::getDistance(cs[0]->getCoordinates(),_midNode);
  _dx(1,1) = MathFunctions::getDistance(cs[1]->getCoordinates(),_midNode);  
  
  State& leftValues  = getValues(LEFT);
  State& rightValues = getValues(RIGHT);
  
  CFuint count = startID;
  for  (CFuint iVar = 0; iVar < nbEqs; ++iVar) {
    const CFreal duL = newLimiter[count++]*_drUlUll[iVar]/_dx(0,0)*_dx(0,1);
    leftValues[iVar]  = (*cs[0])[iVar] + duL;
    
    const CFreal duR = newLimiter[count++]*_drUrrUr[iVar]/_dx(1,0)*_dx(1,1);
    rightValues[iVar] = (*cs[1])[iVar] - duR;
  }
}
      
//////////////////////////////////////////////////////////////////////////////

void WMUSCLPolyRec::extrapolateBoundary(GeometricEntity* const face)
{
  // if we are not perturbing the residual compute and store the 
  // face mid point
  if (!getMethodData().isPerturb()) {
    computeMidFaceNode(face);
  } 
  
  State *const innerState = face->getState(LEFT);
  State *const ghostState = face->getState(RIGHT);
  const Node& inNode = innerState->getCoordinates();
  const Node& ghNode = ghostState->getCoordinates();
  const CFreal drXiXg = MathFunctions::getDistance(inNode, ghNode);
  const CFreal drXgXw = drXiXg - _drXiXw;
  const CFreal invDr = 1./drXiXg;
  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  
  for (CFuint iVar = 0; iVar < nbEqs; ++iVar) {
    if (!(*_zeroGradient)[iVar]) {	   
      // inverse distance based weighted linear reconstruction
      getValues(LEFT)[iVar] = getValues(RIGHT)[iVar] =
	(drXgXw*(*innerState)[iVar] + _drXiXw*(*ghostState)[iVar])*invDr;
    }
    else {
      // ghostState == innerState
      getValues(LEFT)[iVar] = getValues(RIGHT)[iVar] = (*innerState)[iVar];
    }
  }
}
      
//////////////////////////////////////////////////////////////////////////////

void WMUSCLPolyRec::extrapolateImpl(GeometricEntity* const face,
				   CFuint iVar, CFuint leftOrRight)
{
  getMethodData().getUpdateVar()->setExtraData(true);
  
  copyBackupValues();
  if (!isBoundaryFace()) {
    computeReconstruction(face, leftOrRight, iVar);
  }
  else {
    DataHandle<vector<State*> > stencil = socket_stencil.getDataHandle();
//     const CFuint faceID = face->getID();
    cf_assert (face->getID() < stencil.size());
//     const vector<State*>& cs = stencil[faceID];
    
    cf_assert (isBoundaryFace()); 
    State *const innerState = face->getState(LEFT);
    State *const ghostState = face->getState(RIGHT);
    const Node& inNode = innerState->getCoordinates();
    const Node& ghNode = ghostState->getCoordinates();
    const CFreal drXiXg = MathFunctions::getDistance(inNode, ghNode);
    const CFreal drXgXw = drXiXg - _drXiXw;
    cf_assert (drXgXw > 0.);
    const CFreal invDr = 1./drXiXg;
    
    // inverse UNLIMITED distance based weighted linear reconstruction
    const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
    for (CFuint iVar = 0; iVar < nbEqs; ++iVar) {
      if (!(*_zeroGradient)[iVar]) {
	getValues(leftOrRight)[iVar] = 
	  (drXgXw*(*innerState)[iVar] + _drXiXw*(*ghostState)[iVar])*invDr;
      }
      else {
        // ghostState == innerState
	getValues(leftOrRight)[iVar] = (*innerState)[iVar];
      } 
    }
  }
  
  getMethodData().getUpdateVar()->setExtraData(false);
}

//////////////////////////////////////////////////////////////////////////////
      
void WMUSCLPolyRec::updateWeights()
{
  //  FVMCC_PolyRec::updateWeights();
//  cout << "WMUSCLPolyRec::updateWeights()" << endl;
 // throw Common::NotImplementedException (FromHere(),"WMUSCLPolyRec::updateWeights()");
}

//////////////////////////////////////////////////////////////////////////////

void WMUSCLPolyRec::setup()
{
  FVMCC_PolyRec::setup();
  
  // re-setup the variable transformer
  getMethodData().getUpdateToSolutionVecTrans()->setup(4);

  _lStateBkp.resize(PhysicalModelStack::getActive()->getNbEq());
  _rStateBkp.resize(PhysicalModelStack::getActive()->getNbEq());
  _vLimiter.resize(2*PhysicalModelStack::getActive()->getNbEq());
  _drUlUll.resize(PhysicalModelStack::getActive()->getNbEq());
  _drUrUl.resize(PhysicalModelStack::getActive()->getNbEq());
  _drUrrUr.resize(PhysicalModelStack::getActive()->getNbEq());
  _r1.resize(PhysicalModelStack::getActive()->getNbEq());
  _r2.resize(PhysicalModelStack::getActive()->getNbEq());
  _faceNormal.resize(PhysicalModelStack::getActive()->getDim());
  _midNode.resize(PhysicalModelStack::getActive()->getDim());
  
  // distances matrix
  _dx.resize(2,2);
  
  // reconstruction stencil pattern for left(0) and right(1) states
  _stencilPattern.resize(2,2);
  _stencilPattern(0,0) = 0;
  _stencilPattern(0,1) = 2;
  _stencilPattern(1,0) = 1;
  _stencilPattern(1,1) = 3;
}
      
//////////////////////////////////////////////////////////////////////////////

void WMUSCLPolyRec::computeMidFaceNode(GeometricEntity* const face)
{
  const CFuint dim = PhysicalModelStack::getActive()->getDim();
  const CFuint faceID = face->getID();
  const CFuint startID = faceID*dim;
  
  DataHandle<CFreal> normals = socket_normals.getDataHandle();
  
  // set the current normal
  for (CFuint i = 0; i < dim; ++i) {
    _faceNormal[i] = normals[startID + i];
  }

  // compute the original position of the ghost state @see ComputeDummyState
  const CFreal k = - MathFunctions::innerProd(_faceNormal, (*face->getNode(0)));
  const CFreal n2 = MathFunctions::innerProd(_faceNormal, _faceNormal);
  cf_assert (std::abs(n2) > 0.0);
  State *const innerState = face->getState(0);
  const Node& innerNode = innerState->getCoordinates();
  const CFreal t = (MathFunctions::innerProd(_faceNormal,innerNode) + k)/n2;
  _midNode = innerNode - t*_faceNormal;
  
  // first calculate the "unmodified distances" inner-wall, inner-ghost
  _drXiXw = MathFunctions::getDistance(innerNode,_midNode);
}

//////////////////////////////////////////////////////////////////////////////

void WMUSCLPolyRec::computeReconstruction(GeometricEntity* const face, 
					  CFuint leftOrRight, CFuint iVar)
{
  DataHandle<CFreal> newLimiter = socket_limiter.getDataHandle();
  DataHandle<vector<State*> > stencil = socket_stencil.getDataHandle();
  cf_assert (_zeroGradient != CFNULL);
  
  copyBackupValues();
  
  const CFuint faceID = face->getID();
  cf_assert (face->getID() < stencil.size());
  const vector<State*>& cs = stencil[faceID];
  cf_assert (cs.size() == 4);
  
  const CFreal residual = SubSystemStatusStack::getActive()->getResidual();
  const CFuint startID = faceID*PhysicalModelStack::getActive()->getNbEq();
  
  CFreal r = 0.0;
  const CFreal drUrUl = (*cs[1])[iVar] - (*cs[0])[iVar];
  if (leftOrRight == LEFT) {
    r = ((*cs[0])[iVar] - (*cs[2])[iVar])/max(drUrUl,MathTools::MathConsts::CFrealEps());
  }
  else {
    r = drUrUl/max((*cs[3])[iVar] - (*cs[1])[iVar], MathTools::MathConsts::CFrealEps());
  }
  
  SafePtr<Limiter<CellCenterFVMData> > limiter = getMethodData().getLimiter();
  
  CFreal limiterValue = 1.0;
  const CFuint limiterID = startID + iVar*2 + leftOrRight;
  if (residual > _limitRes) {
    limiter->limitScalar(r, limiterValue);
  }
  else {
    if (!_freezeLimiter) { 
      CFreal tmpLimiter = 1.0;
      // historical modification of the limiter
      limiter->limitScalar(r, tmpLimiter);
      limiterValue = min(tmpLimiter,newLimiter[limiterID]);
    }
  }    
  
  const CFreal u0 = (*cs[_stencilPattern(leftOrRight,0)])[iVar];
  const CFreal u1 = (*cs[_stencilPattern(leftOrRight,1)])[iVar];
  const CFreal dr = (u0 - u1)/_dx(leftOrRight,0)*_dx(leftOrRight,1);
  
  getValues(leftOrRight)[iVar]  = u0 + limiterValue*dr;
  
  //  (*getValues(leftOrRight)[0])[iVar]  = u0 + newLimiter[startID + iVar*2 + leftOrRight]*dr;
}
      
//////////////////////////////////////////////////////////////////////////////
      
    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD
