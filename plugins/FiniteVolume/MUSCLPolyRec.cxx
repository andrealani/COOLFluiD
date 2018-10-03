#include "MUSCLPolyRec.hh"
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

MethodStrategyProvider<MUSCLPolyRec, 
		       CellCenterFVMData, 
		       PolyReconstructor<CellCenterFVMData>, 
		       FiniteVolumeModule> 
musclPolyRecProvider("MUSCL");

//////////////////////////////////////////////////////////////////////////////

MUSCLPolyRec::MUSCLPolyRec(const std::string& name) :
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
  _midNode()
{
}

//////////////////////////////////////////////////////////////////////////////

MUSCLPolyRec::~MUSCLPolyRec()
{
}
      
//////////////////////////////////////////////////////////////////////////////

void MUSCLPolyRec::configure ( Config::ConfigArgs& args )
{
  FVMCC_PolyRec::configure(args);
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
MUSCLPolyRec::needsSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSink> > result =
    FVMCC_PolyRec::needsSockets();

  // Add the needed DataSocketSinks
  result.push_back(&socket_stencil);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

void MUSCLPolyRec::computeGradients()
{
}

//////////////////////////////////////////////////////////////////////////////

void MUSCLPolyRec::extrapolateImpl(GeometricEntity* const face)
{ 
  FVMCC_PolyRec::baseExtrapolateImpl(face);
  
  getMethodData().getUpdateVar()->setExtraData(true);
  (!isBoundaryFace()) ? extrapolateInner(face) : extrapolateBoundary(face);
  
  getBackupValues(LEFT)  = getValues(LEFT);
  getBackupValues(RIGHT) = getValues(RIGHT);
    
  getMethodData().getUpdateVar()->setExtraData(false);
}

//////////////////////////////////////////////////////////////////////////////

void MUSCLPolyRec::extrapolateInner(GeometricEntity* const face)
{
  // U_i-1  U_i   U_i+1  U_i+2
  // c[2]   c[0]  c[1]   c[3]
  
  DataHandle<CFreal> newLimiter = socket_limiter.getDataHandle();
  
  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  const CFreal residual = SubSystemStatusStack::getActive()->getResidual();
  const CFuint faceID = face->getID();
  
  const vector<State*>& cs = *getStencil(faceID,false);
  
  SafePtr<Limiter<CellCenterFVMData> > limiter = getMethodData().getLimiter();
  
  cf_assert (_zeroGradient != CFNULL);
  
  cf_assert (cs.size() == 4);
  _drUlUll = (*cs[0]) - (*cs[2]);
  _drUrUl  = (*cs[1]) - (*cs[0]);
  _drUrrUr = (*cs[3]) - (*cs[1]);
  
  // fix for division by zero
  for  (CFuint iVar = 0; iVar < nbEqs; ++iVar) { 
    _r1[iVar] = (std::abs(_drUrUl[iVar]) > 0.) ? _drUlUll[iVar]/_drUrUl[iVar] : 0.0;
    _r2[iVar] = (std::abs(_drUrrUr[iVar]) > 0.) ? _drUrUl[iVar]/_drUrrUr[iVar] : 0.0;
    //  _r2[iVar] = (std::abs(_drUrUl[iVar]) > 0.) ? _drUrrUr[iVar]/_drUrUl[iVar] : 0.0;
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

  if (PhysicalModelStack::getActive()->getDim() > DIM_1D) {
    _midNode = 0.5*(*face->getNode(0) + *face->getNode(1));
  }
  else {
    _midNode = *face->getNode(0);
  }

  State& leftValues  = getValues(LEFT);
  State& rightValues = getValues(RIGHT);
  CFuint count = startID;
  for  (CFuint iVar = 0; iVar < nbEqs; ++iVar) {
    const CFreal duL = newLimiter[count++]*_drUrUl[iVar]*0.5;
    //const CFreal duL = newLimiter[count++]*_drUlUll[iVar]*0.5;
    leftValues[iVar]  = (*cs[0])[iVar] + duL;
    
    const CFreal duR = newLimiter[count++]*_drUrrUr[iVar]*0.5;
    rightValues[iVar] = (*cs[1])[iVar] - duR;
  }
}
      
//////////////////////////////////////////////////////////////////////////////

void MUSCLPolyRec::extrapolateBoundary(GeometricEntity* const face)
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

  // inverse distance based weighted linear reconstruction
  // *getValues(LEFT)[0] = *getValues(RIGHT)[0] = 
  //  (drXgXw*(*innerState) + _drXiXw*(*ghostState))*invDr;

  const vector<State*>& cs = *getStencil(face->getID(),true);
  getValues(RIGHT) = (drXgXw*(*cs[0]) + _drXiXw*(*cs[1]))*invDr;
  
  //   DataHandle<vector<State*> > stencil = socket_stencil.getDataHandle();
  //   const vector<State*>& ss = stencil[face->getID()];
  
  if (PhysicalModelStack::getActive()->getDim() > DIM_1D) {
    _midNode = 0.5*(*face->getNode(0) + *face->getNode(1));
  }
  else {
    _midNode = *face->getNode(0);
  }
  
  getValues(LEFT) = (*cs[0]) + 0.5*((*cs[1]) - (*cs[0]));
  //*getValues(LEFT)[0] = (*cs[0]) + 0.5*((*cs[0]) - (*cs[2]));
}

//////////////////////////////////////////////////////////////////////////////

void MUSCLPolyRec::extrapolateImpl(GeometricEntity* const face,
				   CFuint iVar, CFuint leftOrRight)
{
  getMethodData().getUpdateVar()->setExtraData(true);
  
  copyBackupValues();
  if (!isBoundaryFace()) {
    computeReconstruction(face, leftOrRight, iVar);
  }
  else {

    cf_assert (isBoundaryFace()); 
    const vector<State*>& cs = *getStencil(face->getID(),true);
    // inverse UNLIMITED distance based weighted linear reconstruction
    if (leftOrRight == RIGHT) {
      State *const innerState = face->getState(LEFT);
      State *const ghostState = face->getState(RIGHT);
      const Node& inNode = innerState->getCoordinates();
      const Node& ghNode = ghostState->getCoordinates();
      const CFreal drXiXg = MathFunctions::getDistance(inNode, ghNode);
      const CFreal drXgXw = drXiXg - _drXiXw;
      cf_assert (drXgXw > 0.);
      const CFreal invDr = 1./drXiXg;
      getValues(leftOrRight) = (drXgXw*(*cs[0]) + _drXiXw*(*cs[1]))*invDr;
    }
    else {
      //  getValues(leftOrRight) = (*cs[0]) + 0.5*((*cs[0]) - (*cs[2]));
      getValues(leftOrRight) = (*cs[0]) + 0.5*((*cs[1]) - (*cs[0]));
    }
  }
  
  getMethodData().getUpdateVar()->setExtraData(false);
}

//////////////////////////////////////////////////////////////////////////////
      
void MUSCLPolyRec::updateWeights()
{
  throw Common::NotImplementedException (FromHere(),"MUSCLPolyRec::updateWeights()");
}

//////////////////////////////////////////////////////////////////////////////

void MUSCLPolyRec::setup()
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
}
      
//////////////////////////////////////////////////////////////////////////////

void MUSCLPolyRec::computeMidFaceNode(GeometricEntity* const face)
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

void MUSCLPolyRec::computeReconstruction(GeometricEntity* const face, 
					 CFuint leftOrRight, CFuint iVar)
{
  DataHandle<CFreal> newLimiter = socket_limiter.getDataHandle();
  cf_assert (_zeroGradient != CFNULL);
  
  SafePtr<Limiter<CellCenterFVMData> > limiter = getMethodData().getLimiter();
  
  const CFuint faceID = face->getID();
  const vector<State*>& cs = *getStencil(faceID,false);
  cf_assert (cs.size() == 4);
  
  const CFreal residual = SubSystemStatusStack::getActive()->getResidual();
  const CFuint startID = faceID*PhysicalModelStack::getActive()->getNbEq();
  
  CFreal r = 0.0;
  const CFreal drUrUl = (*cs[1])[iVar] - (*cs[0])[iVar];
  if (leftOrRight == LEFT) {
    r = (std::abs(drUrUl) > 0.0) ? ((*cs[0])[iVar] - (*cs[2])[iVar])/drUrUl : 0.0;
  }
  else {
    const CFreal drUrrUr = (*cs[3])[iVar] - (*cs[1])[iVar];
    r = (std::abs(drUrrUr) > 0.0) ? drUrUl/drUrrUr : 0.0;
    // r = (std::abs(drUrUl > 0.0)) ? drUrrUr/drUrUl : 0.0;
  }
  
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
  
  if (leftOrRight == LEFT) {
    // (*getValues(leftOrRight)[0])[iVar] = 
    // (*cs[0])[iVar] + limiterValue*0.5*((*cs[0])[iVar] - (*cs[2])[iVar]);
    getValues(leftOrRight)[iVar] = 
      (*cs[0])[iVar] + limiterValue*0.5*((*cs[1])[iVar] - (*cs[0])[iVar]); 
  }
  else {
    getValues(leftOrRight)[iVar] = 
      (*cs[1])[iVar] - limiterValue*0.5*((*cs[3])[iVar] - (*cs[1])[iVar]); 
  }
}
      
//////////////////////////////////////////////////////////////////////////////
      
vector<State*>* MUSCLPolyRec::getStencil(CFuint faceID, bool isBoundaryFace)
{  
  DataHandle<vector<State*> > stencil = socket_stencil.getDataHandle();
  cf_assert (faceID < stencil.size());
  if (getMethodData().reconstructSolVars()) {
    const vector<State*>& ss = stencil[faceID];  
    _upStates.resize((!isBoundaryFace) ? 4 : 3);
    _upStates[0] = ss[0];
    _upStates[1] = ss[1];
    _upStates[2] = ss[2];
    if (!isBoundaryFace) {
      _upStates[3] = ss[3];
    }
    return getMethodData().getUpdateToSolutionVecTrans()->transform(&_upStates);
  }
  else {
    return &stencil[faceID];
  }
}

//////////////////////////////////////////////////////////////////////////////
 
   } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD
