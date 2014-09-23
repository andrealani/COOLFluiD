#include "CorrectedDerivative3D.hh"
#include "Framework/MethodStrategyProvider.hh"
#include "Framework/GeometricEntity.hh"
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

MethodStrategyProvider<CorrectedDerivative3D,
                       CellCenterFVMData,
		       DerivativeComputer,
                       FiniteVolumeModule>
correctedDerivative3DProvider("Corrected3D");

//////////////////////////////////////////////////////////////////////////////

CorrectedDerivative3D::CorrectedDerivative3D(const std::string& name) :
  CorrectedDerivative2D(name),
  socket_uZ("uZ"),
  _lf3()
{
}

//////////////////////////////////////////////////////////////////////////////

CorrectedDerivative3D::~CorrectedDerivative3D()
{
}

//////////////////////////////////////////////////////////////////////////////

void CorrectedDerivative3D::setup()
{
  CorrectedDerivative2D::setup();

  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  _lf3.resize(nbEqs);
}

//////////////////////////////////////////////////////////////////////////////

void CorrectedDerivative3D::computeInnerGradients(const RealMatrix& values,
						  vector<RealVector*>& gradients)
{
  const CFuint nbVars = values.nbRows();
  cf_assert(gradients.size() == nbVars);

  State* const stateL = getMethodData().getCurrentFace()->getState(0);
  State* const stateR = getMethodData().getCurrentFace()->getState(1);
  
  DataHandle<CFreal> uX = socket_uX.getDataHandle();
  DataHandle<CFreal> uY = socket_uY.getDataHandle();
  DataHandle<CFreal> uZ = socket_uZ.getDataHandle();
  
  // if the face is a boundary face, consider the stencil made by
  // the ghost state and the inner cells sharing at least a node
  // with the face*
  for (CFuint i = 0; i < nbVars; ++i) {
    RealVector& grad = *gradients[i];
    cf_assert(grad.size() == 3);
   
    const CFreal dv = (values(i,1) - values(i,0))/_dr;
    const CFuint stateIDR = stateR->getLocalID();
    const CFuint stateIDL = stateL->getLocalID();
    
    _gradLR[XX] = 0.5*(uX(stateIDL,i,nbVars) + uX(stateIDR,i,nbVars));
    _gradLR[YY] = 0.5*(uY(stateIDL,i,nbVars) + uY(stateIDR,i,nbVars));
    _gradLR[ZZ] = 0.5*(uZ(stateIDL,i,nbVars) + uZ(stateIDR,i,nbVars));
    correctGradient(dv,_gradLR, grad);
  }
}

//////////////////////////////////////////////////////////////////////////////

void CorrectedDerivative3D::computeBoundaryGradients(const RealMatrix& values,
						     vector<RealVector*>& gradients)
{
  // if the face is a boundary face, consider the stencil made by
  // the ghost state and the inner cells sharing at least a node
  // with the face and apply a LS (least square) algorithm
  
  // left + right + nodal states + neighbor states
  GeometricEntity *const face = getMethodData().getCurrentFace();
  State* const stateL = face->getState(0);
  State* const stateR = face->getState(1);
  const CFuint nbVars = values.nbRows();
  cf_assert(gradients.size() == nbVars);
  
  // storage of states
  const CFuint startID = 2 + face->nbNodes();
  const CFuint stencilSize = _stencil.size();
  const RealVector& nodeR = stateR->getCoordinates();

  CFreal l11 = 0.0;
  CFreal l12 = 0.0;
  CFreal l13 = 0.0;
  CFreal l22 = 0.0;
  CFreal l23 = 0.0;
  CFreal l33 = 0.0;
  _lf1 = 0.0;
  _lf2 = 0.0;
  _lf3 = 0.0;

  // loop over the neighbor cells belonging to the chosen stencil
  for(CFuint in = 0; in < stencilSize; ++in) {
    const RealVector& nodeLast = _stencil[in]->getCoordinates();
    const CFreal deltaR = MathFunctions::getDistance(nodeR,nodeLast);
    const CFreal weight = 1.0/deltaR;
    const CFreal dx = weight*(nodeLast[XX] - nodeR[XX]);
    const CFreal dy = weight*(nodeLast[YY] - nodeR[YY]);
    const CFreal dz = weight*(nodeLast[ZZ] - nodeR[ZZ]);
    
    l11 += dx*dx;
    l12 += dx*dy;
    l13 += dx*dz;
    l22 += dy*dy;
    l23 += dy*dz;
    l33 += dz*dz;
    
    const CFuint nID = startID + in;
    for (CFuint iVar = 0; iVar < nbVars; ++iVar) {
      const CFreal du = weight*(values(iVar,nID) - values(iVar,1));
      _lf1[iVar] += dx*du;
      _lf2[iVar] += dy*du;
      _lf3[iVar] += dz*du;
    }
  }
  
  DataHandle<CFreal> uX = socket_uX.getDataHandle();
  DataHandle<CFreal> uY = socket_uY.getDataHandle();
  DataHandle<CFreal> uZ = socket_uZ.getDataHandle();
  
  const CFreal det = (l11*l22*l33 - l11*l23*l23 - l12*l12*l33 +
		      l12*l13*l23 + l13*l12*l23 - l13*l13*l22);
  cf_assert(std::abs(det) >= 0.);
  const CFreal invDet = 1./det;
  
  const CFreal linv11 = l22*l33 - l23*l23;
  const CFreal linv22 = l11*l33 - l13*l13;
  const CFreal linv33 = l11*l22 - l12*l12;
  const CFreal linv12 = -(l12*l33 - l13*l23);
  const CFreal linv13 = l12*l23 - l13*l22;
  const CFreal linv23 = -(l11*l23 - l13*l12);
  
  for (CFuint i = 0; i < nbVars; ++i) {
    RealVector& grad = *gradients[i];
    const CFreal dv = (values(i,1) - values(i,0))/_dr;
    const CFreal uXR = (linv11*_lf1[i] +
			linv12*_lf2[i] +
			linv13*_lf3[i])*invDet;
    const CFreal uYR = (linv12*_lf1[i] +
			linv22*_lf2[i] +
			linv23*_lf3[i])*invDet;
    const CFreal uZR = (linv13*_lf1[i] +
			linv23*_lf2[i] +
			linv33*_lf3[i])*invDet;
    
    const CFuint stateID = stateL->getLocalID();
    _gradLR[XX] = 0.5*(uX(stateID,i,nbVars) + uXR);
    _gradLR[YY] = 0.5*(uY(stateID,i,nbVars) + uYR);
    _gradLR[ZZ] = 0.5*(uZ(stateID,i,nbVars) + uZR);
    correctGradient(dv,_gradLR, grad);
  }
}

//////////////////////////////////////////////////////////////////////////////

vector<SafePtr<BaseDataSocketSink> > CorrectedDerivative3D::needsSockets()
{
  vector<SafePtr<BaseDataSocketSink> > result = CorrectedDerivative2D::needsSockets();
  result.push_back(&socket_uZ);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

void CorrectedDerivative3D::computeOEminusNN(const RealVector& n)
{
  const CFreal nx = n[XX];
  const CFreal ny = n[YY];
  const CFreal nz = n[ZZ];
  
  const CFreal nxny = nx*ny; 
  const CFreal nxnz = nx*nz;
  const CFreal nynz = ny*nz;
  
  _oEminusNN(0,0) = 1.0 - nx*nx;
  _oEminusNN(0,1) = - nxny;
  _oEminusNN(0,2) = - nxnz;

  _oEminusNN(1,0) = - nxny;
  _oEminusNN(1,1) = 1.0 - ny*ny;
  _oEminusNN(1,2) = - nynz;
  
  _oEminusNN(2,0) = - nxnz;
  _oEminusNN(2,1) = - nynz;
  _oEminusNN(2,2) = 1.0 - nz*nz;
}

//////////////////////////////////////////////////////////////////////////////

} // namespace FiniteVolume

} // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
