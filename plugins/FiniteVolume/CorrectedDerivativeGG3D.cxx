#include "CorrectedDerivativeGG3D.hh"
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

MethodStrategyProvider<CorrectedDerivativeGG3D,
                       CellCenterFVMData,
		       DerivativeComputer,
                       FiniteVolumeModule>
corrrectedDerivativeGG3DProvider("CorrectedGG3D");

//////////////////////////////////////////////////////////////////////////////

CorrectedDerivativeGG3D::CorrectedDerivativeGG3D(const std::string& name) :
  CorrectedDerivativeGG2D(name),
  _lf3()
{
}

//////////////////////////////////////////////////////////////////////////////

CorrectedDerivativeGG3D::~CorrectedDerivativeGG3D()
{
}

//////////////////////////////////////////////////////////////////////////////

void CorrectedDerivativeGG3D::setup()
{
  CorrectedDerivativeGG2D::setup();
  
  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  _lf3.resize(nbEqs);
}

//////////////////////////////////////////////////////////////////////////////

void CorrectedDerivativeGG3D::computeBoundaryGradients(const RealMatrix& values,
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
    const CFreal weight = (!_useWeights) ? 1.0 :
      1.0/MathFunctions::getDistance(nodeR,nodeLast);
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

  const CFreal invDet = 1./(l11*l22*l33
			    - l11*l23*l23
			    - l12*l12*l33
			    + l12*l13*l23
			    + l13*l12*l23
			    - l13*l13*l22);

  const CFreal linv11 = l22*l33 - l23*l23;
  const CFreal linv22 = l11*l33 - l13*l13;
  const CFreal linv33 = l11*l22 - l12*l12;
  const CFreal linv12 = -(l12*l33 - l13*l23);
  const CFreal linv13 = l12*l23 - l13*l22;
  const CFreal linv23 = -(l11*l23 - l13*l12);

  CFuint counter = 2;
  _uxL = 0.0;
  computeGradientsGG(_faceNodesL, _cellFaceIDsL, stateL, values, counter, _uxL);

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

    _gradLR[XX] = 0.5*(_uxL(XX,i) + uXR);
    _gradLR[YY] = 0.5*(_uxL(YY,i) + uYR);
    _gradLR[ZZ] = 0.5*(_uxL(ZZ,i) + uZR);

    correctGradient(dv,_gradLR, grad);
  }
}

//////////////////////////////////////////////////////////////////////////////

void CorrectedDerivativeGG3D::computeOEminusNN(const RealVector& n)
{
  const CFreal nx = n[XX];
  const CFreal ny = n[YY];
  const CFreal nz = n[ZZ];

  _oEminusNN(0,0) = 1.0 - nx*nx;
  _oEminusNN(0,1) = - ny*nx;
  _oEminusNN(0,2) = - nz*nx;

  _oEminusNN(1,0) = - nx*ny;
  _oEminusNN(1,1) = 1.0 - ny*ny;
  _oEminusNN(1,2) = - nz*ny;

  _oEminusNN(2,0) = - nx*nz;
  _oEminusNN(2,1) = - ny*nz;
  _oEminusNN(2,2) = 1.0 - nz*nz;
}

//////////////////////////////////////////////////////////////////////////////

} // namespace FiniteVolume

} // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
