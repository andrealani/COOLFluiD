#include "FluctSplit/FluctSplitNavierStokes.hh"
#include "WeakSlipWallEuler2DImpl.hh"
#include "InwardNormalsData.hh"
#include "Framework/PhysicalModel.hh"
#include "Framework/MethodCommandProvider.hh"
#include "NavierStokes/Euler2DVarSet.hh"
#include "Framework/MeshData.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Physics::NavierStokes;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<WeakSlipWallEuler2DImpl, 
		      FluctuationSplitData, 
		      FluctSplitNavierStokesModule> 
weakSlipWallEuler2DImplProvider("WeakSlipWallEuler2DImpl");

//////////////////////////////////////////////////////////////////////////////

WeakSlipWallEuler2DImpl::WeakSlipWallEuler2DImpl(const std::string& name) :
   WeakSlipWall2DImpl(name),
   _varSet(CFNULL)
{
}

//////////////////////////////////////////////////////////////////////////////

WeakSlipWallEuler2DImpl::~WeakSlipWallEuler2DImpl()
{
}

//////////////////////////////////////////////////////////////////////////////

void WeakSlipWallEuler2DImpl::setup()
{
  WeakSlipWall2DImpl::setup();
  
  _varSet = getMethodData().getUpdateVar().d_castTo<Euler2DVarSet>();
  
  // set up the physical data
  _varSet->getModel()->resizePhysicalData(_physicalData);
}

//////////////////////////////////////////////////////////////////////////////

void WeakSlipWallEuler2DImpl::computeNormalFluxAndJacob
(const State& state,
 const RealVector& normal,
 RealVector& flux,
 RealMatrix& fluxJacob)
{
  State& ss = *(const_cast<State*>(&state));
  _varSet->computePhysicalData(ss, _physicalData); 
  
  const CFreal nx = normal[XX];
  const CFreal ny = normal[YY];
  const CFreal rho = _physicalData[EulerTerm::RHO];
  const CFreal u = _physicalData[EulerTerm::VX];
  const CFreal v = _physicalData[EulerTerm::VY];
  const CFreal un = u*nx + v*ny;
  const CFreal vel2 = u*u + v*v;
  const CFreal H = _physicalData[EulerTerm::H];
  const CFreal gamma = _varSet->getModel()->getGamma();
  const CFreal gammaMinus1 = gamma - 1.;
  const CFreal halfGammaMinus1 = 0.5*gammaMinus1;

   flux[0] = rho*un;
   flux[1] = un*rho*u;
   flux[2] = un*rho*v;
   flux[3] = un*rho*H;

   if (!getMethodData().isResidualTransformationNeeded()) {
     // the analytical jacobian of the normal fluxes
     fluxJacob(0,1) = nx;
     fluxJacob(0,2) = ny;
     fluxJacob(1,0) = -u*un;
     fluxJacob(1,1) = un + u*nx;
     fluxJacob(1,2) = u*ny;
     fluxJacob(2,0) = -v*un;
     fluxJacob(2,1) = v*nx;
     fluxJacob(2,2) = un + v*ny;
     fluxJacob(3,0) = un*(halfGammaMinus1*vel2 - H);
     fluxJacob(3,1) = H*nx - un*u*gammaMinus1;
     fluxJacob(3,2) = H*ny - un*v*gammaMinus1;
     fluxJacob(3,3) = _varSet->getModel()->getGamma()*un;
   }
   else {
     // the analytical jacobian of the normal fluxes
     _tJacob(0,1) = nx;
     _tJacob(0,2) = ny;
     _tJacob(1,0) = -u*un;
     _tJacob(1,1) = un + u*nx;
     _tJacob(1,2) = u*ny;
     _tJacob(2,0) = -v*un;
     _tJacob(2,1) = v*nx;
     _tJacob(2,2) = un + v*ny;
     _tJacob(3,0) = un*(halfGammaMinus1*vel2 - H);
     _tJacob(3,1) = H*nx - un*u*gammaMinus1;
     _tJacob(3,2) = H*ny - un*v*gammaMinus1;
     _tJacob(3,3) = _varSet->getModel()->getGamma()*un;

     // set the transformation from update to solution in update
     SafePtr<VarSetMatrixTransformer> updateToSolInUpdate =
       getMethodData().getUpdateToSolutionInUpdateMatTrans();

     updateToSolInUpdate->setMatrix(state);
     const RealMatrix& tMatrix = *updateToSolInUpdate->getMatrix();

     fluxJacob = _tJacob*tMatrix;
   }
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
