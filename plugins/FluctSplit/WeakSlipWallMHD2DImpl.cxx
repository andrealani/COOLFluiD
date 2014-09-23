#include "FluctSplit/FluctSplitMHD.hh"
#include "WeakSlipWallMHD2DImpl.hh"
#include "InwardNormalsData.hh"
#include "Framework/PhysicalModel.hh"
#include "Framework/MethodCommandProvider.hh"
#include "MHD/MHD2DVarSet.hh"
#include "Framework/MeshData.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Physics::MHD;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<WeakSlipWallMHD2DImpl, 
		      FluctuationSplitData, 
		      FluctSplitMHDModule> 
weakSlipWallMHD2DImplProvider("WeakSlipWallMHD2DImpl");

//////////////////////////////////////////////////////////////////////////////

WeakSlipWallMHD2DImpl::WeakSlipWallMHD2DImpl(const std::string& name) :
  WeakSlipWall2DImpl(name),
  _varSet(CFNULL)
{
}

//////////////////////////////////////////////////////////////////////////////

WeakSlipWallMHD2DImpl::~WeakSlipWallMHD2DImpl()
{
}

//////////////////////////////////////////////////////////////////////////////

void WeakSlipWallMHD2DImpl::setup()
{
  WeakSlipWall2DImpl::setup();
  
  _varSet = getMethodData().getUpdateVar().d_castTo<MHD2DVarSet>();

  // set up the physical data
  _varSet->getModel()->resizePhysicalData(_physicalData);
 }

//////////////////////////////////////////////////////////////////////////////

void WeakSlipWallMHD2DImpl::computeNormalFluxAndJacob
(const State& state,
 const RealVector& normal,
 RealVector& flux,
 RealMatrix& fluxJacob)
{
  State& ss = *(const_cast<State*>(&state));
  _varSet->computePhysicalData(ss, _physicalData);

   const CFreal nx = normal[XX];
   const CFreal ny = normal[YY];
   const CFreal u = _physicalData[MHDTerm::VX];
   const CFreal v = _physicalData[MHDTerm::VY];
   const CFreal w = _physicalData[MHDTerm::VZ];
   const CFreal Bx = _physicalData[MHDTerm::BX];
   const CFreal By = _physicalData[MHDTerm::BY];
   const CFreal Bz = _physicalData[MHDTerm::BZ];
   const CFreal a = _physicalData[MHDTerm::A];

   const CFreal un = u*nx + v*ny;
   // unused // const CFreal Bn = Bx*nx + By*ny;
   const CFreal gamma = _varSet->getModel()->getGamma();
   const CFreal V2 = u*u + v*v + w*w;
   const CFreal B2 = Bx*Bx + By*By + Bz*Bz;
   // unused // const CFreal VB = u*Bx + v*By + w*Bz;
   const CFreal rho = _physicalData[MHDTerm::RHO];
   const CFreal a2 = a*a;

   flux[0] = rho*un;
   flux[1] = un*rho*u;
   flux[2] = un*rho*v;
   flux[3] = un*rho*w;
   flux[4] = un*Bx;
   flux[5] = un*By;
   flux[6] = un*Bz;
   flux[7] = un*(gamma/(gamma - 1.)*_physicalData[MHDTerm::P] + 0.5*rho*V2 + B2);

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

     fluxJacob(3,0) = -w*un;
     fluxJacob(3,1) = w*nx;
     fluxJacob(3,2) = w*ny;
     fluxJacob(3,3) = un;

     fluxJacob(4,0) = -un*Bx/rho;
     fluxJacob(4,1) = Bx/rho*nx;
     fluxJacob(4,2) = Bx/rho*ny;
     fluxJacob(4,4) = un;

     fluxJacob(5,0) = -un*By/rho;
     fluxJacob(5,1) = By/rho*nx;
     fluxJacob(5,2) = By/rho*ny;
     fluxJacob(5,5) = un;

     fluxJacob(6,0) = -un*Bz/rho;
     fluxJacob(6,1) = Bz/rho*nx;
     fluxJacob(6,2) = Bz/rho*ny;
     fluxJacob(6,6) = un;

     fluxJacob(7,0) = -a2*un/(gamma - 1.) + 0.5*(gamma - 2.)*V2*un - un*B2/rho;
     fluxJacob(7,1) = (a2/(gamma - 1.) + V2/2. + B2/rho)*nx - (gamma - 1.)*u*un;
     fluxJacob(7,2) = (a2/(gamma - 1.) + V2/2. + B2/rho)*ny - (gamma - 1.)*v*un;
     fluxJacob(7,3) = -(gamma - 1.)*w*un;
     fluxJacob(7,4) = (2. - gamma)*Bx*un;
     fluxJacob(7,5) = (2. - gamma)*By*un;
     fluxJacob(7,6) = (2. - gamma)*Bz*un;
     fluxJacob(7,7) = un*gamma;
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

     _tJacob(3,0) = -w*un;
     _tJacob(3,1) = w*nx;
     _tJacob(3,2) = w*ny;
     _tJacob(3,3) = un;

     _tJacob(4,0) = -un*Bx/rho;
     _tJacob(4,1) = Bx/rho*nx;
     _tJacob(4,2) = Bx/rho*ny;
     _tJacob(4,4) = un;

     _tJacob(5,0) = -un*By/rho;
     _tJacob(5,1) = By/rho*nx;
     _tJacob(5,2) = By/rho*ny;
     _tJacob(5,5) = un;

     _tJacob(6,0) = -un*Bz/rho;
     _tJacob(6,1) = Bz/rho*nx;
     _tJacob(6,2) = Bz/rho*ny;
     _tJacob(6,6) = un;

     _tJacob(7,0) = -a2*un/(gamma - 1.) - (1. - gamma/2.)*V2*un - un*B2/rho;
     _tJacob(7,1) = (a2/(gamma - 1.) + V2/2. + B2/rho)*nx - (gamma - 1.)*u*un;
     _tJacob(7,2) = (a2/(gamma - 1.) + V2/2. + B2/rho)*ny - (gamma - 1.)*v*un;
     _tJacob(7,3) = -(gamma - 1.)*w*un;
     _tJacob(7,4) = (2. - gamma)*Bx*un;
     _tJacob(7,5) = (2. - gamma)*By*un;
     _tJacob(7,6) = (2. - gamma)*Bz*un;
     _tJacob(7,7) = un*gamma;

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
