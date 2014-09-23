#include "FluctSplit/FluctSplitMHD.hh"
#include "WeakSlipWallMHD2DProjectionImpl.hh"
#include "InwardNormalsData.hh"
#include "Framework/PhysicalModel.hh"
#include "Framework/MethodCommandProvider.hh"
#include "MHD/MHD2DProjectionVarSet.hh"
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

MethodCommandProvider<WeakSlipWallMHD2DProjectionImpl, 
		      FluctuationSplitData, 
		      FluctSplitMHDModule> 
weakSlipWallMHD2DProjectionImplProvider("WeakSlipWallMHD2DProjectionImpl");

//////////////////////////////////////////////////////////////////////////////

void WeakSlipWallMHD2DProjectionImpl::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOption< CFreal >("refPhi","Reference phi value imposed at the wall.");
}

//////////////////////////////////////////////////////////////////////////////

WeakSlipWallMHD2DProjectionImpl::WeakSlipWallMHD2DProjectionImpl(const std::string& name) :
  WeakSlipWall2DImpl(name),
  _varSet(CFNULL)
{
  addConfigOptionsTo(this);
  _refPhi = 0.;
  setParameter("refPhi",&_refPhi);
}

//////////////////////////////////////////////////////////////////////////////

WeakSlipWallMHD2DProjectionImpl::~WeakSlipWallMHD2DProjectionImpl()
{
}

//////////////////////////////////////////////////////////////////////////////

void WeakSlipWallMHD2DProjectionImpl::setup()
{
  WeakSlipWall2DImpl::setup();
  
  _varSet = getMethodData().getUpdateVar().d_castTo<MHD2DProjectionVarSet>();
  
  // set up the physical data
  _varSet->getModel()->resizePhysicalData(_physicalData);
}
      
//////////////////////////////////////////////////////////////////////////////

void WeakSlipWallMHD2DProjectionImpl::computeNormalFluxAndJacob
(const State& state,
 const RealVector& normal,
 RealVector& flux,
 RealMatrix& fluxJacob)
{
  State& ss = *(const_cast<State*>(&state));
  _varSet->computePhysicalData(ss, _physicalData);

   const CFreal nx = normal[XX];
   const CFreal ny = normal[YY];
   const CFreal u = _physicalData[MHDProjectionTerm::VX];
   const CFreal v = _physicalData[MHDProjectionTerm::VY];
   const CFreal w = _physicalData[MHDProjectionTerm::VZ];
   const CFreal Bx = _physicalData[MHDProjectionTerm::BX];
   const CFreal By = _physicalData[MHDProjectionTerm::BY];
   const CFreal Bz = _physicalData[MHDProjectionTerm::BZ];

   const CFreal un = u*nx + v*ny;
   const CFreal Bn = Bx*nx + By*ny;
   const CFreal gamma = _varSet->getModel()->getGamma();
   const CFreal beta = _varSet->getModel()->getRefSpeed();
   const CFreal V2 = u*u + v*v + w*w;
   const CFreal B2 = Bx*Bx + By*By + Bz*Bz;
   const CFreal VB = u*Bx + v*By + w*Bz;
   const CFreal rho = _physicalData[MHDProjectionTerm::RHO];
   const CFreal a = _physicalData[MHDProjectionTerm::A];
   const CFreal a2 = a*a;
//   const CFreal _refPhi = _physicalData[MHDProjectionTerm::PHI];

   flux[0] = rho*un;
   flux[1] = un*rho*u - Bn*Bx;
   flux[2] = un*rho*v - Bn*By;
   flux[3] = un*rho*w - Bn*Bz;
   flux[4] = un*Bx - Bn*u + nx*_refPhi;
   flux[5] = un*By - Bn*v + ny*_refPhi;
   flux[6] = un*Bz - Bn*w;
   flux[7] = un*(gamma/(gamma - 1.)*_physicalData[MHDProjectionTerm::P] + 0.5*rho*V2 + B2) - Bn*VB;
   flux[8] = Bn*beta*beta;

   if (!getMethodData().isResidualTransformationNeeded()) {
     // the analytical jacobian of the normal fluxes
     fluxJacob(0,1) = nx;
     fluxJacob(0,2) = ny;

     fluxJacob(1,0) = -u*un;
     fluxJacob(1,1) = un + u*nx;
     fluxJacob(1,2) = u*ny;
     fluxJacob(1,4) = -Bn - Bx*nx;
     fluxJacob(1,5) = -Bx*ny;

     fluxJacob(2,0) = -v*un;
     fluxJacob(2,1) = v*nx;
     fluxJacob(2,2) = un + v*ny;
     fluxJacob(2,4) = -By*nx;
     fluxJacob(2,5) = -Bn - By*ny;

     fluxJacob(3,0) = -w*un;
     fluxJacob(3,1) = w*nx;
     fluxJacob(3,2) = w*ny;
     fluxJacob(3,3) = un;
     fluxJacob(3,4) = -Bz*nx;
     fluxJacob(3,5) = -Bz*ny;
     fluxJacob(3,6) = -Bn;

     fluxJacob(4,0) = -un*Bx/rho + u*Bn/rho;
     fluxJacob(4,1) = -By/rho*ny;
     fluxJacob(4,2) = Bx/rho*ny;
     fluxJacob(4,4) = v*ny;
     fluxJacob(4,5) = -u*ny;
//     fluxJacob(4,8) = nx;

     fluxJacob(5,0) = -un*By/rho + v*Bn/rho;
     fluxJacob(5,1) = By/rho*nx;
     fluxJacob(5,2) = -Bx/rho*nx;
     fluxJacob(5,4) = -v*nx;
     fluxJacob(5,5) = u*nx;
//     fluxJacob(5,8) = ny;

     fluxJacob(6,0) = -un*Bz/rho + w*Bn/rho;
     fluxJacob(6,1) = Bz/rho*nx;
     fluxJacob(6,2) = Bz/rho*ny;
     fluxJacob(6,3) = -Bn/rho;
     fluxJacob(6,4) = -w*nx;
     fluxJacob(6,5) = -w*ny;
     fluxJacob(6,6) = un;

     fluxJacob(7,0) = -a2*un/(gamma - 1.) - (2. - gamma)/2.*V2*un - 0.5*un*B2/rho + Bn/rho*VB;
     fluxJacob(7,1) = (a2/(gamma - 1.) + V2/2. + B2/rho)*nx - (gamma - 1.)*u*un - Bx/rho*Bn;
     fluxJacob(7,2) = (a2/(gamma - 1.) + V2/2. + B2/rho)*ny - (gamma - 1.)*v*un - By/rho*Bn;
     fluxJacob(7,3) = -(gamma - 1.)*w*un - Bz/rho*Bn;
     fluxJacob(7,4) = (2. - gamma)*Bx*un - VB*nx - u*Bn;
     fluxJacob(7,5) = (2. - gamma)*By*un - VB*ny - v*Bn;
     fluxJacob(7,6) = (2. - gamma)*Bz*un - w*Bn;
     fluxJacob(7,7) = un*gamma;

     fluxJacob(8,4) = beta*beta*nx;
     fluxJacob(8,5) = beta*beta*ny;
   }

   else {
     // the analytical jacobian of the normal fluxes
     _tJacob(0,1) = nx;
     _tJacob(0,2) = ny;

     _tJacob(1,0) = -u*un;
     _tJacob(1,1) = un + u*nx;
     _tJacob(1,2) = u*ny;
     _tJacob(1,4) = -Bn-Bx*nx;
     _tJacob(1,5) = -Bx*ny;

     _tJacob(2,0) = -v*un;
     _tJacob(2,1) = v*nx;
     _tJacob(2,2) = un + v*ny;
     _tJacob(2,4) = -By*nx;
     _tJacob(2,5) = -Bn-By*ny;

     _tJacob(3,0) = -w*un;
     _tJacob(3,1) = w*nx;
     _tJacob(3,2) = w*ny;
     _tJacob(3,3) = un;
     _tJacob(3,4) = -Bz*nx;
     _tJacob(3,5) = -Bz*ny;
     _tJacob(3,6) = -Bn;

     _tJacob(4,0) = -un*Bx/rho + u*Bn/rho;
     _tJacob(4,1) = -By/rho*ny;
     _tJacob(4,2) = Bx/rho*ny;
     _tJacob(4,4) = v*ny;
     _tJacob(4,5) = -u*ny;
//     _tJacob(4,8) = nx;

     _tJacob(5,0) = -un*By/rho + v*Bn/rho;
     _tJacob(5,1) = By/rho*nx;
     _tJacob(5,2) = -Bx/rho*nx;
     _tJacob(5,4) = -v*nx;
     _tJacob(5,5) = u*nx;
//     _tJacob(5,8) = ny;

     _tJacob(6,0) = -un*Bz/rho + w*Bn/rho;
     _tJacob(6,1) = Bz/rho*nx;
     _tJacob(6,2) = Bz/rho*ny;
     _tJacob(6,3) = -Bn/rho;
     _tJacob(6,4) = -w*nx;
     _tJacob(6,5) = -w*ny;
     _tJacob(6,6) = un;

     _tJacob(7,0) = -a2*un/(gamma - 1.) - (2. - gamma)/2.*V2*un - 0.5*un*B2/rho + Bn/rho*VB;
     _tJacob(7,1) = (a2/(gamma - 1.) + V2/2. + B2/rho)*nx - (gamma - 1.)*u*un - Bx/rho*Bn;
     _tJacob(7,2) = (a2/(gamma - 1.) + V2/2. + B2/rho)*ny - (gamma - 1.)*v*un - By/rho*Bn;
     _tJacob(7,3) = -(gamma - 1.)*w*un - Bz/rho*Bn;
     _tJacob(7,4) = (2. - gamma)*Bx*un - VB*nx - u*Bn;
     _tJacob(7,5) = (2. - gamma)*By*un - VB*ny - v*Bn;
     _tJacob(7,6) = (2. - gamma)*Bz*un - w*Bn;
     _tJacob(7,7) = un*gamma;

     _tJacob(8,4) = beta*beta*nx;
     _tJacob(8,5) = beta*beta*ny;

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
