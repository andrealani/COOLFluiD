#include "FluctSplit/FluctSplitLTE.hh"
#include "WeakSlipWallEuler2DPrvtLTEImpl.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/PhysicalChemicalLibrary.hh"
#include "NavierStokes/Euler2DVarSet.hh"

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

MethodCommandProvider<WeakSlipWallEuler2DPrvtLTEImpl, FluctuationSplitData, FluctSplitLTEModule> weakSlipWallEuler2DPrvtLTEImplProvider("WeakSlipWallEuler2DPrvtLTEImpl");

//////////////////////////////////////////////////////////////////////////////

WeakSlipWallEuler2DPrvtLTEImpl::WeakSlipWallEuler2DPrvtLTEImpl
(const std::string& name) :
  WeakSlipWallEuler2DImpl(name),
  _library(CFNULL),
  _tol(10e-4),
  _dhe(3)
{
}

//////////////////////////////////////////////////////////////////////////////

WeakSlipWallEuler2DPrvtLTEImpl::~WeakSlipWallEuler2DPrvtLTEImpl()
{
}

//////////////////////////////////////////////////////////////////////////////

void WeakSlipWallEuler2DPrvtLTEImpl::setup()
{
  _library = PhysicalModelStack::getActive()->getImplementor()->
    getPhysicalPropertyLibrary<PhysicalChemicalLibrary>();
  cf_assert(_library.isNotNull());
}

//////////////////////////////////////////////////////////////////////////////

void WeakSlipWallEuler2DPrvtLTEImpl::computeNormalFluxAndJacob
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
  const CFreal V2 = u*u + v*v;
  const CFreal H = _physicalData[EulerTerm::H];

  flux[0] = rho*un;
  flux[1] = un*rho*u;
  flux[2] = un*rho*v;
  flux[3] = un*rho*H;

  const RealVector& refData = _varSet->getModel()->getReferencePhysicalData();
  const CFreal p = _physicalData[EulerTerm::P];
  const CFreal T = _physicalData[EulerTerm::T];
  CFreal pdim = p*refData[EulerTerm::P];
  CFreal Tdim = T*refData[EulerTerm::T];
  const CFreal h = H - 0.5*V2;

  // perturb the pressure and compute derivatives for rho and e
  const CFreal dp = _tol*p;
  const CFreal p1 = p + dp;
  pdim = p1*refData[EulerTerm::P]; // Tdim unchanges so far, set new pdim
  _library->setComposition(Tdim, pdim);
  _library->setDensityEnthalpyEnergy(Tdim, pdim, _dhe);
  CFreal rho1 = _dhe[0]/refData[EulerTerm::RHO];
  CFreal h1 = _dhe[1]/refData[EulerTerm::H];
  const CFreal drdp = (rho1 - rho)/dp;
  const CFreal dhdp = (h1 - h)/dp;

  pdim = p*refData[EulerTerm::P]; // restore original pdim
  // perturb the temperature and compute derivatives for rho and e
  const CFreal dT = _tol*T;
  const CFreal T1 = T + dT;
  Tdim = T1*refData[EulerTerm::T]; // set new Tdim
  _library->setComposition(Tdim, pdim); 
  _library->setDensityEnthalpyEnergy(Tdim, pdim, _dhe);
  rho1 = _dhe[0]/refData[EulerTerm::RHO];
  h1 = _dhe[1]/refData[EulerTerm::H];
  const CFreal drdT = (rho1 - rho)/dT;
  const CFreal dhdT = (h1 - h)/dT;

  // the analytical jacobian of the normal fluxes
  fluxJacob(0,1) = nx;
  fluxJacob(0,2) = ny;

  fluxJacob(1,0) = -u*un*drdp;
  fluxJacob(1,1) = un + u*nx;
  fluxJacob(1,2) = u*ny;
  fluxJacob(1,3) = -u*un*drdT;

  fluxJacob(2,0) = -v*un*drdp;
  fluxJacob(2,1) = v*nx;
  fluxJacob(2,2) = un + v*ny;
  fluxJacob(2,3) = -v*un*drdT;

  fluxJacob(3,0) = rho*un*(dhdp - V2*drdp/rho);
  fluxJacob(3,1) = H*nx + u*un;
  fluxJacob(3,2) = H*ny + v*un;
  fluxJacob(3,3) = rho*un*(dhdT - V2*drdT/rho);
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
