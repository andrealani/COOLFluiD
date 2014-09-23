#include "LTE.hh"
#include "Euler2DConsToPrvtLTEInRef.hh"
#include "Framework/PhysicalModel.hh"
#include "Environment/ObjectProvider.hh"
#include "Framework/PhysicalChemicalLibrary.hh"
#include "NavierStokes/EulerTerm.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Physics::NavierStokes;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace LTE {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<Euler2DConsToPrvtLTEInRef, VarSetMatrixTransformer, LTEModule, 1> euler2DConsToPrvtLTEInRefProvider("Euler2DConsToPrvtLTEInRef");

//////////////////////////////////////////////////////////////////////////////

Euler2DConsToPrvtLTEInRef::Euler2DConsToPrvtLTEInRef(Common::SafePtr<Framework::PhysicalModelImpl> model) :
  VarSetMatrixTransformer(model),
  _model(model->getConvectiveTerm().d_castTo<EulerTerm>()),
  _tol(10e-4),
  _dhe(3)
{
}

//////////////////////////////////////////////////////////////////////////////

Euler2DConsToPrvtLTEInRef::~Euler2DConsToPrvtLTEInRef()
{
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DConsToPrvtLTEInRef::setMatrixFromRef()
{
  // find a way of storing this pointer (sdd setup() function)
  static Common::SafePtr<PhysicalChemicalLibrary> library =
    PhysicalModelStack::getActive()->getImplementor()->getPhysicalPropertyLibrary<PhysicalChemicalLibrary>();

  const RealVector& phyData = _model->getPhysicalData();
  // set reference density and internal energy
  const CFreal p = phyData[EulerTerm::P];
  const CFreal T = phyData[EulerTerm::T];
  const CFreal u = phyData[EulerTerm::VX];
  const CFreal v = phyData[EulerTerm::VY];
  const CFreal V2 = phyData[EulerTerm::V]*phyData[EulerTerm::V];
  const CFreal rho = phyData[EulerTerm::RHO];
  const CFreal e = phyData[EulerTerm::E] - 0.5*V2;

  const RealVector& refData = _model->getReferencePhysicalData();
  const CFreal pref  = refData[EulerTerm::P];
  const CFreal Tref  = refData[EulerTerm::T];
  CFreal pdim = p*pref;
  CFreal Tdim = T*Tref;

  // perturb the pressure and compute derivatives for rho and e
  const CFreal dp = _tol*p;
  const CFreal p1 = p + dp;
  pdim = p1*pref; // Tdim unchanges so far, set new pdim
  library->setComposition(Tdim, pdim);
  library->setDensityEnthalpyEnergy(Tdim,pdim,_dhe);
  CFreal rho1 = _dhe[0]/refData[EulerTerm::RHO];
  CFreal e1 = _dhe[2]/refData[EulerTerm::H];
  const CFreal drdp = (rho1 - rho)/dp;
  const CFreal dedp = (e1 - e)/dp;
  
  pdim = p*pref; // restore original pdim
  // perturb the temperature and compute derivatives for rho and e
  const CFreal dT = _tol*T;
  const CFreal T1 = T + dT;
  Tdim = T1*Tref; // set new Tdim
  library->setComposition(Tdim, pdim);
  library->setDensityEnthalpyEnergy(Tdim,pdim,_dhe);
  rho1 = _dhe[0]/refData[EulerTerm::RHO];
  e1 = _dhe[2]/refData[EulerTerm::H];
  const CFreal drdT = (rho1 - rho)/dT;
  const CFreal dedT = (e1 - e)/dT;
  const CFreal kE = e - 0.5*V2;
  const CFreal ka = drdp*kE + rho*dedp;
  const CFreal kb = drdT*kE + rho*dedT;
  const CFreal kd = 1./(kb*drdp - ka*drdT);

  _transMatrix(0,0) = kb*kd;
  _transMatrix(0,1) = drdT*u*kd;
  _transMatrix(0,2) = drdT*v*kd;
  _transMatrix(0,3) = -drdT*kd;

  _transMatrix(1,1) = 1.0;

  _transMatrix(2,2) = 1.0;

  _transMatrix(3,0) = -ka*kd;
  _transMatrix(3,1) = -drdp*u*kd;
  _transMatrix(3,2) = -drdp*v*kd;
  _transMatrix(3,3) = drdp*kd;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace LTE

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
