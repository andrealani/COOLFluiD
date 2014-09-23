#include "LTE.hh"
#include "Euler2DPrvtLTEToConsInRef.hh"
#include "NavierStokes/EulerTerm.hh"
#include "Framework/PhysicalModel.hh"
#include "Environment/ObjectProvider.hh"
#include "Framework/PhysicalChemicalLibrary.hh"

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

Environment::ObjectProvider<Euler2DPrvtLTEToConsInRef, VarSetMatrixTransformer, LTEModule, 1> euler2DPrvtLTEToConsInRefProvider("Euler2DPrvtLTEToConsInRef");

//////////////////////////////////////////////////////////////////////////////

Euler2DPrvtLTEToConsInRef::Euler2DPrvtLTEToConsInRef(Common::SafePtr<Framework::PhysicalModelImpl> model) :
  VarSetMatrixTransformer(model),
  _model(model->getConvectiveTerm().d_castTo<EulerTerm>()),
  _tol(10e-4),
  _dhe(3)
{
}

//////////////////////////////////////////////////////////////////////////////

Euler2DPrvtLTEToConsInRef::~Euler2DPrvtLTEToConsInRef()
{
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DPrvtLTEToConsInRef::setMatrixFromRef()
{
  // find a way of storing this pointer (sdd setup() function)
  static Common::SafePtr<PhysicalChemicalLibrary> library =
    PhysicalModelStack::getActive()->getImplementor()->getPhysicalPropertyLibrary<PhysicalChemicalLibrary>();

  const RealVector& phyData = _model->getPhysicalData();
  // set reference density and internal energy
  const CFreal p = _model->getPressureFromState(phyData[EulerTerm::P]);
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

  _transMatrix(0,0) = drdp;
  _transMatrix(0,3) = drdT;

  _transMatrix(1,1) = 1.0;

  _transMatrix(2,2) = 1.0;

  _transMatrix(3,0) = drdp*kE + rho*dedp;
  _transMatrix(3,1) = u;
  _transMatrix(3,2) = v;
  _transMatrix(3,3) = drdT*kE + rho*dedT;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace LTE

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
