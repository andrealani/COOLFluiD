#include "LTE.hh"
#include "Euler2DPrvtLTEToConsInPrvtLTE.hh"
#include "Framework/PhysicalModel.hh"
#include "Environment/ObjectProvider.hh"
#include "NavierStokes/EulerTerm.hh"
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

Environment::ObjectProvider<Euler2DPrvtLTEToConsInPrvtLTE, VarSetMatrixTransformer, LTEModule, 1> euler2DPrvtLTEToConsInPrvtLTEProvider("Euler2DPrvtLTEToConsInPrvtLTE");

//////////////////////////////////////////////////////////////////////////////

Euler2DPrvtLTEToConsInPrvtLTE::Euler2DPrvtLTEToConsInPrvtLTE(Common::SafePtr<Framework::PhysicalModelImpl> model) :
  VarSetMatrixTransformer(model),
  _model(model->getConvectiveTerm().d_castTo<EulerTerm>()),
  _tol(10e-4),
  _dhe(3)
{
}

//////////////////////////////////////////////////////////////////////////////

Euler2DPrvtLTEToConsInPrvtLTE::~Euler2DPrvtLTEToConsInPrvtLTE()
{
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DPrvtLTEToConsInPrvtLTE::setMatrix(const RealVector& state)
{
  // find a way of storing this pointer (sdd setup() function)
  static Common::SafePtr<PhysicalChemicalLibrary> library =
    PhysicalModelStack::getActive()->getImplementor()->getPhysicalPropertyLibrary<PhysicalChemicalLibrary>();

  const RealVector& refData = _model->getReferencePhysicalData();

  // set reference density and internal energy
  const CFreal p = _model->getPressureFromState(state[0]);
  const CFreal T = state[3];
  CFreal pdim = p*refData[EulerTerm::P];
  CFreal Tdim = T*_model->getTempRef();
  library->setComposition(Tdim, pdim);
  library->setDensityEnthalpyEnergy(Tdim,pdim,_dhe);
  const CFreal rho = _dhe[0]/refData[EulerTerm::RHO];
  const CFreal e = _dhe[2]/refData[EulerTerm::H];
  
  // perturb the pressure and compute derivatives for rho and e
  const CFreal dp = _tol*p;
  const CFreal p1 = p + dp;
  pdim = p1*refData[EulerTerm::P]; // Tdim unchanges so far, set new pdim
  // this additional call to set the composition is useless
  // (result is the same even without and run-time speed a bit better)
  library->setComposition(Tdim, pdim);
  library->setDensityEnthalpyEnergy(Tdim,pdim,_dhe);
  CFreal rho1 = _dhe[0]/refData[EulerTerm::RHO];
  CFreal e1 = _dhe[2]/refData[EulerTerm::H];
  const CFreal drdp = (rho1 - rho)/dp;
  const CFreal dedp = (e1 - e)/dp;

  pdim = p*refData[EulerTerm::P]; // restore original pdim
  // perturb the temperature and compute derivatives for rho and e
  const CFreal dT = _tol*T;
  const CFreal T1 = T + dT;
  Tdim = T1*_model->getTempRef(); // set new Tdim
  // this additional call to set the composition is useless
  // (result is the same even without and run-time speed a bit better)
  library->setComposition(Tdim, pdim);
  library->setDensityEnthalpyEnergy(Tdim,pdim,_dhe);
  rho1 = _dhe[0]/refData[EulerTerm::RHO];
  e1 = _dhe[2]/refData[EulerTerm::H];
  const CFreal drdT = (rho1 - rho)/dT;
  const CFreal dedT = (e1 - e)/dT;
  const CFreal u = state[1]/rho;
  const CFreal v = state[2]/rho;
  const CFreal V2 = u*u + v*v;
//unused//  const CFreal E = e + 0.5*V2;

  _transMatrix(0,0) = drdp;
  _transMatrix(0,3) = drdT;

  _transMatrix(1,1) = 1.0;

  _transMatrix(2,2) = 1.0;

  _transMatrix(3,0) = drdp*(e - 0.5*V2) + rho*dedp;
  _transMatrix(3,1) = u;
  _transMatrix(3,2) = v;
  _transMatrix(3,3) = drdT*(e - 0.5*V2) + rho*dedT;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace LTE

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
