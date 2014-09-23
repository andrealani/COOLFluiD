#include "NavierStokes/NavierStokes.hh"
#include "Euler2DConsToCharInRef.hh"
#include "EulerPhysicalModel.hh"
#include "Framework/PhysicalModel.hh"
#include "Environment/ObjectProvider.hh"
#include "EulerPhysicalModel.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace NavierStokes {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<Euler2DConsToCharInRef, VarSetMatrixTransformer, NavierStokesModule, 1> euler2DConsToCharInRefProvider("Euler2DConsToCharInRef");

//////////////////////////////////////////////////////////////////////////////

Euler2DConsToCharInRef::Euler2DConsToCharInRef(Common::SafePtr<Framework::PhysicalModelImpl> model) :
  VarSetMatrixTransformer(model),
  _model(model->getConvectiveTerm().d_castTo<EulerTerm>()),
  _temp(model->getNbEquations(), model->getNbEquations())
{
}

//////////////////////////////////////////////////////////////////////////////

Euler2DConsToCharInRef::~Euler2DConsToCharInRef()
{
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DConsToCharInRef::setMatrixFromRef()
{    
  const RealVector& linearData = _model->getPhysicalData();
  
  const CFreal rho = linearData[EulerTerm::RHO];
  const CFreal u   = linearData[EulerTerm::VX];
  const CFreal v   = linearData[EulerTerm::VY];
  const CFreal speed = linearData[EulerTerm::V];
  const CFreal cosD = u/speed;
  const CFreal sinD = v/speed;
  const CFreal gamma = _model->getGamma();
  const CFreal gammaMinus1 = gamma - 1.;
  const CFreal a = linearData[EulerTerm::A];
  const CFreal rhoDivA = rho/a;
  const CFreal invA2   = 1./(a*a);
  const CFreal M     = speed/a;
  const CFreal M2    = M*M;
  const CFreal invM  = 1./M;
  const CFreal eps   = 0.05;
  const CFreal beta  = sqrt(max(eps*eps,std::abs(M2 - 1.)));
  const CFreal chi   = beta/max(M,1.);
  const CFreal k1  = 0.5*beta/(chi*M2);
  const CFreal k2  = k1*rhoDivA;
  const CFreal k3  = 0.5*rho/(chi*M);
  const CFreal k4  = rhoDivA*speed*speed*0.5 + rho*a/gammaMinus1;
  const CFreal k5  = -u*sinD + v*cosD;
  const CFreal k6  = u*cosD + v*sinD;
  
  _temp(0,0) = k2;
  _temp(0,1) = k2;
  _temp(0,2) = rhoDivA/M2;
  _temp(0,3) = -invA2;
  _temp(1,0) = k2*u - k3*sinD;
  _temp(1,1) = k2*u + k3*sinD;
  _temp(1,2) = rhoDivA*u/M2 + invM*rho*cosD;
  _temp(1,3) = -invA2*u;
  _temp(2,0) = k2*v + k3*cosD;
  _temp(2,1) = k2*v - k3*cosD;
  _temp(2,2) = rhoDivA*v/M2 + invM*rho*sinD;
  _temp(2,3) = -invA2*v;
  _temp(3,0) = k1*k4 + k3*k5;
  _temp(3,1) = k1*k4 - k3*k5;
  _temp(3,2) = k4/M2 + k6*rho*invM;
  _temp(3,3) = -0.5*M2;

  m_inverter4.invert(_temp,_transMatrix);
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace NavierStokes

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
