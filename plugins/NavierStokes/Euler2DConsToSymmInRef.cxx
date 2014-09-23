#include "NavierStokes/NavierStokes.hh"
#include "Euler2DConsToSymmInRef.hh"
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

Environment::ObjectProvider<Euler2DConsToSymmInRef, VarSetMatrixTransformer, NavierStokesModule, 1> euler2DConsToSymmInRefProvider("Euler2DConsToSymmInRef");

//////////////////////////////////////////////////////////////////////////////

Euler2DConsToSymmInRef::Euler2DConsToSymmInRef(Common::SafePtr<Framework::PhysicalModelImpl> model) :
  VarSetMatrixTransformer(model),
  _model(model->getConvectiveTerm().d_castTo<EulerTerm>()),
  _temp(model->getNbEquations(), model->getNbEquations())
{
}

//////////////////////////////////////////////////////////////////////////////

Euler2DConsToSymmInRef::~Euler2DConsToSymmInRef()
{
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DConsToSymmInRef::setMatrixFromRef()
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
  
  _temp(0,0) = rhoDivA;
  _temp(0,3) = -invA2;
  _temp(1,0) = rhoDivA*u;
  _temp(1,1) = rho*cosD;
  _temp(1,2) = -rho*sinD;
  _temp(1,3) = -invA2*u;
  _temp(2,0) = rhoDivA*v;
  _temp(2,1) = rho*sinD;
  _temp(2,2) = rho*cosD;
  _temp(2,3) = -invA2*v;
  _temp(3,0) = rhoDivA*speed*speed*0.5 + rho*a/gammaMinus1;
  _temp(3,1) = rho*(u*cosD + v*sinD);
  _temp(3,2) = -rho*(u*sinD + v*cosD);
  _temp(3,3) = -0.5*speed*speed*invA2;
  
  m_inverter4.invert(_temp,_transMatrix);
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace NavierStokes

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
