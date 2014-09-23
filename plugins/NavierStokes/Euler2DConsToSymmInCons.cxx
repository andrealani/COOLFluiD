#include "NavierStokes/NavierStokes.hh"
#include "Euler2DConsToSymmInCons.hh"
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

Environment::ObjectProvider<Euler2DConsToSymmInCons, VarSetMatrixTransformer, NavierStokesModule, 1> euler2DConsToSymmInConsProvider("Euler2DConsToSymmInCons");

//////////////////////////////////////////////////////////////////////////////

Euler2DConsToSymmInCons::Euler2DConsToSymmInCons(Common::SafePtr<Framework::PhysicalModelImpl> model) :
  VarSetMatrixTransformer(model),
  _model(model->getConvectiveTerm().d_castTo<EulerTerm>()),
  _temp(model->getNbEquations(), model->getNbEquations())
{
}

//////////////////////////////////////////////////////////////////////////////

Euler2DConsToSymmInCons::~Euler2DConsToSymmInCons()
{
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DConsToSymmInCons::setMatrix(const RealVector& state)
{
  const CFreal rho = state[0];
  const CFreal u   = state[1]/rho;
  const CFreal v   = state[2]/rho;
  const CFreal speed = sqrt(u*u + v*v);
  const CFreal cosD = u/speed;
  const CFreal sinD = v/speed;
  const CFreal gamma = _model->getGamma();
  const CFreal gammaMinus1 = gamma - 1.;
  const CFreal p = gammaMinus1*(state[3] - 0.5*rho*speed*speed);
  const CFreal a = sqrt(gamma*p/rho);
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
  _temp(3,1) = rho*u*cosD + rho*v*sinD;
  _temp(3,2) = -rho*u*sinD + rho*v*cosD;
  _temp(3,3) = -0.5*speed*speed*invA2;

  m_inverter4.invert(_temp,_transMatrix);
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace NavierStokes

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
