#include "NavierStokes/NavierStokes.hh"
#include "Euler3DConsToRoeInCons.hh"
#include "EulerPhysicalModel.hh"
#include "Framework/PhysicalModel.hh"
#include "Environment/ObjectProvider.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace NavierStokes {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<Euler3DConsToRoeInCons, VarSetMatrixTransformer, NavierStokesModule, 1> euler3DConsToRoeInConsProvider("Euler3DConsToRoeInCons");

//////////////////////////////////////////////////////////////////////////////

Euler3DConsToRoeInCons::Euler3DConsToRoeInCons(Common::SafePtr<Framework::PhysicalModelImpl> model) :
  VarSetMatrixTransformer(model),
   _model(model->getConvectiveTerm().d_castTo<EulerTerm>())
{
}

//////////////////////////////////////////////////////////////////////////////

Euler3DConsToRoeInCons::~Euler3DConsToRoeInCons()
{
}

//////////////////////////////////////////////////////////////////////////////

void Euler3DConsToRoeInCons::setMatrix(const RealVector& state)
{
  const CFreal rho  = state[0];
  const CFreal rhoU = state[1];
  const CFreal rhoV = state[2];
  const CFreal rhoW = state[3];
  const CFreal rhoE = state[4];

  (*_transState)[0] = sqrt(rho);
  const CFreal invZ0 = 1./(*_transState)[0];
  const CFreal halfInvZ0Z0 = 0.5*invZ0*invZ0;

  const CFreal gamma = _model->getGamma();
  const CFreal gammaMinus1 = gamma - 1.;

  (*_transState)[1] = rhoU*invZ0;
  (*_transState)[2] = rhoV*invZ0;
  (*_transState)[3] = rhoW*invZ0;
  (*_transState)[4] = (gamma*rhoE -
		    0.5*gammaMinus1*
		    (rhoU*rhoU + rhoV*rhoV + rhoW*rhoW)/rho)*invZ0;

  const CFreal z12sq = (*_transState)[1]*(*_transState)[1] +
    (*_transState)[2]*(*_transState)[2] + (*_transState)[3]*(*_transState)[3];

  _transMatrix(0,0) = 0.5*invZ0;
  _transMatrix(1,0) = -halfInvZ0Z0*(*_transState)[1];
  _transMatrix(1,1) = invZ0;
  _transMatrix(2,0) = -halfInvZ0Z0*(*_transState)[2];
  _transMatrix(2,2) = invZ0;
  _transMatrix(3,0) = -halfInvZ0Z0*(*_transState)[3];
  _transMatrix(3,3) = invZ0;
  _transMatrix(4,0) = halfInvZ0Z0*(gammaMinus1*z12sq*invZ0 -(*_transState)[4]);
  _transMatrix(4,1) = -invZ0*invZ0*gammaMinus1*(*_transState)[1];
  _transMatrix(4,2) = -invZ0*invZ0*gammaMinus1*(*_transState)[2];
  _transMatrix(4,3) = -invZ0*invZ0*gammaMinus1*(*_transState)[3];
  _transMatrix(4,4) = gamma*invZ0;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace NavierStokes

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
