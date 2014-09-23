#include "NavierStokes/NavierStokes.hh"
#include "Euler2DConsToRoeInCons.hh"
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

Environment::ObjectProvider<Euler2DConsToRoeInCons, VarSetMatrixTransformer, NavierStokesModule, 1> euler2DConsToRoeInConsProvider("Euler2DConsToRoeInCons");

//////////////////////////////////////////////////////////////////////////////

Euler2DConsToRoeInCons::Euler2DConsToRoeInCons(Common::SafePtr<Framework::PhysicalModelImpl> model) :
  VarSetMatrixTransformer(model),
 _model(model->getConvectiveTerm().d_castTo<EulerTerm>())
{
}

//////////////////////////////////////////////////////////////////////////////

Euler2DConsToRoeInCons::~Euler2DConsToRoeInCons()
{
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DConsToRoeInCons::setMatrix(const RealVector& state)
{
  const CFreal rho  = state[0];
  const CFreal rhoU = state[1];
  const CFreal rhoV = state[2];
  const CFreal rhoE = state[3];
  const CFreal gamma = _model->getGamma();
  const CFreal gammaMinus1 = gamma - 1.;

  (*_transState)[0] = sqrt(rho);
  const CFreal invZ0 = 1./(*_transState)[0];
  const CFreal halfInvZ0Z0 = 0.5*invZ0*invZ0;

  (*_transState)[1] = rhoU*invZ0;
  (*_transState)[2] = rhoV*invZ0;
  (*_transState)[3] = (gamma*rhoE -
                    0.5*gammaMinus1*
                    (rhoU*rhoU+rhoV*rhoV)/rho)*invZ0;

  const CFreal z12sq = (*_transState)[1]*(*_transState)[1] +
    (*_transState)[2]*(*_transState)[2];

  _transMatrix(0,0) = 0.5*invZ0;
  _transMatrix(1,0) = -halfInvZ0Z0*(*_transState)[1];
  _transMatrix(1,1) = invZ0;
  _transMatrix(2,0) = -halfInvZ0Z0*(*_transState)[2];
  _transMatrix(2,2) = invZ0;
  _transMatrix(3,0) = halfInvZ0Z0*(gammaMinus1*z12sq*invZ0 -(*_transState)[3]);
  _transMatrix(3,1) = -invZ0*invZ0*gammaMinus1*(*_transState)[1];
  _transMatrix(3,2) = -invZ0*invZ0*gammaMinus1*(*_transState)[2];
  _transMatrix(3,3) = _model->getGamma()*invZ0;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace NavierStokes

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
