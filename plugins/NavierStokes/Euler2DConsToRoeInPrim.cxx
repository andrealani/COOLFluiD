#include "NavierStokes/NavierStokes.hh"
#include "Euler2DConsToRoeInPrim.hh"
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

Environment::ObjectProvider<Euler2DConsToRoeInPrim, VarSetMatrixTransformer, NavierStokesModule, 1> euler2DConsToRoeInPrimProvider("Euler2DConsToRoeInPrim");

//////////////////////////////////////////////////////////////////////////////

Euler2DConsToRoeInPrim::Euler2DConsToRoeInPrim(Common::SafePtr<Framework::PhysicalModelImpl> model) :
  VarSetMatrixTransformer(model),
  _model(model->getConvectiveTerm().d_castTo<EulerTerm>())
{
}

//////////////////////////////////////////////////////////////////////////////

Euler2DConsToRoeInPrim::~Euler2DConsToRoeInPrim()
{
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DConsToRoeInPrim::setMatrix(const RealVector& state)
{
  // AL: not sure this works for incompressible
  const CFreal gammaMinus1 = _model->getGamma() - 1.;
  const CFreal rho  = state[0];
  const CFreal rhoU = rho*state[1];
  const CFreal rhoV = rho*state[2];
  const CFreal p = _model->getPressureFromState(state[3]);
  const CFreal rhoE = p/gammaMinus1 + 0.5*(rhoU*rhoU + rhoV*rhoV)/rho;
  const CFreal z0 = sqrt(rho);
  const CFreal invZ0 = 1./z0;
  const CFreal halfInvZ0Z0 = 0.5*invZ0*invZ0;
  
  (*_transState)[0] = z0;
  (*_transState)[1] = rhoU*invZ0;
  (*_transState)[2] = rhoV*invZ0;
  (*_transState)[3] = (rhoE + p)*invZ0;
  
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
