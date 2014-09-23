#include "NavierStokes/NavierStokes.hh"
#include "Euler2DRoeToConsInCons.hh"
#include "EulerPhysicalModel.hh"
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

Environment::ObjectProvider<Euler2DRoeToConsInCons, VarSetMatrixTransformer, NavierStokesModule, 1> euler2DRoeToConsInConsProvider("Euler2DRoeToConsInCons");

//////////////////////////////////////////////////////////////////////////////

Euler2DRoeToConsInCons::Euler2DRoeToConsInCons(Common::SafePtr<Framework::PhysicalModelImpl> model) :
  VarSetMatrixTransformer(model),
  _model(model->getConvectiveTerm().d_castTo<EulerTerm>())
{
}

//////////////////////////////////////////////////////////////////////////////

Euler2DRoeToConsInCons::~Euler2DRoeToConsInCons()
{
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DRoeToConsInCons::setMatrix(const RealVector& state)
{
  cf_assert(_model.isNotNull());

  const CFreal gamma = _model->getGamma();
  const CFreal gammaMinus1 = gamma - 1.;
  const CFreal gammaMinus1DivGamma = gammaMinus1/gamma;
  const CFreal rho  = state[0];
  const CFreal rhoU = state[1];
  const CFreal rhoV = state[2];
  const CFreal rhoE = state[3];
  const CFreal rhoV2 = (rhoU*rhoU + rhoV*rhoV)/rho;
  const CFreal z0 = sqrt(rho);
  const CFreal z1 = rhoU/z0;
  const CFreal z2 = rhoV/z0;
  const CFreal z3 = (gamma*rhoE - 0.5*gammaMinus1*
		     rhoV2)/z0;

  _transMatrix(0,0) = 2.0*z0;
  _transMatrix(1,0) = z1;
  _transMatrix(1,1) = z0;
  _transMatrix(2,0) = z2;
  _transMatrix(2,2) = z0;
  _transMatrix(3,0) = z3/gamma;
  _transMatrix(3,1) = z1*gammaMinus1DivGamma;
  _transMatrix(3,2) = z2*gammaMinus1DivGamma;
  _transMatrix(3,3) = z0/gamma;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace NavierStokes

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
