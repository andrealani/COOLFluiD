#include "NavierStokes/NavierStokes.hh"
#include "Euler3DConsToRoe.hh"
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

Environment::ObjectProvider<Euler3DConsToRoe, VarSetTransformer, NavierStokesModule, 1> euler3DConsToRoeProvider("Euler3DConsToRoe");

//////////////////////////////////////////////////////////////////////////////

Euler3DConsToRoe::Euler3DConsToRoe(Common::SafePtr<Framework::PhysicalModelImpl> model) :
  VarSetTransformer(model),
 _model(model->getConvectiveTerm().d_castTo<EulerTerm>())
{
}

//////////////////////////////////////////////////////////////////////////////

Euler3DConsToRoe::~Euler3DConsToRoe()
{
}

//////////////////////////////////////////////////////////////////////////////

void Euler3DConsToRoe::transform(const State& state, State& result)
{
  const CFreal rho  = state[0];
  const CFreal rhoU = state[1];
  const CFreal rhoV = state[2];
  const CFreal rhoW = state[3];
  const CFreal rhoE = state[4];
  const CFreal gamma = _model->getGamma();
  const CFreal gammaMinus1 = gamma - 1.;

  result[0] = sqrt(rho);
  result[1] = rhoU/result[0];
  result[2] = rhoV/result[0];
  result[3] = rhoW/result[0];
  result[4] = (gamma*rhoE - 0.5*gammaMinus1*
		       (rhoU*rhoU + rhoV*rhoV + rhoW*rhoW)/rho)/
    result[0];
}

//////////////////////////////////////////////////////////////////////////////

void Euler3DConsToRoe::transformFromRef(const RealVector& data, State& result)
{
  const CFreal sqRho = sqrt(data[EulerTerm::RHO]);

  result[0] = sqRho;
  result[1] = sqRho*data[EulerTerm::VX];
  result[2] = sqRho*data[EulerTerm::VY];
  result[3] = sqRho*data[EulerTerm::VZ];
  result[4] = sqRho*data[EulerTerm::H];
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace NavierStokes

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
