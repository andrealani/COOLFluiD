#include "NavierStokes/NavierStokes.hh"
#include "Euler3DPrimToRoe.hh"
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

Environment::ObjectProvider<Euler3DPrimToRoe, VarSetTransformer, NavierStokesModule, 1> euler3DPrimToRoeProvider("Euler3DPrimToRoe");

//////////////////////////////////////////////////////////////////////////////

Euler3DPrimToRoe::Euler3DPrimToRoe(Common::SafePtr<Framework::PhysicalModelImpl> model) :
  VarSetTransformer(model),
  _model(model->getConvectiveTerm().d_castTo<EulerTerm>())
{
}

//////////////////////////////////////////////////////////////////////////////

Euler3DPrimToRoe::~Euler3DPrimToRoe()
{
}

//////////////////////////////////////////////////////////////////////////////

void Euler3DPrimToRoe::transform(const State& state, State& result)
{
  const CFreal rho  = state[0];
  const CFreal rhoU = rho*state[1];
  const CFreal rhoV = rho*state[2];
  const CFreal rhoW = rho*state[3];
  const CFreal rhoV2 = rhoU*rhoU + rhoV*rhoV + rhoW*rhoW;
  const CFreal gamma = _model->getGamma();
  const CFreal gammaMinus1 = gamma - 1.;
  const CFreal rhoE = _model->getPressureFromState(state[4])/
    gammaMinus1 + 0.5*rhoV2/rho;
  
  result[0] = sqrt(rho);
  result[1] = rhoU/result[0];
  result[2] = rhoV/result[0];
  result[3] = rhoW/result[0];
  result[4] = (gamma*rhoE - 0.5*gammaMinus1*rhoV2/rho)/result[0];
}

//////////////////////////////////////////////////////////////////////////////

void Euler3DPrimToRoe::transformFromRef(const RealVector& data, State& result)
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
