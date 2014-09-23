#include "NavierStokes/NavierStokes.hh"
#include "Euler2DPrimToRoe.hh"
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

Environment::ObjectProvider<Euler2DPrimToRoe, VarSetTransformer, NavierStokesModule, 1> euler2DPrimToRoeProvider("Euler2DPrimToRoe");

//////////////////////////////////////////////////////////////////////////////

Euler2DPrimToRoe::Euler2DPrimToRoe(Common::SafePtr<Framework::PhysicalModelImpl> model) :
  VarSetTransformer(model),
  _model(model->getConvectiveTerm().d_castTo<EulerTerm>())
{
}

//////////////////////////////////////////////////////////////////////////////

Euler2DPrimToRoe::~Euler2DPrimToRoe()
{
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DPrimToRoe::transform(const State& state, State& result)
{
  const CFreal rho  = state[0];
  const CFreal rhoU = rho*state[1];
  const CFreal rhoV = rho*state[2];
  const CFreal gamma = _model->getGamma();
  const CFreal gammaMinus1 = gamma - 1.;
  const CFreal rhoE = _model->getPressureFromState(state[3])/
    gammaMinus1 + 0.5*(rhoU*rhoU + rhoV*rhoV)/rho;
  
  result[0] = sqrt(rho);
  result[1] = rhoU/result[0];
  result[2] = rhoV/result[0];
  result[3] = (gamma*rhoE - 0.5*gammaMinus1*
	       (rhoU*rhoU+rhoV*rhoV)/rho)/result[0];
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DPrimToRoe::transformFromRef(const RealVector& data, State& result)
{
  const CFreal sqRho = sqrt(data[EulerTerm::RHO]);
  
  result[0] = sqRho;
  result[1] = sqRho*data[EulerTerm::VX];
  result[2] = sqRho*data[EulerTerm::VY];
  result[3] = sqRho*data[EulerTerm::H];
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace NavierStokes

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
