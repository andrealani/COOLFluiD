#include "NavierStokes/NavierStokes.hh"
#include "Euler2DPuvtToRoe.hh"
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

Environment::ObjectProvider<Euler2DPuvtToRoe, VarSetTransformer, NavierStokesModule, 1> euler2DPuvtToRoeProvider("Euler2DPuvtToRoe");

//////////////////////////////////////////////////////////////////////////////

Euler2DPuvtToRoe::Euler2DPuvtToRoe(Common::SafePtr<Framework::PhysicalModelImpl> model) :
  VarSetTransformer(model),
 _model(model->getConvectiveTerm().d_castTo<EulerTerm>())
{
}

//////////////////////////////////////////////////////////////////////////////

Euler2DPuvtToRoe::~Euler2DPuvtToRoe()
{
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DPuvtToRoe::transform(const State& state, State& result)
{
  const CFreal p = _model->getPressureFromState(state[0]);
  const CFreal rho  = _model->getDensity(p,state[3]);
  const CFreal rhoU = rho*state[1];
  const CFreal rhoV = rho*state[2];
  const CFreal gamma = _model->getGamma();
  const CFreal gammaMinus1 = gamma - 1.;
  const CFreal rhoE = p/gammaMinus1 + 0.5*(rhoU*rhoU + rhoV*rhoV)/rho;
  
  result[0] = sqrt(rho);
  result[1] = rhoU/result[0];
  result[2] = rhoV/result[0];
  result[3] = (gamma*rhoE - 0.5*gammaMinus1*(rhoU*rhoU+rhoV*rhoV)/rho)/result[0];
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DPuvtToRoe::transformFromRef(const RealVector& data, State& result)
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
