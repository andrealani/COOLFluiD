#include "NavierStokes/NavierStokes.hh"
#include "Euler1DPvtToRoe.hh"
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

Environment::ObjectProvider<Euler1DPvtToRoe, VarSetTransformer, NavierStokesModule, 1> euler1DPuvtToRoeProvider("Euler1DPvtToRoe");

//////////////////////////////////////////////////////////////////////////////

Euler1DPvtToRoe::Euler1DPvtToRoe(Common::SafePtr<Framework::PhysicalModelImpl> model) :
  VarSetTransformer(model),
 _model(model->getConvectiveTerm().d_castTo<EulerTerm>())
{
}

//////////////////////////////////////////////////////////////////////////////

Euler1DPvtToRoe::~Euler1DPvtToRoe()
{
}

//////////////////////////////////////////////////////////////////////////////

void Euler1DPvtToRoe::transform(const State& state, State& result)
{
  const CFreal p = _model->getPressureFromState(state[0]);
  const CFreal rho  = _model->getDensity(p,state[2]);
  const CFreal rhoU = rho*state[1];
  const CFreal gamma = _model->getGamma();
  const CFreal gammaMinus1 = gamma - 1.;
  const CFreal rhoE = p/gammaMinus1 + 0.5*rhoU*rhoU/rho;
  
  result[0] = sqrt(rho);
  result[1] = rhoU/result[0];
  result[2] = (gamma*rhoE - 0.5*gammaMinus1*rhoU*rhoU/rho)/result[0];

}

//////////////////////////////////////////////////////////////////////////////

void Euler1DPvtToRoe::transformFromRef(const RealVector& data, State& result)
{
  const CFreal sqRho = sqrt(data[EulerTerm::RHO]);

  result[0] = sqRho;
  result[1] = sqRho*data[EulerTerm::VX];
  result[2] = sqRho*data[EulerTerm::H];

}

//////////////////////////////////////////////////////////////////////////////

    } // namespace NavierStokes

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
