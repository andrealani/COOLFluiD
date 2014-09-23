#include "KOmega.hh"
#include "Euler3DKOmegaConsToRoe.hh"
#include "Framework/PhysicalModel.hh"
#include "Environment/ObjectProvider.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Physics::NavierStokes;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace KOmega {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<Euler3DKOmegaConsToRoe, VarSetTransformer, KOmegaModule, 1>
euler3DKOmegaConsToRoeProvider("Euler3DKOmegaConsToRoe");

//////////////////////////////////////////////////////////////////////////////

Euler3DKOmegaConsToRoe::Euler3DKOmegaConsToRoe(Common::SafePtr<Framework::PhysicalModelImpl> model) :
 VarSetTransformer(model),
 _model(model->getConvectiveTerm().d_castTo<EulerKOmegaTerm>())
{
}

//////////////////////////////////////////////////////////////////////////////

Euler3DKOmegaConsToRoe::~Euler3DKOmegaConsToRoe()
{
}

//////////////////////////////////////////////////////////////////////////////

void Euler3DKOmegaConsToRoe::transform(const State& state, State& result)
{
  const CFreal rho  = state[0];
  const CFreal rhoU = state[1];
  const CFreal rhoV = state[2];
  const CFreal rhoW = state[3];
  const CFreal rhoE = state[4];
  const CFreal rhoK = state[5];
  const CFreal rhoOmega = state[6];
  const CFreal gamma = _model->getGamma();
  const CFreal gammaMinus1 = gamma - 1.;

  result[0] = sqrt(rho);
  result[1] = rhoU/result[0];
  result[2] = rhoV/result[0];
  result[3] = rhoW/result[0];
  result[4] = (gamma*rhoE - 0.5*gammaMinus1*
              (rhoU*rhoU + rhoV*rhoV + rhoW*rhoW)/rho)/result[0];
  result[5] = rhoK/result[0];
  result[6] = rhoOmega/result[0];
}

//////////////////////////////////////////////////////////////////////////////

void Euler3DKOmegaConsToRoe::transformFromRef(const RealVector& data, State& result)
{
  const CFreal sqRho = sqrt(data[EulerKOmegaTerm::RHO]);
  
  result[0] = sqRho;
  result[1] = sqRho*data[EulerKOmegaTerm::VX];
  result[2] = sqRho*data[EulerKOmegaTerm::VY];
  result[3] = sqRho*data[EulerKOmegaTerm::VZ];
  result[4] = sqRho*data[EulerKOmegaTerm::H];

  const CFuint iK = _model->getFirstScalarVar(0);
  result[5] = sqRho*data[iK];
  result[6] = sqRho*data[iK+1];
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace KOmega

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
