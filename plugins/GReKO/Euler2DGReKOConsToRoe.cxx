#include "GReKO.hh"
#include "Euler2DGReKOConsToRoe.hh"
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

    namespace GReKO {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<Euler2DGReKOConsToRoe, VarSetTransformer, GReKOModule, 1>
euler2DGReKOConsToRoeProvider("Euler2DGReKOConsToRoe");

//////////////////////////////////////////////////////////////////////////////

Euler2DGReKOConsToRoe::Euler2DGReKOConsToRoe(Common::SafePtr<Framework::PhysicalModelImpl> model) :
 VarSetTransformer(model),
 _model(model->getConvectiveTerm().d_castTo<EulerGReKOTerm>())
{
}

//////////////////////////////////////////////////////////////////////////////

Euler2DGReKOConsToRoe::~Euler2DGReKOConsToRoe()
{
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DGReKOConsToRoe::transform(const State& state, State& result)
{
  const CFreal rho  = state[0];
  const CFreal rhoU = state[1];
  const CFreal rhoV = state[2];
  const CFreal rhoE = state[3];
  const CFreal rhoK = state[4];
  const CFreal rhoOmega = state[5];
  const CFreal rhoGa = state[6];
  const CFreal rhoRe = state[7];
  const CFreal gamma = _model->getGamma();
  const CFreal gammaMinus1 = gamma - 1.;

  result[0] = sqrt(rho);
  result[1] = rhoU/result[0];
  result[2] = rhoV/result[0];
  result[3] = (gamma*rhoE - 0.5*gammaMinus1*
		    (rhoU*rhoU+rhoV*rhoV)/rho)/result[0];

  result[4] = rhoK/result[0];
  result[5] = rhoOmega/result[0];
  result[6] = rhoGa/result[0];
  result[7] = rhoRe/result[0];
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DGReKOConsToRoe::transformFromRef(const RealVector& data, State& result)
{
  const CFreal sqRho = sqrt(data[EulerGReKOTerm::RHO]);

  result[0] = sqRho;
  result[1] = sqRho*data[EulerGReKOTerm::VX];
  result[2] = sqRho*data[EulerGReKOTerm::VY];
  result[3] = sqRho*data[EulerGReKOTerm::H];

  const CFuint iK = _model->getFirstScalarVar(0);
  result[4] = sqRho*data[iK];
  result[5] = sqRho*data[iK+1];
  result[6] = sqRho*data[iK+2];
  result[7] = sqRho*data[iK+3];
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace GReKO

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
