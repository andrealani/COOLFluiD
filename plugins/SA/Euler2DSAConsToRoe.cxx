#include "SA/SA.hh"
#include "Euler2DSAConsToRoe.hh"
#include "Framework/PhysicalModel.hh"
#include "Environment/ObjectProvider.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace SA {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<Euler2DSAConsToRoe, VarSetTransformer, SAModule, 1>
euler2DSAConsToRoeProvider("Euler2DSAConsToRoe");

//////////////////////////////////////////////////////////////////////////////

Euler2DSAConsToRoe::Euler2DSAConsToRoe(Common::SafePtr<Framework::PhysicalModelImpl> model) :
 VarSetTransformer(model),
 _model(model->getConvectiveTerm().d_castTo<EulerSATerm>())
{
}

//////////////////////////////////////////////////////////////////////////////

Euler2DSAConsToRoe::~Euler2DSAConsToRoe()
{
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DSAConsToRoe::transform(const State& state, State& result)
{
  const CFreal rho  = state[0];
  const CFreal rhoU = state[1];
  const CFreal rhoV = state[2];
  const CFreal rhoE = state[3];
  const CFreal rhoK = state[4];
  const CFreal gamma = _model->getGamma();
  const CFreal gammaMinus1 = gamma - 1.;

  result[0] = sqrt(rho);

  const CFreal overSqrtRho = 1./result[0];

  result[1] = rhoU*overSqrtRho;
  result[2] = rhoV*overSqrtRho;
  result[3] = (gamma*rhoE - 0.5*gammaMinus1*
		    (rhoU*rhoU+rhoV*rhoV)/rho)*overSqrtRho;

  result[4] = rhoK*overSqrtRho;
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DSAConsToRoe::transformFromRef(const RealVector& data, State& result)
{
  const CFreal sqRho = sqrt(data[EulerSATerm::RHO]);

  result[0] = sqRho;
  result[1] = sqRho*data[EulerSATerm::VX];
  result[2] = sqRho*data[EulerSATerm::VY];
  result[3] = sqRho*data[EulerSATerm::H];

  const CFuint iK = _model->getFirstScalarVar(0);
  result[4] = sqRho*data[iK];
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace SA

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
