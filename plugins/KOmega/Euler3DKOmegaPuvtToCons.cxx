#include "KOmega.hh"
#include "Euler3DKOmegaPuvtToCons.hh"
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

Environment::ObjectProvider<Euler3DKOmegaPuvtToCons, VarSetTransformer, KOmegaModule, 1>
euler3DKOmegaPuvtToConsProvider("Euler3DKOmegaPuvtToCons");

//////////////////////////////////////////////////////////////////////////////

Euler3DKOmegaPuvtToCons::Euler3DKOmegaPuvtToCons(Common::SafePtr<Framework::PhysicalModelImpl> model) :
  VarSetTransformer(model),
  _model(model->getConvectiveTerm().
	 d_castTo<EulerKOmegaTerm>())
{
}

//////////////////////////////////////////////////////////////////////////////

Euler3DKOmegaPuvtToCons::~Euler3DKOmegaPuvtToCons()
{
}

//////////////////////////////////////////////////////////////////////////////

void Euler3DKOmegaPuvtToCons::transform(const State& state, State& result)
{
  const CFreal R = _model->getR();
  const CFreal p = state[0];
  const CFreal u = state[1];
  const CFreal v = state[2];
  const CFreal w = state[3];
  const CFreal T = state[4];
  const CFreal K = state[5];
  const CFreal Omega = state[6];
  const CFreal rho = p/(R*T);
  const CFreal rhoV2 = rho*(u*u + v*v + w*w);

  result[0] = rho;
  result[1] = rho*u;
  result[2] = rho*v;
  result[3] = rho*w;
  result[4] = p/(_model->getGamma() - 1.) + 0.5*rhoV2;
  result[5] = rho*K;
  result[6] = rho*Omega;
}

//////////////////////////////////////////////////////////////////////////////

void Euler3DKOmegaPuvtToCons::transformFromRef(const RealVector& data, State& result)
{
  const CFreal rho = data[EulerKOmegaTerm::RHO];

  result[0] = rho;
  result[1] = rho*data[EulerKOmegaTerm::VX];
  result[2] = rho*data[EulerKOmegaTerm::VY];
  result[3] = rho*data[EulerKOmegaTerm::VZ];
  result[4] = rho*data[EulerKOmegaTerm::H] - data[EulerKOmegaTerm::P];

  const CFuint iK = _model->getFirstScalarVar(0);
  result[5] = rho*data[iK];
  result[6] = rho*data[iK+1];
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace KOmega

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
