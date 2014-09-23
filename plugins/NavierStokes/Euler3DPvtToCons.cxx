#include "NavierStokes/NavierStokes.hh"
#include "Euler3DPvtToCons.hh"
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

Environment::ObjectProvider<Euler3DPvtToCons, VarSetTransformer, NavierStokesModule, 1>
euler3DPvtToConsProvider("Euler3DPvtToCons");

Environment::ObjectProvider<Euler3DPvtToCons, VarSetTransformer, NavierStokesModule, 1>
euler3DRotationPvtToRotationConsProvider("Euler3DRotationPvtToRotationCons");
        
//////////////////////////////////////////////////////////////////////////////

Euler3DPvtToCons::Euler3DPvtToCons(Common::SafePtr<Framework::PhysicalModelImpl> model) :
  VarSetTransformer(model),
 _model(model->getConvectiveTerm().d_castTo<EulerTerm>())
{
}

//////////////////////////////////////////////////////////////////////////////

Euler3DPvtToCons::~Euler3DPvtToCons()
{
}

//////////////////////////////////////////////////////////////////////////////

void Euler3DPvtToCons::transform(const State& state, State& result)
{
  const CFreal p = _model->getPressureFromState(state[0]);
  const CFreal u = state[1];
  const CFreal v = state[2];
  const CFreal w = state[3];
  const CFreal T = state[4];
  const CFreal rho = _model->getDensity(p,T);
  const CFreal rhoV2 = rho*(u*u + v*v + w*w);

  result[0] = rho;
  result[1] = rho*u;
  result[2] = rho*v;
  result[3] = rho*w;
  result[4] = p/(_model->getGamma() - 1.) + 0.5*rhoV2;
}

//////////////////////////////////////////////////////////////////////////////

void Euler3DPvtToCons::transformFromRef(const RealVector& data, State& result)
{
  const CFreal rho = data[EulerTerm::RHO];

  result[0] = rho;
  result[1] = rho*data[EulerTerm::VX];
  result[2] = rho*data[EulerTerm::VY];
  result[3] = rho*data[EulerTerm::VZ];
  result[4] = rho*data[EulerTerm::H] - _model->getPressureFromState(data[EulerTerm::P]);
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace NavierStokes

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
