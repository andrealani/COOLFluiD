#include "NavierStokes/NavierStokes.hh"
#include "Euler3DRhovtToCons.hh"
#include "Framework/PhysicalModel.hh"
#include "Environment/ObjectProvider.hh"
#include "EulerTerm.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace NavierStokes {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<Euler3DRhovtToCons, VarSetTransformer, NavierStokesModule, 1> euler3DRhovtToConsProvider("Euler3DRhovtToCons");

//////////////////////////////////////////////////////////////////////////////

Euler3DRhovtToCons::Euler3DRhovtToCons(Common::SafePtr<Framework::PhysicalModelImpl> model) :
  VarSetTransformer(model),
  _model(model->getConvectiveTerm().d_castTo<EulerTerm>())
{
}

//////////////////////////////////////////////////////////////////////////////

Euler3DRhovtToCons::~Euler3DRhovtToCons()
{
}

//////////////////////////////////////////////////////////////////////////////

void Euler3DRhovtToCons::transform(const State& state, State& result)
{
  const CFreal R = _model->getR();
  const CFreal rho = state[0];
  const CFreal u = state[1];
  const CFreal v = state[2];
  const CFreal w = state[3];
  const CFreal T = state[4];
  const CFreal p = rho*R*T;
  const CFreal rhoV2 = rho*(u*u + v*v + w*w);
  
  result[0] = rho;
  result[1] = rho*u;
  result[2] = rho*v;
  result[3] = rho*w;
  result[4] = p/(_model->getGamma() - 1.) + 0.5*rhoV2;
}

//////////////////////////////////////////////////////////////////////////////

void Euler3DRhovtToCons::transformFromRef(const RealVector& data, State& result)
{
  const CFreal rho = data[EulerTerm::RHO];
  
  result[0] = rho;
  result[1] = rho*data[EulerTerm::VX];
  result[2] = rho*data[EulerTerm::VY];
  result[3] = rho*data[EulerTerm::VZ];
  result[4] = rho*data[EulerTerm::H] - data[EulerTerm::P];
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace NavierStokes

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
