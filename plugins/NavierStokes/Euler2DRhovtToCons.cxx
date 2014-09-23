#include "NavierStokes/NavierStokes.hh"
#include "Euler2DRhovtToCons.hh"
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

Environment::ObjectProvider<Euler2DRhovtToCons, VarSetTransformer, NavierStokesModule, 1> euler2DRhovtToConsProvider("Euler2DRhovtToCons");

//////////////////////////////////////////////////////////////////////////////

Euler2DRhovtToCons::Euler2DRhovtToCons(Common::SafePtr<Framework::PhysicalModelImpl> model) :
  VarSetTransformer(model),
  _model(model->getConvectiveTerm().d_castTo<EulerTerm>())
{
}

//////////////////////////////////////////////////////////////////////////////

Euler2DRhovtToCons::~Euler2DRhovtToCons()
{
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DRhovtToCons::transform(const State& state, State& result)
{
  const CFreal R = _model->getR();
  const CFreal rho = state[0];
  const CFreal u = state[1];
  const CFreal v = state[2];
  const CFreal T = state[3];
  const CFreal p = rho*R*T;
  const CFreal rhoV2 = rho*(u*u + v*v);
  
  result[0] = rho;
  result[1] = rho*u;
  result[2] = rho*v;
  result[3] = p/(_model->getGamma() - 1.) + 0.5*rhoV2;
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DRhovtToCons::transformFromRef(const RealVector& data, State& result)
{
  const CFreal rho = data[EulerTerm::RHO];
  
  result[0] = rho;
  result[1] = rho*data[EulerTerm::VX];
  result[2] = rho*data[EulerTerm::VY];
  result[3] = rho*data[EulerTerm::H] - data[EulerTerm::P];
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace NavierStokes

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
