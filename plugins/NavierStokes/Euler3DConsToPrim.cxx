#include "NavierStokes/NavierStokes.hh"
#include "NavierStokes/Euler3DConsToPrim.hh"
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

Environment::ObjectProvider<Euler3DConsToPrim, VarSetTransformer, NavierStokesModule, 1> Euler3DConsToPrimProvider("Euler3DConsToPrim");

//////////////////////////////////////////////////////////////////////////////

Euler3DConsToPrim::Euler3DConsToPrim(Common::SafePtr<Framework::PhysicalModelImpl> model) :
  VarSetTransformer(model),
  _model(model->getConvectiveTerm().d_castTo<EulerTerm>())
{
}

//////////////////////////////////////////////////////////////////////////////

Euler3DConsToPrim::~Euler3DConsToPrim()
{
}

//////////////////////////////////////////////////////////////////////////////

void Euler3DConsToPrim::transform(const State& state, State& result)
{
  const CFreal gamma = _model->getGamma();
  const CFreal rho = state[0];
  const CFreal rhoInv = 1./rho;
  const CFreal u = state[1]*rhoInv;
  const CFreal v = state[2]*rhoInv;
  const CFreal w = state[3]*rhoInv;
  const CFreal p = (gamma-1.0)*(state[4] - 0.5*rho*(u*u+v*v+w*w)) - _model->getPressInf(); // p or dp
  
  result[0] = rho;
  result[1] = u;
  result[2] = v;
  result[3] = w;
  result[4] = p;
}

//////////////////////////////////////////////////////////////////////////////

void Euler3DConsToPrim::transformFromRef(const RealVector& data, State& result)
{
  result[0] = data[EulerTerm::RHO];
  result[1] = data[EulerTerm::VX];
  result[2] = data[EulerTerm::VY];
  result[3] = data[EulerTerm::VZ];
  result[4] = data[EulerTerm::P];
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace NavierStokes

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
