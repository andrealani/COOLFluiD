#include "NavierStokes/NavierStokes.hh"
#include "NavierStokes/Euler2DConsToPrim.hh"
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

Environment::ObjectProvider<Euler2DConsToPrim, VarSetTransformer, NavierStokesModule, 1> euler2DConsToPrimProvider("Euler2DConsToPrim");

//////////////////////////////////////////////////////////////////////////////

Euler2DConsToPrim::Euler2DConsToPrim(Common::SafePtr<Framework::PhysicalModelImpl> model) :
  VarSetTransformer(model),
  _model(model->getConvectiveTerm().d_castTo<EulerTerm>())
{
}

//////////////////////////////////////////////////////////////////////////////

Euler2DConsToPrim::~Euler2DConsToPrim()
{
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DConsToPrim::transform(const State& state, State& result)
{
  const CFreal gamma = _model->getGamma();
  const CFreal rho = state[0];
  const CFreal rhoInv = 1./rho;
  const CFreal u = state[1]*rhoInv;
  const CFreal v = state[2]*rhoInv;
  
  result[0] = rho;
  result[1] = u;
  result[2] = v;
  result[3] = (gamma-1.0)*(state[3] - 0.5*rho*(u*u+v*v)) - _model->getPressInf();
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DConsToPrim::transformFromRef(const RealVector& data, State& result)
{
  result[0] = data[EulerTerm::RHO];
  result[1] = data[EulerTerm::VX];
  result[2] = data[EulerTerm::VY];
  result[3] = data[EulerTerm::P];
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace NavierStokes

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
