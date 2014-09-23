#include "NavierStokes/NavierStokes.hh"
#include "Euler2DPrimToCons.hh"
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

Environment::ObjectProvider<Euler2DPrimToCons, VarSetTransformer, NavierStokesModule, 1> euler2DPrimToConsProvider("Euler2DPrimToCons");

//////////////////////////////////////////////////////////////////////////////

Euler2DPrimToCons::Euler2DPrimToCons(Common::SafePtr<Framework::PhysicalModelImpl> model) :
  VarSetTransformer(model),
 _model(model->getConvectiveTerm().d_castTo<EulerTerm>())
{
}

//////////////////////////////////////////////////////////////////////////////

Euler2DPrimToCons::~Euler2DPrimToCons()
{
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DPrimToCons::transform(const State& state, State& result)
{
  const CFreal rho  = state[0];
  const CFreal u = state[1];
  const CFreal v = state[2];
  const CFreal p = _model->getPressureFromState(state[3]);
  const CFreal rhoV2 = rho*(u*u + v*v);
  const CFreal gamma = _model->getGamma();
  const CFreal gammaMinus1 = gamma - 1.;
  
  result[0] = rho;
  result[1] = rho*u;
  result[2] = rho*v;
  result[3] = p/gammaMinus1 + 0.5*rhoV2;
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DPrimToCons::transformFromRef(const RealVector& data, State& result)
{
  const CFreal rho = data[EulerTerm::RHO];
  
  result[0] = rho;
  result[1] = rho*data[EulerTerm::VX];
  result[2] = rho*data[EulerTerm::VY];
  result[3] = rho*data[EulerTerm::H] - _model->getPressureFromState(data[EulerTerm::P]);
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace NavierStokes

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
