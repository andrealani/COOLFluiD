#include "NavierStokes/NavierStokes.hh"
#include "Euler1DPvtToCons.hh"
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

Environment::ObjectProvider<Euler1DPvtToCons, VarSetTransformer, NavierStokesModule, 1> euler1DPvtToConsProvider("Euler1DPvtToCons");

//////////////////////////////////////////////////////////////////////////////

Euler1DPvtToCons::Euler1DPvtToCons(Common::SafePtr<Framework::PhysicalModelImpl> model) :
  VarSetTransformer(model),
  _model(model->getConvectiveTerm().d_castTo<EulerTerm>())
{
}

//////////////////////////////////////////////////////////////////////////////

Euler1DPvtToCons::~Euler1DPvtToCons()
{
}

//////////////////////////////////////////////////////////////////////////////

void Euler1DPvtToCons::transform(const State& state, State& result)
{
  const CFreal p = _model->getPressureFromState(state[0]);
  const CFreal u = state[1];
  const CFreal T = state[2];
  const CFreal rho = _model->getDensity(p,T);
  const CFreal rhoV2 = rho*u*u;
  
  result[0] = rho;
  result[1] = rho*u;
  result[2] = p/(_model->getGamma() - 1.) + 0.5*rhoV2;
}

//////////////////////////////////////////////////////////////////////////////

void Euler1DPvtToCons::transformFromRef(const RealVector& data, State& result)
{
  const CFreal rho = data[EulerTerm::RHO];
  
  result[0] = rho;
  result[1] = rho*data[EulerTerm::VX];
  result[2] = rho*data[EulerTerm::H] - _model->getPressureFromState(data[EulerTerm::P]);
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace NavierStokes

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
