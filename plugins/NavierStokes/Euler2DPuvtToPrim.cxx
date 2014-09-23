#include "NavierStokes/NavierStokes.hh"
#include "Euler2DPuvtToPrim.hh"
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

Environment::ObjectProvider<Euler2DPuvtToPrim, VarSetTransformer, NavierStokesModule, 1> euler2DPuvtToPrimProvider("Euler2DPuvtToPrim");

//////////////////////////////////////////////////////////////////////////////

Euler2DPuvtToPrim::Euler2DPuvtToPrim(Common::SafePtr<Framework::PhysicalModelImpl> model) :
  VarSetTransformer(model),
  _model(model->getConvectiveTerm().d_castTo<EulerTerm>())
{
}

//////////////////////////////////////////////////////////////////////////////

Euler2DPuvtToPrim::~Euler2DPuvtToPrim()
{
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DPuvtToPrim::transform(const State& state, State& result)
{
  const CFreal p = _model->getPressureFromState(state[0]);
  result[0] = _model->getDensity(p, state[3]);
  result[1] = state[1];
  result[2] = state[2];
  result[3] = state[0];
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DPuvtToPrim::transformFromRef(const RealVector& data, State& result)
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
