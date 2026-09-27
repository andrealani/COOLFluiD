#include <cmath>

#include "NEQ.hh"
#include "Euler2DNEQLogRhoivLogTTvToCons.hh"
#include "Environment/ObjectProvider.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace NEQ {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<Euler2DNEQLogRhoivLogTTvToCons, VarSetTransformer, NEQModule, 1>
euler2DNEQLogRhoivLogTTvToConsProvider("Euler2DNEQLogRhoivLogTTvToCons");

//////////////////////////////////////////////////////////////////////////////

Euler2DNEQLogRhoivLogTTvToCons::Euler2DNEQLogRhoivLogTTvToCons
(Common::SafePtr<Framework::PhysicalModelImpl> model) :
  Euler2DNEQRhoivtTvToCons(model),
  m_rvtState()
{
}

//////////////////////////////////////////////////////////////////////////////

Euler2DNEQLogRhoivLogTTvToCons::~Euler2DNEQLogRhoivLogTTvToCons()
{
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DNEQLogRhoivLogTTvToCons::transform(const State& state, State& result)
{
  if (m_rvtState.get() == CFNULL) {
    m_rvtState.reset(new State());
  }

  const CFuint nbSpecies = _model->getNbScalarVars(0);
  const CFuint TID = nbSpecies + 2;
  State& rvt = *m_rvtState;

  // exp of ln rho_i, u and v copied, exp of ln T and ln Tv
  for (CFuint i = 0; i < state.size(); ++i) {
    rvt[i] = (i < nbSpecies || i >= TID) ? std::exp(state[i]) : state[i];
  }

  Euler2DNEQRhoivtTvToCons::transform(rvt, result);
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace NEQ

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
