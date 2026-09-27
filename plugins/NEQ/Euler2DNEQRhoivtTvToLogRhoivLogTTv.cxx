#include <cmath>
#include <algorithm>

#include "NEQ.hh"
#include "Euler2DNEQRhoivtTvToLogRhoivLogTTv.hh"
#include "Environment/ObjectProvider.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace NEQ {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<Euler2DNEQRhoivtTvToLogRhoivLogTTv, VarSetTransformer, NEQModule, 1>
euler2DNEQRhoivtTvToLogRhoivLogTTvProvider("Euler2DNEQRhoivtTvToLogRhoivLogTTv");

// BCDirichlet builds the transformer name from the physical model name
Environment::ObjectProvider<Euler2DNEQRhoivtTvToLogRhoivLogTTv, VarSetTransformer, NEQModule, 1>
ns2DNEQRhoivtTvToLogRhoivLogTTvProvider("NavierStokes2DNEQRhoivtTvToLogRhoivLogTTv");

//////////////////////////////////////////////////////////////////////////////

Euler2DNEQRhoivtTvToLogRhoivLogTTv::Euler2DNEQRhoivtTvToLogRhoivLogTTv
(Common::SafePtr<Framework::PhysicalModelImpl> model) :
  VarSetTransformer(model),
  m_model(model->getConvectiveTerm().d_castTo<NEQTerm>())
{
}

//////////////////////////////////////////////////////////////////////////////

Euler2DNEQRhoivtTvToLogRhoivLogTTv::~Euler2DNEQRhoivtTvToLogRhoivLogTTv()
{
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DNEQRhoivtTvToLogRhoivLogTTv::transform(const State& state, State& result)
{
  const CFuint nbSpecies = m_model->getNbScalarVars(0);
  const CFuint TID = nbSpecies + 2;

  // partial densities below 1e-30 are taken as 1e-30
  for (CFuint ie = 0; ie < nbSpecies; ++ie) {
    result[ie] = std::log(std::max(state[ie], 1e-30));
  }
  for (CFuint i = nbSpecies; i < TID; ++i) {
    result[i] = state[i];
  }
  for (CFuint i = TID; i < state.size(); ++i) {
    result[i] = std::log(state[i]);
  }
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace NEQ

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
