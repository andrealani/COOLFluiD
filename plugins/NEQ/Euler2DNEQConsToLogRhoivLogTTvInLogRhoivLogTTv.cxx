#include <cmath>

#include "NEQ.hh"
#include "Euler2DNEQConsToLogRhoivLogTTvInLogRhoivLogTTv.hh"
#include "Environment/ObjectProvider.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace NEQ {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<Euler2DNEQConsToLogRhoivLogTTvInLogRhoivLogTTv,
                            VarSetMatrixTransformer,
                            NEQModule, 1>
euler2DNEQConsToLogRhoivLogTTvInLogRhoivLogTTvProvider("Euler2DNEQConsToLogRhoivLogTTvInLogRhoivLogTTv");

//////////////////////////////////////////////////////////////////////////////

Euler2DNEQConsToLogRhoivLogTTvInLogRhoivLogTTv::Euler2DNEQConsToLogRhoivLogTTvInLogRhoivLogTTv
(Common::SafePtr<Framework::PhysicalModelImpl> model) :
  Euler2DNEQConsToRhoivtTvInRhoivtTv(model),
  m_nbSpecies(model->getConvectiveTerm().d_castTo<NEQTerm>()->getNbScalarVars(0)),
  m_rvtState()
{
}

//////////////////////////////////////////////////////////////////////////////

Euler2DNEQConsToLogRhoivLogTTvInLogRhoivLogTTv::~Euler2DNEQConsToLogRhoivLogTTvInLogRhoivLogTTv()
{
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DNEQConsToLogRhoivLogTTvInLogRhoivLogTTv::setMatrix(const RealVector& state)
{
  const CFuint TID = m_nbSpecies + 2;

  // exp of ln rho_i, u and v copied, exp of ln T and ln Tv
  m_rvtState.resize(state.size());
  for (CFuint i = 0; i < state.size(); ++i) {
    m_rvtState[i] = (i < m_nbSpecies || i >= TID) ? std::exp(state[i]) : state[i];
  }

  Euler2DNEQConsToRhoivtTvInRhoivtTv::setMatrix(m_rvtState);

  // d(ln x) = dx/x for the partial densities, T and the Tv's
  const CFuint nbCols = _transMatrix.nbCols();
  for (CFuint i = 0; i < state.size(); ++i) {
    if (i < m_nbSpecies || i >= TID) {
      const CFreal ovX = 1./m_rvtState[i];
      for (CFuint j = 0; j < nbCols; ++j) {
        _transMatrix(i,j) *= ovX;
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace NEQ

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
