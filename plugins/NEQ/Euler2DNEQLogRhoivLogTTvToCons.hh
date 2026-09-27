#ifndef COOLFluiD_Physics_NEQ_Euler2DNEQLogRhoivLogTTvToCons_hh
#define COOLFluiD_Physics_NEQ_Euler2DNEQLogRhoivLogTTvToCons_hh

//////////////////////////////////////////////////////////////////////////////

#include <memory>

#include "NEQ/Euler2DNEQRhoivtTvToCons.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace NEQ {

//////////////////////////////////////////////////////////////////////////////

/**
 * Transformer from (ln rho_i, u, v, ln T, ln Tv) to conservative variables:
 * exp of the logarithms, then Euler2DNEQRhoivtTvToCons.
 *
 * @author Rayan Dhib
 */
class Euler2DNEQLogRhoivLogTTvToCons : public Euler2DNEQRhoivtTvToCons {
public:

  /// Constructor
  Euler2DNEQLogRhoivLogTTvToCons(Common::SafePtr<Framework::PhysicalModelImpl> model);

  /// Default destructor
  virtual ~Euler2DNEQLogRhoivLogTTvToCons();

  /// Transform a state into another one
  virtual void transform(const Framework::State& state, Framework::State& result);

private:

  /// the (rho_i, u, v, T, Tv) copy of the state
  std::unique_ptr<Framework::State> m_rvtState;

}; // end of class Euler2DNEQLogRhoivLogTTvToCons

//////////////////////////////////////////////////////////////////////////////

    } // namespace NEQ

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_NEQ_Euler2DNEQLogRhoivLogTTvToCons_hh
