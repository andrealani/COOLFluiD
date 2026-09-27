#ifndef COOLFluiD_Physics_NEQ_Euler2DNEQConsToLogRhoivLogTTvInLogRhoivLogTTv_hh
#define COOLFluiD_Physics_NEQ_Euler2DNEQConsToLogRhoivLogTTvInLogRhoivLogTTv_hh

//////////////////////////////////////////////////////////////////////////////

#include "NEQ/Euler2DNEQConsToRhoivtTvInRhoivtTv.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace NEQ {

//////////////////////////////////////////////////////////////////////////////

/**
 * Matrix transformer from conservative to (ln rho_i, u, v, ln T, ln Tv)
 * variables, evaluated at a (ln rho_i, u, v, ln T, ln Tv) state: the
 * RhoivtTv matrix at the physical state, with the rows of the logarithmic
 * variables divided by the variable (d(ln x) = dx/x).
 *
 * @author Rayan Dhib
 */
class Euler2DNEQConsToLogRhoivLogTTvInLogRhoivLogTTv : public Euler2DNEQConsToRhoivtTvInRhoivtTv {
public:

  /// Constructor
  Euler2DNEQConsToLogRhoivLogTTvInLogRhoivLogTTv(Common::SafePtr<Framework::PhysicalModelImpl> model);

  /// Default destructor
  ~Euler2DNEQConsToLogRhoivLogTTvInLogRhoivLogTTv();

  /// Set the transformation matrix from a (ln rho_i, u, v, ln T, ln Tv) state
  void setMatrix(const RealVector& state);

private:

  /// number of species
  CFuint m_nbSpecies;

  /// the (rho_i, u, v, T, Tv) copy of the state
  RealVector m_rvtState;

}; // end of class Euler2DNEQConsToLogRhoivLogTTvInLogRhoivLogTTv

//////////////////////////////////////////////////////////////////////////////

    } // namespace NEQ

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_NEQ_Euler2DNEQConsToLogRhoivLogTTvInLogRhoivLogTTv_hh
