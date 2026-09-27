#ifndef COOLFluiD_Physics_NEQ_Euler2DNEQLogRhoivLogTTv_hh
#define COOLFluiD_Physics_NEQ_Euler2DNEQLogRhoivLogTTv_hh

//////////////////////////////////////////////////////////////////////////////

#include <memory>

#include "NEQ/Euler2DNEQRhoivtTv.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace NEQ {

//////////////////////////////////////////////////////////////////////////////

/**
 * 2D Euler variable set with chemical and thermal NEQ in the variables
 * (ln rho_i, u, v, ln T, ln Tv): logarithms of the partial densities rho_i,
 * of the temperature T and of every vibrational temperature Tv.
 *
 * The method interpolates the stored variables, so every reconstructed
 * partial density and temperature (exp of a polynomial) is positive. Each
 * function turns the state into its (rho_i, u, v, T, Tv) copy and calls
 * Euler2DNEQRhoivtTv on it.
 *
 * @author Rayan Dhib
 */
class Euler2DNEQLogRhoivLogTTv : public Euler2DNEQRhoivtTv {

public:

  /// Constructor
  Euler2DNEQLogRhoivLogTTv(Common::SafePtr<Framework::BaseTerm> term);

  /// Default destructor
  virtual ~Euler2DNEQLogRhoivLogTTv();

  /// Set up the private data
  virtual void setup();

  /// Names of the extra output variables: the parent ones, rho_i, T and Tv
  virtual std::vector<std::string> getExtraVarNames() const;

  /// Physical data of a (ln rho_i, u, v, ln T, ln Tv) state
  virtual void computePhysicalData(const Framework::State& state,
                                   RealVector& data);

  /**
   * Physical data of a state with the variable iVar perturbed. A perturbed
   * u or v only changes the kinetic energy (parent branch, which reads u and
   * v only), any other variable needs the full computation.
   */
  virtual void computePerturbedPhysicalData(const Framework::State& state,
                                            const RealVector& pdataBkp,
                                            RealVector& pdata,
                                            CFuint iVar);

  /// Dimensional values: ln rho_i + ln rho_ref, u, v times their reference, ln T + ln T_ref
  virtual void setDimensionalValues(const Framework::State& state,
                                    RealVector& result);

  /// Adimensional values: ln rho_i - ln rho_ref, u, v over their reference, ln T - ln T_ref
  virtual void setAdimensionalValues(const Framework::State& state,
                                     RealVector& result);

  /// Dimensional values plus the parent extra values, the partial densities and the temperatures
  virtual void setDimensionalValuesPlusExtraValues(const Framework::State& state,
                                                   RealVector& result,
                                                   RealVector& extra);

  /// Pressure derivatives: dp/d(ln rho_i) = rho_i dp/d(rho_i), dp/d(ln T) = T dp/dT
  virtual void computePressureDerivatives(const Framework::State& state, RealVector& dp);

  /// Checks a (ln rho_i, u, v, ln T, ln Tv) state: every logarithm finite
  virtual bool isValid(const RealVector& data);

private:

  /// fills m_rvtState with (rho_i, u, v, T, Tv) of state
  void setRhoivtTvState(const Framework::State& state);

private:

  /// the (rho_i, u, v, T, Tv) copy of the current state
  std::unique_ptr<Framework::State> m_rvtState;

  /// parent extra values
  RealVector m_parentExtra;

}; // end of class Euler2DNEQLogRhoivLogTTv

//////////////////////////////////////////////////////////////////////////////

    } // namespace NEQ

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_NEQ_Euler2DNEQLogRhoivLogTTv_hh
