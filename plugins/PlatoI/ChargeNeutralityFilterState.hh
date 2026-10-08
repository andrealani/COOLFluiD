#ifndef COOLFluiD_Physics_Plato_ChargeNeutralityFilterState_hh
#define COOLFluiD_Physics_Plato_ChargeNeutralityFilterState_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/FilterState.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace Plato {

//////////////////////////////////////////////////////////////////////////////

/**
 * State filter for ionized PLATO mixtures run with Plato.ChargeNeutrality = true. After every
 * update of the solution it sets the electron partial density of the stored state from the ions,
 * rho_e = sum_k c_k rho_k with c_k = m_e q_k / m_k (m = molar mass, q_k = charge number of species
 * k, c_k = 0 for the non-ions). Without it the electron continuity equation carries a copy of rho_e
 * that PLATO no longer uses, so the written electron density and the electron pressure computed
 * from the state are not the charge-neutral ones.
 * The electron is the first species and the partial densities are the first entries of the state.
 * Use: NewtonIterator.Data.FilterState = ChargeNeutrality with
 * NewtonIterator.Data.ChargeNeutrality.DensityVariables = Log (state stores ln rho_i, update
 * variables LogRhoivLogTTv) or Linear (state stores rho_i, e.g. RhoivtTv).
 *
 * @author Rayan Dhib
 */
class ChargeNeutralityFilterState : public Framework::FilterState {
public:

  /// Defines the Config Option's of this class
  /// @param options a OptionList where to add the Option's
  static void defineConfigOptions(Config::OptionList& options);

  /// Constructor
  ChargeNeutralityFilterState(const std::string& name);

  /// Destructor
  virtual ~ChargeNeutralityFilterState();

  /// Configures the filter and checks DensityVariables
  virtual void configure(Config::ConfigArgs& args);

  /// Sets the electron partial density of the state from the ions
  /// @param state state to be filtered
  virtual void filter(RealVector& state) const;

private:

  /// Reads c_k from the PLATO library (first call, when the library is set up)
  void setCoefficients() const;

private:

  /// "Log" if the state stores ln(rho_i), "Linear" if it stores rho_i
  std::string m_densityVariables;

  /// true if the state stores ln(rho_i)
  bool m_logDensities;

  /// c_k = m_e q_k / m_k for the ions, 0 for the other species
  mutable RealVector m_coef;

  /// true once m_coef is set
  mutable bool m_isSet;

}; // end of class ChargeNeutralityFilterState

//////////////////////////////////////////////////////////////////////////////

    } // namespace Plato

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_Plato_ChargeNeutralityFilterState_hh
