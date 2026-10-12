// Copyright (C) 2016 KU Leuven, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Numerics_NewtonMethod_PhysicalityCFLUpdateSolNEQ_hh
#define COOLFluiD_Numerics_NewtonMethod_PhysicalityCFLUpdateSolNEQ_hh

//////////////////////////////////////////////////////////////////////////////

#include <vector>

#include "NewtonMethod/PhysicalityCFLUpdateSol.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace NewtonMethod {

//////////////////////////////////////////////////////////////////////////////

/// Physicality-limited Newton update for thermochemical nonequilibrium flows
/// (species partial densities, velocity, T and Tv). On top of
/// PhysicalityCFLUpdateSol:
///
///  - PartialDensityVars (rho_i): each change is bounded by PartialDensityEtaMax
///    (default 0.1) times rho_i before the trial states are evaluated; a negative
///    or non-finite density or step rejects the update.
///
///  - LogVars (e.g. ln rho_i): each change is capped at LogVarsMaxChange, so exp()
///    of a trial stays finite. LogBoundedVars (e.g. ln T, ln Tv): each change is
///    capped at ln(1 + BoundedVarsEtaMax).
///
///  - Layout from the model: LogRhoivLogTTv gives the LogVars (ln rho_i) and
///    LogBoundedVars (T, Tv), RhoivtTv / Rhoivt the PartialDensityVars. Lists
///    set in the case win.
///
///  - ChargeNeutralElectrons (with FilterState = ChargeNeutrality): the electron
///    entry (state index 0) is set by the filter after every update, so it
///    leaves the checks above and every trial state is filtered before its rho
///    and p are checked.
///
/// @author Rayan Dhib

class PhysicalityCFLUpdateSolNEQ : public PhysicalityCFLUpdateSol {
public:

  /// Defines the Config Option's of this class
  /// @param options a OptionList where to add the Option's
  static void defineConfigOptions(Config::OptionList& options);

  /// Constructor.
  explicit PhysicalityCFLUpdateSolNEQ(const std::string& name);

  /// Destructor.
  virtual ~PhysicalityCFLUpdateSolNEQ();

  /// Fill the index lists from the model, then set up the base
  virtual void setup();

protected:

  /// PartialDensityVars and LogVars bounds into omegaTrial, LogBoundedVars cap
  /// into omegaLinear
  virtual void boundEntrySteps(const Framework::State& state, const CFreal* dU,
                               CFreal& omegaTrial, CFreal& omegaLinear);

  /// the base test, and every PartialDensityVars entry of the full step positive
  virtual bool fullStepEvaluable(const Framework::State& state,
                                 const Framework::State& trial) const;

  /// with ChargeNeutralElectrons: the trial state gets the ChargeNeutrality filter
  virtual void filterTrialState(Framework::State& trial);

  /// the binding-state print also gives the full step
  virtual bool reportFullStep() const { return true; }

private:

  /// fill the species and temperature index lists from the model
  void fillLayout();

private:

  /// indices of partial densities (rho_i) in the stored update state
  std::vector<CFuint> m_partialDensityVars;

  /// largest fractional increase or decrease allowed for each partial density
  CFreal m_partialDensityEtaMax;

  /// indices of logarithmic variables (e.g. ln rho_i) in the stored update state
  std::vector<CFuint> m_logVars;

  /// largest absolute change allowed for each LogVars entry per update
  CFreal m_logVarsMaxChange;

  /// stored logarithms of variables (e.g. ln T) whose relative change is bounded by BoundedVarsEtaMax
  std::vector<CFuint> m_logBoundedVars;

  /// the electron entry is set by the ChargeNeutrality filter: drop its step and
  /// filter the trial states
  bool m_chargeNeutralElectrons;

  /// number of species
  CFuint m_nbSpecies;

  /// number of Tv energies in the physical data
  CFuint m_nbTvData;

}; // class PhysicalityCFLUpdateSolNEQ

//////////////////////////////////////////////////////////////////////////////

    } // namespace NewtonMethod

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_NewtonMethod_PhysicalityCFLUpdateSolNEQ_hh
