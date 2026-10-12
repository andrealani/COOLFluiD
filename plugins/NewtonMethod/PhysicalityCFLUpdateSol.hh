// Copyright (C) 2016 KU Leuven, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Numerics_NewtonMethod_PhysicalityCFLUpdateSol_hh
#define COOLFluiD_Numerics_NewtonMethod_PhysicalityCFLUpdateSol_hh

//////////////////////////////////////////////////////////////////////////////

#include "NewtonMethod/StdUpdateSol.hh"
#include "Framework/ConvectiveVarSet.hh"
#include "MathTools/RealVector.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework { class State; }

  namespace Numerics {

    namespace NewtonMethod {

//////////////////////////////////////////////////////////////////////////////

/// Physicality-limited Newton update with a CFL schedule driven by the
/// under-relaxation the update needed.
///
/// Before the states are modified, the largest omega in [0,1] such that
/// state + omega*Relaxation*dU keeps density and pressure above
/// (1 - EtaMax) times their current values at every state is computed
/// (global minimum over all ranks). Then:
///  - omega < OmegaMin: the update direction is rejected. The states are left
///    untouched, the CFL is multiplied by CFLBackoff and the same Newton
///    iteration is repeated (iteration counter rewound, update norm forced
///    above any stop target). MaxRejections consecutive rejections abort the
///    run: the step size is no longer the limiting factor.
///  - otherwise the update is applied as state += omega*Relaxation*dU. The
///    CFL grows by CFLGrowth, capped at CFLMax, only when omega == 1 and no
///    rejection happened in the current time step; a limited update or a
///    retried step holds the CFL, so the schedule stops growing before the
///    stability wall instead of after hitting it. CFLMax <= 0 disables growth.
///
/// Density and pressure are read from the update variable set's physical data
/// at RhoIndex and PIndex, so the command works for any variable choice of a
/// physics whose physical data carries rho and p. This command owns the CFL
/// value: run it with Data.CFL.ComputeCFL = Null.
///
/// Optional BoundedVars are indices in the stored update state, not physical
/// data. The change of each one is bounded by BoundedVarsEtaMax times its
/// current magnitude. This scale also bounds the segment on which rho and p
/// are checked, so a variable that must stay positive (p and T of Puvt, with
/// BoundedVarsEtaMax < 1) is never evaluated past its bound. An empty index
/// list disables these checks.
///
/// PhysicalityCFLUpdateSolNEQ extends the limiter through the protected hooks.
///
/// Limiter: Ceze and Fidkowski, Int. J. Numer. Meth. Engng 102 (2015),
/// Algorithm 2. CFL gating on the relaxation factor: Ceze and Fidkowski,
/// AIAA 2013-2686 (exponential progression with under-relaxation).
///
/// @author Rayan Dhib

class PhysicalityCFLUpdateSol : public StdUpdateSol {
public:

  /// Defines the Config Option's of this class
  /// @param options a OptionList where to add the Option's
  static void defineConfigOptions(Config::OptionList& options);

  /// Constructor.
  explicit PhysicalityCFLUpdateSol(const std::string& name);

  /// Destructor.
  virtual ~PhysicalityCFLUpdateSol();

  /// Set up private data and data of the aggregated classes
  /// in this command before processing phase
  virtual void setup();

  /// Unset up private data and data of the aggregated classes
  virtual void unsetup();

  /// Execute Processing actions
  virtual void execute();

protected:

  /// Hook: lower the two step scales of one state with bounds on its single
  /// entries. omegaTrial bounds the segment on which rho and p are checked, so
  /// it must keep the trial states evaluable; omegaLinear is only combined with
  /// the result. omegaTrial = 0 rejects the update.
  virtual void boundEntrySteps(const Framework::State& state, const CFreal* dU,
                               CFreal& omegaTrial, CFreal& omegaLinear) {}

  /// Hook: also print the full step of the binding state
  virtual bool reportFullStep() const { return false; }

  /// Hook: change a trial state the way the update will (e.g. a state filter),
  /// before its rho and p are checked
  virtual void filterTrialState(Framework::State& trial) {}

  /// Hook: true when the physical data of the full-step trial of a state can be
  /// evaluated for the binding-state print (every BoundedVars entry keeps its sign)
  virtual bool fullStepEvaluable(const Framework::State& state,
                                 const Framework::State& trial) const;

  /// rho and p of a state through the update variable set's physical data
  /// (also fills m_pdata)
  void computeRhoP(const Framework::State& state, CFreal& rho, CFreal& p);

private:

  /// Largest omega in [0,1] such that state + omega*Relaxation*dU keeps
  /// rho and p above their floors. Returns 0 when the state itself is not
  /// physical or the update is not finite, which forces a rejection.
  CFreal computeStateOmega(const Framework::State& state, const CFreal* dU);

  /// Fill the trial state with state + omega*Relaxation*dU
  void setTrialState(const Framework::State& state, const CFreal* dU,
                     const CFreal omega);

  /// Reject the current update: cut the CFL and make the Newton loop repeat
  /// this iteration with the states untouched.
  void rejectUpdate(const CFreal omega);

protected:

  /// physical data scratch (filled by computeRhoP)
  RealVector m_pdata;

  /// state variables whose relative change per update is bounded by BoundedVarsEtaMax
  std::vector<CFuint> m_boundedVars;

  /// largest relative change allowed for the BoundedVars
  CFreal m_boundedEtaMax;

private:

  /// update variable set (owned by the space method)
  Common::SafePtr<Framework::ConvectiveVarSet> m_varSet;

  /// trial state scratch
  Framework::State* m_trial;

  /// largest fractional decrease of rho and p allowed per update
  CFreal m_etaMax;

  /// relaxation factor below which the update is rejected
  CFreal m_omegaMin;

  /// CFL growth factor after an unlimited update
  CFreal m_cflGrowth;

  /// CFL cut factor on a rejected update
  CFreal m_cflBackoff;

  /// CFL ceiling for the growth (<= 0: no growth)
  CFreal m_cflMax;

  /// consecutive rejected updates that abort the run
  CFuint m_maxRejections;

  /// position of the density in the physical data
  CFuint m_rhoIndex;

  /// position of the pressure in the physical data
  CFuint m_pIndex;

  /// consecutive rejected updates so far
  CFuint m_nbConsecutiveRejections;

  /// global iteration seen by the last execute (detects a new time step)
  CFuint m_lastGlobalIter;

  /// an update was rejected in the current time step
  bool m_rejectedThisStep;

}; // class PhysicalityCFLUpdateSol

//////////////////////////////////////////////////////////////////////////////

    } // namespace NewtonMethod

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_NewtonMethod_PhysicalityCFLUpdateSol_hh
