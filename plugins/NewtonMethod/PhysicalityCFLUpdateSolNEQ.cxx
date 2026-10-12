// Copyright (C) 2016 KU Leuven, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "NewtonMethod/NewtonMethod.hh"
#include "NewtonMethod/PhysicalityCFLUpdateSolNEQ.hh"

#include "Common/BadValueException.hh"
#include "Common/CFLog.hh"
#include "Framework/BaseTerm.hh"
#include "Framework/FilterState.hh"
#include "Framework/MultiScalarTerm.hh"
#include "Framework/PhysicalModel.hh"
#include "Framework/SpaceMethod.hh"
#include "Framework/SpaceMethodData.hh"
#include "Framework/State.hh"
#include "NavierStokes/EulerTerm.hh"

#include <algorithm>
#include <cmath>
#include <sstream>

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Physics::NavierStokes;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace NewtonMethod {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<PhysicalityCFLUpdateSolNEQ, NewtonIteratorData, NewtonMethodModule>
physicalityCFLUpdateSolNEQProvider("PhysicalityCFLUpdateSolNEQ");

//////////////////////////////////////////////////////////////////////////////

void PhysicalityCFLUpdateSolNEQ::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< std::vector<CFuint> >("PartialDensityVars","Indices of partial densities (rho_i) in the stored update state. Empty disables their check (filled from the model with RhoivtTv / Rhoivt).");
  options.addConfigOption< CFreal, Config::DynamicOption<> >("PartialDensityEtaMax","Largest fractional increase or decrease of each partial density per accepted Newton update (default 0.1).");
  options.addConfigOption< std::vector<CFuint> >("LogVars","Indices of logarithmic variables (e.g. ln rho_i) in the stored update state. Empty disables their check (filled from the model with LogRhoivLogTTv).");
  options.addConfigOption< CFreal, Config::DynamicOption<> >("LogVarsMaxChange","Largest absolute change of each LogVars entry per accepted Newton update (default 2.3, a factor 10 on the variable).");
  options.addConfigOption< std::vector<CFuint> >("LogBoundedVars","Stored logarithms of variables (e.g. ln T, ln Tv) whose relative change per update is bounded by BoundedVarsEtaMax, i.e. |d ln x| <= ln(1 + BoundedVarsEtaMax) (filled from the model with LogRhoivLogTTv).");
  options.addConfigOption< bool >("ChargeNeutralElectrons","With FilterState = ChargeNeutrality: the electron entry (state index 0) is set by the filter after the update, so its own step is dropped before omega is computed and every trial state is filtered before its rho and p are checked (default false).");
}

//////////////////////////////////////////////////////////////////////////////

PhysicalityCFLUpdateSolNEQ::PhysicalityCFLUpdateSolNEQ(const std::string& name) :
  PhysicalityCFLUpdateSol(name),
  m_nbSpecies(0),
  m_nbTvData(0)
{
  addConfigOptionsTo(this);

  m_partialDensityVars = std::vector<CFuint>();
  setParameter("PartialDensityVars",&m_partialDensityVars);

  m_partialDensityEtaMax = 0.1;
  setParameter("PartialDensityEtaMax",&m_partialDensityEtaMax);

  m_logVars = std::vector<CFuint>();
  setParameter("LogVars",&m_logVars);

  m_logVarsMaxChange = 2.3;
  setParameter("LogVarsMaxChange",&m_logVarsMaxChange);

  m_logBoundedVars = std::vector<CFuint>();
  setParameter("LogBoundedVars",&m_logBoundedVars);

  m_chargeNeutralElectrons = false;
  setParameter("ChargeNeutralElectrons",&m_chargeNeutralElectrons);
}

//////////////////////////////////////////////////////////////////////////////

PhysicalityCFLUpdateSolNEQ::~PhysicalityCFLUpdateSolNEQ()
{
}

//////////////////////////////////////////////////////////////////////////////

void PhysicalityCFLUpdateSolNEQ::setup()
{
  fillLayout();

  PhysicalityCFLUpdateSol::setup();

  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  for (CFuint i = 0; i < m_partialDensityVars.size(); ++i) {
    if (m_partialDensityVars[i] >= nbEqs) {
      throw BadValueException(FromHere(), "PhysicalityCFLUpdateSolNEQ: PartialDensityVars index outside the state");
    }
  }
  if (!(m_partialDensityEtaMax > 0. && m_partialDensityEtaMax < 1.)) {
    throw BadValueException(FromHere(), "PhysicalityCFLUpdateSolNEQ: PartialDensityEtaMax must lie in (0,1)");
  }
  for (CFuint i = 0; i < m_logVars.size(); ++i) {
    if (m_logVars[i] >= nbEqs) {
      throw BadValueException(FromHere(), "PhysicalityCFLUpdateSolNEQ: LogVars index outside the state");
    }
  }
  if (!(m_logVarsMaxChange > 0.)) {
    throw BadValueException(FromHere(), "PhysicalityCFLUpdateSolNEQ: LogVarsMaxChange must be > 0");
  }
  for (CFuint i = 0; i < m_logBoundedVars.size(); ++i) {
    if (m_logBoundedVars[i] >= nbEqs) {
      throw BadValueException(FromHere(), "PhysicalityCFLUpdateSolNEQ: LogBoundedVars index outside the state");
    }
  }
  if (m_chargeNeutralElectrons) {
    // the filter writes state[0]: it must be the ChargeNeutrality filter and
    // index 0 must be a species entry
    SafePtr<FilterState> filter = getMethodData().getFilterState();
    if (filter.isNull() || filter->getName() != "ChargeNeutrality") {
      throw BadValueException(FromHere(), "PhysicalityCFLUpdateSolNEQ: ChargeNeutralElectrons needs Data.FilterState = ChargeNeutrality");
    }
    const bool speciesZero =
      std::find(m_logVars.begin(), m_logVars.end(), 0) != m_logVars.end() ||
      std::find(m_partialDensityVars.begin(), m_partialDensityVars.end(), 0) != m_partialDensityVars.end();
    if (!speciesZero) {
      throw BadValueException(FromHere(), "PhysicalityCFLUpdateSolNEQ: ChargeNeutralElectrons needs the electrons at state index 0 (LogVars or PartialDensityVars)");
    }
    // the filter overwrites entry 0 after the update: it leaves the checked lists
    m_logVars.erase(std::remove(m_logVars.begin(), m_logVars.end(), 0), m_logVars.end());
    m_partialDensityVars.erase(std::remove(m_partialDensityVars.begin(), m_partialDensityVars.end(), 0),
                               m_partialDensityVars.end());
    CFLog(INFO, "PhysicalityCFLUpdateSolNEQ: electron entry 0 set by the ChargeNeutrality filter, its step does not limit omega\n");
  }
}

//////////////////////////////////////////////////////////////////////////////

void PhysicalityCFLUpdateSolNEQ::fillLayout()
{
  // the multi-species Euler term: scalar set 0 = species, set 1 = Tv energies
  typedef MultiScalarTerm<EulerTerm> NEQTerm;
  SafePtr<BaseTerm> conv = PhysicalModelStack::getActive()->getImplementor()->getConvectiveTerm();
  const NEQTerm* neq = dynamic_cast<const NEQTerm*>(&*conv);
  if (neq == CFNULL) {
    throw BadValueException(FromHere(),
      "PhysicalityCFLUpdateSolNEQ: the physical model has no multi-species Euler term (use PhysicalityCFLUpdateSol)");
  }
  m_nbSpecies = neq->getNbScalarVars(0);
  m_nbTvData = (neq->getNbScalarVarSets() > 1) ? neq->getNbScalarVars(1) : 0;

  // state layout of RhoivtTv / LogRhoivLogTTv: rho_i (or ln rho_i), velocity, T, Tv
  const CFuint dim = PhysicalModelStack::getActive()->getDim();
  const CFuint TID = m_nbSpecies + dim;
  const std::string upd = getMethodData().getCollaborator<SpaceMethod>()->getSpaceMethodData()->getUpdateVarStr();
  const bool logLayout = (upd.find("LogRhoiv") != std::string::npos);
  const bool rhoLayout = !logLayout && (upd.find("Rhoivt") != std::string::npos);
  std::string filled;
  if (logLayout) {
    if (m_logVars.empty()) {
      for (CFuint i = 0; i < m_nbSpecies; ++i) { m_logVars.push_back(i); }
      filled += " LogVars";
    }
    if (m_logBoundedVars.empty()) {
      for (CFuint i = 0; i <= m_nbTvData; ++i) { m_logBoundedVars.push_back(TID + i); }
      filled += " LogBoundedVars";
    }
  }
  else if (rhoLayout) {
    if (m_partialDensityVars.empty()) {
      for (CFuint i = 0; i < m_nbSpecies; ++i) { m_partialDensityVars.push_back(i); }
      filled += " PartialDensityVars";
    }
  }

  std::ostringstream lists;
  lists << "LogVars =";
  for (CFuint i = 0; i < m_logVars.size(); ++i) { lists << " " << m_logVars[i]; }
  lists << ", LogBoundedVars =";
  for (CFuint i = 0; i < m_logBoundedVars.size(); ++i) { lists << " " << m_logBoundedVars[i]; }
  lists << ", PartialDensityVars =";
  for (CFuint i = 0; i < m_partialDensityVars.size(); ++i) { lists << " " << m_partialDensityVars[i]; }
  CFLog(INFO, "PhysicalityCFLUpdateSolNEQ: " << m_nbSpecies << " species, update variables " << upd
        << (filled.empty() ? std::string(", lists from the case") : ", filled from the model:" + filled)
        << "; " << lists.str() << "\n");
}

//////////////////////////////////////////////////////////////////////////////

void PhysicalityCFLUpdateSolNEQ::boundEntrySteps(const State& state, const CFreal* dU,
                                                 CFreal& omegaTrial, CFreal& omegaLinear)
{
  // partial densities: each change bounded by a fraction of the density (mixture
  // rho and p can stay positive while one species crosses zero)
  for (CFuint i = 0; i < m_partialDensityVars.size(); ++i) {
    const CFuint iEq = m_partialDensityVars[i];
    const CFreal density = state[iEq];
    const CFreal change = std::abs(m_alpha[iEq]*dU[iEq]);
    if (!std::isfinite(density) || density < 0. || !std::isfinite(change)) {
      omegaTrial = 0.;
      return;
    }
    const CFreal allowed = m_partialDensityEtaMax*density;
    if (change > allowed) {
      omegaTrial = std::min(omegaTrial, allowed/change);
    }
  }

  // logarithmic variables: an absolute cap, so exp() of the trial stays finite
  for (CFuint i = 0; i < m_logVars.size(); ++i) {
    const CFuint iEq = m_logVars[i];
    const CFreal change = std::abs(m_alpha[iEq]*dU[iEq]);
    if (!std::isfinite(change)) {
      omegaTrial = 0.;
      return;
    }
    if (change > m_logVarsMaxChange) {
      omegaTrial = std::min(omegaTrial, m_logVarsMaxChange/change);
    }
  }

  // the relative bound of the BoundedVars on a variable stored as its logarithm
  const CFreal logAllowed = std::log(1. + m_boundedEtaMax);
  for (CFuint i = 0; i < m_logBoundedVars.size(); ++i) {
    const CFuint iEq = m_logBoundedVars[i];
    const CFreal change = std::abs(m_alpha[iEq]*dU[iEq]);
    if (change > logAllowed) {
      omegaLinear = std::min(omegaLinear, logAllowed/change);
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

bool PhysicalityCFLUpdateSolNEQ::fullStepEvaluable(const State& state, const State& trial) const
{
  // the thermodynamics of a negative partial density cannot be evaluated
  for (CFuint i = 0; i < m_partialDensityVars.size(); ++i) {
    if (!(trial[m_partialDensityVars[i]] > 0.)) { return false; }
  }
  return PhysicalityCFLUpdateSol::fullStepEvaluable(state, trial);
}

//////////////////////////////////////////////////////////////////////////////

void PhysicalityCFLUpdateSolNEQ::filterTrialState(State& trial)
{
  // the trial gets the electrons the update will get from the filter
  if (m_chargeNeutralElectrons) { getMethodData().getFilterState()->filter(trial); }
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace NewtonMethod

  } // namespace Numerics

} // namespace COOLFluiD
