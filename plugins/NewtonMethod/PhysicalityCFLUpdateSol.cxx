// Copyright (C) 2016 KU Leuven, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "NewtonMethod/NewtonMethod.hh"
#include "NewtonMethod/PhysicalityCFLUpdateSol.hh"

#include "Common/BadValueException.hh"
#include "Common/CFLog.hh"
#include "Common/PE.hh"
#include "Framework/BaseTerm.hh"
#include "Framework/MeshData.hh"
#include "Framework/PhysicalModel.hh"
#include "Framework/SpaceMethod.hh"
#include "Framework/SpaceMethodData.hh"
#include "Framework/State.hh"
#include "Framework/SubSystemStatus.hh"

#include <algorithm>
#include <cmath>
#include <limits>

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace NewtonMethod {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<PhysicalityCFLUpdateSol, NewtonIteratorData, NewtonMethodModule>
physicalityCFLUpdateSolProvider("PhysicalityCFLUpdateSol");

//////////////////////////////////////////////////////////////////////////////

void PhysicalityCFLUpdateSol::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< CFreal, Config::DynamicOption<> >("EtaMax","Largest fractional decrease of density and pressure allowed per update.");
  options.addConfigOption< CFreal, Config::DynamicOption<> >("OmegaMin","Relaxation factor below which the update is rejected and the CFL cut.");
  options.addConfigOption< CFreal, Config::DynamicOption<> >("CFLGrowth","CFL growth factor after an unlimited update.");
  options.addConfigOption< CFreal, Config::DynamicOption<> >("CFLBackoff","CFL cut factor on a rejected update.");
  options.addConfigOption< CFreal, Config::DynamicOption<> >("CFLMax","CFL ceiling for the growth (<= 0: no growth).");
  options.addConfigOption< CFuint >("MaxRejections","Consecutive rejected updates that abort the run.");
  options.addConfigOption< CFuint >("RhoIndex","Position of the density in the physical data.");
  options.addConfigOption< CFuint >("PIndex","Position of the pressure in the physical data.");
  options.addConfigOption< std::vector<CFuint> >("BoundedVars","State variables whose relative change per update is bounded by BoundedVarsEtaMax (e.g. the temperatures, or p and T of Puvt). The bound also limits the segment on which rho and p are checked.");
  options.addConfigOption< CFreal, Config::DynamicOption<> >("BoundedVarsEtaMax","Largest relative change allowed for the BoundedVars.");
}

//////////////////////////////////////////////////////////////////////////////

PhysicalityCFLUpdateSol::PhysicalityCFLUpdateSol(const std::string& name) :
  StdUpdateSol(name),
  m_pdata(),
  m_varSet(CFNULL),
  m_trial(CFNULL),
  m_nbConsecutiveRejections(0),
  m_lastGlobalIter(0),
  m_rejectedThisStep(false)
{
  addConfigOptionsTo(this);

  m_etaMax = 0.1;
  setParameter("EtaMax",&m_etaMax);

  m_omegaMin = 0.01;
  setParameter("OmegaMin",&m_omegaMin);

  m_cflGrowth = 1.5;
  setParameter("CFLGrowth",&m_cflGrowth);

  m_cflBackoff = 0.1;
  setParameter("CFLBackoff",&m_cflBackoff);

  m_cflMax = 0.;
  setParameter("CFLMax",&m_cflMax);

  m_maxRejections = 10;
  setParameter("MaxRejections",&m_maxRejections);

  m_rhoIndex = 0;
  setParameter("RhoIndex",&m_rhoIndex);

  m_pIndex = 1;
  setParameter("PIndex",&m_pIndex);

  m_boundedVars = std::vector<CFuint>();
  setParameter("BoundedVars",&m_boundedVars);

  m_boundedEtaMax = 0.5;
  setParameter("BoundedVarsEtaMax",&m_boundedEtaMax);
}

//////////////////////////////////////////////////////////////////////////////

PhysicalityCFLUpdateSol::~PhysicalityCFLUpdateSol()
{
  if (m_trial != CFNULL) { m_trial->resetSpaceCoordinates(); }
  delete m_trial;
}

//////////////////////////////////////////////////////////////////////////////

void PhysicalityCFLUpdateSol::setup()
{
  StdUpdateSol::setup();

  if (!(m_etaMax > 0. && m_etaMax < 1.)) {
    throw BadValueException(FromHere(), "PhysicalityCFLUpdateSol: EtaMax must lie in (0,1)");
  }
  if (!(m_omegaMin > 0. && m_omegaMin <= 1.)) {
    throw BadValueException(FromHere(), "PhysicalityCFLUpdateSol: OmegaMin must lie in (0,1]");
  }
  if (!(m_cflBackoff > 0. && m_cflBackoff < 1.)) {
    throw BadValueException(FromHere(), "PhysicalityCFLUpdateSol: CFLBackoff must lie in (0,1)");
  }
  if (!(m_cflGrowth >= 1.)) {
    throw BadValueException(FromHere(), "PhysicalityCFLUpdateSol: CFLGrowth must be >= 1");
  }
  if (m_maxRejections == 0) {
    throw BadValueException(FromHere(), "PhysicalityCFLUpdateSol: MaxRejections must be > 0");
  }

  SafePtr<BaseTerm> convTerm =
    PhysicalModelStack::getActive()->getImplementor()->getConvectiveTerm();
  convTerm->resizePhysicalData(m_pdata);
  for (CFuint i = 0; i < m_boundedVars.size(); ++i) {
    if (m_boundedVars[i] >= PhysicalModelStack::getActive()->getNbEq()) {
      throw BadValueException(FromHere(), "PhysicalityCFLUpdateSol: BoundedVars index outside the state");
    }
  }
  if (!(m_boundedEtaMax > 0.)) {
    throw BadValueException(FromHere(), "PhysicalityCFLUpdateSol: BoundedVarsEtaMax must be > 0");
  }
  if (m_rhoIndex >= m_pdata.size() || m_pIndex >= m_pdata.size()) {
    throw BadValueException(FromHere(),
      "PhysicalityCFLUpdateSol: RhoIndex/PIndex outside the physical data of size " +
      StringOps::to_str(m_pdata.size()));
  }

  m_trial = new State();
}

//////////////////////////////////////////////////////////////////////////////

void PhysicalityCFLUpdateSol::unsetup()
{
  if (m_trial != CFNULL) { m_trial->resetSpaceCoordinates(); }
  delete m_trial;
  m_trial = CFNULL;

  StdUpdateSol::unsetup();
}

//////////////////////////////////////////////////////////////////////////////

void PhysicalityCFLUpdateSol::execute()
{
  CFAUTOTRACE;

  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();
  DataHandle<CFreal> dU = socket_rhs.getDataHandle();

  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  const CFuint nbStates = states.size();

  m_varSet = getMethodData().getCollaborator<SpaceMethod>()->getSpaceMethodData()->getUpdateVar();

  // a new time step starts with a clean retry record
  const CFuint globalIter = SubSystemStatusStack::getActive()->getNbIter();
  if (globalIter != m_lastGlobalIter) {
    m_lastGlobalIter = globalIter;
    m_rejectedThisStep = false;
  }

  // 1. relaxation factor: smallest admissible scale over all states and ranks
  CFreal omega = 1.;
  CFuint iBind = 0;
  for (CFuint iState = 0; iState < nbStates; ++iState) {
    const State& state = *states[iState];
    if (state.isParUpdatable()) {
      const CFreal omegaState = computeStateOmega(state, &dU[iState*nbEqs]);
      if (omegaState < omega) { omega = omegaState; iBind = iState; }
    }
  }
  const CFreal omegaLocal = omega;
#ifdef CF_HAVE_MPI
  if (PE::GetPE().IsParallel()) {
    const std::string nsp = getMethodData().getNamespace();
    MPI_Allreduce(&omegaLocal, &omega, 1, MPI_DOUBLE, MPI_MIN,
                  PE::GetPE().GetCommunicator(nsp));
  }
#endif

  // Say which state is throttling the whole field when it bites hard. The
  // owning rank packs the information and rank 0 prints it, since only rank 0
  // writes to stdout by default.
  if (omega < 0.1) {
    const CFuint dim = PhysicalModelStack::getActive()->getDim();
    const CFuint nbBlocks = reportFullStep() ? 2 : 1;
    std::vector<CFreal> info(6 + dim + nbBlocks*nbEqs, 0.);
    if (omegaLocal == omega) {
      const State& state = *states[iBind];
      CFreal rho0 = 0., p0 = 0., rho1 = 0., p1 = 0.;
      computeRhoP(state, rho0, p0);
      setTrialState(state, &dU[iBind*nbEqs], 1.);
      // the full step may leave the admissible set, where the physical data cannot be
      // evaluated (Euler2DPuvt asserts p > 0): report it only when the hook says the
      // full-step state can be evaluated, NaN otherwise
      if (fullStepEvaluable(state, *m_trial)) {
        computeRhoP(*m_trial, rho1, p1);
      }
      else {
        rho1 = std::numeric_limits<CFreal>::quiet_NaN();
        p1   = std::numeric_limits<CFreal>::quiet_NaN();
      }
      info[0] = 1.; info[1] = rho0; info[2] = rho1; info[3] = p0; info[4] = p1;
      info[5] = static_cast<CFreal>(iBind);
      for (CFuint d = 0; d < dim; ++d) { info[6+d] = state.getCoordinates()[d]; }
      for (CFuint iEq = 0; iEq < nbEqs; ++iEq) { info[6+dim+iEq] = state[iEq]; }
      if (reportFullStep()) {
        for (CFuint iEq = 0; iEq < nbEqs; ++iEq) { info[6+dim+nbEqs+iEq] = m_alpha[iEq]*dU[iBind*nbEqs+iEq]; }
      }
    }
#ifdef CF_HAVE_MPI
    if (PE::GetPE().IsParallel()) {
      // one rank has info[0] = 1 (ties are harmless), the others zeros
      const std::string nsp = getMethodData().getNamespace();
      std::vector<CFreal> infoLocal = info;
      MPI_Allreduce(&infoLocal[0], &info[0], info.size(), MPI_DOUBLE, MPI_SUM,
                    PE::GetPE().GetCommunicator(nsp));
    }
#endif
    const CFreal nOwners = std::max(info[0], 1.);
    CFLog(INFO, "PhysicalityCFLUpdateSol: binding state at (");
    for (CFuint d = 0; d < dim; ++d) { CFLog(INFO, info[6+d]/nOwners << (d+1 < dim ? ", " : ")")); }
    CFLog(INFO, ", full step would take rho " << info[1]/nOwners << " -> " << info[2]/nOwners
          << ", p " << info[3]/nOwners << " -> " << info[4]/nOwners << ", state = [");
    for (CFuint iEq = 0; iEq < nbEqs; ++iEq) { CFLog(INFO, info[6+dim+iEq]/nOwners << " "); }
    if (reportFullStep()) {
      CFLog(INFO, "], full step = [");
      for (CFuint iEq = 0; iEq < nbEqs; ++iEq) { CFLog(INFO, info[6+dim+nbEqs+iEq]/nOwners << " "); }
    }
    CFLog(INFO, "]\n");
  }

  // 2. hopeless direction: cut the CFL and redo this iteration
  if (!(omega >= m_omegaMin)) {
    rejectUpdate(omega);
    return;
  }

  // 3. accepted: apply omega*Relaxation*dU through the base class so that
  //    the filters, the validation and the updateCoeff reset stay in one place
  m_nbConsecutiveRejections = 0;
  if (omega < 1.) {
    const std::vector<CFreal> alpha = m_alpha;
    for (CFuint iEq = 0; iEq < nbEqs; ++iEq) { m_alpha[iEq] *= omega; }
    StdUpdateSol::execute();
    m_alpha = alpha;
  }
  else {
    StdUpdateSol::execute();
  }

  // 4. CFL schedule: grow only after a clean update in an unretried step
  SafePtr<CFL> cfl = getMethodData().getCFL();
  const CFreal cflValue = cfl->getCFLValue();
  if (omega < 1. || m_rejectedThisStep) {
    CFLog(INFO, "PhysicalityCFLUpdateSol: omega = " << omega
          << (m_rejectedThisStep ? " after a retry" : "")
          << ", CFL held at " << cflValue << "\n");
  }
  else if (m_cflMax > 0. && cflValue < m_cflMax) {
    const CFreal newCFL = std::min(cflValue*m_cflGrowth, m_cflMax);
    getMethodData().setNextCFL(newCFL);
    CFLog(VERBOSE, "PhysicalityCFLUpdateSol: omega = 1, CFL " << cflValue
          << " -> " << newCFL << "\n");
  }
}

//////////////////////////////////////////////////////////////////////////////

void PhysicalityCFLUpdateSol::rejectUpdate(const CFreal omega)
{
  ++m_nbConsecutiveRejections;
  m_rejectedThisStep = true;

  SafePtr<CFL> cfl = getMethodData().getCFL();
  const CFreal oldCFL = cfl->getCFLValue();
  const CFreal newCFL = oldCFL*m_cflBackoff;

  CFLog(INFO, "PhysicalityCFLUpdateSol: omega = " << omega << " < OmegaMin = "
        << m_omegaMin << ", update rejected, CFL " << oldCFL << " -> " << newCFL
        << " (rejection " << m_nbConsecutiveRejections << "/" << m_maxRejections
        << ")\n");

  if (m_nbConsecutiveRejections >= m_maxRejections) {
    throw BadValueException(FromHere(),
      "PhysicalityCFLUpdateSol: " + StringOps::to_str(m_maxRejections) +
      " consecutive updates rejected (omega < OmegaMin) down to CFL " +
      StringOps::to_str(newCFL) +
      ": the step size is no longer the limiting factor, the discrete state is nonphysical or unreachable");
  }

  cfl->setCFLValue(newCFL);

  // The states are untouched. Keep the update norm above any stop target,
  // drop the accumulated wave speeds and rewind the Newton iteration counter
  // so the loop re-assembles and re-solves this iteration at the new CFL.
  socket_rhs.getDataHandle() = 1.e6;
  socket_updateCoeff.getDataHandle() = 0.;
  getMethodData().getConvergenceStatus().iter -= 1;
}

//////////////////////////////////////////////////////////////////////////////

/// rho and p pass their floors and are finite: an overflowed trial state
/// (rho = inf, p = inf) would otherwise pass the floor test
static bool isAdmissible(const CFreal rho, const CFreal p,
                         const CFreal rhoFloor, const CFreal pFloor)
{
  return std::isfinite(rho) && std::isfinite(p) && rho >= rhoFloor && p >= pFloor;
}

//////////////////////////////////////////////////////////////////////////////

CFreal PhysicalityCFLUpdateSol::computeStateOmega(const State& state, const CFreal* dU)
{
  // Bounded state variables first: the change is linear in omega, so the
  // admissible scale is exact. It also bounds the segment on which rho and p
  // are checked, so p and T of Puvt are never evaluated past their bound.
  CFreal omegaTrial = 1.;
  for (CFuint i = 0; i < m_boundedVars.size(); ++i) {
    const CFuint iEq = m_boundedVars[i];
    const CFreal change = std::abs(m_alpha[iEq]*dU[iEq]);
    if (!std::isfinite(change)) { return 0.; }
    const CFreal allowed = m_boundedEtaMax*std::abs(state[iEq]);
    if (change > allowed) {
      omegaTrial = std::min(omegaTrial, (allowed > 0.) ? allowed/change : 0.);
    }
  }

  // bounds of a subclass on single entries: omegaTrial bounds the segment of
  // the rho and p check, omegaLinear is only combined with the result
  CFreal omegaLinear = 1.;
  boundEntrySteps(state, dU, omegaTrial, omegaLinear);
  if (!(omegaTrial > 0.)) { return 0.; }

  CFreal rho0 = 0., p0 = 0.;
  computeRhoP(state, rho0, p0);
  if (!(rho0 > 0.) || !(p0 > 0.)) {
    return 0.;
  }
  const CFreal rhoFloor = (1. - m_etaMax)*rho0;
  const CFreal pFloor   = (1. - m_etaMax)*p0;

  CFreal rho = 0., p = 0.;
  setTrialState(state, dU, omegaTrial);
  computeRhoP(*m_trial, rho, p);
  if (isAdmissible(rho, p, rhoFloor, pFloor)) {
    return std::min(omegaTrial, omegaLinear);
  }

  // The admissible set along the segment is an interval starting at the
  // state: bisect for its end. A non-finite update fails every test and
  // leaves omega at 0.
  CFreal lo = 0., hi = omegaTrial;
  for (CFuint it = 0; it < 30; ++it) {
    const CFreal mid = 0.5*(lo + hi);
    setTrialState(state, dU, mid);
    computeRhoP(*m_trial, rho, p);
    if (isAdmissible(rho, p, rhoFloor, pFloor)) { lo = mid; }
    else                                        { hi = mid; }
  }
  return std::min(lo, omegaLinear);
}

//////////////////////////////////////////////////////////////////////////////

bool PhysicalityCFLUpdateSol::fullStepEvaluable(const State& state, const State& trial) const
{
  // a full step that flips the sign of a bounded variable (p or T of Puvt)
  // leaves the set where the physical data can be evaluated
  for (CFuint i = 0; i < m_boundedVars.size(); ++i) {
    const CFuint iEq = m_boundedVars[i];
    if (!(trial[iEq]*state[iEq] > 0.)) { return false; }
  }
  return true;
}

//////////////////////////////////////////////////////////////////////////////

void PhysicalityCFLUpdateSol::computeRhoP(const State& state, CFreal& rho, CFreal& p)
{
  m_varSet->computePhysicalData(state, m_pdata);
  rho = m_pdata[m_rhoIndex];
  p   = m_pdata[m_pIndex];
}

//////////////////////////////////////////////////////////////////////////////

void PhysicalityCFLUpdateSol::setTrialState(const State& state, const CFreal* dU,
                                         const CFreal omega)
{
  // MHD read the state coordinates while filling the physical data, 
  // so the trial has to borrow the node of the state it copies.
  // The node is detached again before m_trial is deleted.
  m_trial->setSpaceCoordinates(state.getNodePtr());

  const CFuint nbEqs = state.size();
  for (CFuint iEq = 0; iEq < nbEqs; ++iEq) {
    (*m_trial)[iEq] = state[iEq] + omega*m_alpha[iEq]*dU[iEq];
  }
  filterTrialState(*m_trial);
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace NewtonMethod

  } // namespace Numerics

} // namespace COOLFluiD
