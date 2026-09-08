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
}

//////////////////////////////////////////////////////////////////////////////

PhysicalityCFLUpdateSol::PhysicalityCFLUpdateSol(const std::string& name) :
  StdUpdateSol(name),
  m_varSet(CFNULL),
  m_trial(CFNULL),
  m_pdata(),
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
  for (CFuint iState = 0; iState < nbStates; ++iState) {
    const State& state = *states[iState];
    if (state.isParUpdatable()) {
      omega = std::min(omega, computeStateOmega(state, &dU[iState*nbEqs]));
    }
  }
#ifdef CF_HAVE_MPI
  if (PE::GetPE().IsParallel()) {
    const std::string nsp = getMethodData().getNamespace();
    CFreal omegaLocal = omega;
    MPI_Allreduce(&omegaLocal, &omega, 1, MPI_DOUBLE, MPI_MIN,
                  PE::GetPE().GetCommunicator(nsp));
  }
#endif

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
    cfl->setCFLValue(newCFL);
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

CFreal PhysicalityCFLUpdateSol::computeStateOmega(const State& state, const CFreal* dU)
{
  CFreal rho0 = 0., p0 = 0.;
  computeRhoP(state, rho0, p0);
  if (!(rho0 > 0.) || !(p0 > 0.)) {
    return 0.;
  }
  const CFreal rhoFloor = (1. - m_etaMax)*rho0;
  const CFreal pFloor   = (1. - m_etaMax)*p0;

  CFreal rho = 0., p = 0.;
  setTrialState(state, dU, 1.);
  computeRhoP(*m_trial, rho, p);
  if (rho >= rhoFloor && p >= pFloor) {
    return 1.;
  }

  // The admissible set along the segment is an interval starting at the
  // state: bisect for its end. A non-finite update fails every test and
  // leaves omega at 0.
  CFreal lo = 0., hi = 1.;
  for (CFuint it = 0; it < 30; ++it) {
    const CFreal mid = 0.5*(lo + hi);
    setTrialState(state, dU, mid);
    computeRhoP(*m_trial, rho, p);
    if (rho >= rhoFloor && p >= pFloor) { lo = mid; }
    else                                { hi = mid; }
  }
  return lo;
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
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace NewtonMethod

  } // namespace Numerics

} // namespace COOLFluiD
