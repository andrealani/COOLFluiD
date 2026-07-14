#include "FluxReconstructionEntropyMHD/EntropyMHDSourceTerm.hh"
#include "FluxReconstructionEntropyMHD/FluxReconstructionEntropyMHD.hh"
#include "Framework/MethodCommandProvider.hh"
#include "EntropyMHD/EntropyMHDTerm.hh"
#include "EntropyMHD/EntropyMHD3DVarSet.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<EntropyMHDSourceTerm, FluxReconstructionSolverData, FluxReconstructionEntropyMHDModule>
EntropyMHDSourceTermProvider("EntropyMHDSourceTerm");

//////////////////////////////////////////////////////////////////////////////

void EntropyMHDSourceTerm::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< bool >("Gravity", "Switch on gravity source term.");
  options.addConfigOption< bool >("Rotation", "Switch on rotation (Coriolis + centrifugal) source term.");
  options.addConfigOption< CFreal >("OmegaSun", "Non-dimensional solar rotation rate Omega*l0/V0 (default: 0.0409 for helio normalization).");
  options.addConfigOption< bool >("DednerDamping",
    "Enable Dedner parabolic damping on the GLM phi field. Default false. "
    "Required when the inlet phi BC is set to absorbing (zero-gradient); "
    "without damping the cleaning channel has no volumetric sink.");
  options.addConfigOption< CFreal >("DednerCr",
    "Dedner damping coefficient c_r (dimensionless, ~0.1-0.5). "
    "Default 0.18 per Mignone & Tzeferacos 2010 (JCP 229:5896).");
}

//////////////////////////////////////////////////////////////////////////////

EntropyMHDSourceTerm::EntropyMHDSourceTerm(const std::string& name)
  : StdSourceTerm(name),
    socket_updateCoeff("updateCoeff"),
    m_gravFactor(0.0),
    m_model(CFNULL),
    m_varSet(CFNULL)
{
  addConfigOptionsTo(this);

  m_gravity = true;
  setParameter("Gravity", &m_gravity);

  m_rotation = false;
  setParameter("Rotation", &m_rotation);

  m_omegaSun = 0.0409;
  setParameter("OmegaSun", &m_omegaSun);

  m_dednerDamping = false;
  setParameter("DednerDamping", &m_dednerDamping);

  m_dednerCr = 0.18;
  setParameter("DednerCr", &m_dednerCr);
}

//////////////////////////////////////////////////////////////////////////////

EntropyMHDSourceTerm::~EntropyMHDSourceTerm()
{
}

//////////////////////////////////////////////////////////////////////////////

void EntropyMHDSourceTerm::setup()
{
  StdSourceTerm::setup();

  // Get EntropyMHDTerm for gravity factor and equilibrium
  m_model = PhysicalModelStack::getActive()->getImplementor()
    ->getConvectiveTerm().d_castTo<Physics::EntropyMHD::EntropyMHDTerm>();

  m_gravFactor = m_model->getGravFactorNonDim();

  m_varSet = getMethodData().getUpdateVar();
  m_model->resizePhysicalData(m_solPhysData);

  // Guard against silent refSpeed misconfiguration: Dedner damping rate is
  // c_r * c_h, where c_h = EntropyMHDTerm::refSpeed. The term default is 1.0,
  // which is rarely the intended GLM speed for a real run.
  if (m_dednerDamping && std::abs(m_model->getRefSpeed() - 1.0) < 1.0e-10) {
    CFLog(WARN, "EntropyMHDSourceTerm: DednerDamping=true with refSpeed=1.0 "
                "(default). Set the model refSpeed explicitly so the "
                "Dedner damping rate matches the GLM characteristic speed.\n");
  }
}

//////////////////////////////////////////////////////////////////////////////

void EntropyMHDSourceTerm::unsetup()
{
  StdSourceTerm::unsetup();
}

//////////////////////////////////////////////////////////////////////////////

vector<SafePtr<BaseDataSocketSink> > EntropyMHDSourceTerm::needsSockets()
{
  vector<SafePtr<BaseDataSocketSink> > result = StdSourceTerm::needsSockets();
  result.push_back(&socket_updateCoeff);
  return result;
}

//////////////////////////////////////////////////////////////////////////////

void EntropyMHDSourceTerm::addSourceTerm(RealVector& resUpdates)
{
  const CFuint nbrStates = m_cellStates->size();

  for (CFuint iSol = 0; iSol < nbrStates; ++iSol)
  {
    // Initialize all source terms to zero
    for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
    {
      resUpdates[m_nbrEqs*iSol + iEq] = 0.0;
    }

    // Gravity source (momentum eqs 1-3 only)
    // Eq 7 (entropy sigma): NO source; d(sigma)/dt + div(sigma*v) = 0 is exact
    if (m_gravity)
    {
      m_varSet->computePhysicalData(*((*m_cellStates)[iSol]), m_solPhysData);

      const CFreal rho = m_solPhysData[Physics::EntropyMHD::EntropyMHDTerm::RHO];
      const CFreal x   = m_solPhysData[Physics::EntropyMHD::EntropyMHDTerm::XP];
      const CFreal y   = m_solPhysData[Physics::EntropyMHD::EntropyMHDTerm::YP];
      const CFreal z   = m_solPhysData[Physics::EntropyMHD::EntropyMHDTerm::ZP];

      const CFreal r2 = x*x + y*y + z*z;
      const CFreal r  = sqrt(r2);
      const CFreal gFactor = -m_gravFactor / (r2 * r);

      resUpdates[m_nbrEqs*iSol + 1] += rho * gFactor * x;
      resUpdates[m_nbrEqs*iSol + 2] += rho * gFactor * y;
      resUpdates[m_nbrEqs*iSol + 3] += rho * gFactor * z;
    }

    // Rotation source (Coriolis + centrifugal); momentum only
    // Rotation axis = z, rate = m_omegaSun (non-dimensional)
    // Coriolis:    a_cor  = -2 Omega x v
    // Centrifugal: a_cent = -Omega^2 (x, y, 0)
    // Entropy sigma: NO source; rotation is a conservative force, no entropy production
    if (m_rotation)
    {
      if (!m_gravity) {
        m_varSet->computePhysicalData(*((*m_cellStates)[iSol]), m_solPhysData);
      }

      const CFreal rho = m_solPhysData[Physics::EntropyMHD::EntropyMHDTerm::RHO];
      const CFreal x   = m_solPhysData[Physics::EntropyMHD::EntropyMHDTerm::XP];
      const CFreal y   = m_solPhysData[Physics::EntropyMHD::EntropyMHDTerm::YP];
      const CFreal Vx  = m_solPhysData[Physics::EntropyMHD::EntropyMHDTerm::VX];
      const CFreal Vy  = m_solPhysData[Physics::EntropyMHD::EntropyMHDTerm::VY];

      // Coriolis: -2*rho*(Omega x V) = -2*rho*Omega*(-Vy, Vx, 0) for Omega along z
      const CFreal cor_x =  2.0 * m_omegaSun * Vy;
      const CFreal cor_y = -2.0 * m_omegaSun * Vx;

      // Centrifugal: +rho*Omega^2*(x, y, 0) [outward, = -Omega x (Omega x r)]
      const CFreal O2 = m_omegaSun * m_omegaSun;
      const CFreal cent_x = O2 * x;
      const CFreal cent_y = O2 * y;

      resUpdates[m_nbrEqs*iSol + 1] += rho * (cor_x + cent_x);
      resUpdates[m_nbrEqs*iSol + 2] += rho * (cor_y + cent_y);
    }

    // Dedner parabolic GLM damping (eq 8, phi):
    //   d(phi)/dt + c_h^2 div(B) = -(c_r * c_h) * phi
    // The pure-hyperbolic GLM (c_h^2 div(B) in the flux) has no volumetric
    // sink on its own; without this term, residual div(B) production at
    // helmet streamer cusps drives phi to accumulate. Coupled with the
    // absorbing inlet BC (ghost = +int), this gives Dedner 2002 Sec. 3.3
    // mixed GLM-MHD cleaning. Calibration c_r = 0.18 from
    // Mignone & Tzeferacos 2010 (JCP 229:5896).
    if (m_dednerDamping)
    {
      if (!m_gravity && !m_rotation) {
        m_varSet->computePhysicalData(*((*m_cellStates)[iSol]), m_solPhysData);
      }
      const CFreal phi = m_solPhysData[Physics::EntropyMHD::EntropyMHDTerm::PHI];
      const CFreal c_h = m_model->getRefSpeed();
      resUpdates[m_nbrEqs*iSol + 8] -= m_dednerCr * c_h * phi;
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void EntropyMHDSourceTerm::getSToStateJacobian(const CFuint iState)
{
  // Reset Jacobian
  for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
  {
    m_stateJacobian[iEq] = 0.0;
  }

  if (m_gravity)
  {
    m_varSet->computePhysicalData(*((*m_cellStates)[iState]), m_solPhysData);

    const CFreal x = m_solPhysData[Physics::EntropyMHD::EntropyMHDTerm::XP];
    const CFreal y = m_solPhysData[Physics::EntropyMHD::EntropyMHDTerm::YP];
    const CFreal z = m_solPhysData[Physics::EntropyMHD::EntropyMHDTerm::ZP];

    const CFreal r2 = x*x + y*y + z*z;
    const CFreal r  = sqrt(r2);
    const CFreal gFactor = -m_gravFactor / (r2 * r);

    // Momentum: S_mom_i = rho*g_i -> dS_mom_i/d(rho) = g_i
    m_stateJacobian[0][1] = gFactor * x;
    m_stateJacobian[0][2] = gFactor * y;
    m_stateJacobian[0][3] = gFactor * z;

    // Eq 7 (entropy sigma): no source term, no Jacobian entries
  }

  if (m_rotation)
  {
    if (!m_gravity) {
      m_varSet->computePhysicalData(*((*m_cellStates)[iState]), m_solPhysData);
    }

    const CFreal rho = m_solPhysData[Physics::EntropyMHD::EntropyMHDTerm::RHO];
    const CFreal x   = m_solPhysData[Physics::EntropyMHD::EntropyMHDTerm::XP];
    const CFreal y   = m_solPhysData[Physics::EntropyMHD::EntropyMHDTerm::YP];
    const CFreal Vx  = m_solPhysData[Physics::EntropyMHD::EntropyMHDTerm::VX];
    const CFreal Vy  = m_solPhysData[Physics::EntropyMHD::EntropyMHDTerm::VY];

    const CFreal O2 = m_omegaSun * m_omegaSun;

    // dS_mom_x/d(rho) = cor_x + cent_x = 2*Omega*Vy + Omega^2*x
    // dS_mom_y/d(rho) = cor_y + cent_y = -2*Omega*Vx + Omega^2*y
    m_stateJacobian[0][1] += (2.0 * m_omegaSun * Vy + O2 * x);
    m_stateJacobian[0][2] += (-2.0 * m_omegaSun * Vx + O2 * y);

    // dS_mom_x/d(Vy) = 2*rho*Omega (Coriolis)
    m_stateJacobian[2][1] += 2.0 * rho * m_omegaSun;
    // dS_mom_y/d(Vx) = -2*rho*Omega (Coriolis)
    m_stateJacobian[1][2] += -2.0 * rho * m_omegaSun;
  }

  // Dedner damping Jacobian: dS_phi / dphi = -c_r * c_h
  // Local linear damping on slot 8. c_h is model-level, so no
  // computePhysicalData call is needed here.
  if (m_dednerDamping)
  {
    const CFreal c_h = m_model->getRefSpeed();
    m_stateJacobian[8][8] += -m_dednerCr * c_h;
  }
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace FluxReconstructionMethod

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
