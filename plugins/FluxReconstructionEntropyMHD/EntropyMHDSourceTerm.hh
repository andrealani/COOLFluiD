#ifndef COOLFluiD_FluxReconstructionEntropyMHD_EntropyMHDSourceTerm_hh
#define COOLFluiD_FluxReconstructionEntropyMHD_EntropyMHDSourceTerm_hh

//////////////////////////////////////////////////////////////////////////////

#include "FluxReconstructionMethod/StdSourceTerm.hh"
#include "Framework/ConvectiveVarSet.hh"
#include "EntropyMHD/EntropyMHDTerm.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

/**
 * Source terms for the entropy-variable MHD model (EntropyMHD):
 *  - Gravity: rho*g in momentum (eqs 1-3)
 *  - Eq 7 (entropy sigma = p*rho^(1-gamma)): NO source; d(sigma)/dt + div(sigma*v) = 0
 *
 * Gravity does NOT appear in eq 7. The entropy equation is a pure
 * conservation law derived from dp/dt + gamma*p*div(v) = 0 and continuity.
 *
 * @author Rayan Dhib
 */
class EntropyMHDSourceTerm : public StdSourceTerm {
public:

  static void defineConfigOptions(Config::OptionList& options);

  EntropyMHDSourceTerm(const std::string& name);
  virtual ~EntropyMHDSourceTerm();

  virtual void setup();
  virtual void unsetup();

  virtual void addSourceTerm(RealVector& resUpdates);

  virtual void getSToStateJacobian(const CFuint iState);

  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();

private:

  /// sink socket for updateCoeff
  Framework::DataSocketSink<CFreal> socket_updateCoeff;

  /// switch gravity on/off
  bool m_gravity;

  /// switch rotation (Coriolis + centrifugal) on/off
  bool m_rotation;

  /// non-dimensional solar rotation rate Omega*l0/V0
  CFreal m_omegaSun;

  /// switch Dedner parabolic damping on the GLM phi field on/off
  bool m_dednerDamping;

  /// Dedner damping coefficient c_r (dimensionless, ~0.1-0.5)
  /// Damping rate alpha_p = c_r * c_h / L_ref with L_ref = 1 (non-dim).
  /// Default 0.18 per Mignone & Tzeferacos 2010 (JCP 229:5896) calibration.
  CFreal m_dednerCr;

  /// non-dimensional gravity factor: G*M_sun*mu_0*rho_ref / (B_ref^2 * L_ref)
  CFreal m_gravFactor;

  /// EntropyMHD physical term (for equilibrium density)
  Common::SafePtr<Physics::EntropyMHD::EntropyMHDTerm> m_model;

  /// update variable set (for computePhysicalData)
  Common::SafePtr<Framework::ConvectiveVarSet> m_varSet;

  /// physical data array
  RealVector m_solPhysData;

}; // end of class EntropyMHDSourceTerm

//////////////////////////////////////////////////////////////////////////////

  } // namespace FluxReconstructionMethod

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_FluxReconstructionEntropyMHD_EntropyMHDSourceTerm_hh
