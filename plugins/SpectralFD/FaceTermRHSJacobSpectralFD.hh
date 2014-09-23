#ifndef COOLFluiD_SpectralFD_FaceTermRHSJacobSpectralFD_hh
#define COOLFluiD_SpectralFD_FaceTermRHSJacobSpectralFD_hh

//////////////////////////////////////////////////////////////////////////////

#include "SpectralFD/FaceTermRHSSpectralFD.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace SpectralFD {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represent a basic command that computes the contribution of face terms for the
 * spectral finite difference schemes for both convective and diffusive terms to the RHS and the Jacobian
 * if the full scheme is compact.
 *
 * @author Kris Van den Abeele
 *
 */
class FaceTermRHSJacobSpectralFD : public FaceTermRHSSpectralFD {
public:

  /**
   * Constructor.
   */
  explicit FaceTermRHSJacobSpectralFD(const std::string& name);

  /**
   * Destructor.
   */
  virtual ~FaceTermRHSJacobSpectralFD();

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Set up private data and data of the aggregated classes
   * in this command before processing phase
   */
  virtual void setup();

  /**
   * Unsetup private data
   */
  virtual void unsetup();

  /**
   * Configures the command.
   */
  virtual void configure ( Config::ConfigArgs& args );

  /**
   * Execute Processing actions
   */
  virtual void execute();

protected: // functions

  /// backup and reconstruct a physical variable and its gradient in the face flux points
  void backupAndReconstructFacePhysVarsAndGrad(const CFuint side, const CFuint iVar);

  /// restore a physical variable and its gradient in the face flux points
  void restoreFacePhysVarsAndGrad(const CFuint side, const CFuint iVar);

  /// compute the contribution of the face term to both Jacobians
  void computeBothJacobsFaceTerm();

  /// compute the contribution of the face term to one Jacobians
  void computeOneJacobFaceTerm(const CFuint side);

protected: // data

  /// pointer to the linear system solver
  Common::SafePtr<Framework::LinearSystemSolver> m_lss;

  /// pointer to the numerical Jacobian computer
  Common::SafePtr<Framework::NumericalJacobian> m_numJacob;

  /// accumulator for LSSMatrix
  std::auto_ptr<Framework::BlockAccumulator> m_acc;

  /// perturbed updates to the residuals
  std::vector< RealVector > m_pertResUpdates;

  /// derivative of update to one CV-residual
  RealVector m_derivResUpdates;

}; // class FaceTermRHSJacobSpectralFD

//////////////////////////////////////////////////////////////////////////////

  } // namespace SpectralFD

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_SpectralFD_ComputeConvVolTerms_hh
