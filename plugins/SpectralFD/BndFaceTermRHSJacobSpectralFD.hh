#ifndef COOLFluiD_SpectralFD_BndFaceTermRHSJacobSpectralFD_hh
#define COOLFluiD_SpectralFD_BndFaceTermRHSJacobSpectralFD_hh

//////////////////////////////////////////////////////////////////////////////

#include "SpectralFD/BndFaceTermRHSSpectralFD.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace SpectralFD {

//////////////////////////////////////////////////////////////////////////////

  /**
   * This class represents a command that computes contribution of the boundary face terms for the
   * spectral finite difference schemes for both convective and diffusive terms to the RHS and the
   * Jacobian matrix if the full scheme is compact.
   *
   * @author Kris Van Den Abeele
   *
   */
class BndFaceTermRHSJacobSpectralFD : public BndFaceTermRHSSpectralFD {

public:
  typedef Framework::BaseMethodCommandProvider<
      SpectralFDMethodData,BndFaceTermRHSJacobSpectralFD > PROVIDER;

public:

  /**
   * Constructor
   */
  BndFaceTermRHSJacobSpectralFD(const std::string& name);

  /**
   * Default destructor
   */
  virtual ~BndFaceTermRHSJacobSpectralFD();

  /**
   * Set up private data and data of the aggregated classes
   * in this command before processing phase
   */
  virtual void setup();

  /**
   * Unset up private data and data of the aggregated classes
   * in this command after processing phase
   */
  virtual void unsetup();

  /**
   * Configures the command.
   */
  virtual void configure ( Config::ConfigArgs& args );

protected: // functions

  /**
   * Execute on the current TRS
   */
  virtual void executeOnTrs();

  /// backup and reconstruct a physical variable and its gradient in the face flux points
  void backupAndReconstructFacePhysVarsAndGrad(const CFuint iVar);

  /// restore a physical variable and its gradient in the face flux points
  void restoreFacePhysVarsAndGrad(const CFuint iVar);

  /// compute the contribution of the boundary face term to the Jacobian
  void computeJacobBndFaceTerm();

protected: // data

  /// pointer to the linear system solver
  Common::SafePtr<Framework::LinearSystemSolver> m_lss;

  /// pointer to the numerical Jacobian computer
  Common::SafePtr<Framework::NumericalJacobian> m_numJacob;

  /// accumulator for LSSMatrix
  std::auto_ptr<Framework::BlockAccumulator> m_acc;

  /// perturbed updates to the residuals
  RealVector m_pertResUpdates;

  /// derivative of update to one CV-residual
  RealVector m_derivResUpdates;

}; // end of class BndFaceTermRHSJacobSpectralFD

//////////////////////////////////////////////////////////////////////////////

 } // namespace SpectralFD

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_SpectralFD_BndFaceTermRHSJacobSpectralFD_hh
