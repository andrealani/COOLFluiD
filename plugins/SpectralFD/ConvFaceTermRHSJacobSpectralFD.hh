#ifndef COOLFluiD_SpectralFD_ConvFaceTermRHSJacobSpectralFD_hh
#define COOLFluiD_SpectralFD_ConvFaceTermRHSJacobSpectralFD_hh

//////////////////////////////////////////////////////////////////////////////

#include "SpectralFD/ConvFaceTermRHSSpectralFD.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {
    class BlockAccumulator;
  }

  namespace SpectralFD {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a command that computes contribution of the face terms for the
 * spectral finite difference schemes for convection terms, for both the RHS and the Jacobian
 *
 * @author Kris Van den Abeele
 *
 */
class ConvFaceTermRHSJacobSpectralFD : public ConvFaceTermRHSSpectralFD {
public:

  /**
   * Constructor.
   */
  explicit ConvFaceTermRHSJacobSpectralFD(const std::string& name);

  /**
   * Destructor.
   */
  virtual ~ConvFaceTermRHSJacobSpectralFD();

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

  /**
   * compute the contribution of the convective face term to both Jacobians
   */
  void computeBothJacobsConvFaceTerm();

  /**
   * compute the contribution of the convective face term to one Jacobians
   */
  void computeOneJacobConvFaceTerm(const CFuint side);

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

}; // class ConvFaceTermRHSJacobSpectralFD

//////////////////////////////////////////////////////////////////////////////

  } // namespace SpectralFD

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_SpectralFD_ComputeConvVolTerms_hh
