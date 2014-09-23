#ifndef COOLFluiD_SpectralFD_ConvVolTermRHSJacobSpectralFD_hh
#define COOLFluiD_SpectralFD_ConvVolTermRHSJacobSpectralFD_hh

//////////////////////////////////////////////////////////////////////////////

#include "SpectralFD/ConvVolTermRHSSpectralFD.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {
    class BlockAccumulator;
  }

  namespace SpectralFD {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a command that computes contribution of the volume terms for the
 * spectral finite difference schemes for convection terms, for both the RHS and the Jacobian
 *
 * @author Kris Van den Abeele
 *
 */
class ConvVolTermRHSJacobSpectralFD : public ConvVolTermRHSSpectralFD {
public:

  /**
   * Constructor.
   */
  explicit ConvVolTermRHSJacobSpectralFD(const std::string& name);

  /**
   * Destructor.
   */
  virtual ~ConvVolTermRHSJacobSpectralFD();

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Setup private data and data of the aggregated classes
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
   * set the data required to compute the volume term
   */

  void setVolumeTermData();

  /**
   * resize the variable m_pertResUpdates corresponding to the current element type
   */
  void resizeResAndGradUpdates();

  /**
   * compute the contribution of the convective volume term to the Jacobian
   */
  void computeJacobConvVolTerm();

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

}; // class ConvVolTermRHSJacobSpectralFD

//////////////////////////////////////////////////////////////////////////////

  } // namespace SpectralFD

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_SpectralFD_ConvVolTermRHSJacobSpectralFD_hh
