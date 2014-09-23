#ifndef COOLFluiD_SpectralFD_VolTermRHSJacobSpectralFD_hh
#define COOLFluiD_SpectralFD_VolTermRHSJacobSpectralFD_hh

//////////////////////////////////////////////////////////////////////////////

#include "SpectralFD/VolTermRHSSpectralFD.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace SpectralFD {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a basic command that computes the volume terms for the
 * spectral finite difference schemes for both convective and diffusive terms
 * for the rhs and the Jacobian matrix if the full scheme is compact.
 *
 * @author Kris Van den Abeele
 *
 */
class VolTermRHSJacobSpectralFD : public VolTermRHSSpectralFD {
public:

  /**
   * Constructor.
   */
  explicit VolTermRHSJacobSpectralFD(const std::string& name);

  /**
   * Destructor.
   */
  virtual ~VolTermRHSJacobSpectralFD();

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

  /// resize the residual and gradient updates corresponding to the current element type
  void resizeResAndGradUpdates();

  /// set the required data in the volume and face term computers
  void setVolumeAndFaceTermComputersData();

  /// compute cell gradients minus boundary face terms
  void computeCellGradsMinusBndFaceTerms();

  /// compute and reconstruct perturbed cell physical variable gradient
  void computeAndReconstructPertPhysVarGrad(const CFuint iVar);

  /// back up and reconstruct variables in the volume and face term flux points
  void backupAndReconstructVolumeAndFacePhysVars(const CFuint iVar);

  /// restore variables the volume and face term flux points
  void restoreVolumeAndFacePhysVars(const CFuint iVar);

  /// compute the contribution of the volume term to the Jacobian
  void computeJacobVolTerm();

protected: // data

  /// pointer to the linear system solver
  Common::SafePtr<Framework::LinearSystemSolver> m_lss;

  /// pointer to the numerical Jacobian computer
  Common::SafePtr<Framework::NumericalJacobian> m_numJacob;

  /// accumulator for LSSMatrix
  std::auto_ptr<Framework::BlockAccumulator> m_acc;

  /// perturbed updates to the residuals
  RealVector m_pertResUpdates;

  /// derivative of update to the residuals
  RealVector m_derivResUpdates;

  /// gradients minus boundary face terms
  std::vector< std::vector< RealVector > > m_gradsMBndFaceTerm;

  /// physical variable gradient
  std::vector< std::vector< RealVector > > m_physVarGradUpdates;

}; // class VolTermRHSJacobSpectralFD

//////////////////////////////////////////////////////////////////////////////

  } // namespace SpectralFD

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_SpectralFD_VolTermRHSJacobSpectralFD_hh
