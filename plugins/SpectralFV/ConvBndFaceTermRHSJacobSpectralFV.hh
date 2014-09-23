#ifndef COOLFluiD_SpectralFV_ConvBndFaceTermRHSJacobSpectralFV_hh
#define COOLFluiD_SpectralFV_ConvBndFaceTermRHSJacobSpectralFV_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/DataSocketSink.hh"
#include "SpectralFV/ConvBndFaceTermRHSSpectralFV.hh"
#include "SpectralFV/SpectralFVMethodData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {
    class BlockAccumulator;
  }

    namespace SpectralFV {

//////////////////////////////////////////////////////////////////////////////


/**
  * This class represents a command that computes contribution of the boundary face terms for the
  * spectral finite volume schemes for convection terms, for both the RHS and the Jacobian
  *
  * @author Kris Van Den Abeele
  *
  */
class ConvBndFaceTermRHSJacobSpectralFV : public ConvBndFaceTermRHSSpectralFV {

public:
  typedef Framework::BaseMethodCommandProvider<
      SpectralFVMethodData,ConvBndFaceTermRHSJacobSpectralFV > PROVIDER;

public:

  /**
   * Constructor
   */
  ConvBndFaceTermRHSJacobSpectralFV(const std::string& name);

  /**
   * Default destructor
   */
  virtual ~ConvBndFaceTermRHSJacobSpectralFV();

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

  /**
   * compute the contribution of the convective boundary face term to the Jacobian
   */
  void computeJacobConvBndFaceTerm();

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

}; // end of class ConvBndFaceTermRHSJacobSpectralFV

//////////////////////////////////////////////////////////////////////////////

 } // namespace SpectralFV

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_SpectralFV_ConvBndFaceTermRHSJacobSpectralFV_hh
