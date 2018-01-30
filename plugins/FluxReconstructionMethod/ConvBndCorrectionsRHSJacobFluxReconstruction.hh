#ifndef COOLFluiD_FluxReconstructionMethod_ConvBndCorrectionsRHSJacobFluxReconstruction_hh
#define COOLFluiD_FluxReconstructionMethod_ConvBndCorrectionsRHSJacobFluxReconstruction_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/DataSocketSink.hh"
#include "FluxReconstructionMethod/BCStateComputer.hh"
#include "FluxReconstructionMethod/FluxReconstructionSolverData.hh"
#include "FluxReconstructionMethod/RiemannFlux.hh"
#include "FluxReconstructionMethod/BaseCorrectionFunction.hh"
#include "FluxReconstructionMethod/ConvBndCorrectionsRHSFluxReconstruction.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  
  namespace Framework {
    class BlockAccumulator;
  }

    namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

  /**
   * This class represents a command that computes contribution of the boundary face terms for the
   * Flux Reconstruction schemes for convective terms to the RHS for implicit schemes
   *
   * @author Ray Vandenhoeck
   * @author Alexander Papen
   *
   */
class ConvBndCorrectionsRHSJacobFluxReconstruction : public ConvBndCorrectionsRHSFluxReconstruction {

public:
  typedef Framework::BaseMethodCommandProvider<
      FluxReconstructionSolverData,ConvBndCorrectionsRHSJacobFluxReconstruction > PROVIDER;

public:

  /**
   * Constructor
   */
  ConvBndCorrectionsRHSJacobFluxReconstruction(const std::string& name);

  /**
   * Default destructor
   */
  virtual ~ConvBndCorrectionsRHSJacobFluxReconstruction();

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
   * compute the contribution of the convective boundary flux correction to the Jacobian
   */
  void computeJacobConvBndCorrection();


protected: // data
  
  /// pointer to the linear system solver
  Common::SafePtr<Framework::LinearSystemSolver> m_lss;

  /// pointer to the numerical Jacobian computer
  Common::SafePtr<Framework::NumericalJacobian> m_numJacob;

  /// accumulator for LSSMatrix
  std::auto_ptr<Framework::BlockAccumulator> m_acc;

  /// perturbed updates to the residuals
  RealVector m_pertResUpdates;
  
  /// unperturbed updates to the residuals
  RealVector m_resUpdates;

  /// derivative of update to one CV-residual
  RealVector m_derivResUpdates;
  
  /// perturbed corrections due to the boundary faces for the Jacobian
  std::vector< RealVector> m_pertCorrections;

}; // end of class ConvBndCorrectionsRHSJacobFluxReconstruction

//////////////////////////////////////////////////////////////////////////////

 } // namespace FluxReconstructionMethod

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_FluxReconstructionMethod_ConvBndCorrectionsRHSJacobFluxReconstruction_hh
