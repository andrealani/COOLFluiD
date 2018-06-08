#ifndef COOLFluiD_FluxReconstructionMethod_ConvRHSJacobFluxReconstruction_hh
#define COOLFluiD_FluxReconstructionMethod_ConvRHSJacobFluxReconstruction_hh

//////////////////////////////////////////////////////////////////////////////

#include "FluxReconstructionMethod/ConvRHSFluxReconstruction.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {
    class BlockAccumulator;
  }

  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a command that computes contribution for the
 * FR schemes for convection terms, for both the RHS and the Jacobian
 *
 * @author Ray Vandenhoeck
 * @author Alexander Papen
 *
 */
class ConvRHSJacobFluxReconstruction : public ConvRHSFluxReconstruction {
public:

  /**
   * Constructor.
   */
  explicit ConvRHSJacobFluxReconstruction(const std::string& name);

  /**
   * Destructor.
   */
  virtual ~ConvRHSJacobFluxReconstruction();

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
   * compute the contribution of the convective face term to one Jacobians
   */
  void computeJacobConvCorrection();
  
  /**
   * 
   */
  void computeBothJacobs();
  
  /**
   * 
   */
  void computeOneJacob(const CFuint side);

protected: // data

  /// pointer to the linear system solver
  Common::SafePtr<Framework::LinearSystemSolver> m_lss;

  /// pointer to the numerical Jacobian computer
  Common::SafePtr<Framework::NumericalJacobian> m_numJacob;

  /// accumulator for LSSMatrix
  std::auto_ptr<Framework::BlockAccumulator> m_acc;
  
  /// accumulator for LSSMatrix for the loop over faces
  std::auto_ptr<Framework::BlockAccumulator> m_accFace;
  
  /// perturbed updates to the residuals
  std::vector< RealVector > m_pertResUpdates;
  
  /// unperturbed updates to the residuals
  std::vector< RealVector > m_resUpdates;

  /// derivative of update to one element-residual
  RealVector m_derivResUpdates;
  
  /// perturbed corrections
  std::vector< RealVector> m_pertCorrections;
  
  /// Divergence of the continuous flux at the solution points of the left neighbour
  std::vector< RealVector> m_divContFlxL;
  
  /// Divergence of the continuous flux at the solution points of the right neighbour
  std::vector< RealVector> m_divContFlxR;
  
  /// Perturbed divergence of the continuous flux at the solution points of the neighbours
  std::vector< std::vector< RealVector> > m_pertDivContFlx;

}; // class ConvRHSJacobFluxReconstruction

//////////////////////////////////////////////////////////////////////////////

  } // namespace FluxReconstructionMethod

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_FluxReconstructionMethod_ComputeConvVolTerms_hh
