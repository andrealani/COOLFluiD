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
   * compute the Jacobians of both neighbor cells
   */
  void computeBothJacobs();
  
  /**
   * compute the jacobians of one neighbor cell
   */
  void computeOneJacob(const CFuint side);
  
  /**
   * compute the contribution of the convective boundary flux correction to the Jacobian
   */
  void extrapolatePerturbedState();
  
  /**
   * store backups of values before perturbing the states
   */
  void storeBackups();
  
  /**
   * restore values after perturbing a state
   */
  void restoreFromBackups();
  
  /**
   * store backups of values before perturbing the states for loop over cells
   */
  void storeBackupsCell();
  
  /**
   * restore values after perturbing a state for loop over cells
   */
  void restoreFromBackupsCell();
  
  /// compute the perturbed interface flux
  virtual void computePertInterfaceFlxCorrection();
  
  /// compute the total perturbed correction
  void computePertCorrection(CFuint sided, RealVector& corrections);
  
  /// compute the divergence of the perturbed discontinuous flx (-divFD+divhFD)
  void computePertDivDiscontFlx(RealVector& residuals);

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
  
  /// Perturbed divergence of the continuous flux at the solution points of the neighbours
  std::vector< std::vector< RealVector> > m_pertDivContFlx;
  
  /// backup of extrapolated states in the flux points of the cell
  std::vector< std::vector< RealVector > > m_cellStatesFlxPntBackup;
  
  /// flux pnt idx influenced by the perturbation
  CFuint m_influencedFlxPnt;

  /// influenced flx pnts idx (by perturbation)
  std::vector< CFuint> m_influencedFlxPnts;

  /// Number of influenced flx pnts (by perturbation)
  CFuint m_NbInfluencedFlxPnts;

  /// backup of interface fluxes at the flux points of a face
  std::vector< std::vector< RealVector> > m_flxPntRiemannFluxBackup;
  
  /// index of the perturbed solution point
  CFuint m_pertSol;
  
  /// index of the perturbed variable
  CFuint m_pertVar;
  
  /// index of the side which is being perturbed
  CFuint m_pertSide;

  /// Element shape
  CFGeoShape::Type elemShape;
  
  /// the discontinuous flux extrapolated to the flux points backup
  std::vector< RealVector > m_extrapolatedFluxesBackup;
  
  /// Continuous flux at the solution points backup
  std::vector< std::vector< RealVector> > m_contFlxBackup;
  
  /// flags to track which solution points have been updated (to avoid double-counting)
  std::vector< bool > m_solPntUpdated;

}; // class ConvRHSJacobFluxReconstruction

//////////////////////////////////////////////////////////////////////////////

  } // namespace FluxReconstructionMethod

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_FluxReconstructionMethod_ComputeConvVolTerms_hh
