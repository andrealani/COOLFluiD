#ifndef COOLFluiD_FluxReconstructionMethod_ConvBndCorrectionsRHSJacobFluxReconstructionBlending_hh
#define COOLFluiD_FluxReconstructionMethod_ConvBndCorrectionsRHSJacobFluxReconstructionBlending_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/DataSocketSink.hh"
#include "FluxReconstructionMethod/BCStateComputer.hh"
#include "FluxReconstructionMethod/FluxReconstructionSolverData.hh"
#include "FluxReconstructionMethod/RiemannFlux.hh"
#include "FluxReconstructionMethod/BaseCorrectionFunction.hh"
#include "FluxReconstructionMethod/ConvBndCorrectionsRHSFluxReconstructionBlending.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  
  namespace Framework {
    class BlockAccumulator;
  }

    namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

  /**
   * This class represents a command that computes contribution of the boundary face terms for the
   * Flux Reconstruction schemes for convective terms to the RHS with order blending for implicit schemes
   *
   * @author Rayan Dhib
   * 
   * based on the class ConvBndCorrectionsRHSJacobFluxReconstruction
   * @author Ray Vandenhoeck
   * @author Alexander Papen
   *
   */
class ConvBndCorrectionsRHSJacobFluxReconstructionBlending : public ConvBndCorrectionsRHSFluxReconstructionBlending {

public:
  typedef Framework::BaseMethodCommandProvider<
      FluxReconstructionSolverData,ConvBndCorrectionsRHSJacobFluxReconstructionBlending > PROVIDER;

public:

  /**
   * Constructor
   */
  ConvBndCorrectionsRHSJacobFluxReconstructionBlending(const std::string& name);

  /**
   * Default destructor
   */
  virtual ~ConvBndCorrectionsRHSJacobFluxReconstructionBlending();

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
  
  /// compute the perturbed interface flux
  virtual void computePertInterfaceFlxCorrection();
  
  /// compute the total perturbed correction
  void computePertCorrection(std::vector< RealVector >& corrections);

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
  
  /// index of the perturbed solution point
  CFuint m_pertSol;
  
  /// index of the perturbed variable
  CFuint m_pertVar;
  
  /// dependencies of sol pnts on flx pnts
  Common::SafePtr< std::vector< std::vector< CFuint > > > m_solFlxDep;

  /// nbr of flx pnts on which a sol pnt is dependent
  CFuint m_nbrFlxDep;
  
  /// backup of extrapolated states in the flux points of the cell
  std::vector< RealVector > m_cellStatesFlxPntBackup;
  
  /// influenced flx pnt idx (by perturbation)
  CFuint m_influencedFlxPnt;

  /// influenced flx pnts idx (by perturbation)
  std::vector< CFuint> m_influencedFlxPnts;
  
  /// Number of influenced flx pnts (by perturbation)
  CFuint m_NbInfluencedFlxPnts;
  
  /// Element shape
  CFGeoShape::Type elemShape;

  /// backup of interface fluxes at the flux points of a face
  std::vector< RealVector> m_flxPntRiemannFluxBackup;

  std::vector< RealVector> m_flxPntRiemannFluxBackupP0;
  std::vector< RealVector > m_cellStatesFlxPntBackupP0;

  std::vector< Framework::State* >* m_cellStatesP0;

  Common::SafePtr<Framework::NumericalJacobian> m_numJacob_P0;

}; // end of class ConvBndCorrectionsRHSJacobFluxReconstructionBlending

//////////////////////////////////////////////////////////////////////////////

 } // namespace FluxReconstructionMethod

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_FluxReconstructionMethod_ConvBndCorrectionsRHSJacobFluxReconstructionBlending_hh
