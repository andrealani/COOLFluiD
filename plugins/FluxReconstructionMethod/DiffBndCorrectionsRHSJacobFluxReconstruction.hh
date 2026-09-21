#ifndef COOLFluiD_FluxReconstructionMethod_DiffBndCorrectionsRHSJacobFluxReconstruction_hh
#define COOLFluiD_FluxReconstructionMethod_DiffBndCorrectionsRHSJacobFluxReconstruction_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/DataSocketSink.hh"
#include "FluxReconstructionMethod/BCStateComputer.hh"
#include "FluxReconstructionMethod/FluxReconstructionSolverData.hh"
#include "FluxReconstructionMethod/RiemannFlux.hh"
#include "FluxReconstructionMethod/BaseCorrectionFunction.hh"
#include "FluxReconstructionMethod/DiffBndCorrectionsRHSFluxReconstruction.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

  /**
   * This class represents a command that computes contribution of the boundary faces for the
   * Flux Reconstruction schemes for diffusive terms to the RHS for an implicit scheme
   *
   * @author Ray Vandenhoeck
   * @author Alexander Papen
   * @author Rayan Dhib
   *
   */
class DiffBndCorrectionsRHSJacobFluxReconstruction : public DiffBndCorrectionsRHSFluxReconstruction {

public:
  typedef Framework::BaseMethodCommandProvider<
      FluxReconstructionSolverData,DiffBndCorrectionsRHSJacobFluxReconstruction > PROVIDER;

public:

  /**
   * Constructor
   */
  DiffBndCorrectionsRHSJacobFluxReconstruction(const std::string& name);

  /**
   * Default destructor
   */
  virtual ~DiffBndCorrectionsRHSJacobFluxReconstruction();

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
   * compute the contribution of the diffusive boundary face term to the Jacobian
   */
  void computeJacobDiffBndContribution();
  
  /**
   * compute the terms for the gradient computation for a bnd face
   */
  virtual void computeBndGradTerms(RealMatrix& gradTerm, RealMatrix& ghostGradTerm);
  
  /**
   * compute the terms for the gradient computation for a bnd face
   */
  virtual void computeBndGradTerms2(RealMatrix& gradTerm, RealMatrix& ghostGradTerm);

  /**
   * compute the term for the gradient computation for the cell
   */
  virtual void computeCellGradTerm(RealMatrix& gradTerm);
  
  /**
   * compute the terms for the gradient computation for a face
   */
  virtual void computeFaceGradTerms(RealMatrix& gradTermL, RealMatrix& gradTermR);

  /**
   * compute the terms for the gradient computation for one side of a face,
   * from the states currently held in m_cellStatesFlxPnt2
   */
  virtual void computeFlxPntGradTerm2(RealMatrix& gradTerm);

  /**
   * Rebuild, after a perturbation of the cell states, what the boundary flux
   * depends on: the states and ghost states at the flux points and the compact
   * boundary face gradient.
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

protected: // data
  
  /// builder of cells
  Common::SafePtr< Framework::GeometricEntityPool<CellToFaceGEBuilder> > m_cellBuilder;
  
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
  
  /// perturbed corrections due to the boundary faces for the Jacobian
  std::vector< RealVector> m_pertCorrections;
  
  /// unperturbed updates to the residuals
  RealVector m_resUpdates;
  
  /// extrapolated perturbed states in the flux points of the cell
  std::vector< std::vector< Framework::State* > > m_pertCellStatesFlxPnt;
  
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

  /// influenced flx pnts idx (by perturbation)
  std::vector< CFuint> m_influencedFlxPnts;

  /// Number of influenced flx pnts (by perturbation)
  CFuint m_NbInfluencedFlxPnts;
  
  /// Element shape
  CFGeoShape::Type elemShape;
  
  /// extrapolated states in the flux points of the cell
  std::vector< Framework::State* > m_cellStatesFlxPnt2;

}; // end of class DiffBndCorrectionsRHSJacobFluxReconstruction

//////////////////////////////////////////////////////////////////////////////

 } // namespace FluxReconstructionMethod

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_FluxReconstructionMethod_DiffBndCorrectionsRHSJacobFluxReconstruction_hh
