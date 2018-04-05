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
   * set the local indexes of the other faces (not the current boundary face)
   * @pre m_faces is set
   */
  void setOtherFacesLocalIdxs();

  /**
   * set the face neighbour states
   * @pre setOtherFacesLocalIdxs()
   */
  void setFaceNeighbourStates();
  
  /**
   * compute the contribution of the diffusive boundary face term to the Jacobian
   */
  void computeJacobDiffBndContribution();
  
  /**
   * recompute the cell gradients from the current cell and the neighbouring cells,
   * after perturbation
   */
  void computePerturbedGradients();
  
  /**
   * compute the terms for the gradient computation for a bnd face
   */
  virtual void computeBndGradTerms(RealMatrix& gradTerm, RealMatrix& ghostGradTerm);
  
  /**
   * compute the term for the gradient computation for the cell
   */
  virtual void computeCellGradTerm(RealMatrix& gradTerm);
  
  /**
   * compute the terms for the gradient computation for a face
   */
  virtual void computeFaceGradTerms(RealMatrix& gradTermL, RealMatrix& gradTermR);

protected: // data
  
  /// builder of cells
  Common::SafePtr< Framework::GeometricEntityPool<CellToFaceGEBuilder> > m_cellBuilder;
  
  /// pointer to the linear system solver
  Common::SafePtr<Framework::LinearSystemSolver> m_lss;

  /// pointer to the numerical Jacobian computer
  Common::SafePtr<Framework::NumericalJacobian> m_numJacob;

  /// accumulator for LSSMatrix
  std::auto_ptr<Framework::BlockAccumulator> m_acc;

  /// variable for faces
  const std::vector< Framework::GeometricEntity* >* m_faces;

  /// vector containing pointers to the left and right states with respect to a face
  std::vector< std::vector< std::vector< Framework::State* >* > > m_faceNghbrStates;

  /// perturbed updates to the residuals
  RealVector m_pertResUpdates;

  /// derivative of update to one CV-residual
  RealVector m_derivResUpdates;

  /// updates to the gradients
  std::vector< std::vector< std::vector< RealVector > > > m_gradUpdates;

  /// perturbed boundary cell gradients
  std::vector< std::vector< RealVector >* > m_pertGrads;

  /// Jacobian determinants
  std::valarray< CFreal > m_solJacobDet;

  /// cell local indexes of the other faces (not the boundary face itself)
  std::vector< CFuint > m_otherFaceLocalIdxs;

  /// pointer to booleans telling whether a face is on the boundary
  Common::SafePtr< std::vector< bool > > m_isFaceOnBoundary;

  /// pointer to neighbouring cell side vector
  Common::SafePtr< std::vector< CFuint > > m_nghbrCellSide;

  /// pointer to current cell side vector
  Common::SafePtr< std::vector< CFuint > > m_currCellSide;

  /// pointer to orientation vector
  Common::SafePtr< std::vector< CFuint > > m_faceOrients;

  /// pointer to BC index vector
  Common::SafePtr< std::vector< CFuint > > m_faceBCIdx;

  /// boundary condition state computers
  Common::SafePtr< std::vector< Common::SafePtr< BCStateComputer > > > m_bcStateComputers;
  
  /// perturbed corrections due to the boundary faces for the Jacobian
  std::vector< RealVector> m_pertCorrections;
  
  /// unperturbed updates to the residuals
  RealVector m_resUpdates;
  
  /// backup of the gradients in the neighbouring cell
  std::vector< std::vector< RealVector > > m_cellGradsBackUp;
  
  /// coefs to compute the derivative of the states in the sol pnts
  Common::SafePtr< std::vector< std::vector< std::vector< CFreal > > > > m_solPolyDerivAtSolPnts;
  
  /// flx pnt - face connectivity per orient
  Common::SafePtr< std::vector< std::vector< std::vector< CFuint > > > > m_faceFlxPntConnPerOrient;
  
  /// extrapolated perturbed states in the flux points of the cell
  std::vector< std::vector< Framework::State* > > m_pertCellStatesFlxPnt;
  
  /// local cell face - mapped coordinate direction per orientation
  Common::SafePtr< std::vector< std::vector< CFint > > > m_faceMappedCoordDirPO;

}; // end of class DiffBndCorrectionsRHSJacobFluxReconstruction

//////////////////////////////////////////////////////////////////////////////

 } // namespace FluxReconstructionMethod

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_FluxReconstructionMethod_DiffBndCorrectionsRHSJacobFluxReconstruction_hh
