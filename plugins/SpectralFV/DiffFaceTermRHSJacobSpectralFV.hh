#ifndef COOLFluiD_SpectralFV_DiffFaceTermRHSJacobSpectralFV_hh
#define COOLFluiD_SpectralFV_DiffFaceTermRHSJacobSpectralFV_hh

//////////////////////////////////////////////////////////////////////////////

#include "SpectralFV/BaseBndFaceTermComputer.hh"
#include "SpectralFV/BaseVolTermComputer.hh"
#include "SpectralFV/BCStateComputer.hh"
#include "SpectralFV/DiffFaceTermRHSSpectralFV.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {
    class BlockAccumulator;
  }

  namespace SpectralFV {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a command that computes contribution of the face terms for the
 * spectral finite volume schemes for diffusion terms, for both the RHS and the Jacobian
 *
 * @author Kris Van den Abeele
 *
 */
class DiffFaceTermRHSJacobSpectralFV : public DiffFaceTermRHSSpectralFV {
public:

  /**
   * Constructor.
   */
  explicit DiffFaceTermRHSJacobSpectralFV(const std::string& name);

  /**
   * Destructor.
   */
  virtual ~DiffFaceTermRHSJacobSpectralFV();

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
   * set face term data
   */
  void setFaceTermData();

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
   * set the face neighbour gradients
   * @pre setOtherFacesLocalIdxs()
   * @pre setFaceNeighbourStates()
   */
  void setFaceNeighbourGradients();

  /**
   * set the data related to the neighbouring cells
   * @pre setOtherFacesLocalIdxs()
   * @pre setFaceNeighbourStates()
   * @pre setFaceNeighbourGradients()
   */
  void setCellsData();

  /**
   * compute the unperturbed gradients minus the current face term
   * @pre m_faceTermComputer->reconstructFluxPntsStates
   * @pre setCellsData()
   */
  void computeCellGradsMinusFaceTerm();

  /**
   * compute the unperturbed gradients minus the other face terms
   * @pre m_faceTermComputers->reconstructFluxPntsStates
   * @pre setCellsData()
   */
  void computeCellGradsMinusOtherFaceTerms(const CFuint side);

  /**
   * compute the unperturbed cell diffusive residuals
   * @pre m_faceTermComputers->computeDiffFaceTermAndUpdateCoefContributions
   * @pre setCellsData()
   */
  void computeUnpertCellDiffResiduals();

  /**
   * back up and reconstruct variables the other face flux points and in the cell flux points,
   * for the left or the right side
   */
  void backupAndReconstructOtherFacesAndCellPhysVars(const CFuint side, const CFuint iVar);

  /**
   * restore variables the other face flux points and in the cell flux points,
   * for the left or the right side
   */
  void restoreOtherFacesAndCellPhysVars(const CFuint side, const CFuint iVar);

  /**
   * recompute the cell gradients from the current cell and the neighbouring cells solutions,
   * after perturbation
   * @pre setCellsData()
   * @pre backupAndReconstructOtherFacesAndCellPhysVars()
   */
  void computePerturbedGradients(const CFuint side);

  /**
   * recompute the cell gradients from the current cell and the neighbouring cells,
   * after perturbation
   * @pre computePerturbedGradients() (for the other cell)
   * @pre computeCellGradsMinusFaceTerm()
   */
  void computePertGradsFromFaceTerm(const CFuint side);

  /**
   * recompute the cell gradients, from one of the other face terms
   * @pre computeCellGradsMinusOtherFaceTerm()
   */
  void computePertGradsFromOtherFaceTerm(const CFuint side, const CFuint iFace);

  /**
   * reconstruct the gradients in the flux points,
   * after the perturbed cell gradients have been computed
   * @pre computePerturbedGradients(), computePertGradsFromFaceTerm()
   */
  void reconstructOtherFacesAndCellGradients(const CFuint side);

  /**
   * compute the perturbed cell diffusive residuals for one cell
   * @pre m_faceTermComputer->computeDiffFaceTerm
   * @pre backupAndReconstructOtherFacesAndCellPhysVars(
   * @pre reconstructOtherFacesAndCellGradients()
   */
  void computePertCellDiffResiduals(const CFuint side);

  /**
   * compute the contribution of the diffusive face term to both Jacobians
   */
  void computeBothJacobsDiffFaceTerm();

  /**
   * compute the contribution of the diffusive face term to one Jacobians
   */
  void computeOneJacobDiffFaceTerm(const CFuint side);

protected: // data

  /// builder of cells
  std::vector< Common::SafePtr< Framework::GeometricEntityPool<CellToFaceGEBuilder> > > m_cellBuilders;

  /// pointer to the linear system solver
  Common::SafePtr<Framework::LinearSystemSolver> m_lss;

  /// pointer to the numerical Jacobian computer
  Common::SafePtr<Framework::NumericalJacobian> m_numJacob;

  /// accumulator for LSSMatrix
  std::auto_ptr<Framework::BlockAccumulator> m_acc;

  /// single cell accumulator for LSSMatrix
  std::auto_ptr<Framework::BlockAccumulator> m_accSC;

  /// volume term computers
  std::vector< Common::SafePtr< BaseVolTermComputer > > m_volTermComputers;

  /// face term computers
  std::vector< std::vector< Common::SafePtr< BaseFaceTermComputer    > > > m_faceTermComputers;

  /// boundary face term computers
  std::vector< std::vector< Common::SafePtr< BaseBndFaceTermComputer > > > m_bndFaceTermComputers;

  /// variable for faces
  std::vector< const std::vector< Framework::GeometricEntity* >* > m_faces;

  /// vector containing pointers to the left and right states with respect to a face
  std::vector< std::vector< std::vector< std::vector< Framework::State* >* > > > m_faceNghbrStates;

  /// vector containing pointers to the left and right gradients with respect to a face
  std::vector< std::vector< std::vector< std::vector< std::vector< RealVector >* > > > > m_faceNghbrGrads;

  /// volume fractions of CVs
  Common::SafePtr< std::vector< CFreal > > m_invVolFracCVs;

  /// boundary face CV connectivity on SV face
  Common::SafePtr< std::vector< std::vector< CFuint > > > m_svFaceCVConn;

  /// perturbed updates to the residuals
  RealVector m_pertResUpdates;

  /// derivative of update to one CV-residual
  RealVector m_derivResUpdates;

  /// updates to the gradients
  std::vector< std::vector< RealVector > > m_gradUpdates;

  /// perturbed left and right cell gradients
  std::vector< std::vector< std::vector< RealVector >* > > m_pertGrads;

  /// left and right cell gradients minus current face term
  std::vector< std::vector< std::vector< RealVector > > > m_cellGradsMinusFaceTerm;

  /// left and right cell gradients minus other face terms
  std::vector< std::vector< std::vector< std::vector< RealVector > > > > m_cellGradsMinusOtherFaceTerm;

  /// unperturbed diffusive residuals
  std::vector< RealVector > m_unpertCellDiffRes;

  /// unperturbed diffusive residuals
  RealVector m_pertCellDiffRes;

  /// unperturbed diffusive residuals
  RealVector m_derivCellDiffRes;

  /// CV inverse volumes of left and right cell
  std::vector< std::vector< CFreal > > m_cvInvVols;

  /// cell local indexes of the other faces (not the face itself)
  std::vector< std::vector< CFuint > > m_otherFaceLocalIdxs;

  /// pointer to booleans telling whether a face is on the boundary
  std::vector< Common::SafePtr< std::vector< bool > > > m_isFaceOnBoundary;

  /// pointer to neighbouring cell side vector
  std::vector< Common::SafePtr< std::vector< CFuint > > > m_nghbrCellSide;

  /// pointer to current cell side vector
  std::vector< Common::SafePtr< std::vector< CFuint > > > m_currCellSide;

  /// pointer to orientation vector
  std::vector< Common::SafePtr< std::vector< CFuint > > > m_faceOrients;

  /// pointer to BC index vector
  std::vector< Common::SafePtr< std::vector< CFuint > > > m_faceBCIdx;

  /// boundary condition state computers
  Common::SafePtr< std::vector< Common::SafePtr< BCStateComputer > > > m_bcStateComputers;

}; // class DiffFaceTermRHSJacobSpectralFV

//////////////////////////////////////////////////////////////////////////////

  } // namespace SpectralFV

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_SpectralFV_ComputeConvVolTerms_hh
