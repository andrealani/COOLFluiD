#ifndef COOLFluiD_SpectralFD_DiffVolTermRHSJacobSpectralFD_hh
#define COOLFluiD_SpectralFD_DiffVolTermRHSJacobSpectralFD_hh

//////////////////////////////////////////////////////////////////////////////

#include "SpectralFD/BaseBndFaceTermComputer.hh"
#include "SpectralFD/BaseFaceTermComputer.hh"
#include "SpectralFD/BCStateComputer.hh"
#include "SpectralFD/DiffVolTermRHSSpectralFD.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {
    class BlockAccumulator;
  }

  namespace SpectralFD {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a command that computes contribution of the volume terms for the
 * spectral finite difference schemes for diffusion terms, for both the RHS and the Jacobian
 *
 * @author Kris Van den Abeele
 *
 */
class DiffVolTermRHSJacobSpectralFD : public DiffVolTermRHSSpectralFD {
public:

  /**
   * Constructor.
   */
  explicit DiffVolTermRHSJacobSpectralFD(const std::string& name);

  /**
   * Destructor.
   */
  virtual ~DiffVolTermRHSJacobSpectralFD();

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
   * set the face neighbour states
   * @pre m_faces is set
   */
  void setFaceNeighbourStates();

  /**
   * set data related to the cell (CV volumes and face computer data)
   * @pre setFaceNeighbourStates()
   */
  void setCellData();

  /**
   * back up and reconstruct variables the face flux points
   */
  void backupAndReconstructFacePhysVars(const CFuint iVar);

  /**
   * restore variables the face flux points
   */
  void restoreFacePhysVars(const CFuint iVar);

  /**
   * recompute the cell gradients from the current cell and the neighbouring cells,
   * after perturbation
   * @pre m_volTermComputer->backupAndReconstructPhysVar
   * @pre m_faceTermComputers->backupAndReconstructPhysVar
   * @pre m_bndFaceTermComputers->backupAndReconstructPhysVar
   */
  void computePerturbedGradients();

  /**
   * compute the contribution of the diffusive volume term to the Jacobian
   */
  void computeJacobDiffVolTerm();

protected: // data

  /// builder of cells
  Common::SafePtr< Framework::GeometricEntityPool<CellToFaceGEBuilder> > m_cellBuilder;

  /// face term computers
  std::vector< Common::SafePtr< BaseFaceTermComputer    > > m_faceTermComputers;

  /// boundary face term computers
  std::vector< Common::SafePtr< BaseBndFaceTermComputer > > m_bndFaceTermComputers;

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

  /// derivative of update to the residuals
  RealVector m_derivResUpdates;

  /// updates to the gradients
  std::vector< std::vector< std::vector< RealVector > > > m_gradUpdates;

  /// perturbed gradients
  std::vector< std::vector< RealVector >* > m_pertGrads;

  /// Jacobian determinants
  std::valarray<CFreal> m_solJacobDet;

  /// solution point mapped coordinates
  Common::SafePtr< std::vector< RealVector > > m_solPntsLocalCoords;

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

  /// number of dimensions in the physical model
  CFuint m_dim;

}; // class DiffVolTermRHSJacobSpectralFD

//////////////////////////////////////////////////////////////////////////////

  } // namespace SpectralFD

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_SpectralFD_DiffVolTermRHSJacobSpectralFD_hh
