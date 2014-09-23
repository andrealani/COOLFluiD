#ifndef COOLFluiD_SpectralFV_DiffBndFaceTermRHSJacobSpectralFV_hh
#define COOLFluiD_SpectralFV_DiffBndFaceTermRHSJacobSpectralFV_hh

//////////////////////////////////////////////////////////////////////////////

#include "SpectralFV/BaseFaceTermComputer.hh"
#include "SpectralFV/BaseVolTermComputer.hh"
#include "SpectralFV/BCStateComputer.hh"
#include "SpectralFV/DiffBndFaceTermRHSSpectralFV.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {
    class BlockAccumulator;
  }

    namespace SpectralFV {

//////////////////////////////////////////////////////////////////////////////


/**
  * This class represents a command that computes contribution of the boundary face terms for the
  * spectral finite volume schemes for diffusion terms, for both the RHS and the Jacobian
  *
  * @author Kris Van Den Abeele
  *
  */
class DiffBndFaceTermRHSJacobSpectralFV : public DiffBndFaceTermRHSSpectralFV {

public:
  typedef Framework::BaseMethodCommandProvider<
      SpectralFVMethodData,DiffBndFaceTermRHSJacobSpectralFV > PROVIDER;

public:

  /**
   * Constructor
   */
  DiffBndFaceTermRHSJacobSpectralFV(const std::string& name);

  /**
   * Default destructor
   */
  virtual ~DiffBndFaceTermRHSJacobSpectralFV();

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
   * set the data related to the neighbouring cell
   */
  void setCellData();

  /**
   * back up and reconstruct variables the other face flux points and in the cell flux points
   */
  void backupAndReconstructOtherFacesAndCellPhysVars(const CFuint iVar);

  /**
   * restore variables the other face flux points and in the cell flux points
   */
  void restoreOtherFacesAndCellPhysVars(const CFuint iVar);

  /**
   * recompute the cell gradients from the current cell and the neighbouring cells,
   * after perturbation
   * @pre m_volTermComputer->backupAndReconstructPhysVar
   * @pre m_faceTermComputers->backupAndReconstructPhysVar
   * @pre m_bndFaceTermComputers->backupAndReconstructPhysVar
   */
  void computePerturbedGradients();

  /**
   * compute the contribution of the diffusive boundary face term to the Jacobian
   */
  void computeJacobDiffBndFaceTerm();

protected: // data

  /// builder of cells
  Common::SafePtr< Framework::GeometricEntityPool<CellToFaceGEBuilder> > m_cellBuilder;

  /// volume term computer
  Common::SafePtr< BaseVolTermComputer > m_volTermComputer;

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

  /// volume fractions of CVs
  Common::SafePtr< std::vector< CFreal > > m_invVolFracCVs;

  /// CV-CV connectivity through SV face
  Common::SafePtr< std::vector< std::vector< std::vector< CFuint > > > > m_cvCVConnSVFace;

  /// perturbed updates to the residuals
  RealVector m_pertResUpdates;

  /// derivative of update to one CV-residual
  RealVector m_derivResUpdates;

  /// updates to the gradients
  std::vector< std::vector< RealVector > > m_gradUpdates;

  /// perturbed boundary cell gradients
  std::vector< std::vector< RealVector >* > m_pertGrads;

  /// CV inverse volumes of the boundary cell
  std::vector< CFreal > m_cvInvVols;

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

}; // end of class DiffBndFaceTermRHSJacobSpectralFV

//////////////////////////////////////////////////////////////////////////////

 } // namespace SpectralFV

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_SpectralFV_DiffBndFaceTermRHSJacobSpectralFV_hh
