#ifndef COOLFluiD_SpectralFD_RhsInGivenCellSpectralFD_hh
#define COOLFluiD_SpectralFD_RhsInGivenCellSpectralFD_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/DataSocketSink.hh"

#include "SpectralFD/BaseBndFaceTermComputer.hh"
#include "SpectralFD/BaseFaceTermComputer.hh"
#include "SpectralFD/BaseVolTermComputer.hh"
#include "SpectralFD/BCStateComputer.hh"
#include "SpectralFD/SpectralFDMethodData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace SpectralFD {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represent a command that computes the spectral difference rhs for a given cell.
 * The previous contents of the DataHandle rhsCurrStatesSet is cleared!
 *
 * @author Kris Van den Abeele
 *
 */
class RhsInGivenCellSpectralFD : public SpectralFDMethodCom {
public:

  /**
   * Constructor.
   */
  explicit RhsInGivenCellSpectralFD(const std::string& name);

  /**
   * Destructor.
   */
  virtual ~RhsInGivenCellSpectralFD();

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

  /**
   * Returns the DataSocket's that this command needs as sinks
   * @return a vector of SafePtr with the DataSockets
   */
  std::vector< Common::SafePtr< Framework::BaseDataSocketSink > >
      needsSockets();

protected: // functions

  /**
   * Sets the data for volume term and face term computers.
   */
  void setVolumeAndFaceTermComputersData();

  /**
   * Resizes the residual updates.
   */
  void resizeResUpdates();

  /**
   * clear the residual.
   */
  void clearResidual();

  /**
   * adds the updates to the residual.
   */
  void addUpdatesToResidual();

  /**
   * adds the update to the update coefficients.
   */
  void addUpdateToUpdateCoef();

  /**
   * set and compute volume term computer's data
   * and compute the convective volume terms
   */
  void setVolumeTermDataAndComputeConvVolumeTerm();

  /**
   * set and compute face term computers' data
   * and compute the convective face terms and update coefficient contributions
   */
  void setFaceTermDataAndComputeConvFaceTermsAndUpdateCoefs();

  /**
   * multiply current cell rhs with residual factor
   */
  void multplyRHSWithResFactor();

  /**
   * compute the gradients for the current cell and the diffusive volume terms
   * @pre setVolumeTermDataAndComputeConvVolumeTerm and setFaceTermDataAndComputeConvFaceTermsAndUpdateCoefs
   */
  void computeCurrentCellGradientsAndDiffVolumeTerm();

  /**
   * compute the diffusive face terms and update coefficient contributions
   * @pre computeCurrentCellGradientsAndDiffVolumeTerm
   */
  void computeDiffFaceTermsAndUpdateCoefs();

  /**
   * build the neighbouring cell to the given face
   */
  void buildNeighbourCell(const CFuint currFaceIdx);

  /**
   * set local indexes of the other faces of the current neighbouring cell
   */
  void setOtherFacesLocalIdxs(const CFuint currFaceIdx);

  /**
   * set the face neighbour states for the other faces of the neighbouring cell to the current face
   * @pre setOtherFacesLocalIdxs()
   */
  void setFaceNeighbourStates();

  /**
   * set the data in the neighbouring cell and its faces to compute the gradients
   */
  void setNeighbourCellData();

  /**
   * compute the gradients in the neighbouring cell
   */
  void computeNeighbourCellGradients(const CFuint currFaceIdx);

protected: // data

  /// socket for rhs of current set of states
  Framework::DataSocketSink< CFreal > socket_rhsCurrStatesSet;

  /// socket for update coefficients of current set of states
  Framework::DataSocketSink< CFreal > socket_updateCoeff;

  /// socket for current states set index
  Framework::DataSocketSink< CFint > socket_statesSetIdx;

  /// builder of cells
  Common::SafePtr< Framework::GeometricEntityPool<CellToFaceGEBuilder> > m_cellBuilder;

  /// builder for neighbouring cells
  Common::SafePtr< Framework::GeometricEntityPool<CellToFaceGEBuilder> > m_nghbrCellBuilder;

  /// Strategy that computes the volume terms in a cell
  Common::SafePtr< BaseVolTermComputer > m_volTermComputer;

  /// Strategy that computes the volume terms in a neighbour cell
  Common::SafePtr< BaseVolTermComputer > m_volTermComputerNghbrCell;

  /// face term computers
  std::vector< Common::SafePtr< BaseFaceTermComputer    > > m_faceTermComputers;

  /// boundary face term computers
  std::vector< Common::SafePtr< BaseBndFaceTermComputer > > m_bndFaceTermComputers;

  /// neighbouring cells face term computers
  std::vector< Common::SafePtr< BaseFaceTermComputer    > > m_faceTermComputersNghbrCell;

  /// neighbouring cells boundary face term computers
  std::vector< Common::SafePtr< BaseBndFaceTermComputer > > m_bndFaceTermComputersNghbrCell;

  /// boundary condition state computers
  Common::SafePtr< std::vector< Common::SafePtr< BCStateComputer > > > m_bcStateComputers;

  /// variable for cell
  Framework::GeometricEntity* m_cell;

  /// variable for neighbouring cell
  Framework::GeometricEntity* m_nghbrCell;

  /// vector containing pointers to the states in a cell
  std::vector< Framework::State* >* m_cellStates;

  /// vector containing pointers to the states in the neighbouring cell
  std::vector< Framework::State* >* m_nghbrCellStates;

  /// vector containing pointers to the left and right states with respect to a face
  std::vector< std::vector< std::vector< Framework::State* >* > > m_faceNghbrStates;

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

  /// vector containing pointers to the left and right states with respect to a face of the neighbouring cell
  std::vector< std::vector< std::vector< Framework::State* >* > > m_faceNghbrStatesNghbrCell;

  /// pointer to booleans telling whether a face is on the boundary (neighbour cell builder)
  Common::SafePtr< std::vector< bool > > m_isFaceOnBoundaryNghbrCell;

  /// pointer to neighbouring cell side vector (neighbour cell builder)
  Common::SafePtr< std::vector< CFuint > > m_nghbrCellSideNghbrCell;

  /// pointer to current cell side vector (neighbour cell builder)
  Common::SafePtr< std::vector< CFuint > > m_currCellSideNghbrCell;

  /// pointer to orientation vector (neighbour cell builder)
  Common::SafePtr< std::vector< CFuint > > m_faceOrientsNghbrCell;

  /// pointer to BC index vector (neighbour cell builder)
  Common::SafePtr< std::vector< CFuint > > m_faceBCIdxNghbrCell;

  /// updates to the residuals
  RealVector m_resUpdates;

  /// face updates to the residuals
  std::vector< RealVector > m_faceResUpdates;

  /// update for the wave speed in the neighbouring cell
  CFreal m_updateCoefUpd;

  /// face updates for the wave speed in the neighbouring cells
  std::vector< CFreal > m_updateCoefFaceUpd;

  /// number of equations in the physical model
  CFuint m_nbrEqs;

  /// number of dimensions in the physical model
  CFuint m_dim;

  /// index of element type
  CFuint m_iElemType;

  /// current iteration
  CFuint m_currIter;

  /// boolean telling whether to recompute the update coefficients
  bool m_computeUpdateCoef;

  /// updates to the gradients
  std::vector< std::vector< std::vector< RealVector > > > m_gradUpdates;

  /// current cell gradients
  std::vector< std::vector< RealVector >* > m_currCellGrads;

  /// neighbour cell gradients
  std::vector< std::vector< RealVector >* > m_nghbCellGrads;

  /// Jacobian determinants
  std::valarray<CFreal> m_solJacobDet;

  /// Jacobian determinants of the neighbour cell
  std::valarray<CFreal> m_solJacobDetNghbrCell;

  /// solution point mapped coordinates
  Common::SafePtr< std::vector< RealVector > > m_solPntsLocalCoords;

  /// neighbouring cell local indexes of the other faces (not belonging to the current cell)
  std::vector< CFuint > m_otherFaceLocalIdxs;

}; // class RhsInGivenCellSpectralFD

//////////////////////////////////////////////////////////////////////////////

  } // namespace SpectralFD

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_SpectralFD_RhsInGivenCellSpectralFD_hh
