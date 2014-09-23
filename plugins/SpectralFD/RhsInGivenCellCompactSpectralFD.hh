#ifndef COOLFluiD_SpectralFD_RhsInGivenCellCompactSpectralFD_hh
#define COOLFluiD_SpectralFD_RhsInGivenCellCompactSpectralFD_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/DataSocketSink.hh"

#include "SpectralFD/CompactBndFaceTermComputer.hh"
#include "SpectralFD/CompactFaceTermComputer.hh"
#include "SpectralFD/CompactVolTermComputer.hh"
#include "SpectralFD/BCStateComputer.hh"
#include "SpectralFD/SpectralFDMethodData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace SpectralFD {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represent a command that computes the spectral difference rhs
 * for a given cell with a compact approach.
 * The previous contents of the DataHandle rhsCurrStatesSet is cleared!
 *
 * @author Kris Van den Abeele
 *
 */
class RhsInGivenCellCompactSpectralFD : public SpectralFDMethodCom {
public:

  /**
   * Constructor.
   */
  explicit RhsInGivenCellCompactSpectralFD(const std::string& name);

  /**
   * Destructor.
   */
  virtual ~RhsInGivenCellCompactSpectralFD();

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
   * set and compute face term computers' data
   * and compute the convective face terms and update coefficient contributions
   */
  void setFaceTermDataAndComputeFaceTermsAndUpdateCoefs();

  /**
   * set and compute volume term computer's data
   * and compute the volume terms
   * @pre setFaceTermDataAndComputeFaceTermsAndUpdateCoefs()
   */
  void setVolumeTermDataAndComputeCellGradientsAndVolumeTerm();

  /**
   * multiply current cell rhs with residual factor
   */
  void multplyRHSWithResFactor();

protected: // data

  /// socket for rhs of current set of states
  Framework::DataSocketSink< CFreal > socket_rhsCurrStatesSet;

  /// socket for update coefficients of current set of states
  Framework::DataSocketSink< CFreal > socket_updateCoeff;

  /// socket for current states set index
  Framework::DataSocketSink< CFint > socket_statesSetIdx;

  /// builder of cells
  Common::SafePtr< Framework::GeometricEntityPool<CellToFaceGEBuilder> > m_cellBuilder;

  /// Strategy that computes the volume terms in a cell
  Common::SafePtr< CompactVolTermComputer > m_volTermComputer;

  /// face term computers
  std::vector< Common::SafePtr< CompactFaceTermComputer    > > m_faceTermComputers;

  /// boundary face term computers
  std::vector< Common::SafePtr< CompactBndFaceTermComputer > > m_bndFaceTermComputers;

  /// boundary condition state computers
  Common::SafePtr< std::vector< Common::SafePtr< BCStateComputer > > > m_bcStateComputers;

  /// variable for cell
  Framework::GeometricEntity* m_cell;

  /// variable for faces
  const std::vector< Framework::GeometricEntity* >* m_faces;

  /// vector containing pointers to the states in a cell
  std::vector< Framework::State* >* m_cellStates;

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

  /// updates to the residuals
  RealVector m_resUpdates;

  /// diffusive updates to the residuals
  RealVector m_diffResUpdates;

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

  /// Jacobian determinants
  std::valarray<CFreal> m_solJacobDet;

  /// solution point mapped coordinates
  Common::SafePtr< std::vector< RealVector > > m_solPntsLocalCoords;

  /// neighbouring cell local indexes of the other faces (not belonging to the current cell)
  std::vector< CFuint > m_otherFaceLocalIdxs;

}; // class RhsInGivenCellCompactSpectralFD

//////////////////////////////////////////////////////////////////////////////

  } // namespace SpectralFD

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_SpectralFD_RhsInGivenCellCompactSpectralFD_hh
