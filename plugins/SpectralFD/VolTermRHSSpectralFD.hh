#ifndef COOLFluiD_SpectralFD_VolTermRHSSpectralFD_hh
#define COOLFluiD_SpectralFD_VolTermRHSSpectralFD_hh

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
 * This class represents a basic command that computes the volume terms for the
 * spectral finite difference schemes for both convective and diffusive terms
 * if the full scheme is compact.
 *
 * @author Kris Van den Abeele
 *
 */
class VolTermRHSSpectralFD : public SpectralFDMethodCom {
public:

  /**
   * Constructor.
   */
  explicit VolTermRHSSpectralFD(const std::string& name);

  /**
   * Destructor.
   */
  virtual ~VolTermRHSSpectralFD();

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

  /**
   * Returns the DataSocket's that this command needs as sinks
   * @return a vector of SafePtr with the DataSockets
   */
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();

protected: // functions

  /// set the required data in the volume and face term computers
  void setVolumeAndFaceTermComputersData();

  /// resize the residual and gradient updates corresponding to the current element type
  void resizeResAndGradUpdates();

  /// compute the required data in the volume and face term for the current cell and faces
  void computeCellVolumeAndFaceTermData();

  /// add the residual updates to the RHS
  void addUpdatesToResidual();

  /// compute cell gradients
  void computeAndReconstructGradients();

protected: // data

  /// storage of the rhs
  Framework::DataSocketSink<CFreal> socket_rhs;

  /// builder of cells
  Common::SafePtr<Framework::GeometricEntityPool<CellToFaceGEBuilder> > m_cellBuilder;

  /// Strategy that computes the volume terms in a cell
  Common::SafePtr< CompactVolTermComputer > m_volTermComputer;

  /// face term computers
  std::vector< Common::SafePtr< CompactFaceTermComputer > > m_faceTermComputers;

  /// boundary face term computers
  std::vector< Common::SafePtr< CompactBndFaceTermComputer > > m_bndFaceTermComputers;

  /// boundary condition state computers
  Common::SafePtr< std::vector< Common::SafePtr< BCStateComputer > > > m_bcStateComputers;

  /// index of element type
  CFuint m_iElemType;

  /// variable for cell
  Framework::GeometricEntity* m_cell;

  /// vector containing pointers to the states in a cell
  std::vector< Framework::State* >* m_cellStates;

  /// current cell gradients
  std::vector< std::vector< RealVector >* > m_cellGrads;

  /// vector containing pointers to the left and right states with respect to a face
  std::vector< std::vector< std::vector< Framework::State* >* > > m_faceNghbrStates;

  /// pointer to booleans telling whether a face is on the boundary
  Common::SafePtr< std::vector< bool > > m_isFaceOnBoundary;

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

  /// updates to the gradients
  std::vector< std::vector< std::vector< RealVector > > > m_gradUpdates;

  /// boundary face term updates to the gradients
  std::vector< std::vector< std::vector< RealVector > > > m_bndFaceGradUpdates;

  /// Jacobian determinants
  std::valarray<CFreal> m_solJacobDet;

  /// solution point mapped coordinates
  Common::SafePtr< std::vector< RealVector > > m_solPntsLocalCoords;

  /// number of equations in the physical model
  CFuint m_nbrEqs;

  /// number of dimensions in the physical model
  CFuint m_dim;

}; // class VolTermRHSSpectralFD

//////////////////////////////////////////////////////////////////////////////

  } // namespace SpectralFD

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_SpectralFD_VolTermRHSSpectralFD_hh
