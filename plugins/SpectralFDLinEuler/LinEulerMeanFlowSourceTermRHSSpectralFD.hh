#ifndef COOLFluiD_SpectralFD_LinEulerMeanFlowSourceTermRHSSpectralFD_hh
#define COOLFluiD_SpectralFD_LinEulerMeanFlowSourceTermRHSSpectralFD_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/DataSocketSink.hh"

#include "LinEuler/LinEulerVarSet.hh"

#include "SpectralFD/BaseBndFaceTermComputer.hh"
#include "SpectralFD/BaseFaceTermComputer.hh"
#include "SpectralFD/BaseVolTermComputer.hh"
#include "SpectralFD/BCStateComputer.hh"
#include "SpectralFD/SpectralFDMethodData.hh"
#include "SpectralFD/StdSourceTerm.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace SpectralFD {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a command that computes the contribution to the RHS of
 * source terms associated to the mean flow for the linearized Euler equations.
 *
 * @author Kris Van den Abeele
 *
 */
class LinEulerMeanFlowSourceTermRHSSpectralFD : public StdSourceTerm {
public:

  /**
   * Constructor.
   */
  explicit LinEulerMeanFlowSourceTermRHSSpectralFD(const std::string& name);

  /**
   * Destructor.
   */
  virtual ~LinEulerMeanFlowSourceTermRHSSpectralFD();

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

protected: // functions

  /// get data required for source term computation
  void getSourceTermData();

  /// add the source term
  void addSourceTerm();

  /// set the required data in the volume and face term computers
  void setVolumeAndFaceTermComputersData();

  /// resize the residual and gradient updates corresponding to the current element type
  void resizeGradUpdates();

  /// compute the required data in the volume and face term for the current cell and faces
  void computeCellVolumeAndFaceTermData();

  /// compute cell gradients
  void computeExtraVarsGradients();

protected: // data

  /// Linearized Euler variable set
  Common::SafePtr< Physics::LinearizedEuler::LinEulerVarSet > m_linEulerVarSet;

  /// builder of cells
  Common::SafePtr<Framework::GeometricEntityPool<CellToFaceGEBuilder> > m_cellBuilder;

  /// Strategy that computes the volume terms in a cell
  Common::SafePtr< BaseVolTermComputer > m_volTermComputer;

  /// face term computers
  std::vector< Common::SafePtr< BaseFaceTermComputer > > m_faceTermComputers;

  /// boundary face term computers
  std::vector< Common::SafePtr< BaseBndFaceTermComputer > > m_bndFaceTermComputers;

  /// boundary condition state computers
  Common::SafePtr< std::vector< Common::SafePtr< BCStateComputer > > > m_bcStateComputers;

  /// index of element type
  CFuint m_iElemType;

  /// current cell gradients
  std::vector< std::vector< RealVector >* > m_cellExtraVarGrads;

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

  /// updates to the gradients
  std::vector< std::vector< std::vector< RealVector > > > m_gradUpdates;

  /// number of dimensions in the physical model
  CFuint m_dim;

  /// the source term for one state
  RealVector m_srcTerm;

}; // class LinEulerMeanFlowSourceTermRHSSpectralFD

//////////////////////////////////////////////////////////////////////////////

  } // namespace SpectralFD

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_SpectralFD_LinEulerMeanFlowSourceTermRHSSpectralFD_hh
