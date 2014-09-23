#ifndef COOLFluiD_SpectralFD_ConvFaceTermRHSSpectralFD_hh
#define COOLFluiD_SpectralFD_ConvFaceTermRHSSpectralFD_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/DataSocketSink.hh"

#include "SpectralFD/BaseFaceTermComputer.hh"
#include "SpectralFD/SpectralFDMethodData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace SpectralFD {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represent a basic command that computes the face terms for the
 * spectral finite difference schemes for convective terms
 *
 * @author Kris Van den Abeele
 *
 */
class ConvFaceTermRHSSpectralFD : public SpectralFDMethodCom {
public:

  /**
   * Constructor.
   */
  explicit ConvFaceTermRHSSpectralFD(const std::string& name);

  /**
   * Destructor.
   */
  virtual ~ConvFaceTermRHSSpectralFD();

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

  /**
   * Returns the DataSocket's that this command needs as sinks
   * @return a vector of SafePtr with the DataSockets
   */
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();

protected: // functions

  /// set face term data
  void setFaceTermData();

  /// add the residual updates to the RHS
  void updateRHS();

  /// add the updates to the wave speed
  void updateWaveSpeed();

  /// compute the face term contribution to the gradients
  void computeGradientFaceTerm();

  /// add updates to gradients
  void addGradFaceTerms();

protected: // data

  /// storage of the rhs
  Framework::DataSocketSink<CFreal > socket_rhs;

  /// socket for updateCoeff
  /// denominators of the coefficients for the update
  Framework::DataSocketSink<CFreal > socket_updateCoeff;

  /// socket for gradients
  Framework::DataSocketSink< std::vector< RealVector > > socket_gradients;

  /// builder of faces
  Common::SafePtr<Framework::GeometricEntityPool<Framework::FaceToCellGEBuilder> > m_faceBuilder;

  /// Strategy that computes the face terms for a face
  Common::SafePtr< BaseFaceTermComputer > m_faceTermComputer;

  /// variable for current face
  Framework::GeometricEntity* m_face;

  /// variable for current neighbouring cells
  std::vector< Framework::GeometricEntity* > m_cells;

  /// variable for current face orientation
  CFuint m_orient;

  /// variable for normal vectors
  RealVector m_normal;

  /// variable for the states in the left and right cell
  std::vector< std::vector< Framework::State* >* > m_states;

  /// residual updates
  std::vector< RealVector > m_resUpdates;

  /// updates for the wave speed
  std::vector< CFreal > m_waveSpeedUpd;

  /// gradient updates
  std::vector< std::vector< std::vector< RealVector > > > m_gradUpdates;

  /// number of equations in the physical model
  CFuint m_nbrEqs;

}; // class ConvFaceTermRHSSpectralFD

//////////////////////////////////////////////////////////////////////////////

  } // namespace SpectralFD

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_SpectralFD_ComputeConvVolTerms_hh
