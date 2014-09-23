#ifndef COOLFluiD_SpectralFD_DiffFaceTermRHSSpectralFD_hh
#define COOLFluiD_SpectralFD_DiffFaceTermRHSSpectralFD_hh

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
 * spectral finite difference schemes for diffusive terms
 *
 * @author Kris Van den Abeele
 *
 */
class DiffFaceTermRHSSpectralFD : public SpectralFDMethodCom {
public:

  /**
   * Constructor.
   */
  explicit DiffFaceTermRHSSpectralFD(const std::string& name);

  /**
   * Destructor.
   */
  virtual ~DiffFaceTermRHSSpectralFD();

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

  /// set the gradients of left and right cell
  void setGradients();

  /// add the residual updates to the RHS
  void updateRHS();

  /// add the contributions to the update coefficient
  void addUpdateCoeffContributions();

protected: // data

  /// storage of the rhs
  Framework::DataSocketSink<CFreal > socket_rhs;

  /// socket for updateCoeff
  Framework::DataSocketSink<CFreal > socket_updateCoeff;

  /// socket for gradients
  Framework::DataSocketSink< std::vector< RealVector > > socket_gradients;

  /// Strategy that computes the face terms for a face
  Common::SafePtr< BaseFaceTermComputer > m_faceTermComputer;

  /// builder of faces
  Common::SafePtr<Framework::GeometricEntityPool<Framework::FaceToCellGEBuilder> > m_faceBuilder;

  /// current face
  Framework::GeometricEntity* m_face;

  /// variable for current neighbouring cells
  std::vector< Framework::GeometricEntity* > m_cells;

  /// variable for current face orientation
  CFuint m_orient;

  /// variable for the states in the left and right cell
  std::vector< std::vector< Framework::State* >* > m_states;

  /// variable for the gradients in the left and right cell
  std::vector< std::vector< std::vector< RealVector >* > > m_grads;

  /// residual updates
  std::vector< RealVector > m_resUpdates;

  /// contributions to the update coefficients
  std::vector< CFreal > m_updateCoefContr;

  /// number of equations in the physical model
  CFuint m_nbrEqs;

}; // class DiffFaceTermRHSSpectralFD

//////////////////////////////////////////////////////////////////////////////

  } // namespace SpectralFD

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_SpectralFD_ComputeConvVolTerms_hh
